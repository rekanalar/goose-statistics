// The kernel density estimate of a sample of points.
//
// The estimate is a sum of Gaussian kernels, one for each point of the
// sample. The kernel of a point has the covariance of the sample, scaled
// by the square of a bandwidth factor. The factor of Scott is
// neff^(-1/(d+4)). The factor of Silverman is (neff (d+2)/4)^(-1/(d+4)).
// The neff is the effective number of points. It is one over the sum of
// the squares of the weights of the points.
//
// The data holds the points row by row. Each point holds d values. The
// evaluation points use the same layout.
namespace statistics;

import std;
import math;

// The estimate of a sample of points. The buffer holds the factor of the
// covariance, then the points, then the weights. The factor holds d*d
// values, the points hold n*d values and the weights hold n values.
struct Kde {
    n: i64,             // the number of points
    d: i64,             // the number of dimensions
    factor: f64,        // the bandwidth factor
    neff: f64,          // the effective number of points
    log_det: f64,       // the log of the determinant of the factor
    buf: f64[>..],      // the factor, the points and the weights
}

// An estimate with a fault code. Every field is empty.
fn kde_fault(code: i64) -> Kde {
    return Kde { n: 0, d: 0, factor: nan(), neff: nan(), log_det: nan(), buf: [] };
}

// The value of the factor of the covariance at the row i and the column j.
fn kde_cho(k: Kde&, i: i64, j: i64) -> f64 { return k.buf[i * k.d + j]; }

// The value of the point i at the dimension j.
fn kde_data(k: Kde&, i: i64, j: i64) -> f64 { return k.buf[k.d * k.d + i * k.d + j]; }

// The weight of the point i.
fn kde_weight(k: Kde&, i: i64) -> f64 { return k.buf[k.d * k.d + k.n * k.d + i]; }

// The factor of a matrix, lower triangular. The matrix holds d*d values
// in row order. A matrix that is not positive definite gives the domain
// fault.
fn kde_cholesky(out: f64[>..]&, a: const f64[:], d: i64) -> i64 {
    out.resize(0);
    for k in 0..d * d { out.push(0.0); }
    for i in 0..d {
        for j in 0..i + 1 {
            var s = a[i * d + j];
            for m in 0..j { s = s - out[i * d + m] * out[j * d + m]; }
            if i == j {
                if !(s > 0.0) { return ERR_DOMAIN; }
                out[i * d + i] = sqrt(s);
            } else {
                out[i * d + j] = s / out[j * d + j];
            }
        }
    }
    return OK;
}

// The estimate of a sample of points. The weights may be empty, and then
// every point has the same weight. The bandwidth selects the factor of
// Scott, the factor of Silverman or a factor that the caller gives.
fn kde_fit(data: const f64[:], dims: i64, weights: const f64[:], bw: i64, factor: f64) -> Kde, i64 {
    if dims < 1 { return kde_fault(ERR_DOMAIN), ERR_DOMAIN; }
    if dims > MAX_KDE_DIMS { return kde_fault(ERR_TOO_LONG), ERR_TOO_LONG; }
    if data.len == 0 { return kde_fault(ERR_EMPTY), ERR_EMPTY; }
    if data.len % dims != 0 { return kde_fault(ERR_LENGTH), ERR_LENGTH; }
    let n = data.len / dims;
    if n < 2 { return kde_fault(ERR_TOO_SHORT), ERR_TOO_SHORT; }
    if dims > n { return kde_fault(ERR_DOMAIN), ERR_DOMAIN; }
    for v in data {
        let c = check_scalar(v);
        if c != OK { return kde_fault(c), c; }
    }
    var w: f64[>..] = [];
    var sum: f64 = 0.0;
    if weights.len == 0 {
        for k in 0..n { w.push(1.0 / n as f64); }
        sum = 1.0;
    } else {
        if weights.len != n { return kde_fault(ERR_LENGTH), ERR_LENGTH; }
        for v in weights {
            let c = check_scalar(v);
            if c != OK { return kde_fault(c), c; }
            if !(v > 0.0) { return kde_fault(ERR_DOMAIN), ERR_DOMAIN; }
            sum += v;
        }
        for v in weights { w.push(v / sum); }
    }
    var ss: f64 = 0.0;
    for v in w { ss += v * v; }
    let neff = 1.0 / ss;
    var f: f64 = 0.0;
    let dd = dims as f64;
    if bw == KDE_SCOTT { f = pow(neff, -1.0 / (dd + 4.0)); }
    else if bw == KDE_SILVERMAN { f = pow(neff * (dd + 2.0) / 4.0, -1.0 / (dd + 4.0)); }
    else if bw == KDE_SCALAR { f = factor; }
    else { return kde_fault(ERR_DOMAIN), ERR_DOMAIN; }
    if !(f > 0.0) { return kde_fault(ERR_DOMAIN), ERR_DOMAIN; }
    // the mean of the points
    var mean: f64[>..] = [];
    for j in 0..dims { mean.push(0.0); }
    for i in 0..n {
        for j in 0..dims { mean[j] = mean[j] + w[i] * data[i * dims + j]; }
    }
    // the covariance of the points, with the correction of the weights
    var cov: f64[>..] = [];
    for k in 0..dims * dims { cov.push(0.0); }
    for i in 0..n {
        for a in 0..dims {
            for b in 0..dims {
                let da = data[i * dims + a] - mean[a];
                let db = data[i * dims + b] - mean[b];
                cov[a * dims + b] = cov[a * dims + b] + w[i] * da * db;
            }
        }
    }
    let denom = 1.0 - ss;
    if !(denom > 0.0) { return kde_fault(ERR_DOMAIN), ERR_DOMAIN; }
    for k in 0..dims * dims {
        cov[k] = cov[k] / denom;
        if is_nan(cov[k]) || is_inf(cov[k]) { return kde_fault(ERR_OVERFLOW), ERR_OVERFLOW; }
    }
    var chol: f64[>..] = [];
    let c = kde_cholesky(chol, cov[..], dims);
    if c != OK { return kde_fault(c), c; }
    var cho: f64[>..] = [];
    for k in 0..dims * dims { cho.push(chol[k] * f); }
    var ld: f64 = 0.0;
    for j in 0..dims { ld = ld + log(cho[j * dims + j] * 2.5066282746310002); }
    var k = Kde { n: n, d: dims, factor: f, neff: neff, log_det: 2.0 * ld, buf: [] };
    for v in cho { k.buf.push(v); }
    for v in data { k.buf.push(v); }
    for v in w { k.buf.push(v); }
    return k, OK;
}

// The solution of L y = b, where L is the lower triangular factor of the
// covariance.
fn kde_solve(k: Kde&, b: const f64[:], y: f64[>..]&) {
    y.resize(0);
    for i in 0..k.d {
        var s = b[i];
        for j in 0..i { s = s - kde_cho(k, i, j) * y[j]; }
        y.push(s / kde_cho(k, i, i));
    }
}

// The check of the evaluation points of an estimate.
fn kde_check_points(k: Kde&, xs: const f64[:]) -> i64 {
    if xs.len == 0 { return ERR_EMPTY; }
    if xs.len % k.d != 0 { return ERR_LENGTH; }
    let m = xs.len / k.d;
    if k.n * m > MAX_KDE_WORK { return ERR_TOO_LONG; }
    for v in xs {
        let c = check_scalar(v);
        if c != OK { return c; }
    }
    return OK;
}

// The sum of the kernels of the estimate at the evaluation points. The
// values go into out. The sum of the weights of the estimate is 1, so the
// sum is the density times the norm of the factor.
fn kde_sum(k: Kde&, out: f64[>..]&, xs: const f64[:]) -> i64 {
    out.resize(0);
    let c = kde_check_points(k, xs);
    if c != OK { return c; }
    let m = xs.len / k.d;
    var diff: f64[>..] = [];
    var y: f64[>..] = [];
    for p in 0..m {
        var sum: f64 = 0.0;
        for i in 0..k.n {
            diff.resize(0);
            for j in 0..k.d { diff.push(xs[p * k.d + j] - kde_data(k, i, j)); }
            kde_solve(k, diff[..], y);
            var q: f64 = 0.0;
            for j in 0..k.d { q = q + y[j] * y[j]; }
            sum = sum + exp(-0.5 * q) * kde_weight(k, i);
        }
        out.push(sum);
    }
    return OK;
}

// The density of the estimate at the evaluation points. The values go
// into out.
fn kde_pdf(k: Kde&, out: f64[>..]&, xs: const f64[:]) -> i64 {
    out.resize(0);
    let norm = exp(k.log_det / 2.0);
    var sums: f64[>..] = [];
    let c = kde_sum(k, sums, xs);
    if c != OK { return c; }
    for s in sums { out.push(s / norm); }
    return OK;
}

// The logarithm of the density of the estimate at the evaluation points.
fn kde_logpdf(k: Kde&, out: f64[>..]&, xs: const f64[:]) -> i64 {
    out.resize(0);
    let norm = k.log_det / 2.0;
    var sums: f64[>..] = [];
    let c = kde_sum(k, sums, xs);
    if c != OK { return c; }
    for s in sums {
        if !(s > 0.0) { out.push(-1.0 / 0.0); }
        else { out.push(log(s) - norm); }
    }
    return OK;
}

// The covariance of the estimate, with the bandwidth in it. The values go
// into out, row by row.
fn kde_covariance(k: Kde&, out: f64[>..]&) -> i64 {
    out.resize(0);
    for i in 0..k.d {
        for j in 0..k.d {
            var s: f64 = 0.0;
            let last = min(i, j);
            for m in 0..last + 1 { s = s + kde_cho(k, i, m) * kde_cho(k, j, m); }
            out.push(s);
        }
    }
    return OK;
}

// The share of the estimate between two bounds. The call holds only a
// sample of one dimension.
fn kde_integrate_box_1d(k: Kde&, low: f64, high: f64) -> f64, i64 {
    if k.d != 1 { return nan(), ERR_DOMAIN; }
    if !(low <= high) { return nan(), ERR_DOMAIN; }
    let stdev = kde_cho(k, 0, 0);
    if !(stdev > 0.0) { return nan(), ERR_UNDEFINED; }
    var value: f64 = 0.0;
    for i in 0..k.n {
        let c = kde_data(k, i, 0);
        let lo, clo = norm_cdf((low - c) / stdev);
        let hi, chi = norm_cdf((high - c) / stdev);
        if clo != OK { return nan(), clo; }
        if chi != OK { return nan(), chi; }
        value = value + kde_weight(k, i) * (hi - lo);
    }
    return value, OK;
}

// Draws points from the estimate. The draws go into out, row by row. The
// same seed gives the same draws.
fn kde_resample(k: Kde&, out: f64[>..]&, size: i64, r: rng&) -> i64 {
    out.resize(0);
    if size < 0 { return ERR_DOMAIN; }
    if size > MAX_KDE_POINTS { return ERR_TOO_LONG; }
    var z: f64[>..] = [];
    for s in 0..size {
        let u = rs_unit(r);
        var acc: f64 = 0.0;
        var idx = k.n - 1;
        for i in 0..k.n {
            acc = acc + kde_weight(k, i);
            if u <= acc { idx = i; break; }
        }
        z.resize(0);
        for j in 0..k.d { z.push(rs_normal(r)); }
        for j in 0..k.d {
            var v: f64 = 0.0;
            for i in 0..j + 1 { v = v + kde_cho(k, j, i) * z[i]; }
            out.push(kde_data(k, idx, j) + v);
        }
    }
    return OK;
}
