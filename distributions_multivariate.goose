// The multivariate families of SciPy 1.18.1 that need a matrix of a
// covariance: multivariate_normal, multivariate_t, matrix_normal and
// matrix_t.
//
// Every family writes its result into a flat buffer and gives the count
// of the values and the code. The matrices hold their values in row
// order. The draws take a generator, so one seed gives one answer.
//
// The files of the package hold a helper for the Cholesky factor and a
// helper for a triangular solve. The families below use the helpers of
// this file.
namespace statistics;

import std;
import math;

const MV_MAX_DIM: i64 = 64;        // the largest number of dimensions
const MV_MAX_DRAWS: i64 = 100000;  // the largest count of draws

// The check of a flat buffer of the size n.
fn mv_check_size(n: i64) -> i64 {
    if n < 0 { return ERR_DOMAIN; }
    if n > MV_MAX_DIM { return ERR_DOMAIN; }
    return OK;
}

// The check of a symmetric matrix of the size d, in row order.
fn mv_check_matrix(a: const f64[:], d: i64) -> i64 {
    let c = mv_check_size(d);
    if c != OK { return c; }
    return mv_check_flat(a, d * d);
}

// The check of a vector of the size n.
fn mv_check_vector(v: const f64[:], n: i64) -> i64 {
    let c = mv_check_size(n);
    if c != OK { return c; }
    return mv_check_flat(v, n);
}

// The check of a flat buffer of a known length. The length comes first,
// then a NaN and then an infinity.
fn mv_check_flat(v: const f64[:], need: i64) -> i64 {
    if v.len != need { return ERR_LENGTH; }
    for x in v {
        if is_nan(x) { return ERR_NAN; }
        if is_inf(x) { return ERR_INFINITY; }
    }
    return OK;
}

// The Cholesky factor of a symmetric matrix, lower triangular. The
// matrix holds d*d values in row order.
fn mv_cholesky(out: f64[>..]&, a: const f64[:], d: i64) -> i64 {
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

// The solve of L x = b, where L is the lower triangular factor.
fn mv_solve_lower(l: const f64[:], b: const f64[:], d: i64, x: f64[>..]&) {
    x.resize(0);
    for i in 0..d {
        var s = b[i];
        for j in 0..i { s = s - l[i * d + j] * x[j]; }
        x.push(s / l[i * d + i]);
    }
}

// The solve of L' x = b, where L is the lower triangular factor. The
// factor is read by its transpose.
fn mv_solve_lower_transpose(l: const f64[:], b: const f64[:], d: i64, x: f64[>..]&) {
    x.resize(0);
    // The values keep the row order, so that x[k] holds the value of the
    // row k.
    for i in 0..d { x.push(0.0); }
    var i = d - 1;
    while i >= 0 {
        var s = b[i];
        for k in i + 1..d { s = s - l[k * d + i] * x[k]; }
        x[i] = s / l[i * d + i];
        i -= 1;
    }
}

// The transpose of a matrix, row by row.
fn mv_transpose(out: f64[>..]&, a: const f64[:], r: i64, c: i64) {
    out.resize(0);
    for i in 0..c {
        for j in 0..r { out.push(a[j * c + i]); }
    }
}

// The product of two matrices, AB.
fn mv_matmul(out: f64[>..]&, a: const f64[:], b: const f64[:], r: i64, k: i64, c: i64) {
    out.resize(0);
    for i in 0..r {
        for j in 0..c {
            var s: f64 = 0.0;
            for m in 0..k { s = s + a[i * k + m] * b[m * c + j]; }
            out.push(s);
        }
    }
}

// The trace of a matrix.
fn mv_trace(a: const f64[:], d: i64) -> f64 {
    var s: f64 = 0.0;
    for i in 0..d { s = s + a[i * d + i]; }
    return s;
}

// The log of the determinant of a positive definite matrix. The Cholesky
// factor gives the determinant as the square of the product of its
// diagonal.
fn mv_logdet(a: const f64[:], d: i64) -> f64, i64 {
    var l: f64[>..] = [];
    let c = mv_cholesky(l, a, d);
    if c != OK { return nan(), c; }
    return mv_logdet_factor(l, d), OK;
}

// The log of the determinant from a Cholesky factor that the caller holds.
fn mv_logdet_factor(l: const f64[:], d: i64) -> f64 {
    var s: f64 = 0.0;
    for i in 0..d { s = s + 2.0 * log(l[i * d + i]); }
    return s;
}

// A draw of the chi law over the degrees of freedom. The draw comes from
// the inverse of the cumulative function of the chi law.
fn mv_chisq_draw(df: f64, r: rng&) -> f64 {
    let v, c = chi2_ppf(rs_unit(r), df);
    if c != OK { return 0.0; }
    return v;
}

// ---- multivariate_normal --------------------------------------------------

// The mean of the law, one value per dimension.
fn multivariate_normal_mean(out: f64[>..]&, mean: const f64[:], cov: const f64[:], d: i64) -> i64 {
    out.resize(0);
    let c1 = mv_check_vector(mean, d);
    if c1 != OK { return c1; }
    let c2 = mv_check_matrix(cov, d);
    if c2 != OK { return c2; }
    for v in mean { out.push(v); }
    return OK;
}

// The covariance of the law, row by row.
fn multivariate_normal_cov(out: f64[>..]&, mean: const f64[:], cov: const f64[:], d: i64) -> i64 {
    out.resize(0);
    let c1 = mv_check_vector(mean, d);
    if c1 != OK { return c1; }
    let c2 = mv_check_matrix(cov, d);
    if c2 != OK { return c2; }
    for v in cov { out.push(v); }
    return OK;
}

// The density of the law at a point.
fn multivariate_normal_pdf(x: const f64[:], mean: const f64[:], cov: const f64[:], d: i64) -> f64, i64 {
    let c1 = mv_check_vector(x, d);
    if c1 != OK { return nan(), c1; }
    let c2 = mv_check_vector(mean, d);
    if c2 != OK { return nan(), c2; }
    let c3 = mv_check_matrix(cov, d);
    if c3 != OK { return nan(), c3; }
    var l: f64[>..] = [];
    let cl = mv_cholesky(l, cov, d);
    if cl != OK { return nan(), cl; }
    var diff: f64[>..] = [];
    for i in 0..d { diff.push(x[i] - mean[i]); }
    var y: f64[>..] = [];
    mv_solve_lower(l, diff, d, y);
    var quad: f64 = 0.0;
    for i in 0..d { quad = quad + y[i] * y[i]; }
    let ld = mv_logdet_factor(l, d);
    return exp(-0.5 * quad - 0.5 * ld - 0.5 * d * log(2.0 * PI)), OK;
}

// The draws of the law. The buffer holds size*d values, row by row.
fn multivariate_normal_rvs(out: f64[>..]&, mean: const f64[:], cov: const f64[:], d: i64, size: i64, r: rng&) -> i64 {
    out.resize(0);
    let c1 = mv_check_vector(mean, d);
    if c1 != OK { return c1; }
    let c2 = mv_check_matrix(cov, d);
    if c2 != OK { return c2; }
    if size < 0 || size > MV_MAX_DRAWS { return ERR_DOMAIN; }
    var l: f64[>..] = [];
    let cl = mv_cholesky(l, cov, d);
    if cl != OK { return cl; }
    for s in 0..size {
        var z: f64[>..] = [];
        for i in 0..d { z.push(rs_normal(r)); }
        for i in 0..d {
            var v = mean[i];
            for j in 0..i + 1 { v = v + l[i * d + j] * z[j]; }
            out.push(v);
        }
    }
    return OK;
}

// ---- multivariate_t -------------------------------------------------------

// The mean of the law. It is the mean of the normal law.
fn multivariate_t_mean(out: f64[>..]&, mean: const f64[:], shape: const f64[:], d: i64) -> i64 {
    out.resize(0);
    let c1 = mv_check_vector(mean, d);
    if c1 != OK { return c1; }
    let c2 = mv_check_matrix(shape, d);
    if c2 != OK { return c2; }
    for v in mean { out.push(v); }
    return OK;
}

// The covariance of the law. It is the covariance of the normal law.
fn multivariate_t_cov(out: f64[>..]&, mean: const f64[:], shape: const f64[:], d: i64) -> i64 {
    out.resize(0);
    let c1 = mv_check_vector(mean, d);
    if c1 != OK { return c1; }
    let c2 = mv_check_matrix(shape, d);
    if c2 != OK { return c2; }
    for v in shape { out.push(v); }
    return OK;
}

// The degrees of freedom of the law.
fn multivariate_t_df(df: f64) -> f64, i64 {
    let c = check_scalar(df);
    if c != OK { return nan(), c; }
    if df <= 0.0 { return nan(), ERR_DOMAIN; }
    return df, OK;
}

// The draws of the law. The buffer holds size*d values, row by row. The
// law is a normal draw with the scale of a chi law over the degrees of
// freedom.
fn multivariate_t_rvs(out: f64[>..]&, mean: const f64[:], shape: const f64[:], df: f64, d: i64, size: i64, r: rng&) -> i64 {
    out.resize(0);
    let c1 = mv_check_vector(mean, d);
    if c1 != OK { return c1; }
    let c2 = mv_check_matrix(shape, d);
    if c2 != OK { return c2; }
    let df0, c3 = multivariate_t_df(df);
    if c3 != OK { return c3; }
    if size < 0 || size > MV_MAX_DRAWS { return ERR_DOMAIN; }
    var l: f64[>..] = [];
    let cl = mv_cholesky(l, shape, d);
    if cl != OK { return cl; }
    for s in 0..size {
        var z: f64[>..] = [];
        for i in 0..d { z.push(rs_normal(r)); }
        var v: f64[>..] = [];
        for i in 0..d {
            var t = 0.0;
            for j in 0..i + 1 { t = t + l[i * d + j] * z[j]; }
            v.push(t);
        }
        // The scale is the root of a chi law over the degrees of freedom.
        let g = mv_chisq_draw(df, r);
        let scale = sqrt(df / g);
        for i in 0..d { out.push(mean[i] + v[i] / scale); }
    }
    return OK;
}

// ---- matrix_normal --------------------------------------------------------

// The mean of the law, row by row. The mean holds r*c values.
fn matrix_normal_mean(out: f64[>..]&, mean: const f64[:], cov: const f64[:], d: i64) -> i64 {
    out.resize(0);
    let c1 = mv_check_vector(mean, d * d);
    if c1 != OK { return c1; }
    let c2 = mv_check_matrix(cov, d);
    if c2 != OK { return c2; }
    for v in mean { out.push(v); }
    return OK;
}

// The covariance of the law, row by row. It holds d*d*d*d values.
fn matrix_normal_cov(out: f64[>..]&, mean: const f64[:], cov: const f64[:], d: i64) -> i64 {
    out.resize(0);
    let c1 = mv_check_vector(mean, d * d);
    if c1 != OK { return c1; }
    let c2 = mv_check_matrix(cov, d);
    if c2 != OK { return c2; }
    // The covariance of the rows is the covariance of the law. The
    // covariance of the columns is the same matrix.
    for i in 0..d {
        for j in 0..d {
            for k in 0..d {
                for l in 0..d {
                    if i == k && j == l { out.push(cov[i * d + j]); }
                    else { out.push(0.0); }
                }
            }
        }
    }
    return OK;
}

// The draws of the law. The buffer holds size*d*d values, matrix by
// matrix, row by row.
fn matrix_normal_rvs(out: f64[>..]&, mean: const f64[:], cov: const f64[:], d: i64, size: i64, r: rng&) -> i64 {
    out.resize(0);
    let c1 = mv_check_vector(mean, d * d);
    if c1 != OK { return c1; }
    let c2 = mv_check_matrix(cov, d);
    if c2 != OK { return c2; }
    if size < 0 || size > MV_MAX_DRAWS { return ERR_DOMAIN; }
    var l: f64[>..] = [];
    let cl = mv_cholesky(l, cov, d);
    if cl != OK { return cl; }
    for s in 0..size {
        // Each row of the matrix is a draw of the law.
        for i in 0..d {
            var z: f64[>..] = [];
            for j in 0..d { z.push(rs_normal(r)); }
            for j in 0..d {
                var v = mean[i * d + j];
                for k in 0..j + 1 { v = v + l[j * d + k] * z[k]; }
                out.push(v);
            }
        }
    }
    return OK;
}

// ---- matrix_t -------------------------------------------------------------

// The mean of the law, row by row.
fn matrix_t_mean(out: f64[>..]&, mean: const f64[:], shape: const f64[:], d: i64) -> i64 {
    out.resize(0);
    let c1 = mv_check_vector(mean, d * d);
    if c1 != OK { return c1; }
    let c2 = mv_check_matrix(shape, d);
    if c2 != OK { return c2; }
    for v in mean { out.push(v); }
    return OK;
}

// The draws of the law. The buffer holds size*d*d values, matrix by
// matrix, row by row.
fn matrix_t_rvs(out: f64[>..]&, mean: const f64[:], shape: const f64[:], df: f64, d: i64, size: i64, r: rng&) -> i64 {
    out.resize(0);
    let c1 = mv_check_vector(mean, d * d);
    if c1 != OK { return c1; }
    let c2 = mv_check_matrix(shape, d);
    if c2 != OK { return c2; }
    let df0, c3 = multivariate_t_df(df);
    if c3 != OK { return c3; }
    if size < 0 || size > MV_MAX_DRAWS { return ERR_DOMAIN; }
    var l: f64[>..] = [];
    let cl = mv_cholesky(l, shape, d);
    if cl != OK { return cl; }
    for s in 0..size {
        // The matrix of the draw is the product of a normal draw and the
        // root of the scale.
        var z: f64[>..] = [];
        for i in 0..d * d { z.push(rs_normal(r)); }
        var g = mv_chisq_draw(df, r);
        let scale = sqrt(df / g);
        for i in 0..d {
            for j in 0..d {
                var v = 0.0;
                for k in 0..j + 1 { v = v + l[j * d + k] * z[i * d + k]; }
                out.push(mean[i * d + j] + v / scale);
            }
        }
    }
    return OK;
}
