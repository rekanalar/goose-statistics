// The Box-Cox transform and the shared helpers of the transforms.
//
// Every call returns the result first and a fault code last. A fault
// gives NaN as the result.
namespace statistics;

import std;
import math;

// The criterion of the search for the Box-Cox lambda.
const NORMMAX_PEARSONR: i64 = 0;   // the best fit of the normal probability plot
const NORMMAX_MLE: i64 = 1;        // the maximum of the log-likelihood

// The largest product of rows and values in a normal plot table.
const MAX_PLOT_WORK: i64 = 20000000;

const TF_TINY: f64 = 2.2250738585072014e-308;
const TF_LOG_MAX: f64 = 709.782712893384;       // the log of the largest 64 bit float
const TF_YMAX: f64 = 1.7976931348623157e304;    // the largest transformed value that a lambda may give



// The value of e to the power t, minus 1. A large t gives infinity.
fn tf_expm1(t: f64) -> f64 {
    if t > 709.0 { return inf(); }
    if t < -40.0 { return -1.0; }
    return expm1_core(t);
}

// The population variance of the values F(0) to F(n-1), in one pass.
fn tf_welford<F>(n: i64) -> f64 {
    var mean: f64 = 0.0;
    var m2: f64 = 0.0;
    for i in 0..n {
        let y = F(i);
        let k: f64 = i + 1;
        let d = y - mean;
        mean += d / k;
        m2 += d * (y - mean);
    }
    let nf: f64 = n;
    return m2 / nf;
}

// The log of the sum of the exponentials of the values. The largest value
// sets the scale. A list of minus infinity gives minus infinity.
fn tf_logsumexp(xs: const f64[:]) -> f64 {
    var top = -inf();
    for x in xs {
        if x > top { top = x; }
    }
    if top == -inf() { return -inf(); }
    var s: f64 = 0.0;
    for x in xs { s += exp(x - top); }
    return top + log(s);
}

// The log of the variance of the values exp(logs), in the log space. The
// form never exponentiates a large value, so it does not overflow. It
// follows the scaled form of SciPy: the mean of the logs, then the log of
// each distance to the mean, then the log of the mean of the squares.
//
// The caller passes the lambda. The result is the log of the variance of
// the transformed values. A constant list gives ERR_UNDEFINED.
fn tf_logvar_exp(logs: const f64[:], l: f64) -> f64, i64 {
    let n = logs.len;
    let nf: f64 = n;
    var lx: f64[>..] = [];
    for g in logs { lx.push(l * g); }
    let s = tf_logsumexp(lx[..]);
    if is_nan(s) { return nan(), ERR_OVERFLOW; }
    let logmean = s - log(nf);
    // The list holds the values, then the deviations, then the squares.
    for i in 0..n {
        let d = lx[i] - logmean;
        if d > 0.0 {
            lx[i] = lx[i] + log(-expm1_core(-d));
        } else if d < 0.0 {
            lx[i] = logmean + log(-expm1_core(d));
        } else {
            lx[i] = -inf();
        }
    }
    for i in 0..n { lx[i] = 2.0 * lx[i]; }
    let s2 = tf_logsumexp(lx[..]);
    if is_nan(s2) { return nan(), ERR_OVERFLOW; }
    let r = s2 - log(nf);
    if is_inf(r) { return nan(), ERR_UNDEFINED; }
    return r, OK;
}

// The log of the variance of the values exp(l * logs) / l, in the log
// space. The scale of l does not change the variance of the values, so
// the log of the variance loses the factor 2 log of the size of l. A tiny
// l gives the variance of the logs.
fn tf_logvar(logs: const f64[:], l: f64) -> f64, i64 {
    let n = logs.len;
    if abs(l) < 1e-19 {
        let v = tf_welford(n) { i => logs[i] };
        if v <= 0.0 { return nan(), ERR_UNDEFINED; }
        return log(v), OK;
    }
    let v, c = tf_logvar_exp(logs, l);
    if c != OK { return nan(), c; }
    let r = v - 2.0 * log(abs(l));
    if is_bad(r) { return nan(), ERR_OVERFLOW; }
    return r, OK;
}


// ---- Box-Cox -------------------------------------------------------------

// The Box-Cox value of a log g = log(x).
fn tf_bc_z(g: f64, l: f64) -> f64 {
    if abs(l) < 1e-19 { return g; }
    return tf_expm1(l * g) / l;
}

// Checks a sample for a fit and fills logs with the sorted logs of the
// values. The result is the sum of the logs. The check order is the
// shared one, then the sign of the values, then a constant sample.
fn tf_bc_prepare(xs: const f64[:], logs: f64[>..]&) -> f64, i64 {
    logs.resize(0);
    let c = check_sample(xs, 2);
    if c != OK { return nan(), c; }
    for x in xs {
        if x <= 0.0 { return nan(), ERR_DOMAIN; }
    }
    for x in xs { logs.push(x); }
    sort(logs) { a, b => a < b };
    if logs[0] == logs[logs.len - 1] { logs.resize(0); return nan(), ERR_UNDEFINED; }
    for i in 0..logs.len { logs[i] = log(logs[i]); }
    let s, sc = exact_sum_div(logs[..], 1.0);
    if sc != OK { return nan(), sc; }
    return s, OK;
}

// The log-likelihood of the Box-Cox transform for the sorted logs.
fn tf_bc_llf(l: f64, logs: const f64[:], sumlog: f64) -> f64, i64 {
    let lv, c = tf_logvar(logs, l);
    if c != OK { return nan(), c; }
    let n: f64 = logs.len;
    let r = (l - 1.0) * sumlog - 0.5 * n * lv;
    if is_bad(r) { return nan(), ERR_OVERFLOW; }
    return r, OK;
}

// The value that a search minimizes for the maximum likelihood.
fn tf_bc_neg(l: f64, logs: const f64[:], sumlog: f64) -> f64 {
    let r, c = tf_bc_llf(l, logs, sumlog);
    if c != OK { return inf(); }
    return -r;
}

// The quantiles of the normal law at the medians of the uniform order
// statistics (Filliben), n of them.
fn tf_osm(n: i64, osm: f64[>..]&) {
    osm.resize(n, 0.0);
    let nf: f64 = n;
    let last = pow(0.5, 1.0 / nf);
    osm[n - 1] = norm_std_ppf(last);
    osm[0] = norm_std_ppf(1.0 - last);
    for k in 1..n - 1 {
        let i: f64 = k + 1;
        osm[k] = norm_std_ppf((i - 0.3175) / (nf + 0.365));
    }
}

// The Pearson correlation of the quantiles osm with the values F(0) to
// F(n-1). The values must run in increasing order. A value that is not
// finite gives ERR_OVERFLOW. A constant list gives ERR_UNDEFINED.
fn tf_pearson<F>(osm: const f64[:]) -> f64, i64 {
    let n = osm.len;
    let nf: f64 = n;
    // The values of the block are kept, so the block runs once per element.
    var ys: f64[>..] = [];
    var sx: f64 = 0.0;
    var sy: f64 = 0.0;
    for i in 0..n {
        let y = F(i);
        if is_bad(y) { return nan(), ERR_OVERFLOW; }
        ys.push(y);
        sx += osm[i];
        sy += y;
    }
    let mx = sx / nf;
    let my = sy / nf;
    var sxx: f64 = 0.0;
    var syy: f64 = 0.0;
    var sxy: f64 = 0.0;
    for i in 0..n {
        let dx = osm[i] - mx;
        let dy = ys[i] - my;
        sxx += dx * dx;
        syy += dy * dy;
        sxy += dx * dy;
    }
    if sxx <= 0.0 || syy <= 0.0 { return nan(), ERR_UNDEFINED; }
    let r = sxy / sqrt(sxx * syy);
    if is_bad(r) { return nan(), ERR_OVERFLOW; }
    if r > 1.0 { return 1.0, OK; }
    if r < -1.0 { return -1.0, OK; }
    return r, OK;
}

// The value that a search minimizes for the best normal probability plot.
// A lambda with a value that is not finite counts as the worst fit.
fn tf_bc_pearson_obj(l: f64, logs: const f64[:], osm: const f64[:]) -> f64 {
    let r, c = tf_pearson(osm) { i => tf_bc_z(logs[i], l) };
    if c != OK { return 2.0; }
    return 1.0 - r;
}

// The best lambda for sorted logs. The search is the bracket and the
// Brent minimizer of SciPy, from -2 and 2 with the tolerance 1.48e-8.
// The result is refused when the largest transformed value exceeds
// 1/10000 of the largest 64 bit float. SciPy replaces such a lambda by a
// smaller one that it finds with the Lambert W function. This package
// refuses the data instead.
fn tf_bc_normmax(logs: const f64[:], sumlog: f64, method: i64) -> f64, i64 {
    var lam: f64 = 0.0;
    var c: i64 = OK;
    if method == NORMMAX_MLE {
        lam, c = minimize_scalar(-2.0, 2.0, 1.48e-8) { l => tf_bc_neg(l, logs, sumlog) };
    } else {
        var osm: f64[>..] = [];
        tf_osm(logs.len, osm);
        lam, c = minimize_scalar(-2.0, 2.0, 1.48e-8) { l => tf_bc_pearson_obj(l, logs, osm[..]) };
    }
    if c != OK { return nan(), ERR_FIT; }
    if is_bad(lam) { return nan(), ERR_FIT; }
    let top = abs(tf_bc_z(logs[logs.len - 1], lam));
    let bot = abs(tf_bc_z(logs[0], lam));
    if is_bad(top) || is_bad(bot) || top > TF_YMAX || bot > TF_YMAX { return nan(), ERR_OVERFLOW; }
    return lam, OK;
}

// Signature: boxcox_llf(lmb, xs) -> llf, code
// The log-likelihood of the Box-Cox transform of xs for the lambda lmb.
// A value at or below 0 gives ERR_DOMAIN. A constant sample gives
// ERR_UNDEFINED. The sample needs two values.
fn boxcox_llf(lmb: f64, xs: const f64[:]) -> f64, i64 {
    var logs: f64[>..] = [];
    let sumlog, c = tf_bc_prepare(xs, logs);
    if c != OK { return nan(), c; }
    let cl = check_scalar(lmb);
    if cl != OK { return nan(), cl; }
    return tf_bc_llf(lmb, logs[..], sumlog);
}

// Signature: boxcox(out, xs, lmbda) -> lambda, code
// The Box-Cox transform of xs for the given lambda. The values go into
// out in the order of xs, and the result is lmbda. A value at or below 0
// gives ERR_DOMAIN, where SciPy gives NaN or -inf. A value that the
// transform cannot hold gives ERR_OVERFLOW.
fn boxcox(out: f64[>..]&, xs: const f64[:], lmbda: f64) -> f64, i64 {
    out.resize(0);
    let c = check_sample(xs, 1);
    if c != OK { return nan(), c; }
    let cl = check_scalar(lmbda);
    if cl != OK { return nan(), cl; }
    for x in xs {
        if x <= 0.0 { return nan(), ERR_DOMAIN; }
    }
    for x in xs {
        let z = tf_bc_z(log(x), lmbda);
        if is_bad(z) { out.resize(0); return nan(), ERR_OVERFLOW; }
        out.push(z);
    }
    return lmbda, OK;
}

// Signature: boxcox_normmax(xs, method = NORMMAX_PEARSONR) -> lambda, code
// The lambda that makes the transform most normal. The method is
// NORMMAX_PEARSONR, as in SciPy, or NORMMAX_MLE. The values must be above
// 0. SciPy accepts 0. A lambda that overflows the transform gives
// ERR_OVERFLOW. A search that does not end gives ERR_FIT.
fn boxcox_normmax(xs: const f64[:], method: i64 = NORMMAX_PEARSONR) -> f64, i64 {
    var logs: f64[>..] = [];
    let sumlog, c = tf_bc_prepare(xs, logs);
    if c != OK { return nan(), c; }
    if method != NORMMAX_PEARSONR && method != NORMMAX_MLE { return nan(), ERR_DOMAIN; }
    return tf_bc_normmax(logs[..], sumlog, method);
}

// Signature: boxcox_fit(out, xs) -> lambda, code
// The Box-Cox transform with the lambda of maximum likelihood, as
// boxcox of SciPy does with no lambda. The result is the lambda used.
fn boxcox_fit(out: f64[>..]&, xs: const f64[:]) -> f64, i64 {
    out.resize(0);
    let lam, c = boxcox_normmax(xs, NORMMAX_MLE);
    if c != OK { return nan(), c; }
    let _l, c2 = boxcox(out, xs, lam);
    if c2 != OK { return nan(), c2; }
    return lam, OK;
}

// Signature: boxcox_interval(xs, alpha) -> low, high, lambda, code
// The confidence interval of the lambda of maximum likelihood at the
// level 1 - alpha, by the profile likelihood. alpha lies in 0 to 1. The
// endpoint search steps at most 500 times on each side, and then gives
// ERR_FIT.
fn boxcox_interval(xs: const f64[:], alpha: f64) -> f64, f64, f64, i64 {
    var logs: f64[>..] = [];
    let sumlog, c = tf_bc_prepare(xs, logs);
    if c != OK { return nan(), nan(), nan(), c; }
    if is_nan(alpha) { return nan(), nan(), nan(), ERR_NAN; }
    if alpha <= 0.0 || alpha >= 1.0 { return nan(), nan(), nan(), ERR_DOMAIN; }
    let lmax, c1 = tf_bc_normmax(logs[..], sumlog, NORMMAX_MLE);
    if c1 != OK { return nan(), nan(), nan(), c1; }
    let q, c2 = chi2_ppf(1.0 - alpha, 1.0);
    if c2 != OK { return nan(), nan(), nan(), c2; }
    let top, c3 = tf_bc_llf(lmax, logs[..], sumlog);
    if c3 != OK { return nan(), nan(), nan(), c3; }
    let target = top - 0.5 * q;
    var hi: f64 = 0.0;
    var lo: f64 = 0.0;
    for side in 0..2 {
        var step = 0.1;
        if side == 1 { step = -0.1; }
        var edge = lmax + 5.0 * step;
        var count: i64 = 0;
        var f = tf_bc_neg(edge, logs[..], sumlog);
        while -f - target > 0.0 && count < 500 {
            edge += step;
            count += 1;
            f = tf_bc_neg(edge, logs[..], sumlog);
        }
        if count >= 500 || is_bad(f) { return nan(), nan(), nan(), ERR_FIT; }
        var root: f64 = 0.0;
        var cr: i64 = OK;
        if side == 0 {
            root, cr = brent_root(lmax, edge, 2e-12, BRENT_RTOL) { l => -tf_bc_neg(l, logs[..], sumlog) - target };
            hi = root;
        } else {
            root, cr = brent_root(edge, lmax, 2e-12, BRENT_RTOL) { l => -tf_bc_neg(l, logs[..], sumlog) - target };
            lo = root;
        }
        if cr != OK { return nan(), nan(), nan(), ERR_FIT; }
    }
    return lo, hi, lmax, OK;
}
