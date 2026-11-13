// The Yeo-Johnson transform and the normal plot tables.
//
// Every call returns the result first and a fault code last. A fault
// gives NaN as the result.
namespace statistics;

import std;
import math;

fn tf_yj_y(lg: f64, neg: bool, l: f64) -> f64 {
    if !neg {
        if abs(l) < EPS { return lg; }
        return tf_expm1(l * lg) / l;
    }
    let l2 = 2.0 - l;
    if abs(l2) <= EPS { return -lg; }
    return -tf_expm1(l2 * lg) / l2;
}

// Checks a sample for a fit and fills lg with log1p(|x|) in the sorted
// order of x. The results are the sum of sign(x) * lg, the count of
// negative values, the largest |x| and the code.
fn tf_yj_prepare(xs: const f64[:], lg: f64[>..]&) -> f64, i64, f64, i64 {
    lg.resize(0);
    let c = check_sample(xs, 2);
    if c != OK { return nan(), 0, nan(), c; }
    for x in xs { lg.push(x); }
    sort(lg) { a, b => a < b };
    let n = lg.len;
    if lg[0] == lg[n - 1] { lg.resize(0); return nan(), 0, nan(), ERR_UNDEFINED; }
    let maxabs = max(abs(lg[0]), abs(lg[n - 1]));
    var nneg: i64 = 0;
    var s: f64 = 0.0;
    for i in 0..n {
        let x = lg[i];
        let g = log1p_core(abs(x));
        if x < 0.0 { nneg += 1; s -= g; } else { s += g; }
        lg[i] = g;
    }
    return s, nneg, maxabs, OK;
}

// The log-likelihood of the Yeo-Johnson transform. A sample of one sign
// uses the scaled log-variance. A mixed sample uses the variance of the
// transformed values.
fn tf_yj_llf(l: f64, lg: const f64[:], nneg: i64, sumterm: f64) -> f64, i64 {
    let n = lg.len;
    var lv: f64 = 0.0;
    var c: i64 = OK;
    if nneg == 0 {
        lv, c = tf_logvar(lg, l);
    } else if nneg == n {
        lv, c = tf_logvar(lg, 2.0 - l);
    } else {
        let v = tf_welford(n) { i => tf_yj_y(lg[i], i < nneg, l) };
        if is_bad(v) { return nan(), ERR_OVERFLOW; }
        if v < TF_TINY { return nan(), ERR_UNDEFINED; }
        lv = log(v);
    }
    if c != OK { return nan(), c; }
    let nf: f64 = n;
    let r = (l - 1.0) * sumterm - 0.5 * nf * lv;
    if is_bad(r) { return nan(), ERR_OVERFLOW; }
    return r, OK;
}

fn tf_yj_neg(l: f64, lg: const f64[:], nneg: i64, sumterm: f64) -> f64 {
    let r, c = tf_yj_llf(l, lg, nneg, sumterm);
    if c != OK { return inf(); }
    return -r;
}

// The minimum of F on the interval from x1 to x2 by the bounded Brent
// method of SciPy (fminbound). The search stops after 400 evaluations
// with ERR_FIT.
fn tf_fminbound<F>(x1: f64, x2: f64, xatol: f64) -> f64, i64 {
    let sqrt_eps = sqrt(2.2e-16);
    let golden_mean = 0.5 * (3.0 - sqrt(5.0));
    var a = x1;
    var b = x2;
    var fulc = a + golden_mean * (b - a);
    var nfc = fulc;
    var xf = fulc;
    var rat: f64 = 0.0;
    var e: f64 = 0.0;
    var fx = F(xf);
    var ffulc = fx;
    var fnfc = fx;
    var num: i64 = 1;
    var xm = 0.5 * (a + b);
    var tol1 = sqrt_eps * abs(xf) + xatol / 3.0;
    var tol2 = 2.0 * tol1;
    while abs(xf - xm) > tol2 - 0.5 * (b - a) {
        var golden = true;
        if abs(e) > tol1 {
            golden = false;
            var r = (xf - nfc) * (fx - ffulc);
            var q = (xf - fulc) * (fx - fnfc);
            var p = (xf - fulc) * q - (xf - nfc) * r;
            q = 2.0 * (q - r);
            if q > 0.0 { p = -p; }
            q = abs(q);
            r = e;
            e = rat;
            if abs(p) < abs(0.5 * q * r) && p > q * (a - xf) && p < q * (b - xf) {
                rat = p / q;
                let xt = xf + rat;
                if xt - a < tol2 || b - xt < tol2 {
                    if xm - xf >= 0.0 { rat = tol1; } else { rat = -tol1; }
                }
            } else {
                golden = true;
            }
        }
        if golden {
            if xf >= xm { e = a - xf; } else { e = b - xf; }
            rat = golden_mean * e;
        }
        var x = xf;
        if rat >= 0.0 { x = xf + max(abs(rat), tol1); } else { x = xf - max(abs(rat), tol1); }
        let fu = F(x);
        if is_nan(fu) { return nan(), ERR_FIT; }
        num += 1;
        if fu <= fx {
            if x >= xf { a = xf; } else { b = xf; }
            fulc = nfc; ffulc = fnfc;
            nfc = xf; fnfc = fx;
            xf = x; fx = fu;
        } else {
            if x < xf { a = x; } else { b = x; }
            if fu <= fnfc || nfc == xf {
                fulc = nfc; ffulc = fnfc;
                nfc = x; fnfc = fu;
            } else if fu <= ffulc || fulc == xf || fulc == nfc {
                fulc = x; ffulc = fu;
            }
        }
        xm = 0.5 * (a + b);
        tol1 = sqrt_eps * abs(xf) + xatol / 3.0;
        tol2 = 2.0 * tol1;
        if num >= MAX_FIT_STEPS * 2 { return nan(), ERR_FIT; }
    }
    return xf, OK;
}

// The best Yeo-Johnson lambda. The search is the bounded one of SciPy.
// The bounds keep every transformed value below the largest 64 bit float.
fn tf_yj_normmax(lg: const f64[:], nneg: i64, sumterm: f64, maxabs: f64) -> f64, i64 {
    var lm: f64 = 0.0;
    if maxabs > 1e300 { lm = log(20.0) + log(maxabs); } else { lm = log1p_core(20.0 * maxabs); }
    let log_eps = log(EPS);
    let log_tiny = (log(TF_TINY) - log_eps) / 2.0;
    let log_big = (TF_LOG_MAX + log_eps) / 2.0;
    var lb = log_tiny / lm;
    var ub = log_big / lm;
    if nneg == lg.len {
        let t = lb;
        lb = 2.0 - ub;
        ub = 2.0 - t;
    } else if nneg > 0 {
        let t = lb;
        lb = max(2.0 - ub, lb);
        ub = min(2.0 - t, ub);
    }
    let x, c = tf_fminbound(lb, ub, 1.48e-8) { l => tf_yj_neg(l, lg, nneg, sumterm) };
    if c != OK || is_bad(x) { return nan(), ERR_FIT; }
    return x, OK;
}

// Signature: yeojohnson_llf(lmb, xs) -> llf, code
// The log-likelihood of the Yeo-Johnson transform of xs for the lambda
// lmb. A constant sample gives ERR_UNDEFINED. The sample needs two values.
fn yeojohnson_llf(lmb: f64, xs: const f64[:]) -> f64, i64 {
    var lg: f64[>..] = [];
    let s, nneg, _m, c = tf_yj_prepare(xs, lg);
    if c != OK { return nan(), c; }
    let cl = check_scalar(lmb);
    if cl != OK { return nan(), cl; }
    return tf_yj_llf(lmb, lg[..], nneg, s);
}

// Signature: yeojohnson(out, xs, lmbda) -> lambda, code
// The Yeo-Johnson transform of xs for the given lambda. The values go
// into out in the order of xs, and the result is lmbda. A value that the
// transform cannot hold gives ERR_OVERFLOW.
fn yeojohnson(out: f64[>..]&, xs: const f64[:], lmbda: f64) -> f64, i64 {
    out.resize(0);
    let c = check_sample(xs, 1);
    if c != OK { return nan(), c; }
    let cl = check_scalar(lmbda);
    if cl != OK { return nan(), cl; }
    for x in xs {
        let z = tf_yj_y(log1p_core(abs(x)), x < 0.0, lmbda);
        if is_bad(z) { out.resize(0); return nan(), ERR_OVERFLOW; }
        out.push(z);
    }
    return lmbda, OK;
}

// Signature: yeojohnson_normmax(xs) -> lambda, code
// The lambda of maximum likelihood, from a bounded search. A constant
// sample gives ERR_UNDEFINED, where SciPy gives 1 for a sample of zeros.
// A search that does not end gives ERR_FIT.
fn yeojohnson_normmax(xs: const f64[:]) -> f64, i64 {
    var lg: f64[>..] = [];
    let s, nneg, maxabs, c = tf_yj_prepare(xs, lg);
    if c != OK { return nan(), c; }
    return tf_yj_normmax(lg[..], nneg, s, maxabs);
}

// Signature: yeojohnson_fit(out, xs) -> lambda, code
// The Yeo-Johnson transform with the lambda of maximum likelihood, as
// yeojohnson of SciPy does with no lambda. The result is the lambda used.
fn yeojohnson_fit(out: f64[>..]&, xs: const f64[:]) -> f64, i64 {
    out.resize(0);
    let lam, c = yeojohnson_normmax(xs);
    if c != OK { return nan(), c; }
    let _l, c2 = yeojohnson(out, xs, lam);
    if c2 != OK { return nan(), c2; }
    return lam, OK;
}

// ---- normal plot tables --------------------------------------------------

// The checks of the range and the count of a normal plot table.
fn tf_plot_args(n_values: i64, la: f64, lb: f64, n: i64) -> i64 {
    if is_nan(la) || is_nan(lb) { return ERR_NAN; }
    if is_inf(la) || is_inf(lb) { return ERR_INFINITY; }
    if la >= lb || n < 2 { return ERR_DOMAIN; }
    if n > MAX_INPUT / 2 || n * n_values > MAX_PLOT_WORK { return ERR_TOO_LONG; }
    return OK;
}

// The lambda of row i of n rows, as the linspace of NumPy gives it.
fn tf_plot_lambda(la: f64, lb: f64, n: i64, i: i64) -> f64 {
    if i == n - 1 { return lb; }
    let step = (lb - la) / (n - 1);
    return la + step * i;
}

// Signature: boxcox_normplot(out, xs, la, lb, n) -> rows, code
// The table of the Box-Cox normal plot. Row i holds a lambda and the
// correlation r of the normal probability plot of the transformed sample.
// The table has n evenly spaced lambdas from la to lb, and out holds
// each row as the pair lambda, r. la must be below lb, and n lies from 2
// to 2,097,152 with at most 20,000,000 values in all.
fn boxcox_normplot(out: f64[>..]&, xs: const f64[:], la: f64, lb: f64, n: i64) -> i64, i64 {
    out.resize(0);
    var logs: f64[>..] = [];
    let _s, c = tf_bc_prepare(xs, logs);
    if c != OK { return 0, c; }
    let ca = tf_plot_args(xs.len, la, lb, n);
    if ca != OK { return 0, ca; }
    var osm: f64[>..] = [];
    tf_osm(logs.len, osm);
    for i in 0..n {
        let l = tf_plot_lambda(la, lb, n, i);
        let r, cr = tf_pearson(osm[..]) { k => tf_bc_z(logs[k], l) };
        if cr != OK { out.resize(0); return 0, cr; }
        out.push(l);
        out.push(r);
    }
    return n, OK;
}

// Signature: yeojohnson_normplot(out, xs, la, lb, n) -> rows, code
// The table of the Yeo-Johnson normal plot, in the form of
// boxcox_normplot. The sample may hold values of any sign.
fn yeojohnson_normplot(out: f64[>..]&, xs: const f64[:], la: f64, lb: f64, n: i64) -> i64, i64 {
    out.resize(0);
    var lg: f64[>..] = [];
    let _s, nneg, _m, c = tf_yj_prepare(xs, lg);
    if c != OK { return 0, c; }
    let ca = tf_plot_args(xs.len, la, lb, n);
    if ca != OK { return 0, ca; }
    var osm: f64[>..] = [];
    tf_osm(lg.len, osm);
    for i in 0..n {
        let l = tf_plot_lambda(la, lb, n, i);
        let r, cr = tf_pearson(osm[..]) { k => tf_yj_y(lg[k], k < nneg, l) };
        if cr != OK { out.resize(0); return 0, cr; }
        out.push(l);
        out.push(r);
    }
    return n, OK;
}
