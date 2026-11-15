// Six continuous families whose quantile has no closed form: argus,
// dpareto_lognorm, erlang, exponnorm, foldcauchy and foldnorm.
//
// The names and the parameters follow SciPy 1.18.1.
//
// Every family has three calls: pdf, cdf and ppf. The quantile comes
// from a bracketed root of the cumulative function. Every call returns
// the result first and a fault code last. A value outside the support
// gives a density of 0 and a cumulative function of 0. The quantile of
// 0 or 1 gives the bound of the support.
//
// The laws of this file need the gamma law, the normal law and the
// incomplete gamma function of the package.
namespace statistics;

import std;
import math;

// ---- argus: the Argus law -------------------------------------------------

// The value of the helper of the Argus law. The helper is the regularized
// incomplete gamma of 1.5 at half the square of the argument.
fn argus_phi(chi: f64) -> f64 {
    let g, c = gammainc(1.5, 0.5 * chi * chi);
    if c != OK { return nan(); }
    return 0.5 * g;
}

// The density of the Argus law, without the checks.
fn argus_dens(x: f64, chi: f64) -> f64 {
    let phi = argus_phi(chi);
    if phi <= 0.0 { return 0.0; }
    let y = 1.0 - x * x;
    if y <= 0.0 { return 0.0; }
    let a = 3.0 * log(chi) - 0.5 * log(2.0 * PI) - log(phi);
    return exp(a + log(x) + 0.5 * log(y) - 0.5 * chi * chi * y);
}

fn argus_pdf(x: f64, chi: f64) -> f64, i64 {
    let k = check_shape(chi);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 || x >= 1.0 { return 0.0, OK; }
    return argus_dens(x, chi), OK;
}

fn argus_cdf(x: f64, chi: f64) -> f64, i64 {
    let k = check_shape(chi);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    if x >= 1.0 { return 1.0, OK; }
    // The survival function is the ratio of the helper at the point and
    // the helper at the shape.
    let num = argus_phi(chi * sqrt(1.0 - x * x));
    let den = argus_phi(chi);
    if den <= 0.0 { return nan(), ERR_FIT; }
    let sf = num / den;
    if sf < 0.0 { return 0.0, OK; }
    if sf > 1.0 { return 1.0, OK; }
    return 1.0 - sf, OK;
}

fn argus_ppf(p: f64, chi: f64) -> f64, i64 {
    let k = check_shape(chi);
    if k != OK { return nan(), k; }
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    if p == 0.0 { return 0.0, OK; }
    if p == 1.0 { return 1.0, OK; }
    return fam_search_ppf(p, 0.0, 0.5) { v => argus_cdf_value(v, chi) };
}

// ---- dpareto_lognorm: the double Pareto log-normal law --------------------

// The logarithm of the density of the normal law.
fn dpl_log_phi(z: f64) -> f64 {
    return -0.5 * z * z - 0.5 * log(2.0 * PI);
}

// The logarithm of the survival function of the normal law.
fn dpl_log_phic(z: f64) -> f64 {
    if z >= 0.0 {
        let e, c = erfc(z * INV_SQRT2);
        if c != OK { return nan(); }
        if e <= 0.0 { return -inf(); }
        return log(0.5 * e);
    }
    let lp = fam_log_ndtr(z);
    if is_nan(lp) { return nan(); }
    let e, c = expm1(lp);
    if c != OK { return nan(); }
    return log(-e);
}

// The logarithm of the Mills ratio of the normal law. The ratio is the
// survival function divided by the density.
fn dpl_log_r(z: f64) -> f64 {
    let lc = dpl_log_phic(z);
    if is_nan(lc) { return nan(); }
    return lc - dpl_log_phi(z);
}

// The logarithm of the difference of two terms. The first term must be
// the larger one.
fn dpl_log_sub(a: f64, b: f64) -> f64 {
    if a <= b { return -inf(); }
    if is_inf(b) { return a; }
    let e, c = expm1(b - a);
    if c != OK { return nan(); }
    return a + log(-e);
}

// The logarithm of the sum of two terms.
fn dpl_log_add(a: f64, b: f64) -> f64 {
    let m = max(a, b);
    if is_inf(m) { return m; }
    let d = min(a, b) - m;
    let e = exp(d);
    if is_nan(e) { return nan(); }
    return m + log(1.0 + e);
}

// The density of the double Pareto log-normal law, without the checks.
fn dpl_dens(x: f64, u: f64, s: f64, a: f64, b: f64) -> f64 {
    if x <= 0.0 { return 0.0; }
    let ly = log(x);
    let z = (ly - u) / s;
    let r1 = dpl_log_r(a * s - z);
    let r2 = dpl_log_r(b * s + z);
    if is_nan(r1) || is_nan(r2) { return 0.0; }
    // The sum of the two terms of the density.
    let sum = dpl_log_add(r1, r2);
    if is_nan(sum) { return 0.0; }
    let l = log(a) + log(b) - log(a + b) - ly + dpl_log_phi(z) + sum;
    if is_inf(l) { return 0.0; }
    return exp(l);
}

fn dpareto_lognorm_pdf(x: f64, u: f64, s: f64, a: f64, b: f64) -> f64, i64 {
    let c1 = check_scalar(u);
    if c1 != OK { return nan(), c1; }
    let c2 = check_shape(s);
    if c2 != OK { return nan(), c2; }
    let c3 = check_shape(a);
    if c3 != OK { return nan(), c3; }
    let c4 = check_shape(b);
    if c4 != OK { return nan(), c4; }
    let c5 = fam_x(x);
    if c5 != OK { return nan(), c5; }
    if x <= 0.0 { return 0.0, OK; }
    return dpl_dens(x, u, s, a, b), OK;
}

// The cumulative function of the double Pareto log-normal law. The sum
// of the two terms of the density is a difference of two large values,
// so the form with the logarithm of each part keeps the digits.
fn dpareto_lognorm_cdf(x: f64, u: f64, s: f64, a: f64, b: f64) -> f64, i64 {
    let c1 = check_scalar(u);
    if c1 != OK { return nan(), c1; }
    let c2 = check_shape(s);
    if c2 != OK { return nan(), c2; }
    let c3 = check_shape(a);
    if c3 != OK { return nan(), c3; }
    let c4 = check_shape(b);
    if c4 != OK { return nan(), c4; }
    let c5 = fam_x(x);
    if c5 != OK { return nan(), c5; }
    if x <= 0.0 { return 0.0, OK; }
    let ly = log(x);
    let z = (ly - u) / s;
    let t3 = log(b) + dpl_log_r(a * s - z);
    let t4 = log(a) + dpl_log_r(b * s + z);
    if is_nan(t3) || is_nan(t4) { return nan(), ERR_FIT; }
    // The difference of the two terms of the density. The sign of the
    // difference tells which term is the larger one.
    var t5: f64 = 0.0;
    var sign: f64 = 1.0;
    if t3 >= t4 {
        t5 = dpl_log_sub(t3, t4);
        sign = 1.0;
    } else {
        t5 = dpl_log_sub(t4, t3);
        sign = -1.0;
    }
    if is_nan(t5) { return nan(), ERR_FIT; }
    let t1 = fam_log_ndtr(z);
    if is_nan(t1) { return nan(), ERR_FIT; }
    let t2 = dpl_log_phi(z) + t5 - log(a + b);
    // The cumulative function is the normal cumulative function less the
    // term of the density.
    var l: f64 = 0.0;
    if sign > 0.0 {
        l = dpl_log_sub(t1, t2);
    } else {
        l = dpl_log_add(t1, t2);
    }
    if is_nan(l) { return nan(), ERR_FIT; }
    if is_inf(l) { return 0.0, OK; }
    let v = exp(l);
    if v < 0.0 { return 0.0, OK; }
    if v > 1.0 { return 1.0, OK; }
    return v, OK;
}

// The cumulative function of the law, without the code. The search of
// the quantile takes this value.
fn dpl_cdf_value(x: f64, u: f64, s: f64, a: f64, b: f64) -> f64 {
    let v, c = dpareto_lognorm_cdf(x, u, s, a, b);
    if c != OK { return nan(); }
    return v;
}

fn dpareto_lognorm_ppf(p: f64, u: f64, s: f64, a: f64, b: f64) -> f64, i64 {
    let c1 = check_scalar(u);
    if c1 != OK { return nan(), c1; }
    let c2 = check_shape(s);
    if c2 != OK { return nan(), c2; }
    let c3 = check_shape(a);
    if c3 != OK { return nan(), c3; }
    let c4 = check_shape(b);
    if c4 != OK { return nan(), c4; }
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    if p == 0.0 { return 0.0, OK; }
    return fam_search_ppf(p, 0.0, exp(u)) { v => dpl_cdf_value(v, u, s, a, b) };
}

// The cumulative function of the law, without the code. The search of
// the quantile takes this value.
fn argus_cdf_value(x: f64, chi: f64) -> f64 {
    let v, c = argus_cdf(x, chi);
    if c != OK { return nan(); }
    return v;
}

// The cumulative function of the law, without the code. The search of
// the quantile takes this value.
fn erlang_cdf_value(x: f64, a: f64) -> f64 {
    let v, c = erlang_cdf(x, a);
    if c != OK { return nan(); }
    return v;
}

// The cumulative function of the law, without the code. The search of
// the quantile takes this value.
fn exponnorm_cdf_value(x: f64, k: f64) -> f64 {
    let v, c = exponnorm_cdf(x, k);
    if c != OK { return nan(); }
    return v;
}

// The cumulative function of the law, without the code. The search of
// the quantile takes this value.
fn foldcauchy_cdf_value(x: f64, c: f64) -> f64 {
    let v, c = foldcauchy_cdf(x, c);
    if c != OK { return nan(); }
    return v;
}

// The cumulative function of the law, without the code. The search of
// the quantile takes this value.
fn foldnorm_cdf_value(x: f64, c: f64) -> f64 {
    let v, c = foldnorm_cdf(x, c);
    if c != OK { return nan(); }
    return v;
}

// ---- erlang: the Erlang law -----------------------------------------------

fn erlang_pdf(x: f64, a: f64) -> f64, i64 {
    let c1 = fam_x(x);
    if c1 != OK { return nan(), c1; }
    let c2 = check_whole_shape(a);
    if c2 != OK { return nan(), c2; }
    if x <= 0.0 { return 0.0, OK; }
    let l = (a - 1.0) * log(x) - x - fam_lngamma(a);
    return exp(l), OK;
}

fn erlang_cdf(x: f64, a: f64) -> f64, i64 {
    let c1 = fam_x(x);
    if c1 != OK { return nan(), c1; }
    let c2 = check_whole_shape(a);
    if c2 != OK { return nan(), c2; }
    if x <= 0.0 { return 0.0, OK; }
    let g, c = gammainc(a, x);
    if c != OK { return nan(), c; }
    return g, OK;
}

fn erlang_ppf(p: f64, a: f64) -> f64, i64 {
    let c1 = check_whole_shape(a);
    if c1 != OK { return nan(), c1; }
    let c2 = check_prob(p);
    if c2 != OK { return nan(), c2; }
    if p == 0.0 { return 0.0, OK; }
    return fam_search_ppf(p, 0.0, a * 4.0) { v => erlang_cdf_value(v, a) };
}

// ---- exponnorm: the sum of a normal and an exponential law ----------------

// The scaled complementary error function, which the law needs for a
// shape below one.
fn exponnorm_erfcx(u: f64) -> f64 {
    if u < 0.0 { return exp(u * u) * 2.0 - exponnorm_erfcx_pos(-u); }
    return exponnorm_erfcx_pos(u);
}

// The scaled complementary error function for an argument of 0 or more.
fn exponnorm_erfcx_pos(u: f64) -> f64 {
    // exp(u*u) overflows above 27. The scaled form holds the same value
    // without the exponential, so the far tail keeps its answer.
    if u >= 27.0 { return scaled_q(u); }
    if u < 0.5 { return exp(u * u) * (1.0 - erf_series(u)); }
    return exp(u * u) * erfc_pos(u);
}

fn exponnorm_pdf(x: f64, k: f64) -> f64, i64 {
    let c1 = fam_x(x);
    if c1 != OK { return nan(), c1; }
    let c2 = check_shape(k);
    if c2 != OK { return nan(), c2; }
    let u = (-x + 1.0 / k) * INV_SQRT2;
    if u >= 0.0 && k < 1.0 {
        let l = -log(2.0) - 0.5 * x * x - log(k) + log(exponnorm_erfcx(u));
        return exp(l), OK;
    }
    let inv = 1.0 / k;
    let l = inv * (-x + 0.5 * inv) + fam_log_ndtr(x - inv) - log(k);
    if is_nan(l) { return nan(), ERR_FIT; }
    return exp(l), OK;
}

fn exponnorm_cdf(x: f64, k: f64) -> f64, i64 {
    let c1 = fam_x(x);
    if c1 != OK { return nan(), c1; }
    let c2 = check_shape(k);
    if c2 != OK { return nan(), c2; }
    let u = (-x + 1.0 / k) * INV_SQRT2;
    let n, cn = norm_cdf(x);
    if cn != OK { return nan(), cn; }
    if u >= 0.0 && k < 1.0 {
        let v = n - 0.5 * exp(-0.5 * x * x) * exponnorm_erfcx(u);
        if v < 0.0 { return 0.0, OK; }
        if v > 1.0 { return 1.0, OK; }
        return v, OK;
    }
    let inv = 1.0 / k;
    let l = inv * (0.5 * inv - x) + fam_log_ndtr(x - inv);
    if is_nan(l) { return nan(), ERR_FIT; }
    let v = n - exp(l);
    if v < 0.0 { return 0.0, OK; }
    if v > 1.0 { return 1.0, OK; }
    return v, OK;
}

fn exponnorm_ppf(p: f64, k: f64) -> f64, i64 {
    let c1 = check_shape(k);
    if c1 != OK { return nan(), c1; }
    let c2 = check_prob(p);
    if c2 != OK { return nan(), c2; }
    if p == 0.0 { return -inf(), OK; }
    if p == 1.0 { return inf(), OK; }
    return fam_search_ppf(p, -40.0, 0.0) { v => exponnorm_cdf_value(v, k) };
}

// ---- foldcauchy: the Cauchy law of an absolute value ----------------------

fn foldcauchy_pdf(x: f64, c: f64) -> f64, i64 {
    let c1 = fam_x(x);
    if c1 != OK { return nan(), c1; }
    let c2 = check_shape_nonneg(c);
    if c2 != OK { return nan(), c2; }
    if x <= 0.0 { return 0.0, OK; }
    let d1 = x - c;
    let d2 = x + c;
    return (1.0 / (1.0 + d1 * d1) + 1.0 / (1.0 + d2 * d2)) / PI, OK;
}

fn foldcauchy_cdf(x: f64, c: f64) -> f64, i64 {
    let c1 = fam_x(x);
    if c1 != OK { return nan(), c1; }
    let c2 = check_shape_nonneg(c);
    if c2 != OK { return nan(), c2; }
    if x <= 0.0 { return 0.0, OK; }
    let d1 = x - c;
    let d2 = x + c;
    return (atan(d1) + atan(d2)) / PI, OK;
}

fn foldcauchy_ppf(p: f64, c: f64) -> f64, i64 {
    let c1 = check_shape_nonneg(c);
    if c1 != OK { return nan(), c1; }
    let c2 = check_prob(p);
    if c2 != OK { return nan(), c2; }
    if p == 0.0 { return 0.0, OK; }
    return fam_search_ppf(p, 0.0, c + 1.0) { v => foldcauchy_cdf_value(v, c) };
}

// ---- foldnorm: the normal law of an absolute value ------------------------

fn foldnorm_pdf(x: f64, c: f64) -> f64, i64 {
    let c1 = fam_x(x);
    if c1 != OK { return nan(), c1; }
    let c2 = check_shape_nonneg(c);
    if c2 != OK { return nan(), c2; }
    if x <= 0.0 { return 0.0, OK; }
    let p1, e1 = norm_pdf(x + c);
    if e1 != OK { return nan(), e1; }
    let p2, e2 = norm_pdf(x - c);
    if e2 != OK { return nan(), e2; }
    return p1 + p2, OK;
}

fn foldnorm_cdf(x: f64, c: f64) -> f64, i64 {
    let c1 = fam_x(x);
    if c1 != OK { return nan(), c1; }
    let c2 = check_shape_nonneg(c);
    if c2 != OK { return nan(), c2; }
    if x <= 0.0 { return 0.0, OK; }
    let e1, ce1 = erf((x - c) * INV_SQRT2);
    if ce1 != OK { return nan(), ce1; }
    let e2, ce2 = erf((x + c) * INV_SQRT2);
    if ce2 != OK { return nan(), ce2; }
    return 0.5 * (e1 + e2), OK;
}

fn foldnorm_ppf(p: f64, c: f64) -> f64, i64 {
    let c1 = check_shape_nonneg(c);
    if c1 != OK { return nan(), c1; }
    let c2 = check_prob(p);
    if c2 != OK { return nan(), c2; }
    if p == 0.0 { return 0.0, OK; }
    return fam_search_ppf(p, 0.0, c + 1.0) { v => foldnorm_cdf_value(v, c) };
}
