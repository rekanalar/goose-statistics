// The Weibull, the gamma, the beta and the chi-square families of the
// package. Each family has a density, a cumulative function, a survival
// function, a quantile function and the four moments.
//
// Every call returns the result first and a fault code last. The checks
// run in one order: NaN, infinity, domain, pole.
namespace statistics;

import std;
import math;

fn weibull_check(shape: f64, scale: f64) -> i64 {
    if is_nan(shape) || is_nan(scale) { return ERR_NAN; }
    return first_code(check_shape(shape), check_shape(scale));
}

// The log of x / scale as a pair. The ratio underflows or overflows for an
// extreme scale, and then the difference of the two logs is exact enough.
fn log_ratio(x: f64, scale: f64) -> DD, DD {
    let z = dd_div(dd_from(x), dd_from(scale));
    if z.hi == 0.0 || is_inf(z.hi) { return z, dd_sub(dd_log_d(x), dd_log_d(scale)); }
    return z, dd_log(z);
}

// The log of z = x / scale and the power z^shape as a pair. The power
// overflows far above 1, so the callers test it.
fn weibull_pow(x: f64, shape: f64, scale: f64) -> DD, DD {
    let z, lz = log_ratio(x, scale);
    let l = dd_mul_d(lz, shape);
    if l.hi > 709.0 { return lz, dd_from(1.0 / 0.0); }
    if l.hi < -745.0 { return lz, dd_from(0.0); }
    return lz, dd_exp(l);
}

fn weibull_pdf(x: f64, shape: f64 = 1.0, scale: f64 = 1.0) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    let c = weibull_check(shape, scale);
    if c == ERR_NAN { return nan(), c; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    if c != OK { return nan(), c; }
    if x < 0.0 { return 0.0, OK; }
    if x == 0.0 {
        if shape < 1.0 { return nan(), ERR_OVERFLOW; }
        if shape == 1.0 { return 1.0 / scale, OK; }
        return 0.0, OK;
    }
    let z, t = weibull_pow(x, shape, scale);
    if is_inf(t.hi) { return 0.0, OK; }
    var l = dd_sub(dd_sub(dd_log_d(shape), dd_log_d(scale)), t);
    l = dd_add(l, dd_mul(z, dd_add_d(dd_from(shape), -1.0)));
    return exp_or_overflow(l);
}

fn weibull_cdf(x: f64, shape: f64 = 1.0, scale: f64 = 1.0) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    let c = weibull_check(shape, scale);
    if c == ERR_NAN { return nan(), c; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    if c != OK { return nan(), c; }
    if x <= 0.0 { return 0.0, OK; }
    let z, t = weibull_pow(x, shape, scale);
    if is_inf(t.hi) { return 1.0, OK; }
    if t.hi > 40.0 { return 1.0 - expf_dd(dd_neg(t)), OK; }
    let v, c2 = expm1(-dd_to_f64(t));
    return -v, OK;
}

fn weibull_sf(x: f64, shape: f64 = 1.0, scale: f64 = 1.0) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    let c = weibull_check(shape, scale);
    if c == ERR_NAN { return nan(), c; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    if c != OK { return nan(), c; }
    if x <= 0.0 { return 1.0, OK; }
    let z, t = weibull_pow(x, shape, scale);
    if is_inf(t.hi) { return 0.0, OK; }
    return expf_dd(dd_neg(t)), OK;
}

fn weibull_ppf(p: f64, shape: f64 = 1.0, scale: f64 = 1.0) -> f64, i64 {
    if is_nan(p) { return nan(), ERR_NAN; }
    let c = weibull_check(shape, scale);
    if c != OK { return nan(), c; }
    let c2 = chk_open(p);
    if c2 != OK { return nan(), c2; }
    let v, c3 = log1p(-p);
    let l = dd_add(dd_mul_d(dd_log_d(-v), 1.0 / shape), dd_log_d(scale));
    return exp_or_overflow(l);
}

// The raw moment g_k = Gamma(1 + k / shape) of the Weibull family as a
// pair. The central moments cancel, and the pair keeps their digits.
fn weibull_g(shape: f64, k: f64) -> DD, i64 {
    let z = dd_add_d(dd_div(dd_from(k), dd_from(shape)), 1.0);
    if z.hi > 170.0 { return dd_from(nan()), ERR_OVERFLOW; }
    return dd_exp(dd_lngamma(z)), OK;
}

fn weibull_mean(shape: f64 = 1.0, scale: f64 = 1.0) -> f64, i64 {
    let c = weibull_check(shape, scale);
    if c != OK { return nan(), c; }
    let g1, c1 = weibull_g(shape, 1.0);
    if c1 != OK { return nan(), c1; }
    let v = scale * dd_to_f64(g1);
    if is_inf(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

// The first raw moment and the variance of the standard form, g2 - g1^2.
fn weibull_v(shape: f64) -> DD, DD, DD, i64 {
    let g1, c1 = weibull_g(shape, 1.0);
    let g2, c2 = weibull_g(shape, 2.0);
    let cc = first_code(c1, c2);
    if cc != OK { return dd_from(nan()), dd_from(nan()), dd_from(nan()), cc; }
    return g1, g2, dd_sub(g2, dd_mul(g1, g1)), OK;
}

fn weibull_var(shape: f64 = 1.0, scale: f64 = 1.0) -> f64, i64 {
    let c = weibull_check(shape, scale);
    if c != OK { return nan(), c; }
    let g1, g2, v, c1 = weibull_v(shape);
    if c1 != OK { return nan(), c1; }
    let r = scale * scale * dd_to_f64(v);
    if is_inf(r) || is_nan(r) { return nan(), ERR_OVERFLOW; }
    return r, OK;
}

fn weibull_skew(shape: f64 = 1.0, scale: f64 = 1.0) -> f64, i64 {
    let c = weibull_check(shape, scale);
    if c != OK { return nan(), c; }
    let g1, g2, v, c1 = weibull_v(shape);
    let g3, c3 = weibull_g(shape, 3.0);
    let cc = first_code(c1, c3);
    if cc != OK { return nan(), cc; }
    let g11 = dd_mul(g1, g1);
    var m3 = dd_sub(g3, dd_mul_d(dd_mul(g2, g1), 3.0));
    m3 = dd_add(m3, dd_mul_d(dd_mul(g11, g1), 2.0));
    let r = dd_to_f64(dd_div(m3, dd_mul(v, dd_sqrt(v))));
    if is_inf(r) || is_nan(r) { return nan(), ERR_OVERFLOW; }
    return r, OK;
}

fn weibull_kurt(shape: f64 = 1.0, scale: f64 = 1.0) -> f64, i64 {
    let c = weibull_check(shape, scale);
    if c != OK { return nan(), c; }
    let g1, g2, v, c1 = weibull_v(shape);
    let g3, c3 = weibull_g(shape, 3.0);
    let g4, c4 = weibull_g(shape, 4.0);
    let cc = first_code(c1, first_code(c3, c4));
    if cc != OK { return nan(), cc; }
    let g11 = dd_mul(g1, g1);
    var m4 = dd_sub(g4, dd_mul_d(dd_mul(g3, g1), 4.0));
    m4 = dd_add(m4, dd_mul_d(dd_mul(g2, g11), 6.0));
    m4 = dd_sub(m4, dd_mul_d(dd_mul(g11, g11), 3.0));
    let v2 = dd_mul(v, v);
    let r = dd_to_f64(dd_div(dd_sub(m4, dd_mul_d(v2, 3.0)), v2));
    if is_inf(r) || is_nan(r) { return nan(), ERR_OVERFLOW; }
    return r, OK;
}

// ---- the gamma family ---------------------------------------------------------------------

// The log of the density of the standard gamma at z with the shape a,
// without the scale: (a - 1) log z - z - log gamma(a).
fn gamma_logpdf(z: DD, lz: DD, a: f64) -> DD {
    var l = dd_mul(lz, dd_add_d(dd_from(a), -1.0));
    l = dd_sub(l, z);
    return dd_sub(l, dd_lgam(a));
}

fn gammadist_check(shape: f64, scale: f64) -> i64 {
    if is_nan(shape) || is_nan(scale) { return ERR_NAN; }
    return first_code(check_shape(shape), check_shape(scale));
}

// The density of the gamma family for the shape and the scale.
fn gamma_density(x: f64, a: f64, scale: f64) -> f64, i64 {
    if is_inf(x) { return nan(), ERR_INFINITY; }
    if x < 0.0 { return 0.0, OK; }
    if x == 0.0 {
        if a < 1.0 { return nan(), ERR_OVERFLOW; }
        if a == 1.0 { return 1.0 / scale, OK; }
        return 0.0, OK;
    }
    let z, lz = log_ratio(x, scale);
    let l = dd_sub(gamma_logpdf(z, lz, a), dd_log_d(scale));
    return exp_or_overflow(l);
}

fn gammadist_pdf(x: f64, shape: f64 = 1.0, scale: f64 = 1.0) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    let c = gammadist_check(shape, scale);
    if c == ERR_NAN { return nan(), c; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    if c != OK { return nan(), c; }
    return gamma_density(x, shape, scale);
}

fn gammadist_cdf(x: f64, shape: f64 = 1.0, scale: f64 = 1.0) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    let c = gammadist_check(shape, scale);
    if c == ERR_NAN { return nan(), c; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    if c != OK { return nan(), c; }
    if x <= 0.0 { return 0.0, OK; }
    return gammainc(shape, x / scale);
}

fn gammadist_sf(x: f64, shape: f64 = 1.0, scale: f64 = 1.0) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    let c = gammadist_check(shape, scale);
    if c == ERR_NAN { return nan(), c; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    if c != OK { return nan(), c; }
    if x <= 0.0 { return 1.0, OK; }
    return gammaincc(shape, x / scale);
}

fn gammadist_ppf(p: f64, shape: f64 = 1.0, scale: f64 = 1.0) -> f64, i64 {
    if is_nan(p) { return nan(), ERR_NAN; }
    let c = gammadist_check(shape, scale);
    if c != OK { return nan(), c; }
    let c2 = chk_open(p);
    if c2 != OK { return nan(), c2; }
    let v, c3 = gammaincinv(shape, p);
    if c3 != OK { return nan(), c3; }
    let r = v * scale;
    if is_inf(r) { return nan(), ERR_OVERFLOW; }
    return r, OK;
}

fn gammadist_mean(shape: f64 = 1.0, scale: f64 = 1.0) -> f64, i64 {
    let c = gammadist_check(shape, scale);
    if c != OK { return nan(), c; }
    let v = shape * scale;
    if is_inf(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

fn gammadist_var(shape: f64 = 1.0, scale: f64 = 1.0) -> f64, i64 {
    let c = gammadist_check(shape, scale);
    if c != OK { return nan(), c; }
    let v = shape * scale * scale;
    if is_inf(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

fn gammadist_skew(shape: f64 = 1.0, scale: f64 = 1.0) -> f64, i64 {
    let c = gammadist_check(shape, scale);
    if c != OK { return nan(), c; }
    return 2.0 / sqrt(shape), OK;
}

fn gammadist_kurt(shape: f64 = 1.0, scale: f64 = 1.0) -> f64, i64 {
    let c = gammadist_check(shape, scale);
    if c != OK { return nan(), c; }
    return 6.0 / shape, OK;
}

// ---- the beta family -----------------------------------------------------------------------

fn betadist_check(a: f64, b: f64) -> i64 {
    if is_nan(a) || is_nan(b) { return ERR_NAN; }
    return first_code(check_shape(a), check_shape(b));
}

// The value of a beta tail at the pair u = hi + lo, from the tail at hi and
// a first order step: the density times lo. The argument of a derived
// family (Student t, F) is a ratio, and its rounding would otherwise cost
// as many units as the shape.
fn beta_step(base: f64, a: f64, b: f64, u: DD, sign: f64) -> f64 {
    let pdf, c = betadist_pdf(u.hi, a, b);
    if c != OK || is_inf(pdf) || is_nan(pdf) { return base; }
    return base + sign * pdf * u.lo;
}

fn betadist_pdf(x: f64, a: f64 = 1.0, b: f64 = 1.0) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    let c = betadist_check(a, b);
    if c == ERR_NAN { return nan(), c; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    if c != OK { return nan(), c; }
    if x < 0.0 || x > 1.0 { return 0.0, OK; }
    if x == 0.0 {
        if a < 1.0 { return nan(), ERR_OVERFLOW; }
        if a == 1.0 { return b, OK; }
        return 0.0, OK;
    }
    if x == 1.0 {
        if b < 1.0 { return nan(), ERR_OVERFLOW; }
        if b == 1.0 { return a, OK; }
        return 0.0, OK;
    }
    var l = dd_mul(dd_log_d(x), dd_add_d(dd_from(a), -1.0));
    l = dd_add(l, dd_mul(dd_log_1m(x), dd_add_d(dd_from(b), -1.0)));
    l = dd_sub(l, dd_lnbeta(a, b));
    return exp_or_overflow(l);
}

fn betadist_cdf(x: f64, a: f64 = 1.0, b: f64 = 1.0) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    let c = betadist_check(a, b);
    if c == ERR_NAN { return nan(), c; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    if c != OK { return nan(), c; }
    if x <= 0.0 { return 0.0, OK; }
    if x >= 1.0 { return 1.0, OK; }
    return betainc(a, b, x);
}

fn betadist_sf(x: f64, a: f64 = 1.0, b: f64 = 1.0) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    let c = betadist_check(a, b);
    if c == ERR_NAN { return nan(), c; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    if c != OK { return nan(), c; }
    if x <= 0.0 { return 1.0, OK; }
    if x >= 1.0 { return 0.0, OK; }
    return betainc_upper(a, b, x);
}

fn betadist_ppf(p: f64, a: f64 = 1.0, b: f64 = 1.0) -> f64, i64 {
    if is_nan(p) { return nan(), ERR_NAN; }
    let c = betadist_check(a, b);
    if c != OK { return nan(), c; }
    let c2 = chk_open(p);
    if c2 != OK { return nan(), c2; }
    return betaincinv(a, b, p);
}

fn betadist_mean(a: f64 = 1.0, b: f64 = 1.0) -> f64, i64 {
    let c = betadist_check(a, b);
    if c != OK { return nan(), c; }
    return a / (a + b), OK;
}

fn betadist_var(a: f64 = 1.0, b: f64 = 1.0) -> f64, i64 {
    let c = betadist_check(a, b);
    if c != OK { return nan(), c; }
    let s = a + b;
    let v = (a / s) * (b / s) / (s + 1.0);
    if is_inf(v) || is_nan(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

fn betadist_skew(a: f64 = 1.0, b: f64 = 1.0) -> f64, i64 {
    let c = betadist_check(a, b);
    if c != OK { return nan(), c; }
    let s = a + b;
    let v = 2.0 * (b - a) / (s + 2.0) * sqrt((s + 1.0) / a / b);
    if is_inf(v) || is_nan(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

fn betadist_kurt(a: f64 = 1.0, b: f64 = 1.0) -> f64, i64 {
    let c = betadist_check(a, b);
    if c != OK { return nan(), c; }
    let s = a + b;
    let d = a - b;
    let v = 6.0 * ((d / a) * (d / b) * (s + 1.0) / ((s + 2.0) * (s + 3.0)) - 1.0 / (s + 3.0));
    if is_inf(v) || is_nan(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

// ---- the chi-square family ---------------------------------------------------------

fn chi2_pdf(x: f64, df: f64 = 1.0) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    let c = chk_df(df);
    if c == ERR_NAN { return nan(), c; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    if c != OK { return nan(), c; }
    return gamma_density(x, 0.5 * df, 2.0);
}

fn chi2_cdf(x: f64, df: f64 = 1.0) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    let c = chk_df(df);
    if c == ERR_NAN { return nan(), c; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    if c != OK { return nan(), c; }
    if x <= 0.0 { return 0.0, OK; }
    return gammainc(0.5 * df, 0.5 * x);
}

fn chi2_sf(x: f64, df: f64 = 1.0) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    let c = chk_df(df);
    if c == ERR_NAN { return nan(), c; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    if c != OK { return nan(), c; }
    if x <= 0.0 { return 1.0, OK; }
    return gammaincc(0.5 * df, 0.5 * x);
}

fn chi2_ppf(p: f64, df: f64 = 1.0) -> f64, i64 {
    if is_nan(p) { return nan(), ERR_NAN; }
    let c = chk_df(df);
    if c != OK { return nan(), c; }
    let c2 = chk_open(p);
    if c2 != OK { return nan(), c2; }
    let v, c3 = gammaincinv(0.5 * df, p);
    if c3 != OK { return nan(), c3; }
    return 2.0 * v, OK;
}

fn chi2_mean(df: f64 = 1.0) -> f64, i64 {
    let c = chk_df(df);
    if c != OK { return nan(), c; }
    return df, OK;
}

fn chi2_var(df: f64 = 1.0) -> f64, i64 {
    let c = chk_df(df);
    if c != OK { return nan(), c; }
    return 2.0 * df, OK;
}

fn chi2_skew(df: f64 = 1.0) -> f64, i64 {
    let c = chk_df(df);
    if c != OK { return nan(), c; }
    return sqrt(8.0 / df), OK;
}

fn chi2_kurt(df: f64 = 1.0) -> f64, i64 {
    let c = chk_df(df);
    if c != OK { return nan(), c; }
    return 12.0 / df, OK;
}
