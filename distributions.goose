// The continuous distributions of the package: the normal, the uniform,
// the exponential and the lognormal families. Each family has a density
// or a mass, a cumulative function, a survival function, a quantile
// function and the four moments.
//
// Every call returns the result first and a fault code last. A fault gives
// NaN as the result. The checks run in one order: NaN, infinity, domain,
// pole. The survival function is the tail itself and never one minus the
// cumulative function. Every loop has a fixed bound.
namespace statistics;

import std;
import math;

const SQRT2: f64 = 1.4142135623730951;
const INV_SQRT2: f64 = 0.7071067811865476;
const SQRT_2PI: f64 = 2.5066282746310002;

// ---- checks ------------------------------------------------------------

// The code of a count of degrees of freedom: at or above 1.
fn chk_df(x: f64) -> i64 {
    let c = check_scalar(x);
    if c != OK { return c; }
    if x < 1.0 { return ERR_DOMAIN; }
    return OK;
}

// The code of a whole number parameter at or above 0.
fn chk_whole(x: f64) -> i64 {
    let c = check_scalar(x);
    if c != OK { return c; }
    if x < 0.0 || floor(x) != x { return ERR_DOMAIN; }
    return OK;
}

// The code of a probability parameter in [0, 1].
fn chk_unit(p: f64) -> i64 {
    let c = check_scalar(p);
    if c != OK { return c; }
    if p < 0.0 || p > 1.0 { return ERR_DOMAIN; }
    return OK;
}

// The code of the probability of a continuous quantile: in (0, 1).
fn chk_open(p: f64) -> i64 {
    let c = check_scalar(p);
    if c != OK { return c; }
    if p <= 0.0 || p >= 1.0 { return ERR_DOMAIN; }
    return OK;
}

// The code of a probability of a success: in (0, 1].
fn chk_succ(p: f64) -> i64 {
    let c = check_scalar(p);
    if c != OK { return c; }
    if p <= 0.0 || p > 1.0 { return ERR_DOMAIN; }
    return OK;
}

// The first code that is not OK.
fn first_code(a: i64, b: i64) -> i64 {
    if a != OK { return a; }
    return b;
}

// ---- helpers -------------------------------------------------------------

// The exponential of a pair as one value. Below the smallest subnormal the
// true value rounds to 0. Above the largest double the call gives infinity.
fn expf_dd(l: DD) -> f64 {
    if is_nan(l.hi) { return nan(); }
    if l.hi < -745.3 { return 0.0; }
    if l.hi > 709.9 { return 1.0 / 0.0; }
    let m, k = dd_exp_parts(l);
    return dd_scale(m, k);
}

// The value or the overflow fault of the exponential of a pair.
fn exp_or_overflow(l: DD) -> f64, i64 {
    let v = expf_dd(l);
    if is_inf(v) || is_nan(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

// The log of 1 - x for x in [0, 1) as a pair.
fn dd_log_1m(x: f64) -> DD {
    let s, e = two_sum(1.0, -x);
    return dd_log(DD { hi: s, lo: e });
}

// The log gamma of a positive shape as a pair. The tiny shapes use the
// series of log gamma(1 + a) minus log a.
fn dd_lgam(a: f64) -> DD {
    if a < 1e-3 { return dd_sub(dd_lngamma_tiny(a), dd_log_d(a)); }
    return dd_lngamma(dd_from(a));
}

// ---- the normal family -----------------------------------------------------

fn norm_z(x: f64, mu: f64, sigma: f64) -> f64, i64 {
    if is_nan(x) || is_nan(mu) || is_nan(sigma) { return nan(), ERR_NAN; }
    if is_inf(x) || is_inf(mu) || is_inf(sigma) { return nan(), ERR_INFINITY; }
    if sigma <= 0.0 { return nan(), ERR_DOMAIN; }
    return (x - mu) / sigma, OK;
}

// The scaled difference (x - mu) / sigma as a pair, the difference exact.
fn norm_zdd(x: f64, mu: f64, sigma: f64) -> DD {
    let s, e = two_sum(x, -mu);
    return dd_div(DD { hi: s, lo: e }, dd_from(sigma));
}

// The upper tail of the standard normal at the pair z. Above 0.7 the tail
// is the scaled function times exp(-z*z/2), and the square of z in pairs
// keeps the digits of the far tail. Below 0.7 the tail is 0.5 erfc.
fn norm_upper(zd: DD) -> f64 {
    let z = zd.hi;
    if z > 38.6 { return 0.0; }
    if z < 0.7 {
        let v, c = erfc(z * INV_SQRT2);
        return 0.5 * v;
    }
    let q = scaled_q(z * INV_SQRT2);
    return 0.5 * q * expf_dd(dd_mul_d(dd_mul(zd, zd), -0.5));
}

fn norm_pdf(x: f64, mu: f64 = 0.0, sigma: f64 = 1.0) -> f64, i64 {
    if is_nan(x) || is_nan(mu) || is_nan(sigma) { return nan(), ERR_NAN; }
    if is_inf(mu) { return nan(), ERR_INFINITY; }
    let z, c = norm_z(x, mu, sigma);
    if c != OK { return nan(), c; }
    if is_inf(z) { return nan(), ERR_INFINITY; }
    let zd = norm_zdd(x, mu, sigma);
    return expf_dd(dd_mul_d(dd_mul(zd, zd), -0.5)) / (sigma * SQRT_2PI), OK;
}

fn norm_cdf(x: f64, mu: f64 = 0.0, sigma: f64 = 1.0) -> f64, i64 {
    let z, c = norm_z(x, mu, sigma);
    if c != OK { return nan(), c; }
    if is_inf(z) { return nan(), ERR_INFINITY; }
    return norm_upper(dd_neg(norm_zdd(x, mu, sigma))), OK;
}

fn norm_sf(x: f64, mu: f64 = 0.0, sigma: f64 = 1.0) -> f64, i64 {
    let z, c = norm_z(x, mu, sigma);
    if c != OK { return nan(), c; }
    if is_inf(z) { return nan(), ERR_INFINITY; }
    return norm_upper(norm_zdd(x, mu, sigma)), OK;
}

// The quantile of the standard normal from the inverse of erfc.
fn norm_std_ppf(p: f64) -> f64 {
    if p == 0.5 { return 0.0; }
    if p < 0.5 {
        let v, c = erfcinv(2.0 * p);
        return -SQRT2 * v;
    }
    let v, c = erfcinv(2.0 * (1.0 - p));
    return SQRT2 * v;
}

// The quantile of the standard normal as a pair. One Newton step of the
// tail, in the form that never subtracts from one, gives the low part. The
// step helps in the tails, where the quantile feeds a product.
fn norm_std_ppf_dd(p: f64) -> f64, f64 {
    let z = norm_std_ppf(p);
    if abs(z) > 37.0 || z == 0.0 { return z, 0.0; }
    let hi2, lo2 = two_prod(z, z);
    let pdf = exp(-0.5 * hi2) * (1.0 - 0.5 * lo2) / SQRT_2PI;
    if pdf <= 0.0 { return z, 0.0; }
    var lo: f64 = 0.0;
    if p < 0.5 {
        lo = -(norm_upper(DD { hi: -z, lo: 0.0 }) - p) / pdf;
    } else {
        lo = (norm_upper(DD { hi: z, lo: 0.0 }) - (1.0 - p)) / pdf;
    }
    if abs(lo) > 1e-9 * abs(z) { return z, 0.0; }
    return z, lo;
}

fn norm_ppf(p: f64, mu: f64 = 0.0, sigma: f64 = 1.0) -> f64, i64 {
    let c1 = first_code(check_scalar(mu), check_shape(sigma));
    let c2 = chk_open(p);
    if is_nan(p) || is_nan(mu) || is_nan(sigma) { return nan(), ERR_NAN; }
    if c1 != OK { return nan(), c1; }
    if c2 != OK { return nan(), c2; }
    let zh, zl = norm_std_ppf_dd(p);
    let ph, pl = two_prod(sigma, zh);
    let s, e = two_sum(mu, ph);
    return s + (e + pl + sigma * zl), OK;
}

fn norm_mean(mu: f64 = 0.0, sigma: f64 = 1.0) -> f64, i64 {
    let c = first_code(check_scalar(mu), check_shape(sigma));
    if c != OK { return nan(), c; }
    return mu, OK;
}

fn norm_var(mu: f64 = 0.0, sigma: f64 = 1.0) -> f64, i64 {
    let c = first_code(check_scalar(mu), check_shape(sigma));
    if c != OK { return nan(), c; }
    let v = sigma * sigma;
    if is_inf(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

fn norm_skew(mu: f64 = 0.0, sigma: f64 = 1.0) -> f64, i64 {
    let c = first_code(check_scalar(mu), check_shape(sigma));
    if c != OK { return nan(), c; }
    return 0.0, OK;
}

fn norm_kurt(mu: f64 = 0.0, sigma: f64 = 1.0) -> f64, i64 {
    let c = first_code(check_scalar(mu), check_shape(sigma));
    if c != OK { return nan(), c; }
    return 0.0, OK;
}

// ---- the uniform family ------------------------------------------------------

fn unif_check(lo: f64, hi: f64) -> f64, i64 {
    let c = first_code(check_scalar(lo), check_scalar(hi));
    if is_nan(lo) || is_nan(hi) { return nan(), ERR_NAN; }
    if c != OK { return nan(), c; }
    if hi <= lo { return nan(), ERR_DOMAIN; }
    let w = hi - lo;
    if is_inf(w) { return nan(), ERR_OVERFLOW; }
    return w, OK;
}

fn unif_pdf(x: f64, lo: f64 = 0.0, hi: f64 = 1.0) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    let w, c = unif_check(lo, hi);
    if c != OK { return nan(), c; }
    if x < lo || x > hi { return 0.0, OK; }
    return 1.0 / w, OK;
}

fn unif_cdf(x: f64, lo: f64 = 0.0, hi: f64 = 1.0) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    let w, c = unif_check(lo, hi);
    if c != OK { return nan(), c; }
    if x <= lo { return 0.0, OK; }
    if x >= hi { return 1.0, OK; }
    return (x - lo) / w, OK;
}

fn unif_sf(x: f64, lo: f64 = 0.0, hi: f64 = 1.0) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    let w, c = unif_check(lo, hi);
    if c != OK { return nan(), c; }
    if x <= lo { return 1.0, OK; }
    if x >= hi { return 0.0, OK; }
    return (hi - x) / w, OK;
}

fn unif_ppf(p: f64, lo: f64 = 0.0, hi: f64 = 1.0) -> f64, i64 {
    if is_nan(p) { return nan(), ERR_NAN; }
    let w, c = unif_check(lo, hi);
    if c != OK { return nan(), c; }
    let c2 = chk_open(p);
    if c2 != OK { return nan(), c2; }
    // lo + p (hi - lo) with the product and the difference in pairs
    let wh, wl = two_sum(hi, -lo);
    let ph, pl = two_prod(p, wh);
    let s, e = two_sum(lo, ph);
    return s + (e + pl + p * wl), OK;
}

fn unif_mean(lo: f64 = 0.0, hi: f64 = 1.0) -> f64, i64 {
    let w, c = unif_check(lo, hi);
    if c != OK { return nan(), c; }
    return lo + 0.5 * w, OK;
}

fn unif_var(lo: f64 = 0.0, hi: f64 = 1.0) -> f64, i64 {
    let w, c = unif_check(lo, hi);
    if c != OK { return nan(), c; }
    let v = w * w / 12.0;
    if is_inf(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

fn unif_skew(lo: f64 = 0.0, hi: f64 = 1.0) -> f64, i64 {
    let w, c = unif_check(lo, hi);
    if c != OK { return nan(), c; }
    return 0.0, OK;
}

fn unif_kurt(lo: f64 = 0.0, hi: f64 = 1.0) -> f64, i64 {
    let w, c = unif_check(lo, hi);
    if c != OK { return nan(), c; }
    return -1.2, OK;
}

// ---- the exponential family -----------------------------------------------------

// The product rate * x as a pair, so the exponential keeps its digits.
fn expon_arg(rate: f64, x: f64) -> DD {
    let h, l = two_prod(rate, x);
    return DD { hi: h, lo: l };
}

fn expon_pdf(x: f64, rate: f64 = 1.0) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    let c = check_shape(rate);
    if c == ERR_NAN { return nan(), c; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    if c != OK { return nan(), c; }
    if x < 0.0 { return 0.0, OK; }
    let l = dd_neg(expon_arg(rate, x));
    return rate * expf_dd(l), OK;
}

fn expon_cdf(x: f64, rate: f64 = 1.0) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    let c = check_shape(rate);
    if c == ERR_NAN { return nan(), c; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    if c != OK { return nan(), c; }
    if x <= 0.0 { return 0.0, OK; }
    let t = rate * x;
    if t > 40.0 { return 1.0 - expf_dd(dd_neg(expon_arg(rate, x))), OK; }
    let v, c2 = expm1(-t);
    return -v, OK;
}

fn expon_sf(x: f64, rate: f64 = 1.0) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    let c = check_shape(rate);
    if c == ERR_NAN { return nan(), c; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    if c != OK { return nan(), c; }
    if x <= 0.0 { return 1.0, OK; }
    return expf_dd(dd_neg(expon_arg(rate, x))), OK;
}

fn expon_ppf(p: f64, rate: f64 = 1.0) -> f64, i64 {
    if is_nan(p) { return nan(), ERR_NAN; }
    let c = check_shape(rate);
    if c != OK { return nan(), c; }
    let c2 = chk_open(p);
    if c2 != OK { return nan(), c2; }
    let v, c3 = log1p(-p);
    return -v / rate, OK;
}

fn expon_mean(rate: f64 = 1.0) -> f64, i64 {
    let c = check_shape(rate);
    if c != OK { return nan(), c; }
    let v = 1.0 / rate;
    if is_inf(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

fn expon_var(rate: f64 = 1.0) -> f64, i64 {
    let c = check_shape(rate);
    if c != OK { return nan(), c; }
    let v = 1.0 / (rate * rate);
    if is_inf(v) || is_nan(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

fn expon_skew(rate: f64 = 1.0) -> f64, i64 {
    let c = check_shape(rate);
    if c != OK { return nan(), c; }
    return 2.0, OK;
}

fn expon_kurt(rate: f64 = 1.0) -> f64, i64 {
    let c = check_shape(rate);
    if c != OK { return nan(), c; }
    return 6.0, OK;
}

// ---- the lognormal family ----------------------------------------------------------

fn lognorm_pdf(x: f64, mu: f64 = 0.0, sigma: f64 = 1.0) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    let c = first_code(check_scalar(mu), check_shape(sigma));
    if is_nan(mu) || is_nan(sigma) { return nan(), ERR_NAN; }
    if c == ERR_NAN { return nan(), c; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    if c != OK { return nan(), c; }
    if x <= 0.0 { return 0.0, OK; }
    let lx = dd_log_d(x);
    let w = dd_div(dd_sub(lx, dd_from(mu)), dd_from(sigma));
    var l = dd_mul_d(dd_mul(w, w), -0.5);
    l = dd_sub(l, lx);
    l = dd_sub(l, dd_log_d(sigma));
    l = dd_sub(l, dd_mul_d(DD { hi: LN2PI_H, lo: LN2PI_L }, 0.5));
    return exp_or_overflow(l);
}

// The scaled log (log x - mu) / sigma as a pair.
fn lognorm_z(x: f64, mu: f64, sigma: f64) -> DD {
    return dd_div(dd_sub(dd_log_d(x), dd_from(mu)), dd_from(sigma));
}

fn lognorm_cdf(x: f64, mu: f64 = 0.0, sigma: f64 = 1.0) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    let c = first_code(check_scalar(mu), check_shape(sigma));
    if is_nan(mu) || is_nan(sigma) { return nan(), ERR_NAN; }
    if c == ERR_NAN { return nan(), c; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    if c != OK { return nan(), c; }
    if x <= 0.0 { return 0.0, OK; }
    return norm_upper(dd_neg(lognorm_z(x, mu, sigma))), OK;
}

fn lognorm_sf(x: f64, mu: f64 = 0.0, sigma: f64 = 1.0) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    let c = first_code(check_scalar(mu), check_shape(sigma));
    if is_nan(mu) || is_nan(sigma) { return nan(), ERR_NAN; }
    if c == ERR_NAN { return nan(), c; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    if c != OK { return nan(), c; }
    if x <= 0.0 { return 1.0, OK; }
    return norm_upper(lognorm_z(x, mu, sigma)), OK;
}

fn lognorm_ppf(p: f64, mu: f64 = 0.0, sigma: f64 = 1.0) -> f64, i64 {
    if is_nan(p) || is_nan(mu) || is_nan(sigma) { return nan(), ERR_NAN; }
    let c = first_code(check_scalar(mu), check_shape(sigma));
    if c != OK { return nan(), c; }
    let c2 = chk_open(p);
    if c2 != OK { return nan(), c2; }
    // the exponent mu + sigma z with the product in a pair
    let zh, zl = norm_std_ppf_dd(p);
    let ph, pl = two_prod(sigma, zh);
    let s, e = two_sum(mu, ph);
    return exp_or_overflow(DD { hi: s, lo: e + pl + sigma * zl });
}

fn lognorm_moment_args(mu: f64, sigma: f64) -> i64 {
    if is_nan(mu) || is_nan(sigma) { return ERR_NAN; }
    return first_code(check_scalar(mu), check_shape(sigma));
}

fn lognorm_mean(mu: f64 = 0.0, sigma: f64 = 1.0) -> f64, i64 {
    let c = lognorm_moment_args(mu, sigma);
    if c != OK { return nan(), c; }
    return exp_or_overflow(dd_from(mu + 0.5 * sigma * sigma));
}

fn lognorm_var(mu: f64 = 0.0, sigma: f64 = 1.0) -> f64, i64 {
    let c = lognorm_moment_args(mu, sigma);
    if c != OK { return nan(), c; }
    let s2 = sigma * sigma;
    let em, c2 = expm1(s2);
    if is_inf(em) || c2 != OK { return nan(), ERR_OVERFLOW; }
    let e, c3 = exp_or_overflow(dd_from(2.0 * mu + s2));
    if c3 != OK { return nan(), c3; }
    let v = em * e;
    if is_inf(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

fn lognorm_skew(mu: f64 = 0.0, sigma: f64 = 1.0) -> f64, i64 {
    let c = lognorm_moment_args(mu, sigma);
    if c != OK { return nan(), c; }
    let s2 = sigma * sigma;
    let em, c2 = expm1(s2);
    if is_inf(em) || c2 != OK { return nan(), ERR_OVERFLOW; }
    let v = (em + 3.0) * sqrt(em);
    if is_inf(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

fn lognorm_kurt(mu: f64 = 0.0, sigma: f64 = 1.0) -> f64, i64 {
    let c = lognorm_moment_args(mu, sigma);
    if c != OK { return nan(), c; }
    let s2 = sigma * sigma;
    if 4.0 * s2 > 709.0 { return nan(), ERR_OVERFLOW; }
    // with y = e^(s^2) - 1 the kurtosis is y (16 + y (15 + y (6 + y)))
    let y, c2 = expm1(s2);
    let r = y * (16.0 + y * (15.0 + y * (6.0 + y)));
    if is_inf(r) || is_nan(r) { return nan(), ERR_OVERFLOW; }
    return r, OK;
}

// ---- the Weibull family ---------------------------------------------------------------
