// The probability distributions of the package: ten continuous families
// and five discrete families. Each family has a density or a mass, a
// cumulative function, a survival function, a quantile function and the
// four moments (mean, variance, skewness and excess kurtosis).
//
// Every call returns the result first and a fault code last. A fault gives
// NaN as the result. The checks run in one order: NaN, infinity, domain,
// pole. The survival function is the tail itself and never one minus the
// cumulative function. Every loop has a fixed bound.
//
// The names start with the family: norm_, unif_, expon_, lognorm_, weibull_,
// gammadist_, betadist_, chi2_, tdist_, fdist_, binom_, poisson_, geom_,
// nbinom_ and hypergeom_. The names gamma and beta belong to the special
// functions. The trailing parameters of a family have the values of the
// standard family as defaults.
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

// ---- the Student t family -----------------------------------------------------------

fn tdist_pdf(x: f64, df: f64 = 1.0) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    let c = chk_df(df);
    if c == ERR_NAN { return nan(), c; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    if c != OK { return nan(), c; }
    // 1 / (sqrt(df) B(df/2, 1/2)) (1 + x^2/df)^(-(df+1)/2)
    let q = dd_div(dd_mul(dd_from(x), dd_from(x)), dd_from(df));
    // x^2 overflows above 1e154, and the density falls below the smallest subnormal
    if is_inf(q.hi) { return 0.0, OK; }
    let l1p = dd_log(dd_add_d(q, 1.0));
    var l = dd_neg(dd_mul_d(l1p, 0.5 * (df + 1.0)));
    l = dd_sub(l, dd_lnbeta(0.5 * df, 0.5));
    l = dd_sub(l, dd_mul_d(dd_log_d(df), 0.5));
    return exp_or_overflow(l);
}

// The pair (lower, upper) of the Student t tails at x. The tail beyond
// |x| is half the incomplete beta at df/(df+x^2) for |x| above the root of
// df, and the complement of the central part below it. Both forms keep
// the digits of the small tail, because each passes the exact argument.
fn tdist_tails(x: f64, df: f64) -> f64, f64, i64 {
    let ax = abs(x);
    let a2 = dd_mul(dd_from(ax), dd_from(ax));
    let den = dd_add_d(a2, df);
    var tail: f64 = 0.0;
    if ax * ax < df {
        let y = dd_div(a2, den);
        let t, c = betainc_upper(0.5, 0.5 * df, y.hi);
        if c != OK { return nan(), nan(), c; }
        tail = 0.5 * beta_step(t, 0.5, 0.5 * df, y, -1.0);
    } else {
        let w = dd_div(dd_from(df), den);
        let t, c = betainc(0.5 * df, 0.5, w.hi);
        if c != OK { return nan(), nan(), c; }
        tail = 0.5 * beta_step(t, 0.5 * df, 0.5, w, 1.0);
    }
    if x >= 0.0 { return 1.0 - tail, tail, OK; }
    return tail, 1.0 - tail, OK;
}

fn tdist_cdf(x: f64, df: f64 = 1.0) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    let c = chk_df(df);
    if c == ERR_NAN { return nan(), c; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    if c != OK { return nan(), c; }
    let lo, up, c2 = tdist_tails(x, df);
    return lo, c2;
}

fn tdist_sf(x: f64, df: f64 = 1.0) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    let c = chk_df(df);
    if c == ERR_NAN { return nan(), c; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    if c != OK { return nan(), c; }
    let lo, up, c2 = tdist_tails(x, df);
    return up, c2;
}

fn tdist_ppf(p: f64, df: f64 = 1.0) -> f64, i64 {
    if is_nan(p) { return nan(), ERR_NAN; }
    let c = chk_df(df);
    if c != OK { return nan(), c; }
    let c2 = chk_open(p);
    if c2 != OK { return nan(), c2; }
    if p == 0.5 { return 0.0, OK; }
    let tail = min(p, 1.0 - p);
    let d = abs(2.0 * p - 1.0);
    var t: f64 = 0.0;
    if d < 0.5 {
        let y, c3 = betaincinv(0.5, 0.5 * df, d);
        if c3 != OK { return nan(), c3; }
        t = sqrt(df * y / (1.0 - y));
    } else {
        let w, c3 = betaincinv(0.5 * df, 0.5, 2.0 * tail);
        if c3 != OK { return nan(), c3; }
        t = sqrt(df * (1.0 - w) / w);
    }
    if is_inf(t) || is_nan(t) { return nan(), ERR_OVERFLOW; }
    if p < 0.5 { t = -t; }
    return tdist_polish(t, p, df), OK;
}

// One Newton step on the cumulative function removes the error of the
// inverse of the incomplete beta. The step needs a density above 0 and it
// never moves the value by more than a part in a million.
fn tdist_polish(t: f64, p: f64, df: f64) -> f64 {
    let pdf, c1 = tdist_pdf(t, df);
    if c1 != OK || pdf <= 0.0 { return t; }
    var step: f64 = 0.0;
    if p < 0.5 {
        let f, c2 = tdist_cdf(t, df);
        if c2 != OK { return t; }
        step = (f - p) / pdf;
    } else {
        let f, c2 = tdist_sf(t, df);
        if c2 != OK { return t; }
        step = -(f - (1.0 - p)) / pdf;
    }
    if abs(step) > 1e-6 * abs(t) { return t; }
    return t - step;
}

fn tdist_mean(df: f64 = 1.0) -> f64, i64 {
    let c = chk_df(df);
    if c != OK { return nan(), c; }
    if df <= 1.0 { return nan(), ERR_UNDEFINED; }
    return 0.0, OK;
}

fn tdist_var(df: f64 = 1.0) -> f64, i64 {
    let c = chk_df(df);
    if c != OK { return nan(), c; }
    if df <= 2.0 { return nan(), ERR_UNDEFINED; }
    return df / (df - 2.0), OK;
}

fn tdist_skew(df: f64 = 1.0) -> f64, i64 {
    let c = chk_df(df);
    if c != OK { return nan(), c; }
    if df <= 3.0 { return nan(), ERR_UNDEFINED; }
    return 0.0, OK;
}

fn tdist_kurt(df: f64 = 1.0) -> f64, i64 {
    let c = chk_df(df);
    if c != OK { return nan(), c; }
    if df <= 4.0 { return nan(), ERR_UNDEFINED; }
    return 6.0 / (df - 4.0), OK;
}

// ---- the F family ----------------------------------------------------------------------

fn fdist_check(d1: f64, d2: f64) -> i64 {
    if is_nan(d1) || is_nan(d2) { return ERR_NAN; }
    return first_code(chk_df(d1), chk_df(d2));
}

fn fdist_pdf(x: f64, d1: f64 = 1.0, d2: f64 = 1.0) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    let c = fdist_check(d1, d2);
    if c == ERR_NAN { return nan(), c; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    if c != OK { return nan(), c; }
    if x < 0.0 { return 0.0, OK; }
    if x == 0.0 {
        if d1 < 2.0 { return nan(), ERR_OVERFLOW; }
        if d1 == 2.0 { return 1.0, OK; }
        return 0.0, OK;
    }
    // x^(d1/2-1) (d1/d2)^(d1/2) (1 + d1 x/d2)^(-(d1+d2)/2) / B(d1/2, d2/2)
    let r = dd_div(dd_mul(dd_from(d1), dd_from(x)), dd_from(d2));
    var l = dd_mul(dd_log_d(x), dd_from(0.5 * d1 - 1.0));
    l = dd_add(l, dd_mul_d(dd_sub(dd_log_d(d1), dd_log_d(d2)), 0.5 * d1));
    l = dd_sub(l, dd_mul_d(dd_log(dd_add_d(r, 1.0)), 0.5 * (d1 + d2)));
    l = dd_sub(l, dd_lnbeta(0.5 * d1, 0.5 * d2));
    return exp_or_overflow(l);
}

// The pair u = d1 x / (d1 x + d2) and the pair w = d2 / (d1 x + d2).
fn fdist_ratios(x: f64, d1: f64, d2: f64) -> DD, DD {
    let h, l = two_prod(d1, x);
    let t = DD { hi: h, lo: l };
    let s = dd_add_d(t, d2);
    return dd_div(t, s), dd_div(dd_from(d2), s);
}

fn fdist_cdf(x: f64, d1: f64 = 1.0, d2: f64 = 1.0) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    let c = fdist_check(d1, d2);
    if c == ERR_NAN { return nan(), c; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    if c != OK { return nan(), c; }
    if x <= 0.0 { return 0.0, OK; }
    let t = d1 * x;
    if is_inf(t) { return 1.0, OK; }
    let u, w = fdist_ratios(x, d1, d2);
    if t <= d2 {
        let v, c2 = betainc(0.5 * d1, 0.5 * d2, u.hi);
        if c2 != OK { return nan(), c2; }
        return beta_step(v, 0.5 * d1, 0.5 * d2, u, 1.0), OK;
    }
    let v, c2 = betainc_upper(0.5 * d2, 0.5 * d1, w.hi);
    if c2 != OK { return nan(), c2; }
    return beta_step(v, 0.5 * d2, 0.5 * d1, w, -1.0), OK;
}

fn fdist_sf(x: f64, d1: f64 = 1.0, d2: f64 = 1.0) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    let c = fdist_check(d1, d2);
    if c == ERR_NAN { return nan(), c; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    if c != OK { return nan(), c; }
    if x <= 0.0 { return 1.0, OK; }
    let t = d1 * x;
    if is_inf(t) { return 0.0, OK; }
    let u, w = fdist_ratios(x, d1, d2);
    if t >= d2 {
        let v, c2 = betainc(0.5 * d2, 0.5 * d1, w.hi);
        if c2 != OK { return nan(), c2; }
        return beta_step(v, 0.5 * d2, 0.5 * d1, w, 1.0), OK;
    }
    let v, c2 = betainc_upper(0.5 * d1, 0.5 * d2, u.hi);
    if c2 != OK { return nan(), c2; }
    return beta_step(v, 0.5 * d1, 0.5 * d2, u, -1.0), OK;
}

fn fdist_ppf(p: f64, d1: f64 = 1.0, d2: f64 = 1.0) -> f64, i64 {
    if is_nan(p) { return nan(), ERR_NAN; }
    let c = fdist_check(d1, d2);
    if c != OK { return nan(), c; }
    let c2 = chk_open(p);
    if c2 != OK { return nan(), c2; }
    var r: f64 = 0.0;
    if p <= 0.5 {
        let u, c3 = betaincinv(0.5 * d1, 0.5 * d2, p);
        if c3 != OK { return nan(), c3; }
        r = d2 * u / (d1 * (1.0 - u));
    } else {
        let v, c3 = betaincinv(0.5 * d2, 0.5 * d1, 1.0 - p);
        if c3 != OK { return nan(), c3; }
        r = d2 * (1.0 - v) / (d1 * v);
    }
    if is_inf(r) { return nan(), ERR_OVERFLOW; }
    return fdist_polish(r, p, d1, d2), OK;
}

// One Newton step on the cumulative function, as for the Student t.
fn fdist_polish(x: f64, p: f64, d1: f64, d2: f64) -> f64 {
    let pdf, c1 = fdist_pdf(x, d1, d2);
    if c1 != OK || pdf <= 0.0 || is_inf(pdf) { return x; }
    var step: f64 = 0.0;
    if p <= 0.5 {
        let f, c2 = fdist_cdf(x, d1, d2);
        if c2 != OK { return x; }
        step = (f - p) / pdf;
    } else {
        let f, c2 = fdist_sf(x, d1, d2);
        if c2 != OK { return x; }
        step = -(f - (1.0 - p)) / pdf;
    }
    if abs(step) > 1e-6 * abs(x) { return x; }
    return x - step;
}

fn fdist_mean(d1: f64 = 1.0, d2: f64 = 1.0) -> f64, i64 {
    let c = fdist_check(d1, d2);
    if c != OK { return nan(), c; }
    if d2 <= 2.0 { return nan(), ERR_UNDEFINED; }
    return d2 / (d2 - 2.0), OK;
}

fn fdist_var(d1: f64 = 1.0, d2: f64 = 1.0) -> f64, i64 {
    let c = fdist_check(d1, d2);
    if c != OK { return nan(), c; }
    if d2 <= 4.0 { return nan(), ERR_UNDEFINED; }
    let r = d2 / (d2 - 2.0);
    let v = 2.0 * r * r * ((d1 + d2 - 2.0) / d1) / (d2 - 4.0);
    if is_inf(v) || is_nan(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

fn fdist_skew(d1: f64 = 1.0, d2: f64 = 1.0) -> f64, i64 {
    let c = fdist_check(d1, d2);
    if c != OK { return nan(), c; }
    if d2 <= 6.0 { return nan(), ERR_UNDEFINED; }
    let v = (2.0 * d1 + d2 - 2.0) / (d2 - 6.0) * sqrt(8.0 * (d2 - 4.0) / d1 / (d1 + d2 - 2.0));
    if is_inf(v) || is_nan(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

fn fdist_kurt(d1: f64 = 1.0, d2: f64 = 1.0) -> f64, i64 {
    let c = fdist_check(d1, d2);
    if c != OK { return nan(), c; }
    if d2 <= 8.0 { return nan(), ERR_UNDEFINED; }
    let a = (5.0 * d2 - 22.0) / (d2 - 8.0);
    let b = ((d2 - 4.0) / d1) * ((d2 - 2.0) / (d2 - 8.0)) * ((d2 - 2.0) / (d1 + d2 - 2.0));
    let v = 12.0 / (d2 - 6.0) * (a + b);
    if is_inf(v) || is_nan(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

// ---- the discrete families -------------------------------------------------------------------

// The code of a count argument of a discrete family. A count must be a
// whole number at or above 0. A fraction and a negative count lie outside
// the range, and an infinity is outside it too.
fn chk_count(k: f64) -> i64 {
    if is_nan(k) { return ERR_NAN; }
    if is_inf(k) { return ERR_INFINITY; }
    if k < 0.0 || floor(k) != k { return ERR_DOMAIN; }
    if k > MAX_COUNT { return ERR_DOMAIN; }
    return OK;
}

// ---- the binomial family ----

fn binom_check(n: f64, p: f64) -> i64 {
    if is_nan(n) || is_nan(p) { return ERR_NAN; }
    return first_code(chk_whole(n), chk_unit(p));
}

fn binom_pmf(k: f64, n: f64, p: f64 = 0.5) -> f64, i64 {
    let ck = chk_count(k);
    if ck != OK { return nan(), ck; }
    let c = binom_check(n, p);
    if c != OK { return nan(), c; }
    if k > n { return 0.0, OK; }
    if p == 0.0 { if k == 0.0 { return 1.0, OK; } return 0.0, OK; }
    if p == 1.0 { if k == n { return 1.0, OK; } return 0.0, OK; }
    var l = dd_lchoose(n, k);
    l = dd_add(l, dd_mul_d(dd_log_d(p), k));
    l = dd_add(l, dd_mul_d(dd_log_1m(p), n - k));
    return exp_or_overflow(l);
}

fn binom_cdf(k: f64, n: f64, p: f64 = 0.5) -> f64, i64 {
    let ck = chk_count(k);
    if ck != OK { return nan(), ck; }
    let c = binom_check(n, p);
    if c != OK { return nan(), c; }
    if k >= n { return 1.0, OK; }
    if p == 0.0 { return 1.0, OK; }
    if p == 1.0 { return 0.0, OK; }
    let kf = floor(k);
    return betainc_upper(kf + 1.0, n - kf, p);
}

fn binom_sf(k: f64, n: f64, p: f64 = 0.5) -> f64, i64 {
    let ck = chk_count(k);
    if ck != OK { return nan(), ck; }
    let c = binom_check(n, p);
    if c != OK { return nan(), c; }
    if k >= n { return 0.0, OK; }
    if p == 0.0 { return 0.0, OK; }
    if p == 1.0 { return 1.0, OK; }
    let kf = floor(k);
    return betainc(kf + 1.0, n - kf, p);
}

// ---- the Poisson family ----

fn poisson_pmf(k: f64, rate: f64 = 1.0) -> f64, i64 {
    let ck = chk_count(k);
    if ck != OK { return nan(), ck; }
    let c = check_shape(rate);
    if c != OK { return nan(), c; }
    let l = dd_sub(dd_sub(dd_mul_d(dd_log_d(rate), k), dd_from(rate)), dd_lngamma(dd_from(k + 1.0)));
    return exp_or_overflow(l);
}

fn poisson_cdf(k: f64, rate: f64 = 1.0) -> f64, i64 {
    let ck = chk_count(k);
    if ck != OK { return nan(), ck; }
    let c = check_shape(rate);
    if c != OK { return nan(), c; }
    return gammaincc(floor(k) + 1.0, rate);
}

fn poisson_sf(k: f64, rate: f64 = 1.0) -> f64, i64 {
    let ck = chk_count(k);
    if ck != OK { return nan(), ck; }
    let c = check_shape(rate);
    if c != OK { return nan(), c; }
    return gammainc(floor(k) + 1.0, rate);
}

// ---- the geometric family (the trials up to the first success) ----

fn geom_pmf(k: f64, p: f64 = 0.5) -> f64, i64 {
    let ck = chk_count(k);
    if ck != OK { return nan(), ck; }
    let c = chk_succ(p);
    if c != OK { return nan(), c; }
    if k < 1.0 { return 0.0, OK; }
    if p == 1.0 { if k == 1.0 { return 1.0, OK; } return 0.0, OK; }
    let l = dd_add(dd_mul_d(dd_log_1m(p), k - 1.0), dd_log_d(p));
    return exp_or_overflow(l);
}

fn geom_cdf(k: f64, p: f64 = 0.5) -> f64, i64 {
    let ck = chk_count(k);
    if ck != OK { return nan(), ck; }
    let c = chk_succ(p);
    if c != OK { return nan(), c; }
    if k < 1.0 { return 0.0, OK; }
    if p == 1.0 { return 1.0, OK; }
    let l = dd_mul_d(dd_log_1m(p), floor(k));
    if l.hi < -0.5 { return 1.0 - expf_dd(l), OK; }
    let v, c2 = expm1(dd_to_f64(l));
    return -v, OK;
}

fn geom_sf(k: f64, p: f64 = 0.5) -> f64, i64 {
    let ck = chk_count(k);
    if ck != OK { return nan(), ck; }
    let c = chk_succ(p);
    if c != OK { return nan(), c; }
    if k < 1.0 { return 1.0, OK; }
    if p == 1.0 { return 0.0, OK; }
    return expf_dd(dd_mul_d(dd_log_1m(p), floor(k))), OK;
}

// ---- the negative binomial family (the failures before the r-th success) ----

fn nbinom_check(r: f64, p: f64) -> i64 {
    if is_nan(r) || is_nan(p) { return ERR_NAN; }
    return first_code(check_shape(r), chk_succ(p));
}

fn nbinom_pmf(k: f64, r: f64, p: f64 = 0.5) -> f64, i64 {
    let ck = chk_count(k);
    if ck != OK { return nan(), ck; }
    let c = nbinom_check(r, p);
    if c != OK { return nan(), c; }
    if p == 1.0 { if k == 0.0 { return 1.0, OK; } return 0.0, OK; }
    // Gamma(k + r) / (Gamma(r) k!) = 1 / ((k + r) B(r, k + 1))
    var l = dd_neg(dd_log(dd_add_d(dd_from(r), k)));
    l = dd_sub(l, dd_lnbeta(r, k + 1.0));
    l = dd_add(l, dd_mul_d(dd_log_d(p), r));
    l = dd_add(l, dd_mul_d(dd_log_1m(p), k));
    return exp_or_overflow(l);
}

fn nbinom_cdf(k: f64, r: f64, p: f64 = 0.5) -> f64, i64 {
    let ck = chk_count(k);
    if ck != OK { return nan(), ck; }
    let c = nbinom_check(r, p);
    if c != OK { return nan(), c; }
    if p == 1.0 { return 1.0, OK; }
    return betainc(r, floor(k) + 1.0, p);
}

fn nbinom_sf(k: f64, r: f64, p: f64 = 0.5) -> f64, i64 {
    let ck = chk_count(k);
    if ck != OK { return nan(), ck; }
    let c = nbinom_check(r, p);
    if c != OK { return nan(), c; }
    if p == 1.0 { return 0.0, OK; }
    return betainc_upper(r, floor(k) + 1.0, p);
}

// ---- the hypergeometric family (population m, n successes, d draws) ----

fn hyper_check(m: f64, n: f64, d: f64) -> i64 {
    if is_nan(m) || is_nan(n) || is_nan(d) { return ERR_NAN; }
    let c = first_code(first_code(chk_whole(m), chk_whole(n)), chk_whole(d));
    if c != OK { return c; }
    if m < 1.0 || n > m || d > m || m > MAX_SHAPE { return ERR_DOMAIN; }
    return OK;
}

// The log of the mass at a whole k inside the support, as a pair.
fn hyper_logpmf(k: f64, m: f64, n: f64, d: f64) -> DD {
    let a = dd_lchoose(n, k);
    let b = dd_lchoose(m - n, d - k);
    let t = dd_lchoose(m, d);
    return dd_sub(dd_add(a, b), t);
}

fn hyper_kmin(m: f64, n: f64, d: f64) -> f64 { return max(0.0, d - (m - n)); }
fn hyper_kmax(n: f64, d: f64) -> f64 { return min(n, d); }

fn hypergeom_pmf(k: f64, m: f64, n: f64, d: f64) -> f64, i64 {
    let ck = chk_count(k);
    if ck != OK { return nan(), ck; }
    let c = hyper_check(m, n, d);
    if c != OK { return nan(), c; }
    if k < hyper_kmin(m, n, d) || k > hyper_kmax(n, d) { return 0.0, OK; }
    return expf_dd(hyper_logpmf(k, m, n, d)), OK;
}

// The mode of the hypergeometric family, the count with the largest
// mass.
fn hyper_mode(m: f64, n: f64, d: f64) -> f64 {
    let md = floor((d + 1.0) * (n + 1.0) / (m + 2.0));
    return min(max(md, hyper_kmin(m, n, d)), hyper_kmax(n, d));
}

// The mass at the mode as a pair. It is never below the reciprocal of the
// length of the support, so the pair holds no underflow.
fn hyper_top(m: f64, n: f64, d: f64, mode: f64) -> DD {
    let m2, k2 = dd_exp_parts(hyper_logpmf(mode, m, n, d));
    return dd_scale_dd(m2, k2);
}

// The ratio of the mass at k + 1 to the mass at k as a pair.
fn hyper_up(t: DD, k: f64, m: f64, n: f64, d: f64) -> DD {
    var r = dd_mul_d(dd_mul_d(t, n - k), d - k);
    r = dd_div_d(r, k + 1.0);
    return dd_div_d(r, m - n - d + k + 1.0);
}

// The ratio of the mass at k - 1 to the mass at k as a pair.
fn hyper_down(t: DD, k: f64, m: f64, n: f64, d: f64) -> DD {
    var r = dd_mul_d(dd_mul_d(t, k), m - n - d + k);
    r = dd_div_d(r, n - k + 1.0);
    return dd_div_d(r, d - k + 1.0);
}

// The sums of the masses at and below k, and above k. Both come straight
// from the masses, so a small tail keeps its digits. The support length
// has a bound.
fn hyper_tails(k: f64, m: f64, n: f64, d: f64) -> f64, f64, i64 {
    let kmin = hyper_kmin(m, n, d);
    let kmax = hyper_kmax(n, d);
    if k < kmin { return 0.0, 1.0, OK; }
    if k >= kmax { return 1.0, 0.0, OK; }
    if kmax - kmin + 1.0 > MAX_DIST_SUPPORT { return nan(), nan(), ERR_TOO_LONG; }
    let mode = hyper_mode(m, n, d);
    var t = hyper_top(m, n, d, mode);
    var low = dd_from(0.0);
    var up = dd_from(0.0);
    var j = mode;
    var cur = t;
    while j >= kmin {
        if j <= k { low = dd_add(low, cur); } else { up = dd_add(up, cur); }
        if j == kmin || cur.hi == 0.0 { break; }
        cur = hyper_down(cur, j, m, n, d);
        j -= 1.0;
    }
    j = mode;
    cur = t;
    while j < kmax {
        cur = hyper_up(cur, j, m, n, d);
        j += 1.0;
        if cur.hi == 0.0 { break; }
        if j <= k { low = dd_add(low, cur); } else { up = dd_add(up, cur); }
    }
    return dd_to_f64(low), dd_to_f64(up), OK;
}

fn hypergeom_cdf(k: f64, m: f64, n: f64, d: f64) -> f64, i64 {
    let ck = chk_count(k);
    if ck != OK { return nan(), ck; }
    let c = hyper_check(m, n, d);
    if c != OK { return nan(), c; }
    let lo, up, c2 = hyper_tails(floor(k), m, n, d);
    return lo, c2;
}

fn hypergeom_sf(k: f64, m: f64, n: f64, d: f64) -> f64, i64 {
    let ck = chk_count(k);
    if ck != OK { return nan(), ck; }
    let c = hyper_check(m, n, d);
    if c != OK { return nan(), c; }
    let lo, up, c2 = hyper_tails(floor(k), m, n, d);
    return up, c2;
}

// ---- the quantile of a discrete family ----------------------------------------------------

// The cumulative function of the family number fam at the count k.
// The numbers of the discrete families in the shared helpers.
const D_BINOM: i64 = 0;
const D_POISSON: i64 = 1;
const D_GEOM: i64 = 2;
const D_NBINOM: i64 = 3;
const D_HYPERGEOM: i64 = 4;

fn disc_cdf(fam: i64, k: f64, a: f64, b: f64, c: f64) -> f64, i64 {
    if fam == D_BINOM { return binom_cdf(k, a, b); }
    if fam == D_POISSON { return poisson_cdf(k, a); }
    if fam == D_GEOM { return geom_cdf(k, a); }
    if fam == D_NBINOM { return nbinom_cdf(k, a, b); }
    return hypergeom_cdf(k, a, b, c);
}

// The smallest count whose cumulative function reaches p. The search
// starts at lo and doubles its upper bound until the cumulative function
// reaches p, so a family without an end needs no bound of its own. The
// count of the steps has a bound, and a search that meets its bound gives
// a fault. A probability of 0 gives the smallest count of the support.
fn disc_ppf(fam: i64, p: f64, a: f64, b: f64, c: f64, lo: f64, hi: f64) -> f64, i64 {
    let cp = chk_unit(p);
    if cp != OK { return nan(), cp; }
    if p == 0.0 { return lo, OK; }
    var top = hi;
    if is_inf(hi) {
        top = lo + 1.0;
        var i: i64 = 0;
        var found = false;
        while i < MAX_ITER_SEARCH {
            let v, cc = disc_cdf(fam, top, a, b, c);
            if cc != OK { return nan(), cc; }
            if v >= p { found = true; break; }
            if top > MAX_COUNT * 0.5 { return nan(), ERR_OVERFLOW; }
            top = top * 2.0;
            i += 1;
        }
        if !found { return nan(), ERR_OVERFLOW; }
    } else {
        let v, cc = disc_cdf(fam, top, a, b, c);
        if cc != OK { return nan(), cc; }
        if v < p { return nan(), ERR_OVERFLOW; }
    }
    var l = lo;
    var r = top;
    var steps: i64 = 0;
    while l < r && steps < MAX_ITER_SEARCH {
        let mid = floor(0.5 * (l + r));
        let v, cc = disc_cdf(fam, mid, a, b, c);
        if cc != OK { return nan(), cc; }
        if v >= p { r = mid; } else { l = mid + 1.0; }
        steps += 1;
    }
    if l < r { return nan(), ERR_OVERFLOW; }
    return l, OK;
}

fn binom_ppf(p: f64, n: f64, prob: f64 = 0.5) -> f64, i64 {
    if is_nan(p) { return nan(), ERR_NAN; }
    let c = binom_check(n, prob);
    if c != OK { return nan(), c; }
    return disc_ppf(D_BINOM, p, n, prob, 0.0, 0.0, n);
}

fn poisson_ppf(p: f64, rate: f64 = 1.0) -> f64, i64 {
    if is_nan(p) || is_nan(rate) { return nan(), ERR_NAN; }
    let c = check_shape(rate);
    if c != OK { return nan(), c; }
    return disc_ppf(D_POISSON, p, rate, 0.0, 0.0, 0.0, 1.0 / 0.0);
}

fn geom_ppf(p: f64, prob: f64 = 0.5) -> f64, i64 {
    if is_nan(p) || is_nan(prob) { return nan(), ERR_NAN; }
    let c = chk_succ(prob);
    if c != OK { return nan(), c; }
    return disc_ppf(D_GEOM, p, prob, 0.0, 0.0, 1.0, 1.0 / 0.0);
}

fn nbinom_ppf(p: f64, r: f64, prob: f64 = 0.5) -> f64, i64 {
    if is_nan(p) { return nan(), ERR_NAN; }
    let c = nbinom_check(r, prob);
    if c != OK { return nan(), c; }
    if prob == 1.0 { return 0.0, OK; }
    return disc_ppf(D_NBINOM, p, r, prob, 0.0, 0.0, 1.0 / 0.0);
}

fn hypergeom_ppf(p: f64, m: f64, n: f64, d: f64) -> f64, i64 {
    if is_nan(p) { return nan(), ERR_NAN; }
    let c = hyper_check(m, n, d);
    if c != OK { return nan(), c; }
    return disc_ppf(D_HYPERGEOM, p, m, n, d, hyper_kmin(m, n, d), hyper_kmax(n, d));
}

// ---- the moments of the discrete families --------------------------------------------------

fn binom_mean(n: f64, p: f64 = 0.5) -> f64, i64 {
    let c = binom_check(n, p);
    if c != OK { return nan(), c; }
    return n * p, OK;
}

fn binom_var(n: f64, p: f64 = 0.5) -> f64, i64 {
    let c = binom_check(n, p);
    if c != OK { return nan(), c; }
    return n * p * (1.0 - p), OK;
}

fn binom_skew(n: f64, p: f64 = 0.5) -> f64, i64 {
    let c = binom_check(n, p);
    if c != OK { return nan(), c; }
    let v = n * p * (1.0 - p);
    if v == 0.0 { return nan(), ERR_UNDEFINED; }
    return (1.0 - 2.0 * p) / sqrt(v), OK;
}

fn binom_kurt(n: f64, p: f64 = 0.5) -> f64, i64 {
    let c = binom_check(n, p);
    if c != OK { return nan(), c; }
    let v = n * p * (1.0 - p);
    if v == 0.0 { return nan(), ERR_UNDEFINED; }
    return (1.0 - 6.0 * p * (1.0 - p)) / v, OK;
}

fn poisson_mean(rate: f64 = 1.0) -> f64, i64 {
    let c = check_shape(rate);
    if c != OK { return nan(), c; }
    return rate, OK;
}

fn poisson_var(rate: f64 = 1.0) -> f64, i64 {
    let c = check_shape(rate);
    if c != OK { return nan(), c; }
    return rate, OK;
}

fn poisson_skew(rate: f64 = 1.0) -> f64, i64 {
    let c = check_shape(rate);
    if c != OK { return nan(), c; }
    return 1.0 / sqrt(rate), OK;
}

fn poisson_kurt(rate: f64 = 1.0) -> f64, i64 {
    let c = check_shape(rate);
    if c != OK { return nan(), c; }
    let v = 1.0 / rate;
    if is_inf(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

fn geom_mean(p: f64 = 0.5) -> f64, i64 {
    let c = chk_succ(p);
    if c != OK { return nan(), c; }
    let v = 1.0 / p;
    if is_inf(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

fn geom_var(p: f64 = 0.5) -> f64, i64 {
    let c = chk_succ(p);
    if c != OK { return nan(), c; }
    let v = (1.0 - p) / (p * p);
    if is_inf(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

fn geom_skew(p: f64 = 0.5) -> f64, i64 {
    let c = chk_succ(p);
    if c != OK { return nan(), c; }
    if p == 1.0 { return nan(), ERR_UNDEFINED; }
    return (2.0 - p) / sqrt(1.0 - p), OK;
}

fn geom_kurt(p: f64 = 0.5) -> f64, i64 {
    let c = chk_succ(p);
    if c != OK { return nan(), c; }
    if p == 1.0 { return nan(), ERR_UNDEFINED; }
    return 6.0 + p * p / (1.0 - p), OK;
}

fn nbinom_mean(r: f64, p: f64 = 0.5) -> f64, i64 {
    let c = nbinom_check(r, p);
    if c != OK { return nan(), c; }
    let v = r * (1.0 - p) / p;
    if is_inf(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

fn nbinom_var(r: f64, p: f64 = 0.5) -> f64, i64 {
    let c = nbinom_check(r, p);
    if c != OK { return nan(), c; }
    let v = r * (1.0 - p) / (p * p);
    if is_inf(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

fn nbinom_skew(r: f64, p: f64 = 0.5) -> f64, i64 {
    let c = nbinom_check(r, p);
    if c != OK { return nan(), c; }
    if p == 1.0 { return nan(), ERR_UNDEFINED; }
    return (2.0 - p) / sqrt(r * (1.0 - p)), OK;
}

fn nbinom_kurt(r: f64, p: f64 = 0.5) -> f64, i64 {
    let c = nbinom_check(r, p);
    if c != OK { return nan(), c; }
    if p == 1.0 { return nan(), ERR_UNDEFINED; }
    return 6.0 / r + p * p / (r * (1.0 - p)), OK;
}

fn hypergeom_mean(m: f64, n: f64, d: f64) -> f64, i64 {
    let c = hyper_check(m, n, d);
    if c != OK { return nan(), c; }
    return n * d / m, OK;
}

// The variance d n (m - n) (m - d) / (m^2 (m - 1)) as a pair.
fn hyper_var_dd(m: f64, n: f64, d: f64) -> DD {
    if m == 1.0 { return dd_from(0.0); }
    let num = dd_mul_d(dd_mul_d(dd_mul_d(dd_from(d), n), m - n), m - d);
    let den = dd_mul_d(dd_mul_d(dd_from(m), m), m - 1.0);
    return dd_div(num, den);
}

fn hypergeom_var(m: f64, n: f64, d: f64) -> f64, i64 {
    let c = hyper_check(m, n, d);
    if c != OK { return nan(), c; }
    return dd_to_f64(hyper_var_dd(m, n, d)), OK;
}

// The third and the fourth central moment as sums of the masses of the
// support, in pairs. The sums run outward from the mode, so the large
// terms come first and the loop ends on the bound of the support.
fn hyper_central(m: f64, n: f64, d: f64) -> DD, DD, i64 {
    let kmin = hyper_kmin(m, n, d);
    let kmax = hyper_kmax(n, d);
    if kmax - kmin + 1.0 > MAX_DIST_SUPPORT { return dd_from(nan()), dd_from(nan()), ERR_TOO_LONG; }
    let mean = dd_div_d(dd_mul_d(dd_from(n), d), m);
    let mode = hyper_mode(m, n, d);
    let t = hyper_top(m, n, d, mode);
    var m3 = dd_from(0.0);
    var m4 = dd_from(0.0);
    var j = mode;
    var cur = t;
    while j >= kmin {
        let e = dd_sub(dd_from(j), mean);
        let e2 = dd_mul(e, e);
        m3 = dd_add(m3, dd_mul(cur, dd_mul(e, e2)));
        m4 = dd_add(m4, dd_mul(cur, dd_mul(e2, e2)));
        if j == kmin || cur.hi == 0.0 { break; }
        cur = hyper_down(cur, j, m, n, d);
        j -= 1.0;
    }
    j = mode;
    cur = t;
    while j < kmax {
        cur = hyper_up(cur, j, m, n, d);
        j += 1.0;
        if cur.hi == 0.0 { break; }
        let e = dd_sub(dd_from(j), mean);
        let e2 = dd_mul(e, e);
        m3 = dd_add(m3, dd_mul(cur, dd_mul(e, e2)));
        m4 = dd_add(m4, dd_mul(cur, dd_mul(e2, e2)));
    }
    return m3, m4, OK;
}

fn hypergeom_skew(m: f64, n: f64, d: f64) -> f64, i64 {
    let c = hyper_check(m, n, d);
    if c != OK { return nan(), c; }
    let v = hyper_var_dd(m, n, d);
    if v.hi == 0.0 { return nan(), ERR_UNDEFINED; }
    let m3, m4, c2 = hyper_central(m, n, d);
    if c2 != OK { return nan(), c2; }
    return dd_to_f64(dd_div(m3, dd_mul(v, dd_sqrt(v)))), OK;
}

fn hypergeom_kurt(m: f64, n: f64, d: f64) -> f64, i64 {
    let c = hyper_check(m, n, d);
    if c != OK { return nan(), c; }
    let v = hyper_var_dd(m, n, d);
    if v.hi == 0.0 { return nan(), ERR_UNDEFINED; }
    let m3, m4, c2 = hyper_central(m, n, d);
    if c2 != OK { return nan(), c2; }
    let v2 = dd_mul(v, v);
    return dd_to_f64(dd_div(dd_sub(m4, dd_mul_d(v2, 3.0)), v2)), OK;
}
