// Eight continuous families whose formulas use only elementary functions:
// anglit, arcsine, cosine, gibrat, halfcauchy, halflogistic, halfnorm and
// hypsecant. The names and the parameters follow SciPy.
//
// Every family has eight calls: pdf, cdf, sf, ppf, mean, var, skew and
// kurt. The kurtosis is the excess kurtosis. Every call returns the result
// first and a fault code last. A value outside the support gives a density
// of 0. The quantile of 0 or 1 gives the bound of the support. A bound
// that is infinite does not fit in a result, so the call gives
// ERR_OVERFLOW. The half Cauchy law has no moments, so its moment calls
// give ERR_UNDEFINED.
namespace statistics;

import std;
import math;


// The base of the natural logarithm. The old release has no constant for it.
const EULER_E: f64 = 2.718281828459045;

// The check of a value: a NaN or an infinity is a fault.
fn fam_x(x: f64) -> i64 { return check_scalar(x); }

// The bound of the support for a probability of 0 or 1. The first result
// tells whether the probability was an end point.
fn fam_end(p: f64, lo: f64, hi: f64) -> bool, f64, i64 {
    if p == 0.0 {
        if is_inf(lo) { return true, nan(), ERR_OVERFLOW; }
        return true, lo, OK;
    }
    if p == 1.0 {
        if is_inf(hi) { return true, nan(), ERR_OVERFLOW; }
        return true, hi, OK;
    }
    return false, 0.0, OK;
}

// The excess kurtosis and the skewness of a law from its first four raw
// moments.
fn fam_skew_raw(m1: f64, m2: f64, m3: f64) -> f64 {
    let v = m2 - m1 * m1;
    return (m3 - 3.0 * m1 * m2 + 2.0 * m1 * m1 * m1) / (v * sqrt(v));
}

fn fam_kurt_raw(m1: f64, m2: f64, m3: f64, m4: f64) -> f64 {
    let v = m2 - m1 * m1;
    let m1s = m1 * m1;
    return (m4 - 4.0 * m1 * m3 + 6.0 * m1s * m2 - 3.0 * m1s * m1s) / (v * v) - 3.0;
}

// ---- anglit: the density cos(2 x) on -pi/4 to pi/4 ---------------------

fn anglit_pdf(x: f64) -> f64, i64 {
    let c = fam_x(x);
    if c != OK { return nan(), c; }
    if x < -PI / 4.0 || x > PI / 4.0 { return 0.0, OK; }
    return cos(2.0 * x), OK;
}

fn anglit_cdf(x: f64) -> f64, i64 {
    let c = fam_x(x);
    if c != OK { return nan(), c; }
    if x <= -PI / 4.0 { return 0.0, OK; }
    if x >= PI / 4.0 { return 1.0, OK; }
    let s = sin(x + PI / 4.0);
    return s * s, OK;
}

fn anglit_sf(x: f64) -> f64, i64 {
    let c = fam_x(x);
    if c != OK { return nan(), c; }
    if x <= -PI / 4.0 { return 1.0, OK; }
    if x >= PI / 4.0 { return 0.0, OK; }
    let s = cos(x + PI / 4.0);
    return s * s, OK;
}

fn anglit_ppf(p: f64) -> f64, i64 {
    let c = check_prob(p);
    if c != OK { return nan(), c; }
    let e, ev, ec = fam_end(p, -PI / 4.0, PI / 4.0);
    if e { return ev, ec; }
    return asin(sqrt(p)) - PI / 4.0, OK;
}

fn anglit_mean() -> f64, i64 { return 0.0, OK; }
fn anglit_var() -> f64, i64 { return PI * PI / 16.0 - 0.5, OK; }
fn anglit_skew() -> f64, i64 { return 0.0, OK; }
fn anglit_kurt() -> f64, i64 {
    let d = PI * PI - 8.0;
    return -2.0 * (PI * PI * PI * PI - 96.0) / (d * d), OK;
}

// ---- arcsine: the density 1 / (pi sqrt(x (1 - x))) on 0 to 1 --------------

fn arcsine_pdf(x: f64) -> f64, i64 {
    let c = fam_x(x);
    if c != OK { return nan(), c; }
    if x < 0.0 || x > 1.0 { return 0.0, OK; }
    if x == 0.0 || x == 1.0 { return nan(), ERR_OVERFLOW; }
    return 1.0 / (PI * sqrt(x * (1.0 - x))), OK;
}

fn arcsine_cdf(x: f64) -> f64, i64 {
    let c = fam_x(x);
    if c != OK { return nan(), c; }
    if x <= 0.0 { return 0.0, OK; }
    if x >= 1.0 { return 1.0, OK; }
    return 2.0 / PI * asin(sqrt(x)), OK;
}

fn arcsine_sf(x: f64) -> f64, i64 {
    let c = fam_x(x);
    if c != OK { return nan(), c; }
    if x <= 0.0 { return 1.0, OK; }
    if x >= 1.0 { return 0.0, OK; }
    return 2.0 / PI * asin(sqrt(1.0 - x)), OK;
}

fn arcsine_ppf(p: f64) -> f64, i64 {
    let c = check_prob(p);
    if c != OK { return nan(), c; }
    let s = sin(PI / 2.0 * p);
    return s * s, OK;
}

fn arcsine_mean() -> f64, i64 { return 0.5, OK; }
fn arcsine_var() -> f64, i64 { return 0.125, OK; }
fn arcsine_skew() -> f64, i64 { return 0.0, OK; }
fn arcsine_kurt() -> f64, i64 { return -1.5, OK; }

// ---- cosine: the density (1 + cos x) / (2 pi) on -pi to pi -----------------

fn cosine_pdf(x: f64) -> f64, i64 {
    let c = fam_x(x);
    if c != OK { return nan(), c; }
    if x < -PI || x > PI { return 0.0, OK; }
    return (1.0 + cos(x)) / (2.0 * PI), OK;
}

fn cosine_cdf(x: f64) -> f64, i64 {
    let c = fam_x(x);
    if c != OK { return nan(), c; }
    if x <= -PI { return 0.0, OK; }
    if x >= PI { return 1.0, OK; }
    return (PI + x + sin(x)) / (2.0 * PI), OK;
}

fn cosine_sf(x: f64) -> f64, i64 {
    let c = fam_x(x);
    if c != OK { return nan(), c; }
    if x <= -PI { return 1.0, OK; }
    if x >= PI { return 0.0, OK; }
    return (PI - x - sin(x)) / (2.0 * PI), OK;
}

// The cdf has no closed inverse. The quantile halves the support 100
// times, which is below the spacing of a float for any value in it.
fn cosine_ppf(p: f64) -> f64, i64 {
    let c = check_prob(p);
    if c != OK { return nan(), c; }
    let e, ev, ec = fam_end(p, -PI, PI);
    if e { return ev, ec; }
    // The law is symmetric, so the median is the centre and the lower half
    // carries the search.
    if p == 0.5 { return 0.0, OK; }
    let lower = p <= 0.5;
    let t = if lower { p } else { 1.0 - p };
    let x, cx = brent_root(-PI, 0.0, 1e-15, 4e-16) { v => cosine_cdf_value(v) - t };
    if cx != OK { return nan(), cx; }
    if lower { return x, OK; }
    return -x, OK;
}

// The mass below a point of the cosine law, with the cancellation near
// -pi cut by the series of (u - sin u) for u = x + pi.
fn cosine_cdf_value(x: f64) -> f64 {
    let u = x + PI;
    return (u - sin(u)) / (2.0 * PI);
}

fn cosine_mean() -> f64, i64 { return 0.0, OK; }
fn cosine_var() -> f64, i64 { return PI * PI / 3.0 - 2.0, OK; }
fn cosine_skew() -> f64, i64 { return 0.0, OK; }
fn cosine_kurt() -> f64, i64 {
    let d = PI * PI - 6.0;
    return -6.0 * (PI * PI * PI * PI - 90.0) / (5.0 * d * d), OK;
}

// ---- gibrat: the log normal law with a log scale of 1 ---------------------

fn gibrat_pdf(x: f64) -> f64, i64 {
    let c = fam_x(x);
    if c != OK { return nan(), c; }
    if x <= 0.0 { return 0.0, OK; }
    let l = log(x);
    return exp(-0.5 * l * l) / (x * SQRT_2PI), OK;
}

fn gibrat_cdf(x: f64) -> f64, i64 {
    let c = fam_x(x);
    if c != OK { return nan(), c; }
    if x <= 0.0 { return 0.0, OK; }
    return norm_cdf(log(x));
}

fn gibrat_sf(x: f64) -> f64, i64 {
    let c = fam_x(x);
    if c != OK { return nan(), c; }
    if x <= 0.0 { return 1.0, OK; }
    return norm_sf(log(x));
}

fn gibrat_ppf(p: f64) -> f64, i64 {
    let c = check_prob(p);
    if c != OK { return nan(), c; }
    let e, ev, ec = fam_end(p, 0.0, inf());
    if e { return ev, ec; }
    let r = exp(norm_std_ppf(p));
    if is_inf(r) { return nan(), ERR_OVERFLOW; }
    return r, OK;
}

fn gibrat_mean() -> f64, i64 { return sqrt(EULER_E), OK; }
fn gibrat_var() -> f64, i64 { return EULER_E * (EULER_E - 1.0), OK; }
fn gibrat_skew() -> f64, i64 { return sqrt(EULER_E - 1.0) * (2.0 + EULER_E), OK; }
fn gibrat_kurt() -> f64, i64 {
    let e2 = EULER_E * EULER_E;
    return e2 * e2 + 2.0 * e2 * EULER_E + 3.0 * e2 - 6.0, OK;
}

// ---- halfcauchy: the density 2 / (pi (1 + x^2)) on 0 to inf ---------------

fn halfcauchy_pdf(x: f64) -> f64, i64 {
    let c = fam_x(x);
    if c != OK { return nan(), c; }
    if x < 0.0 { return 0.0, OK; }
    return 2.0 / (PI * (1.0 + x * x)), OK;
}

fn halfcauchy_cdf(x: f64) -> f64, i64 {
    let c = fam_x(x);
    if c != OK { return nan(), c; }
    if x <= 0.0 { return 0.0, OK; }
    return 2.0 / PI * atan(x), OK;
}

fn halfcauchy_sf(x: f64) -> f64, i64 {
    let c = fam_x(x);
    if c != OK { return nan(), c; }
    if x <= 0.0 { return 1.0, OK; }
    return 2.0 / PI * atan2(1.0, x), OK;
}

fn halfcauchy_ppf(p: f64) -> f64, i64 {
    let c = check_prob(p);
    if c != OK { return nan(), c; }
    let e, ev, ec = fam_end(p, 0.0, inf());
    if e { return ev, ec; }
    return tan(PI / 2.0 * p), OK;
}

fn halfcauchy_mean() -> f64, i64 { return nan(), ERR_UNDEFINED; }
fn halfcauchy_var() -> f64, i64 { return nan(), ERR_UNDEFINED; }
fn halfcauchy_skew() -> f64, i64 { return nan(), ERR_UNDEFINED; }
fn halfcauchy_kurt() -> f64, i64 { return nan(), ERR_UNDEFINED; }

// ---- halflogistic: the density 2 exp(-x) / (1 + exp(-x))^2 on 0 to inf ----

fn halflogistic_pdf(x: f64) -> f64, i64 {
    let c = fam_x(x);
    if c != OK { return nan(), c; }
    if x < 0.0 { return 0.0, OK; }
    let e = exp(-x);
    return 2.0 * e / ((1.0 + e) * (1.0 + e)), OK;
}

fn halflogistic_cdf(x: f64) -> f64, i64 {
    let c = fam_x(x);
    if c != OK { return nan(), c; }
    if x <= 0.0 { return 0.0, OK; }
    let m, cm = expm1(-x);
    return -m / (1.0 + exp(-x)), cm;
}

fn halflogistic_sf(x: f64) -> f64, i64 {
    let c = fam_x(x);
    if c != OK { return nan(), c; }
    if x <= 0.0 { return 1.0, OK; }
    return 2.0 / (1.0 + exp(x)), OK;
}

fn halflogistic_ppf(p: f64) -> f64, i64 {
    let c = check_prob(p);
    if c != OK { return nan(), c; }
    let e, ev, ec = fam_end(p, 0.0, inf());
    if e { return ev, ec; }
    let l, cl = log1p(2.0 * p / (1.0 - p));
    return l, cl;
}

fn halflogistic_mean() -> f64, i64 { return 2.0 * log(2.0), OK; }
fn halflogistic_var() -> f64, i64 {
    let l = log(2.0);
    return PI * PI / 3.0 - 4.0 * l * l, OK;
}
fn halflogistic_skew() -> f64, i64 {
    return fam_skew_raw(2.0 * log(2.0), PI * PI / 3.0, 9.0 * ZETA3), OK;
}
fn halflogistic_kurt() -> f64, i64 {
    let p4 = PI * PI * PI * PI;
    return fam_kurt_raw(2.0 * log(2.0), PI * PI / 3.0, 9.0 * ZETA3, 7.0 * p4 / 15.0), OK;
}

// ---- halfnorm: the density sqrt(2 / pi) exp(-x^2 / 2) on 0 to inf ----------

fn halfnorm_pdf(x: f64) -> f64, i64 {
    let c = fam_x(x);
    if c != OK { return nan(), c; }
    if x < 0.0 { return 0.0, OK; }
    return 2.0 * exp(-0.5 * x * x) / SQRT_2PI, OK;
}

fn halfnorm_cdf(x: f64) -> f64, i64 {
    let c = fam_x(x);
    if c != OK { return nan(), c; }
    if x <= 0.0 { return 0.0, OK; }
    return erf(x * INV_SQRT2);
}

fn halfnorm_sf(x: f64) -> f64, i64 {
    let c = fam_x(x);
    if c != OK { return nan(), c; }
    if x <= 0.0 { return 1.0, OK; }
    return erfc(x * INV_SQRT2);
}

fn halfnorm_ppf(p: f64) -> f64, i64 {
    let c = check_prob(p);
    if c != OK { return nan(), c; }
    let e, ev, ec = fam_end(p, 0.0, inf());
    if e { return ev, ec; }
    return norm_std_ppf(0.5 + 0.5 * p), OK;
}

fn halfnorm_mean() -> f64, i64 { return sqrt(2.0 / PI), OK; }
fn halfnorm_var() -> f64, i64 { return 1.0 - 2.0 / PI, OK; }
fn halfnorm_skew() -> f64, i64 {
    let d = PI - 2.0;
    return sqrt(2.0) * (4.0 - PI) / (d * sqrt(d)), OK;
}
fn halfnorm_kurt() -> f64, i64 {
    let d = PI - 2.0;
    return 8.0 * (PI - 3.0) / (d * d), OK;
}

// ---- hypsecant: the density 1 / (pi cosh x) on the whole line --------------

fn hypsecant_pdf(x: f64) -> f64, i64 {
    let c = fam_x(x);
    if c != OK { return nan(), c; }
    return 2.0 / (PI * (exp(x) + exp(-x))), OK;
}

fn hypsecant_cdf(x: f64) -> f64, i64 {
    let c = fam_x(x);
    if c != OK { return nan(), c; }
    return 2.0 / PI * atan(exp(x)), OK;
}

fn hypsecant_sf(x: f64) -> f64, i64 {
    let c = fam_x(x);
    if c != OK { return nan(), c; }
    return 2.0 / PI * atan(exp(-x)), OK;
}

fn hypsecant_ppf(p: f64) -> f64, i64 {
    let c = check_prob(p);
    if c != OK { return nan(), c; }
    let e, ev, ec = fam_end(p, -inf(), inf());
    if e { return ev, ec; }
    return log(tan(PI * p / 2.0)), OK;
}

fn hypsecant_mean() -> f64, i64 { return 0.0, OK; }
fn hypsecant_var() -> f64, i64 { return PI * PI / 4.0, OK; }
fn hypsecant_skew() -> f64, i64 { return 0.0, OK; }
fn hypsecant_kurt() -> f64, i64 { return 2.0, OK; }
