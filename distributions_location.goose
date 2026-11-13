// The Cauchy, Laplace, Rayleigh and Gumbel families of the package. Each
// family has a density, a cumulative function, a survival function, a
// quantile function and the four moments. All take a location and a scale.
//
// Every call returns the result first and a fault code last. The checks
// run in one order: NaN, infinity, domain, pole. A tail never comes from
// one minus the other tail. A quantile at 0 or 1 that has no finite value
// gives the domain fault. The Cauchy moments do not exist and give the
// undefined fault.
namespace statistics;

import std;
import math;

const ZETA3: f64 = 1.2020569031595942;      // the Apery constant
const EULER_GAMMA: f64 = 0.5772156649015329; // the Euler-Mascheroni constant

// The code of a point and of its location and scale.
fn loc_code(x: f64, loc: f64, scale: f64) -> i64 {
    if is_nan(x) || is_nan(loc) || is_nan(scale) { return ERR_NAN; }
    if is_inf(x) || is_inf(loc) || is_inf(scale) { return ERR_INFINITY; }
    if scale <= 0.0 { return ERR_DOMAIN; }
    return OK;
}

// The code of a probability, with the end points as valid input.
fn loc_code_p(p: f64, loc: f64, scale: f64) -> i64 {
    let c = loc_code(p, loc, scale);
    if c != OK { return c; }
    if p < 0.0 || p > 1.0 { return ERR_DOMAIN; }
    return OK;
}

// The standard score of x. The difference of two finite numbers can
// overflow, and then the score is infinite and the callers handle it.
fn loc_z(x: f64, loc: f64, scale: f64) -> f64 {
    return (x - loc) / scale;
}

// A result that is finite, or the overflow fault.
fn loc_result(v: f64) -> f64, i64 {
    if is_nan(v) { return nan(), ERR_UNDEFINED; }
    if is_inf(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

// One minus exp(-a) for a >= 0, exact for a small a and 1 for an a that
// overflows.
fn one_minus_exp(a: f64) -> f64 {
    if is_inf(a) { return 1.0; }
    let v, c = expm1(-a);
    return -v;
}

// ---- Cauchy ---------------------------------------------------------

fn cauchy_pdf(x: f64, loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    let c = loc_code(x, loc, scale);
    if c != OK { return nan(), c; }
    let z = loc_z(x, loc, scale);
    return loc_result(1.0 / (PI * scale * (1.0 + z * z)));
}

// The lower tail from the arc tangent of the inverse score far from 0.
fn cauchy_lower(z: f64) -> f64 {
    if z < -1.0 { return atan(-1.0 / z) / PI; }
    if z > 1.0 { return 1.0 - atan(1.0 / z) / PI; }
    return 0.5 + atan(z) / PI;
}

fn cauchy_cdf(x: f64, loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    let c = loc_code(x, loc, scale);
    if c != OK { return nan(), c; }
    return loc_result(cauchy_lower(loc_z(x, loc, scale)));
}

fn cauchy_sf(x: f64, loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    let c = loc_code(x, loc, scale);
    if c != OK { return nan(), c; }
    return loc_result(cauchy_lower(-loc_z(x, loc, scale)));
}

// The quantile is -cot(pi p). The form keeps the digits of a small p.
fn cauchy_ppf(p: f64, loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    let c = loc_code_p(p, loc, scale);
    if c != OK { return nan(), c; }
    if p == 0.0 || p == 1.0 { return nan(), ERR_DOMAIN; }
    var q: f64 = 0.0;
    if p < 0.5 { q = -1.0 / tan(PI * p); }
    else if p > 0.5 { q = 1.0 / tan(PI * (1.0 - p)); }
    let v = loc + scale * q;
    if is_inf(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

fn cauchy_mean(loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    let c = loc_code(0.0, loc, scale);
    if c != OK { return nan(), c; }
    return nan(), ERR_UNDEFINED;
}

fn cauchy_var(loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    return cauchy_mean(loc, scale);
}

fn cauchy_skew(loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    return cauchy_mean(loc, scale);
}

fn cauchy_kurt(loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    return cauchy_mean(loc, scale);
}

// ---- Laplace --------------------------------------------------------

fn laplace_pdf(x: f64, loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    let c = loc_code(x, loc, scale);
    if c != OK { return nan(), c; }
    return loc_result(exp(-abs(loc_z(x, loc, scale))) / (2.0 * scale));
}

fn laplace_lower(z: f64) -> f64 {
    if z < 0.0 { return 0.5 * exp(z); }
    return 1.0 - 0.5 * exp(-z);
}

fn laplace_cdf(x: f64, loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    let c = loc_code(x, loc, scale);
    if c != OK { return nan(), c; }
    return loc_result(laplace_lower(loc_z(x, loc, scale)));
}

fn laplace_sf(x: f64, loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    let c = loc_code(x, loc, scale);
    if c != OK { return nan(), c; }
    return loc_result(laplace_lower(-loc_z(x, loc, scale)));
}

fn laplace_ppf(p: f64, loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    let c = loc_code_p(p, loc, scale);
    if c != OK { return nan(), c; }
    if p == 0.0 || p == 1.0 { return nan(), ERR_DOMAIN; }
    var q: f64 = 0.0;
    if p < 0.5 { q = log(2.0 * p); }
    else if p > 0.5 { q = -log(2.0 * (1.0 - p)); }
    return loc_result(loc + scale * q);
}

fn laplace_mean(loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    let c = loc_code(0.0, loc, scale);
    if c != OK { return nan(), c; }
    return loc, OK;
}

fn laplace_var(loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    let c = loc_code(0.0, loc, scale);
    if c != OK { return nan(), c; }
    let v = 2.0 * scale * scale;
    if is_inf(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

fn laplace_skew(loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    let c = loc_code(0.0, loc, scale);
    if c != OK { return nan(), c; }
    return 0.0, OK;
}

fn laplace_kurt(loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    let c = loc_code(0.0, loc, scale);
    if c != OK { return nan(), c; }
    return 3.0, OK;
}

// ---- Rayleigh -------------------------------------------------------

fn rayleigh_pdf(x: f64, loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    let c = loc_code(x, loc, scale);
    if c != OK { return nan(), c; }
    let z = loc_z(x, loc, scale);
    if z <= 0.0 { return 0.0, OK; }
    if is_inf(z) { return 0.0, OK; }
    return loc_result(z * exp(-0.5 * z * z) / scale);
}

fn rayleigh_cdf(x: f64, loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    let c = loc_code(x, loc, scale);
    if c != OK { return nan(), c; }
    let z = loc_z(x, loc, scale);
    if z <= 0.0 { return 0.0, OK; }
    return one_minus_exp(0.5 * z * z), OK;
}

fn rayleigh_sf(x: f64, loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    let c = loc_code(x, loc, scale);
    if c != OK { return nan(), c; }
    let z = loc_z(x, loc, scale);
    if z <= 0.0 { return 1.0, OK; }
    return loc_result(exp(-0.5 * z * z));
}

fn rayleigh_ppf(p: f64, loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    let c = loc_code_p(p, loc, scale);
    if c != OK { return nan(), c; }
    if p == 1.0 { return nan(), ERR_DOMAIN; }
    if p == 0.0 { return loc, OK; }
    let l, c2 = log1p(-p);
    return loc_result(loc + scale * sqrt(-2.0 * l));
}

fn rayleigh_mean(loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    let c = loc_code(0.0, loc, scale);
    if c != OK { return nan(), c; }
    return loc_result(loc + scale * sqrt(0.5 * PI));
}

fn rayleigh_var(loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    let c = loc_code(0.0, loc, scale);
    if c != OK { return nan(), c; }
    let v = scale * scale * (4.0 - PI) / 2.0;
    if is_inf(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

fn rayleigh_skew(loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    let c = loc_code(0.0, loc, scale);
    if c != OK { return nan(), c; }
    let d = 4.0 - PI;
    return 2.0 * SQRT_PI * (PI - 3.0) / (d * sqrt(d)), OK;
}

fn rayleigh_kurt(loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    let c = loc_code(0.0, loc, scale);
    if c != OK { return nan(), c; }
    let d = 4.0 - PI;
    return -(6.0 * PI * PI - 24.0 * PI + 16.0) / (d * d), OK;
}

// ---- Gumbel, right-skewed -------------------------------------------

fn gumbel_r_pdf(x: f64, loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    let c = loc_code(x, loc, scale);
    if c != OK { return nan(), c; }
    let z = loc_z(x, loc, scale);
    // The score overflows when the point and the location have opposite
    // signs of a large size. The density then has no defined value.
    if is_inf(z) {
        if z < 0.0 { return nan(), ERR_UNDEFINED; }
        return 0.0, OK;
    }
    return loc_result(exp(-(z + exp(-z))) / scale);
}

fn gumbel_r_cdf(x: f64, loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    let c = loc_code(x, loc, scale);
    if c != OK { return nan(), c; }
    return loc_result(exp(-exp(-loc_z(x, loc, scale))));
}

fn gumbel_r_sf(x: f64, loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    let c = loc_code(x, loc, scale);
    if c != OK { return nan(), c; }
    return one_minus_exp(exp(-loc_z(x, loc, scale))), OK;
}

fn gumbel_r_ppf(p: f64, loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    let c = loc_code_p(p, loc, scale);
    if c != OK { return nan(), c; }
    if p == 0.0 || p == 1.0 { return nan(), ERR_DOMAIN; }
    var l: f64 = log(p);
    if p >= 0.5 { let m, c2 = log1p(p - 1.0); l = m; }
    return loc_result(loc - scale * log(-l));
}

fn gumbel_r_mean(loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    let c = loc_code(0.0, loc, scale);
    if c != OK { return nan(), c; }
    return loc_result(loc + scale * EULER_GAMMA);
}

fn gumbel_r_var(loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    let c = loc_code(0.0, loc, scale);
    if c != OK { return nan(), c; }
    let v = scale * scale * PI * PI / 6.0;
    if is_inf(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

fn gumbel_r_skew(loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    let c = loc_code(0.0, loc, scale);
    if c != OK { return nan(), c; }
    return 12.0 * sqrt(6.0) * ZETA3 / (PI * PI * PI), OK;
}

fn gumbel_r_kurt(loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    let c = loc_code(0.0, loc, scale);
    if c != OK { return nan(), c; }
    return 2.4, OK;
}

// ---- Gumbel, left-skewed --------------------------------------------

fn gumbel_l_pdf(x: f64, loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    let c = loc_code(x, loc, scale);
    if c != OK { return nan(), c; }
    let z = loc_z(x, loc, scale);
    // The score overflows when the point and the location have opposite
    // signs of a large size. The density then has no defined value.
    if is_inf(z) {
        if z > 0.0 { return nan(), ERR_UNDEFINED; }
        return 0.0, OK;
    }
    return loc_result(exp(z - exp(z)) / scale);
}

fn gumbel_l_cdf(x: f64, loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    let c = loc_code(x, loc, scale);
    if c != OK { return nan(), c; }
    return one_minus_exp(exp(loc_z(x, loc, scale))), OK;
}

fn gumbel_l_sf(x: f64, loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    let c = loc_code(x, loc, scale);
    if c != OK { return nan(), c; }
    return loc_result(exp(-exp(loc_z(x, loc, scale))));
}

fn gumbel_l_ppf(p: f64, loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    let c = loc_code_p(p, loc, scale);
    if c != OK { return nan(), c; }
    if p == 0.0 || p == 1.0 { return nan(), ERR_DOMAIN; }
    let l, c2 = log1p(-p);
    return loc_result(loc + scale * log(-l));
}

fn gumbel_l_mean(loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    let c = loc_code(0.0, loc, scale);
    if c != OK { return nan(), c; }
    return loc_result(loc - scale * EULER_GAMMA);
}

fn gumbel_l_var(loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    return gumbel_r_var(loc, scale);
}

fn gumbel_l_skew(loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    let v, c = gumbel_r_skew(loc, scale);
    return -v, c;
}

fn gumbel_l_kurt(loc: f64 = 0.0, scale: f64 = 1.0) -> f64, i64 {
    return gumbel_r_kurt(loc, scale);
}
