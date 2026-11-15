// Seven continuous families with a bounded or a half-line support:
// semicircular, triang, trapezoid, bradford, wald, maxwell and moyal.
// The names and the parameters follow SciPy.
//
// Every family has eight calls: pdf, cdf, sf, ppf, mean, var, skew and
// kurt. The kurtosis is the excess kurtosis. Every call returns the result
// first and a fault code last. A value outside the support gives a density
// of 0. The quantile of 0 or 1 gives the bound of the support. A bound
// that is infinite does not fit in a result, so the call gives
// ERR_OVERFLOW.
namespace statistics;

import std;
import math;


// The check of a shape that lies in [0, 1].
fn fam_unit_shape(x: f64) -> i64 {
    let c = check_scalar(x);
    if c != OK { return c; }
    if x < 0.0 || x > 1.0 { return ERR_DOMAIN; }
    return OK;
}

// ---- semicircular: the density 2 sqrt(1 - x^2) / pi on -1 to 1 -------------

fn semicircular_pdf(x: f64) -> f64, i64 {
    let c = fam_x(x);
    if c != OK { return nan(), c; }
    if x < -1.0 || x > 1.0 { return 0.0, OK; }
    return 2.0 / PI * sqrt(1.0 - x * x), OK;
}

fn semicircular_cdf(x: f64) -> f64, i64 {
    let c = fam_x(x);
    if c != OK { return nan(), c; }
    if x <= -1.0 { return 0.0, OK; }
    if x >= 1.0 { return 1.0, OK; }
    return 0.5 + (x * sqrt(1.0 - x * x) + asin(x)) / PI, OK;
}

fn semicircular_sf(x: f64) -> f64, i64 {
    let c = fam_x(x);
    if c != OK { return nan(), c; }
    if x <= -1.0 { return 1.0, OK; }
    if x >= 1.0 { return 0.0, OK; }
    return 0.5 - (x * sqrt(1.0 - x * x) + asin(x)) / PI, OK;
}

// The quantile is the quantile of the R law with a shape of 3.
fn semicircular_ppf(p: f64) -> f64, i64 {
    let c = check_prob(p);
    if c != OK { return nan(), c; }
    let e, ev, ec = fam_end(p, -1.0, 1.0);
    if e { return ev, ec; }
    let b, cb = betaincinv(1.5, 1.5, p);
    if cb != OK { return nan(), cb; }
    return 2.0 * b - 1.0, OK;
}

fn semicircular_mean() -> f64, i64 { return 0.0, OK; }
fn semicircular_var() -> f64, i64 { return 0.25, OK; }
fn semicircular_skew() -> f64, i64 { return 0.0, OK; }
fn semicircular_kurt() -> f64, i64 { return -1.0, OK; }

// ---- triang: the triangular law with the mode at c -----------------------

fn triang_pdf(x: f64, c: f64) -> f64, i64 {
    let k = fam_unit_shape(c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x < 0.0 || x > 1.0 { return 0.0, OK; }
    if c == 0.0 { return 2.0 - 2.0 * x, OK; }
    if c == 1.0 { return 2.0 * x, OK; }
    if x < c { return 2.0 * x / c, OK; }
    return 2.0 * (1.0 - x) / (1.0 - c), OK;
}

fn triang_cdf(x: f64, c: f64) -> f64, i64 {
    let k = fam_unit_shape(c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    if x >= 1.0 { return 1.0, OK; }
    if c == 0.0 { return 2.0 * x - x * x, OK; }
    if c == 1.0 { return x * x, OK; }
    if x < c { return x * x / c, OK; }
    return (x * x - 2.0 * x + c) / (c - 1.0), OK;
}

fn triang_sf(x: f64, c: f64) -> f64, i64 {
    let v, k = triang_cdf(x, c);
    if k != OK { return nan(), k; }
    return 1.0 - v, OK;
}

fn triang_ppf(p: f64, c: f64) -> f64, i64 {
    let k = fam_unit_shape(c);
    if k != OK { return nan(), k; }
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(p, 0.0, 1.0);
    if e { return ev, ec; }
    if p < c { return sqrt(c * p), OK; }
    return 1.0 - sqrt((1.0 - c) * (1.0 - p)), OK;
}

fn triang_mean(c: f64) -> f64, i64 {
    let k = fam_unit_shape(c);
    if k != OK { return nan(), k; }
    return (c + 1.0) / 3.0, OK;
}

fn triang_var(c: f64) -> f64, i64 {
    let k = fam_unit_shape(c);
    if k != OK { return nan(), k; }
    return (1.0 - c + c * c) / 18.0, OK;
}

fn triang_skew(c: f64) -> f64, i64 {
    let k = fam_unit_shape(c);
    if k != OK { return nan(), k; }
    let d = 1.0 - c + c * c;
    return sqrt(2.0) * (2.0 * c - 1.0) * (c + 1.0) * (c - 2.0) / (5.0 * d * sqrt(d)), OK;
}

fn triang_kurt(c: f64) -> f64, i64 {
    let k = fam_unit_shape(c);
    if k != OK { return nan(), k; }
    return -0.6, OK;
}

// ---- trapezoid: the trapezoid law with the plateau from c to d ------------

fn trapezoid_check(c: f64, d: f64) -> i64 {
    let k = fam_unit_shape(c);
    if k != OK { return k; }
    let k2 = fam_unit_shape(d);
    if k2 != OK { return k2; }
    if d < c { return ERR_DOMAIN; }
    return OK;
}

fn trapezoid_pdf(x: f64, c: f64, d: f64) -> f64, i64 {
    let k = trapezoid_check(c, d);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x < 0.0 || x > 1.0 { return 0.0, OK; }
    let u = 2.0 / (d - c + 1.0);
    if x < c { return u * x / c, OK; }
    if x <= d { return u, OK; }
    return u * (1.0 - x) / (1.0 - d), OK;
}

fn trapezoid_cdf(x: f64, c: f64, d: f64) -> f64, i64 {
    let k = trapezoid_check(c, d);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    if x >= 1.0 { return 1.0, OK; }
    if x < c { return x * x / c / (d - c + 1.0), OK; }
    if x <= d { return (c + 2.0 * (x - c)) / (d - c + 1.0), OK; }
    return 1.0 - (1.0 - x) * (1.0 - x) / (d - c + 1.0) / (1.0 - d), OK;
}

fn trapezoid_sf(x: f64, c: f64, d: f64) -> f64, i64 {
    let v, k = trapezoid_cdf(x, c, d);
    if k != OK { return nan(), k; }
    return 1.0 - v, OK;
}

fn trapezoid_ppf(p: f64, c: f64, d: f64) -> f64, i64 {
    let k = trapezoid_check(c, d);
    if k != OK { return nan(), k; }
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(p, 0.0, 1.0);
    if e { return ev, ec; }
    let qc, kc = trapezoid_cdf(c, c, d);
    if kc != OK { return nan(), kc; }
    let qd, kd = trapezoid_cdf(d, c, d);
    if kd != OK { return nan(), kd; }
    if p < qc { return sqrt(p * c * (1.0 + d - c)), OK; }
    if p <= qd { return 0.5 * p * (1.0 + d - c) + 0.5 * c, OK; }
    return 1.0 - sqrt((1.0 - p) * (d - c + 1.0) * (1.0 - d)), OK;
}

// The raw moment of the trapezoid law of order n. The density has three
// pieces, so the moment is the sum of three integrals.
fn trapezoid_raw(n: i64, c: f64, d: f64) -> f64 {
    let u = 2.0 / (d - c + 1.0);
    var t: f64 = 0.0;
    // The rising part, from 0 to c.
    if c > 0.0 { t = t + u * pow(c, n + 1) / (n + 2.0); }
    // The plateau, from c to d.
    t = t + u * (pow(d, n + 1) - pow(c, n + 1)) / (n + 1.0);
    // The falling part, from d to 1.
    if d < 1.0 {
        let hi = 1.0 / (n + 1.0) - 1.0 / (n + 2.0);
        let lo = pow(d, n + 1.0) / (n + 1.0) - pow(d, n + 2.0) / (n + 2.0);
        t = t + u * (hi - lo) / (1.0 - d);
    }
    return t;
}

fn trapezoid_mean(c: f64, d: f64) -> f64, i64 {
    let k = trapezoid_check(c, d);
    if k != OK { return nan(), k; }
    return trapezoid_raw(1, c, d), OK;
}

fn trapezoid_var(c: f64, d: f64) -> f64, i64 {
    let k = trapezoid_check(c, d);
    if k != OK { return nan(), k; }
    let m1 = trapezoid_raw(1, c, d);
    let m2 = trapezoid_raw(2, c, d);
    return m2 - m1 * m1, OK;
}

fn trapezoid_skew(c: f64, d: f64) -> f64, i64 {
    let k = trapezoid_check(c, d);
    if k != OK { return nan(), k; }
    let m1 = trapezoid_raw(1, c, d);
    let m2 = trapezoid_raw(2, c, d);
    let m3 = trapezoid_raw(3, c, d);
    let v = m2 - m1 * m1;
    return (m3 - 3.0 * m1 * m2 + 2.0 * m1 * m1 * m1) / (v * sqrt(v)), OK;
}

fn trapezoid_kurt(c: f64, d: f64) -> f64, i64 {
    let k = trapezoid_check(c, d);
    if k != OK { return nan(), k; }
    let m1 = trapezoid_raw(1, c, d);
    let m2 = trapezoid_raw(2, c, d);
    let m3 = trapezoid_raw(3, c, d);
    let m4 = trapezoid_raw(4, c, d);
    let v = m2 - m1 * m1;
    let m1s = m1 * m1;
    return (m4 - 4.0 * m1 * m3 + 6.0 * m1s * m2 - 3.0 * m1s * m1s) / (v * v) - 3.0, OK;
}

// ---- bradford: the density c / (1 + c x) / log(1 + c) on 0 to 1 -----------

fn bradford_pdf(x: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x < 0.0 || x > 1.0 { return 0.0, OK; }
    let l, cl = log1p(c);
    if cl != OK { return nan(), cl; }
    return c / (c * x + 1.0) / l, OK;
}

fn bradford_cdf(x: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    if x >= 1.0 { return 1.0, OK; }
    let l, cl = log1p(c);
    if cl != OK { return nan(), cl; }
    let m, cm = log1p(c * x);
    if cm != OK { return nan(), cm; }
    return m / l, OK;
}

fn bradford_sf(x: f64, c: f64) -> f64, i64 {
    let v, k = bradford_cdf(x, c);
    if k != OK { return nan(), k; }
    return 1.0 - v, OK;
}

fn bradford_ppf(p: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(p, 0.0, 1.0);
    if e { return ev, ec; }
    let l, cl = log1p(c);
    if cl != OK { return nan(), cl; }
    let m, cm = expm1(p * l);
    if cm != OK { return nan(), cm; }
    return m / c, OK;
}

fn bradford_mean(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let l, cl = log1p(c);
    if cl != OK { return nan(), cl; }
    return (c - l) / (c * l), OK;
}

fn bradford_var(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let l, cl = log1p(c);
    if cl != OK { return nan(), cl; }
    return ((c + 2.0) * l - 2.0 * c) / (2.0 * c * l * l), OK;
}

fn bradford_skew(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let l, cl = log1p(c);
    if cl != OK { return nan(), cl; }
    let num = sqrt(2.0) * (12.0 * c * c - 9.0 * c * l * (c + 2.0) + 2.0 * l * l * (c * (c + 3.0) + 3.0));
    let den = sqrt(c * (c * (l - 2.0) + 2.0 * l)) * (3.0 * c * (l - 2.0) + 6.0 * l);
    return num / den, OK;
}

fn bradford_kurt(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let l, cl = log1p(c);
    if cl != OK { return nan(), cl; }
    let num = c * c * c * (l - 3.0) * (l * (3.0 * l - 16.0) + 24.0)
        + 12.0 * l * c * c * (l - 4.0) * (l - 3.0)
        + 6.0 * c * l * l * (3.0 * l - 14.0) + 12.0 * l * l * l;
    let den = 3.0 * c * (c * (l - 2.0) + 2.0 * l) * (c * (l - 2.0) + 2.0 * l);
    return num / den, OK;
}

// ---- wald: the inverse Gaussian law with a mean of 1 ----------------------

fn wald_pdf(x: f64) -> f64, i64 {
    let c = fam_x(x);
    if c != OK { return nan(), c; }
    if x <= 0.0 { return 0.0, OK; }
    let d = x - 1.0;
    return exp(-0.5 * d * d / x) / (SQRT_2PI * x * sqrt(x)), OK;
}

fn wald_cdf(x: f64) -> f64, i64 {
    let c = fam_x(x);
    if c != OK { return nan(), c; }
    if x <= 0.0 { return 0.0, OK; }
    // The cumulative function of the inverse Gaussian law with a mean of 1
    // is the sum of two normal functions.
    let r = sqrt(1.0 / x);
    let a = norm_cdf(r * (x - 1.0));
    let b = norm_cdf(-r * (x + 1.0));
    let e = exp(2.0);
    return a + e * b, OK;
}

fn wald_sf(x: f64) -> f64, i64 {
    let v, k = wald_cdf(x);
    if k != OK { return nan(), k; }
    return 1.0 - v, OK;
}

fn wald_ppf(p: f64) -> f64, i64 {
    let c = check_prob(p);
    if c != OK { return nan(), c; }
    let e, ev, ec = fam_end(p, 0.0, inf());
    if e { return ev, ec; }
    // The quantile has no closed form, so the shared search finds it.
    return fam_search_ppf(p, 0.0, 1.0) { v => wald_cdf_value(v) };
}

// The cumulative function of the law, without the code. The search of
// the quantile takes this value.
fn wald_cdf_value(x: f64) -> f64 {
    let v, c = wald_cdf(x);
    if c != OK { return nan(); }
    return v;
}

fn wald_mean() -> f64, i64 { return 1.0, OK; }
fn wald_var() -> f64, i64 { return 1.0, OK; }
fn wald_skew() -> f64, i64 { return 3.0, OK; }
fn wald_kurt() -> f64, i64 { return 15.0, OK; }

// ---- maxwell: the density sqrt(2 / pi) x^2 exp(-x^2 / 2) on 0 to inf -------

fn maxwell_pdf(x: f64) -> f64, i64 {
    let c = fam_x(x);
    if c != OK { return nan(), c; }
    if x < 0.0 { return 0.0, OK; }
    return sqrt(2.0 / PI) * x * x * exp(-0.5 * x * x), OK;
}

fn maxwell_cdf(x: f64) -> f64, i64 {
    let c = fam_x(x);
    if c != OK { return nan(), c; }
    if x <= 0.0 { return 0.0, OK; }
    return gammainc(1.5, 0.5 * x * x);
}

fn maxwell_sf(x: f64) -> f64, i64 {
    let c = fam_x(x);
    if c != OK { return nan(), c; }
    if x <= 0.0 { return 1.0, OK; }
    return gammaincc(1.5, 0.5 * x * x);
}

fn maxwell_ppf(p: f64) -> f64, i64 {
    let c = check_prob(p);
    if c != OK { return nan(), c; }
    let e, ev, ec = fam_end(p, 0.0, inf());
    if e { return ev, ec; }
    let v, k = gammaincinv(1.5, p);
    if k != OK { return nan(), k; }
    return sqrt(2.0 * v), OK;
}

fn maxwell_mean() -> f64, i64 { return 2.0 * sqrt(2.0 / PI), OK; }
fn maxwell_var() -> f64, i64 { return 3.0 - 8.0 / PI, OK; }
fn maxwell_skew() -> f64, i64 {
    let v = 3.0 * PI - 8.0;
    return sqrt(2.0) * (32.0 - 10.0 * PI) / (v * sqrt(v)), OK;
}
fn maxwell_kurt() -> f64, i64 {
    let v = 3.0 * PI - 8.0;
    return (-12.0 * PI * PI + 160.0 * PI - 384.0) / (v * v), OK;
}

// ---- moyal: the density exp(-(x + exp(-x)) / 2) / sqrt(2 pi) -------------

fn moyal_pdf(x: f64) -> f64, i64 {
    let c = fam_x(x);
    if c != OK { return nan(), c; }
    let e = exp(-x);
    return exp(-0.5 * (x + e)) / SQRT_2PI, OK;
}

fn moyal_cdf(x: f64) -> f64, i64 {
    let c = fam_x(x);
    if c != OK { return nan(), c; }
    return erfc(exp(-0.5 * x) * INV_SQRT2);
}

fn moyal_sf(x: f64) -> f64, i64 {
    let c = fam_x(x);
    if c != OK { return nan(), c; }
    return erf(exp(-0.5 * x) * INV_SQRT2);
}

fn moyal_ppf(p: f64) -> f64, i64 {
    let c = check_prob(p);
    if c != OK { return nan(), c; }
    let e, ev, ec = fam_end(p, -inf(), inf());
    if e { return ev, ec; }
    let v, k = erfcinv(p);
    if k != OK { return nan(), k; }
    return -log(2.0 * v * v), OK;
}

fn moyal_mean() -> f64, i64 { return log(2.0) + EULER_GAMMA, OK; }
fn moyal_var() -> f64, i64 { return PI * PI / 2.0, OK; }
fn moyal_skew() -> f64, i64 {
    return 28.0 * sqrt(2.0) * ZETA3 / (PI * PI * PI), OK;
}
fn moyal_kurt() -> f64, i64 { return 4.0, OK; }
