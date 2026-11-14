// Three continuous families of a growth or a logistic kind: genlogistic,
// genhalflogistic and gompertz. The names and the parameters follow
// SciPy.
//
// Every family has eight calls: pdf, cdf, sf, ppf, mean, var, skew and
// kurt. The kurtosis is the excess kurtosis. Every call returns the result
// first and a fault code last. A value outside the support gives a density
// of 0. The quantile of 0 or 1 gives the bound of the support.
//
// The moments of the generalized logistic law use the digamma function and
// the Hurwitz zeta function. The moments of the two other laws come from a
// Gauss-Legendre rule on a finite range.
namespace statistics;

import std;
import math;

// ---- genlogistic: the generalized logistic law --------------------------

fn genlogistic_pdf(x: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    // The form with the negative point loses no digits on either side.
    let a = abs(x);
    let m, cm = log1p(exp(-a));
    if cm != OK { return nan(), cm; }
    let mult = if x < 0.0 { c } else { 1.0 };
    return exp(log(c) - mult * a - (c + 1.0) * m), OK;
}

fn genlogistic_cdf(x: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    // The logarithm of one plus the exponential of the point stays finite
    // on both sides of the line.
    var m = 0.0;
    if x < 0.0 {
        let e, ce = log1p(exp(x));
        if ce != OK { return nan(), ce; }
        m = -x + e;
    } else {
        let e, ce = log1p(exp(-x));
        if ce != OK { return nan(), ce; }
        m = e;
    }
    return exp(-c * m), OK;
}

fn genlogistic_sf(x: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    var m = 0.0;
    if x < 0.0 {
        let e, ce = log1p(exp(x));
        if ce != OK { return nan(), ce; }
        m = -x + e;
    } else {
        let e, ce = log1p(exp(-x));
        if ce != OK { return nan(), ce; }
        m = e;
    }
    // The survival function keeps its digits in the far tail.
    return -expm1(-c * m), OK;
}

fn genlogistic_ppf(p: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(p, -inf(), inf());
    if e { return ev, ec; }
    // The quantile is the negative of the logarithm of the excess of the
    // point over the point raised to the power of the shape. The form
    // keeps its digits at both ends of the range.
    let t = -log(p) / c;
    let m, cm = expm1(t);
    if cm != OK { return nan(), cm; }
    return -log(m), OK;
}

// The Hurwitz zeta of the order s and the shift c. The order runs from 2
// to 4. The sum of 20 terms follows the rule of Euler and Maclaurin with
// six Bernoulli numbers. The relative error stays below 1e-15.
fn genlogistic_zeta(s: i64, c: f64) -> f64 {
    var total: f64 = 0.0;
    for k in 0..20 {
        total += pow(c + k as f64, -(s as f64));
    }
    let q = c + 20.0;
    total += pow(q, 1.0 - s as f64) / (s as f64 - 1.0);
    total += 0.5 * pow(q, -(s as f64));
    // The Bernoulli numbers B2, B4, B6, B8, B10 and B12.
    let bern: f64[6] = [1.0 / 6.0, -1.0 / 30.0, 1.0 / 42.0, -1.0 / 30.0, 5.0 / 66.0, -691.0 / 2730.0];
    for m in 1..6 + 1 {
        var r: f64 = 1.0;
        for i in 0..2 * m - 1 {
            r *= s as f64 + i as f64;
        }
        var f: f64 = 1.0;
        for i in 2..2 * m + 1 {
            f *= i as f64;
        }
        total += bern[m - 1] / f * r * pow(q, -(s as f64) - 2.0 * m as f64 + 1.0);
    }
    return total;
}

fn genlogistic_mean(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let p, cp = digamma(c);
    if cp != OK { return nan(), cp; }
    return EULER_GAMMA + p, OK;
}

fn genlogistic_var(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    return PI * PI / 6.0 + genlogistic_zeta(2, c), OK;
}

fn genlogistic_skew(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let v = PI * PI / 6.0 + genlogistic_zeta(2, c);
    return (-2.0 * genlogistic_zeta(3, c) + 2.0 * ZETA3) / (v * sqrt(v)), OK;
}

fn genlogistic_kurt(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let v = PI * PI / 6.0 + genlogistic_zeta(2, c);
    return (PI * PI * PI * PI / 15.0 + 6.0 * genlogistic_zeta(4, c)) / (v * v), OK;
}

// ---- genhalflogistic: the generalized half logistic law ------------------

fn genhalflogistic_pdf(x: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x < 0.0 || x > 1.0 / c { return 0.0, OK; }
    let t = pow(1.0 - c * x, 1.0 / c - 1.0);
    let t2 = t * (1.0 - c * x);
    return 2.0 * t / ((1.0 + t2) * (1.0 + t2)), OK;
}

fn genhalflogistic_cdf(x: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    if x >= 1.0 / c { return 1.0, OK; }
    let t2 = pow(1.0 - c * x, 1.0 / c);
    return (1.0 - t2) / (1.0 + t2), OK;
}

fn genhalflogistic_sf(x: f64, c: f64) -> f64, i64 {
    let v, k = genhalflogistic_cdf(x, c);
    if k != OK { return nan(), k; }
    return 1.0 - v, OK;
}

fn genhalflogistic_ppf(p: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(p, 0.0, 1.0 / c);
    if e { return ev, ec; }
    let r = pow((1.0 - p) / (1.0 + p), c);
    return (1.0 - r) / c, OK;
}

// The raw moment of order n of the law of the generalized half logistic.
// The law is the law of (1 - y^c) / c, where the law of y has the density
// 2 / (1 + y)^2 on 0 to 1. The change of variable u = (1 - y) / (1 + y)
// puts the integral on 0 to 1.
fn genhalflogistic_raw(n: i64, c: f64) -> f64 {
    var gx: f64[>..] = [];
    var gw: f64[>..] = [];
    gauss_legendre(gx, gw, FAM_NODES);
    var total: f64 = 0.0;
    let panels: i64 = 50;
    let width = 1.0 / panels as f64;
    for k in 0..panels {
        let left = width * k as f64;
        for j in 0..FAM_NODES {
            let u = left + width * (gx[j] + 1.0) / 2.0;
            let y = (1.0 - u) / (1.0 + u);
            // The change of variable gives the uniform law on the point u.
            total += pow((1.0 - pow(y, c)) / c, n as f64) * gw[j] * width / 2.0;
        }
    }
    return total;
}

fn genhalflogistic_mean(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    return genhalflogistic_raw(1, c), OK;
}

fn genhalflogistic_var(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let m1 = genhalflogistic_raw(1, c);
    let m2 = genhalflogistic_raw(2, c);
    return m2 - m1 * m1, OK;
}

fn genhalflogistic_skew(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let m1 = genhalflogistic_raw(1, c);
    let m2 = genhalflogistic_raw(2, c);
    let m3 = genhalflogistic_raw(3, c);
    let v = m2 - m1 * m1;
    return (m3 - 3.0 * m1 * m2 + 2.0 * m1 * m1 * m1) / (v * sqrt(v)), OK;
}

fn genhalflogistic_kurt(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let m1 = genhalflogistic_raw(1, c);
    let m2 = genhalflogistic_raw(2, c);
    let m3 = genhalflogistic_raw(3, c);
    let m4 = genhalflogistic_raw(4, c);
    let v = m2 - m1 * m1;
    let m1s = m1 * m1;
    return (m4 - 4.0 * m1 * m3 + 6.0 * m1s * m2 - 3.0 * m1s * m1s) / (v * v) - 3.0, OK;
}

// ---- gompertz: the law of Gompertz --------------------------------------

fn gompertz_pdf(x: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x < 0.0 { return 0.0, OK; }
    let m, cm = expm1(x);
    if cm != OK { return nan(), cm; }
    return c * exp(x - c * m), OK;
}

fn gompertz_cdf(x: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x < 0.0 { return 0.0, OK; }
    let m, cm = expm1(x);
    if cm != OK { return nan(), cm; }
    let t, ct = expm1(-c * m);
    if ct != OK { return nan(), ct; }
    return -t, OK;
}

fn gompertz_sf(x: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x < 0.0 { return 1.0, OK; }
    let m, cm = expm1(x);
    if cm != OK { return nan(), cm; }
    return exp(-c * m), OK;
}

fn gompertz_ppf(p: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(p, 0.0, inf());
    if e { return ev, ec; }
    let m, cm = log1p(-p);
    if cm != OK { return nan(), cm; }
    let t, ct = log1p(-m / c);
    if ct != OK { return nan(), ct; }
    return t, OK;
}

// The raw moment of order n of the law of Gompertz. The law is the law of
// log(1 + z / c), where the law of z is the exponential law of the rate 1.
// The change of variable y = z / c puts the integral on 0 to inf in the
// form c log(1 + y)^n exp(-c y).
fn gompertz_raw(n: i64, c: f64) -> f64 {
    var gx: f64[>..] = [];
    var gw: f64[>..] = [];
    gauss_legendre(gx, gw, FAM_NODES);
    var total: f64 = 0.0;
    let panels: i64 = 100;
    // The range 0 to 400 / c holds the whole mass of the exponential law.
    let hi = 400.0 / c;
    let width = hi / panels as f64;
    for k in 0..panels {
        let left = width * k as f64;
        for j in 0..FAM_NODES {
            let y = left + width * (gx[j] + 1.0) / 2.0;
            let l = log(1.0 + y);
            total += c * pow(l, n as f64) * exp(-c * y) * gw[j] * width / 2.0;
        }
    }
    return total;
}

fn gompertz_mean(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let v = gompertz_raw(1, c);
    if is_nan(v) || is_inf(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

fn gompertz_var(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let m1 = gompertz_raw(1, c);
    let m2 = gompertz_raw(2, c);
    let v = m2 - m1 * m1;
    if is_nan(v) || is_inf(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

fn gompertz_skew(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let m1 = gompertz_raw(1, c);
    let m2 = gompertz_raw(2, c);
    let m3 = gompertz_raw(3, c);
    let v = m2 - m1 * m1;
    let r = (m3 - 3.0 * m1 * m2 + 2.0 * m1 * m1 * m1) / (v * sqrt(v));
    if is_nan(r) || is_inf(r) { return nan(), ERR_OVERFLOW; }
    return r, OK;
}

fn gompertz_kurt(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let m1 = gompertz_raw(1, c);
    let m2 = gompertz_raw(2, c);
    let m3 = gompertz_raw(3, c);
    let m4 = gompertz_raw(4, c);
    let v = m2 - m1 * m1;
    let m1s = m1 * m1;
    let r = (m4 - 4.0 * m1 * m3 + 6.0 * m1s * m2 - 3.0 * m1s * m1s) / (v * v) - 3.0;
    if is_nan(r) || is_inf(r) { return nan(), ERR_OVERFLOW; }
    return r, OK;
}
