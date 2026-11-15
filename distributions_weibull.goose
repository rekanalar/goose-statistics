// Four continuous families of the Weibull kind: dweibull, invweibull,
// weibull_max and truncweibull_min. The names and the parameters follow
// SciPy.
//
// Every family has eight calls: pdf, cdf, sf, ppf, mean, var, skew and
// kurt. The kurtosis is the excess kurtosis. Every call returns the result
// first and a fault code last. A value outside the support gives a density
// of 0. The quantile of 0 or 1 gives the bound of the support.
//
// The moments of these families have a closed form. The Weibull law of the
// whole line is symmetric, so its mean and its skewness are 0.
namespace statistics;

import std;
import math;

// ---- dweibull: the Weibull law of the whole line ------------------------

fn dweibull_pdf(x: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    let a = abs(x);
    if a == 0.0 { return 0.0, OK; }
    // The log form keeps a huge point from the difference of an infinity
    // and a zero.
    let lp = log(c / 2.0) + (c - 1.0) * log(a) - pow(a, c);
    return exp(lp), OK;
}

fn dweibull_cdf(x: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    let h = 0.5 * exp(-pow(abs(x), c));
    if x > 0.0 { return 1.0 - h, OK; }
    return h, OK;
}

fn dweibull_sf(x: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    let h = 0.5 * exp(-pow(abs(x), c));
    if x > 0.0 { return h, OK; }
    return 1.0 - h, OK;
}

fn dweibull_ppf(p: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(p, -inf(), inf());
    if e { return ev, ec; }
    let f = if p <= 0.5 { 2.0 * p } else { 2.0 * (1.0 - p) };
    let r = pow(-log(f), 1.0 / c);
    if p > 0.5 { return r, OK; }
    return -r, OK;
}

fn dweibull_mean(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    return 0.0, OK;
}

fn dweibull_var(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let g, cg = gamma(1.0 + 2.0 / c);
    if cg != OK { return nan(), cg; }
    return g, OK;
}

fn dweibull_skew(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    return 0.0, OK;
}

fn dweibull_kurt(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let g2, cg = gamma(1.0 + 2.0 / c);
    if cg != OK { return nan(), cg; }
    let g4, cg4 = gamma(1.0 + 4.0 / c);
    if cg4 != OK { return nan(), cg4; }
    return g4 / (g2 * g2) - 3.0, OK;
}

// ---- invweibull: the inverse Weibull law on 0 to inf --------------------

fn invweibull_pdf(x: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    let lp = log(c) + (-c - 1.0) * log(x) - pow(x, -c);
    return exp(lp), OK;
}

fn invweibull_cdf(x: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    return exp(-pow(x, -c)), OK;
}

fn invweibull_sf(x: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 1.0, OK; }
    let m, cm = expm1(-pow(x, -c));
    if cm != OK { return nan(), cm; }
    return -m, OK;
}

fn invweibull_ppf(p: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(p, 0.0, inf());
    if e { return ev, ec; }
    return pow(-log(p), -1.0 / c), OK;
}

// The moment of order n of the inverse Weibull law is gamma(1 - n / c).
fn invweibull_mean(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let g, cg = gamma(1.0 - 1.0 / c);
    if cg != OK { return nan(), cg; }
    return g, OK;
}

fn invweibull_var(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let g1, cg = gamma(1.0 - 1.0 / c);
    if cg != OK { return nan(), cg; }
    let g2, cg2 = gamma(1.0 - 2.0 / c);
    if cg2 != OK { return nan(), cg2; }
    return g2 - g1 * g1, OK;
}

fn invweibull_skew(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let g1, cg = gamma(1.0 - 1.0 / c);
    if cg != OK { return nan(), cg; }
    let g2, cg2 = gamma(1.0 - 2.0 / c);
    if cg2 != OK { return nan(), cg2; }
    let g3, cg3 = gamma(1.0 - 3.0 / c);
    if cg3 != OK { return nan(), cg3; }
    let v = g2 - g1 * g1;
    return (g3 - 3.0 * g1 * g2 + 2.0 * g1 * g1 * g1) / (v * sqrt(v)), OK;
}

fn invweibull_kurt(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let g1, cg = gamma(1.0 - 1.0 / c);
    if cg != OK { return nan(), cg; }
    let g2, cg2 = gamma(1.0 - 2.0 / c);
    if cg2 != OK { return nan(), cg2; }
    let g3, cg3 = gamma(1.0 - 3.0 / c);
    if cg3 != OK { return nan(), cg3; }
    let g4, cg4 = gamma(1.0 - 4.0 / c);
    if cg4 != OK { return nan(), cg4; }
    let v = g2 - g1 * g1;
    let m1s = g1 * g1;
    return (g4 - 4.0 * g1 * g3 + 6.0 * m1s * g2 - 3.0 * m1s * m1s) / (v * v) - 3.0, OK;
}

// ---- weibull_max: the Weibull law of the negative line ------------------

fn weibull_max_pdf(x: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x > 0.0 { return 0.0, OK; }
    let a = -x;
    let lp = log(c) + (c - 1.0) * log(a) - pow(a, c);
    return exp(lp), OK;
}

fn weibull_max_cdf(x: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x >= 0.0 { return 1.0, OK; }
    return exp(-pow(-x, c)), OK;
}

fn weibull_max_sf(x: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x >= 0.0 { return 0.0, OK; }
    let m, cm = expm1(-pow(-x, c));
    if cm != OK { return nan(), cm; }
    return -m, OK;
}

fn weibull_max_ppf(p: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(p, -inf(), 0.0);
    if e { return ev, ec; }
    return -pow(-log(p), 1.0 / c), OK;
}

// The raw moments of the Weibull law of the negative line. The moment of
// order n has the sign of n.
fn weibull_max_raw(n: i64, c: f64) -> f64, i64 {
    let g, cg = gamma(1.0 + n / c);
    if cg != OK { return nan(), cg; }
    if n % 2 == 1 { return -g, OK; }
    return g, OK;
}

fn weibull_max_mean(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    return weibull_max_raw(1, c);
}

fn weibull_max_var(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let m1, c1 = weibull_max_raw(1, c);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = weibull_max_raw(2, c);
    if c2 != OK { return nan(), c2; }
    return m2 - m1 * m1, OK;
}

fn weibull_max_skew(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let m1, c1 = weibull_max_raw(1, c);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = weibull_max_raw(2, c);
    if c2 != OK { return nan(), c2; }
    let m3, c3 = weibull_max_raw(3, c);
    if c3 != OK { return nan(), c3; }
    let v = m2 - m1 * m1;
    return (m3 - 3.0 * m1 * m2 + 2.0 * m1 * m1 * m1) / (v * sqrt(v)), OK;
}

fn weibull_max_kurt(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let m1, c1 = weibull_max_raw(1, c);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = weibull_max_raw(2, c);
    if c2 != OK { return nan(), c2; }
    let m3, c3 = weibull_max_raw(3, c);
    if c3 != OK { return nan(), c3; }
    let m4, c4 = weibull_max_raw(4, c);
    if c4 != OK { return nan(), c4; }
    let v = m2 - m1 * m1;
    let m1s = m1 * m1;
    return (m4 - 4.0 * m1 * m3 + 6.0 * m1s * m2 - 3.0 * m1s * m1s) / (v * v) - 3.0, OK;
}

// ---- truncweibull_min: the Weibull law on a and b -----------------------

fn truncweibull_check(c: f64, a: f64, b: f64) -> i64 {
    let k = check_shape(c);
    if k != OK { return k; }
    let ka = check_scalar(a);
    if ka != OK { return ka; }
    if a < 0.0 { return ERR_DOMAIN; }
    let kb = check_shape(b);
    if kb != OK { return kb; }
    if b <= a { return ERR_DOMAIN; }
    return OK;
}

// The mass of the law between a and b. The difference of two exponentials
// loses its digits when both are small, so the mass comes from the
// exponential of the lower bound times a difference that keeps its
// digits. A mass below the precision of a double number gives the fit
// fault.
fn truncweibull_mass(c: f64, a: f64, b: f64) -> f64, i64 {
    let ea = exp(-pow(a, c));
    let eb = exp(-pow(b, c));
    let d = ea - eb;
    if d <= 0.0 { return nan(), ERR_FIT; }
    return d, OK;
}

fn truncweibull_min_pdf(x: f64, c: f64, a: f64, b: f64) -> f64, i64 {
    let k = truncweibull_check(c, a, b);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x < a || x > b { return 0.0, OK; }
    let d, cd = truncweibull_mass(c, a, b);
    if cd != OK { return nan(), cd; }
    return c * pow(x, c - 1.0) * exp(-pow(x, c)) / d, OK;
}

fn truncweibull_min_cdf(x: f64, c: f64, a: f64, b: f64) -> f64, i64 {
    let k = truncweibull_check(c, a, b);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= a { return 0.0, OK; }
    if x >= b { return 1.0, OK; }
    let d, cd = truncweibull_mass(c, a, b);
    if cd != OK { return nan(), cd; }
    let e = exp(-pow(x, c));
    let ea = exp(-pow(a, c));
    return (ea - e) / d, OK;
}

fn truncweibull_min_sf(x: f64, c: f64, a: f64, b: f64) -> f64, i64 {
    let k = truncweibull_check(c, a, b);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= a { return 1.0, OK; }
    if x >= b { return 0.0, OK; }
    let d, cd = truncweibull_mass(c, a, b);
    if cd != OK { return nan(), cd; }
    let e = exp(-pow(x, c));
    let eb = exp(-pow(b, c));
    return (e - eb) / d, OK;
}

fn truncweibull_min_ppf(p: f64, c: f64, a: f64, b: f64) -> f64, i64 {
    let k = truncweibull_check(c, a, b);
    if k != OK { return nan(), k; }
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(p, a, b);
    if e { return ev, ec; }
    let d, cd = truncweibull_mass(c, a, b);
    if cd != OK { return nan(), cd; }
    let ea = exp(-pow(a, c));
    let eb = exp(-pow(b, c));
    // The mass below x is (exp(-a^c) - exp(-x^c)) / d, so the quantile
    // solves exp(-x^c) = (1 - p) exp(-a^c) + p exp(-b^c).
    let target = (1.0 - p) * ea + p * eb;
    return pow(-log(target), 1.0 / c), OK;
}

// The raw moment of order n of the truncated Weibull law. The change of
// variable u = x^c gives the incomplete gamma.
fn truncweibull_raw(n: i64, c: f64, a: f64, b: f64) -> f64, i64 {
    let d, cd = truncweibull_mass(c, a, b);
    if cd != OK { return nan(), cd; }
    let s = n / c + 1.0;
    let ga, cg = gammaincc(s, pow(a, c));
    if cg != OK { return nan(), cg; }
    let gb, cgb = gammaincc(s, pow(b, c));
    if cgb != OK { return nan(), cgb; }
    let g, cgm = gamma(s);
    if cgm != OK { return nan(), cgm; }
    return (ga - gb) * g / d, OK;
}

fn truncweibull_min_mean(c: f64, a: f64, b: f64) -> f64, i64 {
    let k = truncweibull_check(c, a, b);
    if k != OK { return nan(), k; }
    return truncweibull_raw(1, c, a, b);
}

fn truncweibull_min_var(c: f64, a: f64, b: f64) -> f64, i64 {
    let k = truncweibull_check(c, a, b);
    if k != OK { return nan(), k; }
    let m1, c1 = truncweibull_raw(1, c, a, b);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = truncweibull_raw(2, c, a, b);
    if c2 != OK { return nan(), c2; }
    return m2 - m1 * m1, OK;
}

fn truncweibull_min_skew(c: f64, a: f64, b: f64) -> f64, i64 {
    let k = truncweibull_check(c, a, b);
    if k != OK { return nan(), k; }
    let m1, c1 = truncweibull_raw(1, c, a, b);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = truncweibull_raw(2, c, a, b);
    if c2 != OK { return nan(), c2; }
    let m3, c3 = truncweibull_raw(3, c, a, b);
    if c3 != OK { return nan(), c3; }
    let v = m2 - m1 * m1;
    return (m3 - 3.0 * m1 * m2 + 2.0 * m1 * m1 * m1) / (v * sqrt(v)), OK;
}

fn truncweibull_min_kurt(c: f64, a: f64, b: f64) -> f64, i64 {
    let k = truncweibull_check(c, a, b);
    if k != OK { return nan(), k; }
    let m1, c1 = truncweibull_raw(1, c, a, b);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = truncweibull_raw(2, c, a, b);
    if c2 != OK { return nan(), c2; }
    let m3, c3 = truncweibull_raw(3, c, a, b);
    if c3 != OK { return nan(), c3; }
    let m4, c4 = truncweibull_raw(4, c, a, b);
    if c4 != OK { return nan(), c4; }
    let v = m2 - m1 * m1;
    let m1s = m1 * m1;
    return (m4 - 4.0 * m1 * m3 + 6.0 * m1s * m2 - 3.0 * m1s * m1s) / (v * v) - 3.0, OK;
}
