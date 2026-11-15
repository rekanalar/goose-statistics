// Eight continuous families of a power: fisk, burr, burr12, betaprime,
// genpareto, genlogistic, genhalflogistic and gompertz. The names and
// the parameters follow SciPy.
//
// Every family has eight calls: pdf, cdf, sf, ppf, mean, var, skew and
// kurt. The kurtosis is the excess kurtosis. Every call returns the result
// first and a fault code last. A value outside the support gives a density
// of 0. The quantile of 0 or 1 gives the bound of the support.
//
// The moments of these families have a closed form. The law of the Burr
// family has no mean below a shape of 1. It has no variance below a shape
// of 2. The law of the generalized Pareto family has no mean below a shape
// of 1. The calls give the undefined fault there.
namespace statistics;

import std;
import math;

// ---- fisk: the log-logistic law -----------------------------------------

fn fisk_pdf(x: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    // The form with the logarithm keeps the digits of a large shape. The
    // direct power overflows there.
    let lx = log(x);
    // The log of one plus the power of x. The helper keeps the digits on
    // each side.
    let lp, clp = fam_log1pexp(-c * lx);
    if clp != OK { return nan(), clp; }
    let l = log(c) - c * lx - lx - 2.0 * lp;
    let r = exp(l);
    return fam_finite(r);
}

fn fisk_cdf(x: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    return 1.0 / (1.0 + pow(x, -c)), OK;
}

fn fisk_sf(x: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 1.0, OK; }
    let p = pow(x, -c);
    return p / (1.0 + p), OK;
}

fn fisk_ppf(p: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(p, 0.0, inf());
    if e { return ev, ec; }
    let r = pow(p / (1.0 - p), 1.0 / c);
    return fam_finite(r);
}

// The raw moment of order n of the log-logistic law. It is the beta
// function of d + n / c and 1 - n / c.
fn fisk_raw(n: i64, c: f64) -> f64, i64 {
    let nc = n / c;
    let b, cb = beta(1.0 + nc, 1.0 - nc);
    if cb != OK { return nan(), cb; }
    return b, OK;
}

fn fisk_mean(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    if c <= 1.0 { return nan(), ERR_UNDEFINED; }
    return fisk_raw(1, c);
}

fn fisk_var(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    if c <= 2.0 { return nan(), ERR_UNDEFINED; }
    let m1, c1 = fisk_raw(1, c);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = fisk_raw(2, c);
    if c2 != OK { return nan(), c2; }
    return m2 - m1 * m1, OK;
}

fn fisk_skew(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    if c <= 3.0 { return nan(), ERR_UNDEFINED; }
    let m1, c1 = fisk_raw(1, c);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = fisk_raw(2, c);
    if c2 != OK { return nan(), c2; }
    let m3, c3 = fisk_raw(3, c);
    if c3 != OK { return nan(), c3; }
    let v = m2 - m1 * m1;
    return (m3 - 3.0 * m1 * m2 + 2.0 * m1 * m1 * m1) / (v * sqrt(v)), OK;
}

fn fisk_kurt(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    if c <= 4.0 { return nan(), ERR_UNDEFINED; }
    let m1, c1 = fisk_raw(1, c);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = fisk_raw(2, c);
    if c2 != OK { return nan(), c2; }
    let m3, c3 = fisk_raw(3, c);
    if c3 != OK { return nan(), c3; }
    let m4, c4 = fisk_raw(4, c);
    if c4 != OK { return nan(), c4; }
    let v = m2 - m1 * m1;
    let m1s = m1 * m1;
    return (m4 - 4.0 * m1 * m3 + 6.0 * m1s * m2 - 3.0 * m1s * m1s) / (v * v) - 3.0, OK;
}

// ---- burr: the law of Burr of the third kind ----------------------------

fn burr_pdf(x: f64, c: f64, d: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let kd = check_shape(d);
    if kd != OK { return nan(), kd; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    // The form with the logarithm keeps the digits of a large shape. The
    // direct power overflows there.
    let lx = log(x);
    let lp, clp = fam_log1pexp(-c * lx);
    if clp != OK { return nan(), clp; }
    let l = log(c) + log(d) - c * lx - lx - (d + 1.0) * lp;
    let r = exp(l);
    return fam_finite(r);
}

fn burr_cdf(x: f64, c: f64, d: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let kd = check_shape(d);
    if kd != OK { return nan(), kd; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    return pow(1.0 + pow(x, -c), -d), OK;
}

fn burr_sf(x: f64, c: f64, d: f64) -> f64, i64 {
    let v, k = burr_cdf(x, c, d);
    if k != OK { return nan(), k; }
    return 1.0 - v, OK;
}

fn burr_ppf(p: f64, c: f64, d: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let kd = check_shape(d);
    if kd != OK { return nan(), kd; }
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(p, 0.0, inf());
    if e { return ev, ec; }
    return pow(pow(p, -1.0 / d) - 1.0, -1.0 / c), OK;
}

// The raw moment of order n of the law of Burr. It is the beta function
// of d + n / c and 1 - n / c.
fn burr_raw(n: i64, c: f64, d: f64) -> f64, i64 {
    let nc = n / c;
    let b, cb = beta(d + nc, 1.0 - nc);
    if cb != OK { return nan(), cb; }
    return d * b, OK;
}

fn burr_mean(c: f64, d: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let kd = check_shape(d);
    if kd != OK { return nan(), kd; }
    if c <= 1.0 { return nan(), ERR_UNDEFINED; }
    return burr_raw(1, c, d);
}

fn burr_var(c: f64, d: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let kd = check_shape(d);
    if kd != OK { return nan(), kd; }
    if c <= 2.0 { return nan(), ERR_UNDEFINED; }
    let m1, c1 = burr_raw(1, c, d);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = burr_raw(2, c, d);
    if c2 != OK { return nan(), c2; }
    return m2 - m1 * m1, OK;
}

fn burr_skew(c: f64, d: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let kd = check_shape(d);
    if kd != OK { return nan(), kd; }
    if c <= 3.0 { return nan(), ERR_UNDEFINED; }
    let m1, c1 = burr_raw(1, c, d);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = burr_raw(2, c, d);
    if c2 != OK { return nan(), c2; }
    let m3, c3 = burr_raw(3, c, d);
    if c3 != OK { return nan(), c3; }
    let v = m2 - m1 * m1;
    return (m3 - 3.0 * m1 * m2 + 2.0 * m1 * m1 * m1) / (v * sqrt(v)), OK;
}

fn burr_kurt(c: f64, d: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let kd = check_shape(d);
    if kd != OK { return nan(), kd; }
    if c <= 4.0 { return nan(), ERR_UNDEFINED; }
    let m1, c1 = burr_raw(1, c, d);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = burr_raw(2, c, d);
    if c2 != OK { return nan(), c2; }
    let m3, c3 = burr_raw(3, c, d);
    if c3 != OK { return nan(), c3; }
    let m4, c4 = burr_raw(4, c, d);
    if c4 != OK { return nan(), c4; }
    let v = m2 - m1 * m1;
    let m1s = m1 * m1;
    return (m4 - 4.0 * m1 * m3 + 6.0 * m1s * m2 - 3.0 * m1s * m1s) / (v * v) - 3.0, OK;
}

// ---- burr12: the law of Burr of the twelfth kind ------------------------

fn burr12_pdf(x: f64, c: f64, d: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let kd = check_shape(d);
    if kd != OK { return nan(), kd; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    return c * d * pow(x, c - 1.0) / pow(1.0 + pow(x, c), d + 1.0), OK;
}

fn burr12_cdf(x: f64, c: f64, d: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let kd = check_shape(d);
    if kd != OK { return nan(), kd; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    return 1.0 - pow(1.0 + pow(x, c), -d), OK;
}

fn burr12_sf(x: f64, c: f64, d: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let kd = check_shape(d);
    if kd != OK { return nan(), kd; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 1.0, OK; }
    return pow(1.0 + pow(x, c), -d), OK;
}

fn burr12_ppf(p: f64, c: f64, d: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let kd = check_shape(d);
    if kd != OK { return nan(), kd; }
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(p, 0.0, inf());
    if e { return ev, ec; }
    let m, cm = log1p(-p);
    if cm != OK { return nan(), cm; }
    let t = -m / d;
    return pow(expm1(t), 1.0 / c), OK;
}

// The raw moment of order n of the law of Burr of the twelfth kind. It
// is the beta function of 1 + n / c and d - n / c.
fn burr12_raw(n: i64, c: f64, d: f64) -> f64, i64 {
    let nc = n / c;
    let b, cb = beta(1.0 + nc, d - nc);
    if cb != OK { return nan(), cb; }
    return d * b, OK;
}

fn burr12_mean(c: f64, d: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let kd = check_shape(d);
    if kd != OK { return nan(), kd; }
    if c * d <= 1.0 { return nan(), ERR_UNDEFINED; }
    return burr12_raw(1, c, d);
}

fn burr12_var(c: f64, d: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let kd = check_shape(d);
    if kd != OK { return nan(), kd; }
    if c * d <= 2.0 { return nan(), ERR_UNDEFINED; }
    let m1, c1 = burr12_raw(1, c, d);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = burr12_raw(2, c, d);
    if c2 != OK { return nan(), c2; }
    return m2 - m1 * m1, OK;
}

fn burr12_skew(c: f64, d: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let kd = check_shape(d);
    if kd != OK { return nan(), kd; }
    if c * d <= 3.0 { return nan(), ERR_UNDEFINED; }
    let m1, c1 = burr12_raw(1, c, d);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = burr12_raw(2, c, d);
    if c2 != OK { return nan(), c2; }
    let m3, c3 = burr12_raw(3, c, d);
    if c3 != OK { return nan(), c3; }
    let v = m2 - m1 * m1;
    return (m3 - 3.0 * m1 * m2 + 2.0 * m1 * m1 * m1) / (v * sqrt(v)), OK;
}

fn burr12_kurt(c: f64, d: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let kd = check_shape(d);
    if kd != OK { return nan(), kd; }
    if c * d <= 4.0 { return nan(), ERR_UNDEFINED; }
    let m1, c1 = burr12_raw(1, c, d);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = burr12_raw(2, c, d);
    if c2 != OK { return nan(), c2; }
    let m3, c3 = burr12_raw(3, c, d);
    if c3 != OK { return nan(), c3; }
    let m4, c4 = burr12_raw(4, c, d);
    if c4 != OK { return nan(), c4; }
    let v = m2 - m1 * m1;
    let m1s = m1 * m1;
    return (m4 - 4.0 * m1 * m3 + 6.0 * m1s * m2 - 3.0 * m1s * m1s) / (v * v) - 3.0, OK;
}

// ---- betaprime: the beta prime law --------------------------------------

fn betaprime_pdf(x: f64, a: f64, b: f64) -> f64, i64 {
    let k = check_shape(a);
    if k != OK { return nan(), k; }
    let kb = check_shape(b);
    if kb != OK { return nan(), kb; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    let bt, cb = beta(a, b);
    if cb != OK { return nan(), cb; }
    // The form with the logarithm keeps the digits of a large point.
    let lx = log(x);
    let l = (a - 1.0) * lx - (a + b) * log1p(x) - log(bt);
    let r = exp(l);
    return fam_finite(r);
}

fn betaprime_cdf(x: f64, a: f64, b: f64) -> f64, i64 {
    let k = check_shape(a);
    if k != OK { return nan(), k; }
    let kb = check_shape(b);
    if kb != OK { return nan(), kb; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    // The law of the beta prime is the law of the beta with the point
    // x / (1 + x). The tail of the law of the beta keeps its digits.
    if x > 1.0 {
        let v, cv = betainc(b, a, 1.0 / (1.0 + x));
        if cv != OK { return nan(), cv; }
        return 1.0 - v, OK;
    }
    let v, cv = betainc(a, b, x / (1.0 + x));
    if cv != OK { return nan(), cv; }
    return v, OK;
}

fn betaprime_sf(x: f64, a: f64, b: f64) -> f64, i64 {
    let v, k = betaprime_cdf(x, a, b);
    if k != OK { return nan(), k; }
    return 1.0 - v, OK;
}

fn betaprime_ppf(p: f64, a: f64, b: f64) -> f64, i64 {
    let k = check_shape(a);
    if k != OK { return nan(), k; }
    let kb = check_shape(b);
    if kb != OK { return nan(), kb; }
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(p, 0.0, inf());
    if e { return ev, ec; }
    // The quantile of the beta prime is the quantile of the beta, and
    // the point x / (1 + x) changes it back.
    if p > 0.5 {
        // The quantile of the beta law with the shapes swapped is the
        // complement of the quantile of the beta law.
        let v, cv = betaincinv(b, a, 1.0 - p);
        if cv != OK { return nan(), cv; }
        return (1.0 - v) / v, OK;
    }
    let v, cv = betaincinv(a, b, p);
    if cv != OK { return nan(), cv; }
    return v / (1.0 - v), OK;
}

// The raw moment of order n of the beta prime law. It is the beta
// function of a + n and b - n over the beta function of a and b. The
// moment exists when b is above n.
fn betaprime_raw(n: i64, a: f64, b: f64) -> f64, i64 {
    let top, ct = beta(a + n, b - n);
    if ct != OK { return nan(), ct; }
    let bottom, cb = beta(a, b);
    if cb != OK { return nan(), cb; }
    return top / bottom, OK;
}

fn betaprime_mean(a: f64, b: f64) -> f64, i64 {
    let k = check_shape(a);
    if k != OK { return nan(), k; }
    let kb = check_shape(b);
    if kb != OK { return nan(), kb; }
    if b <= 1.0 { return nan(), ERR_UNDEFINED; }
    return betaprime_raw(1, a, b);
}

fn betaprime_var(a: f64, b: f64) -> f64, i64 {
    let k = check_shape(a);
    if k != OK { return nan(), k; }
    let kb = check_shape(b);
    if kb != OK { return nan(), kb; }
    if b <= 2.0 { return nan(), ERR_UNDEFINED; }
    let m1, c1 = betaprime_raw(1, a, b);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = betaprime_raw(2, a, b);
    if c2 != OK { return nan(), c2; }
    return m2 - m1 * m1, OK;
}

fn betaprime_skew(a: f64, b: f64) -> f64, i64 {
    let k = check_shape(a);
    if k != OK { return nan(), k; }
    let kb = check_shape(b);
    if kb != OK { return nan(), kb; }
    if b <= 3.0 { return nan(), ERR_UNDEFINED; }
    let m1, c1 = betaprime_raw(1, a, b);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = betaprime_raw(2, a, b);
    if c2 != OK { return nan(), c2; }
    let m3, c3 = betaprime_raw(3, a, b);
    if c3 != OK { return nan(), c3; }
    let v = m2 - m1 * m1;
    return (m3 - 3.0 * m1 * m2 + 2.0 * m1 * m1 * m1) / (v * sqrt(v)), OK;
}

fn betaprime_kurt(a: f64, b: f64) -> f64, i64 {
    let k = check_shape(a);
    if k != OK { return nan(), k; }
    let kb = check_shape(b);
    if kb != OK { return nan(), kb; }
    if b <= 4.0 { return nan(), ERR_UNDEFINED; }
    let m1, c1 = betaprime_raw(1, a, b);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = betaprime_raw(2, a, b);
    if c2 != OK { return nan(), c2; }
    let m3, c3 = betaprime_raw(3, a, b);
    if c3 != OK { return nan(), c3; }
    let m4, c4 = betaprime_raw(4, a, b);
    if c4 != OK { return nan(), c4; }
    let v = m2 - m1 * m1;
    let m1s = m1 * m1;
    return (m4 - 4.0 * m1 * m3 + 6.0 * m1s * m2 - 3.0 * m1s * m1s) / (v * v) - 3.0, OK;
}

// ---- genpareto: the generalized Pareto law ------------------------------

fn genpareto_pdf(x: f64, c: f64) -> f64, i64 {
    let k = fam_x(c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x < 0.0 { return 0.0, OK; }
    if c < 0.0 && x > -1.0 / c { return 0.0, OK; }
    if c == 0.0 { return exp(-x), OK; }
    // The form with the logarithm keeps the digits of a small shape. The
    // direct power loses the term c*x when that term is below the
    // precision of the sum.
    let l, cl = log1p(c * x);
    if cl != OK { return nan(), cl; }
    return exp((-1.0 - 1.0 / c) * l), OK;
}

fn genpareto_cdf(x: f64, c: f64) -> f64, i64 {
    let k = fam_x(c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x < 0.0 { return 0.0, OK; }
    if c < 0.0 && x > -1.0 / c { return 1.0, OK; }
    if c == 0.0 { return -expm1(-x), OK; }
    let t = -log1p(c * x) / c;
    return -expm1(t), OK;
}

fn genpareto_sf(x: f64, c: f64) -> f64, i64 {
    let k = fam_x(c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x < 0.0 { return 1.0, OK; }
    if c < 0.0 && x > -1.0 / c { return 0.0, OK; }
    if c == 0.0 { return exp(-x), OK; }
    return exp(-log1p(c * x) / c), OK;
}

fn genpareto_ppf(p: f64, c: f64) -> f64, i64 {
    let k = fam_x(c);
    if k != OK { return nan(), k; }
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    // A negative shape bounds the support above.
    var hi = inf();
    if c < 0.0 { hi = -1.0 / c; }
    let e, ev, ec = fam_end(p, 0.0, hi);
    if e { return ev, ec; }
    if c == 0.0 { return -log1p(-p), OK; }
    return (pow(1.0 - p, -c) - 1.0) / c, OK;
}

fn genpareto_mean(c: f64) -> f64, i64 {
    let k = fam_x(c);
    if k != OK { return nan(), k; }
    if c >= 1.0 { return nan(), ERR_UNDEFINED; }
    return 1.0 / (1.0 - c), OK;
}

fn genpareto_var(c: f64) -> f64, i64 {
    let k = fam_x(c);
    if k != OK { return nan(), k; }
    if c >= 0.5 { return nan(), ERR_UNDEFINED; }
    return 1.0 / (1.0 - c) / (1.0 - c) / (1.0 - 2.0 * c), OK;
}

fn genpareto_skew(c: f64) -> f64, i64 {
    let k = fam_x(c);
    if k != OK { return nan(), k; }
    if c >= 1.0 / 3.0 { return nan(), ERR_UNDEFINED; }
    return 2.0 * (1.0 + c) * sqrt(1.0 - 2.0 * c) / (1.0 - 3.0 * c), OK;
}

fn genpareto_kurt(c: f64) -> f64, i64 {
    let k = fam_x(c);
    if k != OK { return nan(), k; }
    if c >= 0.25 { return nan(), ERR_UNDEFINED; }
    return 3.0 * (1.0 - 2.0 * c) * (2.0 * c * c + c + 3.0) / (1.0 - 3.0 * c) / (1.0 - 4.0 * c) - 3.0, OK;
}
