// Four continuous families of a truncation or a special shape:
// truncexpon, truncpareto, truncnorm and tukeylambda. The names and the
// parameters follow SciPy.
//
// Every family has eight calls: pdf, cdf, sf, ppf, mean, var, skew and
// kurt. The kurtosis is the excess kurtosis. Every call returns the result
// first and a fault code last. A value outside the support gives a density
// of 0. The quantile of 0 or 1 gives the bound of the support.
//
// The support of the truncated Pareto law is bounded, so every moment of
// that law exists.
namespace statistics;

import std;
import math;

// ---- truncexpon: the truncated exponential law --------------------------

fn truncexpon_pdf(x: f64, b: f64) -> f64, i64 {
    let k = check_shape(b);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x < 0.0 || x > b { return 0.0, OK; }
    let d, cd = expm1(-b);
    if cd != OK { return nan(), cd; }
    return exp(-x) / (-d), OK;
}

fn truncexpon_cdf(x: f64, b: f64) -> f64, i64 {
    let k = check_shape(b);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    if x >= b { return 1.0, OK; }
    let d, cd = expm1(-b);
    if cd != OK { return nan(), cd; }
    let e, ce = expm1(-x);
    if ce != OK { return nan(), ce; }
    return e / d, OK;
}

fn truncexpon_sf(x: f64, b: f64) -> f64, i64 {
    let v, k = truncexpon_cdf(x, b);
    if k != OK { return nan(), k; }
    return 1.0 - v, OK;
}

fn truncexpon_ppf(p: f64, b: f64) -> f64, i64 {
    let k = check_shape(b);
    if k != OK { return nan(), k; }
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(p, 0.0, b);
    if e { return ev, ec; }
    let d, cd = expm1(-b);
    if cd != OK { return nan(), cd; }
    let t, ct = log1p(p * d);
    if ct != OK { return nan(), ct; }
    return -t, OK;
}

// The raw moment of order n of the truncated exponential law. It is the
// integral of the point to the n over the exponential law on 0 to b. The
// integral is the factorial of n times one minus the exponential of minus b
// times the sum of the powers of b over the factorials.
fn truncexpon_raw(n: i64, b: f64) -> f64, i64 {
    let d, cd = expm1(-b);
    if cd != OK { return nan(), cd; }
    // The sum of the powers of b over the factorials, from 1 to n.
    var t: f64 = 0.0;
    var term: f64 = 1.0;
    for k in 1..n + 1 {
        term *= b / k;
        t += term;
    }
    // The factorial of n.
    var fact: f64 = 1.0;
    for k in 2..n + 1 { fact *= k; }
    let e = exp(-b);
    let num = -d - e * t;
    return fact * num / -d, OK;
}

fn truncexpon_mean(b: f64) -> f64, i64 {
    let k = check_shape(b);
    if k != OK { return nan(), k; }
    let v, c = truncexpon_raw(1, b);
    if c != OK { return nan(), c; }
    return v, OK;
}

fn truncexpon_var(b: f64) -> f64, i64 {
    let k = check_shape(b);
    if k != OK { return nan(), k; }
    let m1, c1 = truncexpon_raw(1, b);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = truncexpon_raw(2, b);
    if c2 != OK { return nan(), c2; }
    return m2 - m1 * m1, OK;
}

fn truncexpon_skew(b: f64) -> f64, i64 {
    let k = check_shape(b);
    if k != OK { return nan(), k; }
    let m1, c1 = truncexpon_raw(1, b);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = truncexpon_raw(2, b);
    if c2 != OK { return nan(), c2; }
    let m3, c3 = truncexpon_raw(3, b);
    if c3 != OK { return nan(), c3; }
    return fam_skew_raw(m1, m2, m3), OK;
}

fn truncexpon_kurt(b: f64) -> f64, i64 {
    let k = check_shape(b);
    if k != OK { return nan(), k; }
    let m1, c1 = truncexpon_raw(1, b);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = truncexpon_raw(2, b);
    if c2 != OK { return nan(), c2; }
    let m3, c3 = truncexpon_raw(3, b);
    if c3 != OK { return nan(), c3; }
    let m4, c4 = truncexpon_raw(4, b);
    if c4 != OK { return nan(), c4; }
    return fam_kurt_raw(m1, m2, m3, m4), OK;
}

// ---- truncpareto: the truncated Pareto law ------------------------------

fn truncpareto_pdf(x: f64, b: f64, c: f64) -> f64, i64 {
    let kb = truncpareto_check(b, c);
    if kb != OK { return nan(), kb; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x < 1.0 || x > c { return 0.0, OK; }
    let d, cd = truncpareto_den(b, c);
    if cd != OK { return nan(), cd; }
    return b * pow(x, -b - 1.0) / d, OK;
}

fn truncpareto_cdf(x: f64, b: f64, c: f64) -> f64, i64 {
    let kb = truncpareto_check(b, c);
    if kb != OK { return nan(), kb; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 1.0 { return 0.0, OK; }
    if x >= c { return 1.0, OK; }
    let d, cd = truncpareto_den(b, c);
    if cd != OK { return nan(), cd; }
    return (1.0 - pow(x, -b)) / d, OK;
}

fn truncpareto_sf(x: f64, b: f64, c: f64) -> f64, i64 {
    let v, k = truncpareto_cdf(x, b, c);
    if k != OK { return nan(), k; }
    return 1.0 - v, OK;
}

fn truncpareto_ppf(p: f64, b: f64, c: f64) -> f64, i64 {
    let kb = truncpareto_check(b, c);
    if kb != OK { return nan(), kb; }
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(p, 1.0, c);
    if e { return ev, ec; }
    let d, cd = truncpareto_den(b, c);
    if cd != OK { return nan(), cd; }
    return pow(1.0 - d * p, -1.0 / b), OK;
}

// The check of the shapes of the truncated Pareto law. The shape b may be
// negative but not 0. The shape c must be above 1.
fn truncpareto_check(b: f64, c: f64) -> i64 {
    let kb = check_scalar(b);
    if kb != OK { return kb; }
    let kc = check_scalar(c);
    if kc != OK { return kc; }
    if b == 0.0 { return ERR_DOMAIN; }
    if c <= 1.0 { return ERR_DOMAIN; }
    return OK;
}

// The denominator of the law: 1 - 1 / c to the power b. It stays in the
// form of the logarithm, so a shape of a large size keeps its digits.
fn truncpareto_den(b: f64, c: f64) -> f64, i64 {
    let t = b * log(c);
    let e, ce = expm1(-t);
    if ce != OK { return nan(), ce; }
    return -e, OK;
}

// The raw moment of order n of the truncated Pareto law. It is the shape b
// over the shape b minus n, times the end of the support to the n, times
// the power of the end of the support minus the shape b minus one, over the
// end of the support to the shape b minus one. The support of the law is
// bounded, so every moment exists.
fn truncpareto_raw(n: f64, b: f64, c: f64) -> f64, i64 {
    let lc = log(c);
    // The denominator: the end of the support to the shape b minus one.
    let d, cd = expm1(b * lc);
    if cd != OK { return nan(), cd; }
    if b == n {
        // The power of the end of the support times its logarithm.
        return b * exp(b * lc) * lc / d, OK;
    }
    let e, ce = expm1((b - n) * lc);
    if ce != OK { return nan(), ce; }
    return b / (b - n) * pow(c, n) * e / d, OK;
}

fn truncpareto_mean(b: f64, c: f64) -> f64, i64 {
    let kb = truncpareto_check(b, c);
    if kb != OK { return nan(), kb; }
    return truncpareto_raw(1.0, b, c);
}

fn truncpareto_var(b: f64, c: f64) -> f64, i64 {
    let kb = truncpareto_check(b, c);
    if kb != OK { return nan(), kb; }
    let m1, c1 = truncpareto_raw(1.0, b, c);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = truncpareto_raw(2.0, b, c);
    if c2 != OK { return nan(), c2; }
    return m2 - m1 * m1, OK;
}

fn truncpareto_skew(b: f64, c: f64) -> f64, i64 {
    let kb = truncpareto_check(b, c);
    if kb != OK { return nan(), kb; }
    let m1, c1 = truncpareto_raw(1.0, b, c);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = truncpareto_raw(2.0, b, c);
    if c2 != OK { return nan(), c2; }
    let m3, c3 = truncpareto_raw(3.0, b, c);
    if c3 != OK { return nan(), c3; }
    return fam_skew_raw(m1, m2, m3), OK;
}

fn truncpareto_kurt(b: f64, c: f64) -> f64, i64 {
    let kb = truncpareto_check(b, c);
    if kb != OK { return nan(), kb; }
    let m1, c1 = truncpareto_raw(1.0, b, c);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = truncpareto_raw(2.0, b, c);
    if c2 != OK { return nan(), c2; }
    let m3, c3 = truncpareto_raw(3.0, b, c);
    if c3 != OK { return nan(), c3; }
    let m4, c4 = truncpareto_raw(4.0, b, c);
    if c4 != OK { return nan(), c4; }
    return fam_kurt_raw(m1, m2, m3, m4), OK;
}

// ---- truncnorm: the truncated normal law --------------------------------

// The check of the bounds of the truncated normal law. The lower bound must
// be below the upper bound.
fn truncnorm_check(a: f64, b: f64) -> i64 {
    let ka = check_scalar(a);
    if ka != OK { return ka; }
    let kb = check_scalar(b);
    if kb != OK { return kb; }
    if a >= b { return ERR_DOMAIN; }
    return OK;
}

// The mass of the law between the two bounds. The difference of two
// cumulative values loses its digits for bounds far apart, so the mass
// comes from the tail of the upper bound. The bounds far in one tail
// give a mass below the precision of a double number. That mass gives
// the fit fault.
fn truncnorm_mass(a: f64, b: f64) -> f64, i64 {
    let lo, clo = norm_cdf(a);
    if clo != OK { return nan(), clo; }
    let hi, chi = norm_cdf(b);
    if chi != OK { return nan(), chi; }
    var m: f64 = 0.0;
    if lo < 0.5 {
        m = hi - lo;
    } else {
        // The mass is the complement of the lower tail. The tail of the
        // upper bound is below one half there.
        let sl, csl = norm_sf(a);
        if csl != OK { return nan(), csl; }
        let sh, csh = norm_sf(b);
        if csh != OK { return nan(), csh; }
        m = sl - sh;
    }
    if m <= 0.0 { return nan(), ERR_FIT; }
    return m, OK;
}

fn truncnorm_pdf(x: f64, a: f64, b: f64) -> f64, i64 {
    let k = truncnorm_check(a, b);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x < a || x > b { return 0.0, OK; }
    let p, cp = norm_pdf(x);
    if cp != OK { return nan(), cp; }
    let m, cm = truncnorm_mass(a, b);
    if cm != OK { return nan(), cm; }
    return p / m, OK;
}

fn truncnorm_cdf(x: f64, a: f64, b: f64) -> f64, i64 {
    let k = truncnorm_check(a, b);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= a { return 0.0, OK; }
    if x >= b { return 1.0, OK; }
    let p, cp = norm_cdf(x);
    if cp != OK { return nan(), cp; }
    let lo, clo = norm_cdf(a);
    if clo != OK { return nan(), clo; }
    let m, cm = truncnorm_mass(a, b);
    if cm != OK { return nan(), cm; }
    return (p - lo) / m, OK;
}

fn truncnorm_sf(x: f64, a: f64, b: f64) -> f64, i64 {
    let v, k = truncnorm_cdf(x, a, b);
    if k != OK { return nan(), k; }
    return 1.0 - v, OK;
}

fn truncnorm_ppf(p: f64, a: f64, b: f64) -> f64, i64 {
    let k = truncnorm_check(a, b);
    if k != OK { return nan(), k; }
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(p, a, b);
    if e { return ev, ec; }
    let m, cm = truncnorm_mass(a, b);
    if cm != OK { return nan(), cm; }
    // The quantile is the normal quantile of the mass of the lower bound
    // plus the share of the probability.
    let lo, clo = norm_cdf(a);
    if clo != OK { return nan(), clo; }
    let q = lo + p * m;
    if q >= 1.0 {
        // The point is in the far tail. The survival function keeps its
        // digits there.
        let s, cs = norm_sf(a);
        if cs != OK { return nan(), cs; }
        let t = s - p * m;
        if t <= 0.0 { return b, OK; }
        let v, cv = erfcinv(2.0 * t);
        if cv != OK { return nan(), cv; }
        return -SQRT2 * v, OK;
    }
    return norm_std_ppf(q), OK;
}

// The edge of the truncated normal law: the point to the power k times the
// density of the normal law there. An infinite bound gives 0.
fn truncnorm_edge(x: f64, k: i64) -> f64, i64 {
    if is_inf(x) { return 0.0, OK; }
    let p, cp = norm_pdf(x);
    if cp != OK { return nan(), cp; }
    if k == 0 { return p, OK; }
    return p * pow(x, k), OK;
}

// The integral of the point to the n times the density of the normal law
// between the two bounds. The recursion of the normal law gives it. The
// integral of order 0 is the mass of the law.
recursive fn truncnorm_integral(n: i64, a: f64, b: f64) -> f64, i64 {
    if n == 0 {
        let m, cm = truncnorm_mass(a, b);
        if cm != OK { return nan(), cm; }
        return m, OK;
    }
    let ea, cea = truncnorm_edge(a, n - 1);
    if cea != OK { return nan(), cea; }
    let eb, ceb = truncnorm_edge(b, n - 1);
    if ceb != OK { return nan(), ceb; }
    if n == 1 { return ea - eb, OK; }
    let prev, cp = truncnorm_integral(n - 2, a, b);
    if cp != OK { return nan(), cp; }
    return ea - eb + (n - 1.0) * prev, OK;
}

// The raw moment of order n of the truncated normal law.
fn truncnorm_raw(n: i64, a: f64, b: f64) -> f64, i64 {
    let m, cm = truncnorm_mass(a, b);
    if cm != OK { return nan(), cm; }
    let i, ci = truncnorm_integral(n, a, b);
    if ci != OK { return nan(), ci; }
    return i / m, OK;
}

fn truncnorm_mean(a: f64, b: f64) -> f64, i64 {
    let k = truncnorm_check(a, b);
    if k != OK { return nan(), k; }
    return truncnorm_raw(1, a, b);
}

fn truncnorm_var(a: f64, b: f64) -> f64, i64 {
    let k = truncnorm_check(a, b);
    if k != OK { return nan(), k; }
    let m1, c1 = truncnorm_raw(1, a, b);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = truncnorm_raw(2, a, b);
    if c2 != OK { return nan(), c2; }
    return m2 - m1 * m1, OK;
}

fn truncnorm_skew(a: f64, b: f64) -> f64, i64 {
    let k = truncnorm_check(a, b);
    if k != OK { return nan(), k; }
    let m1, c1 = truncnorm_raw(1, a, b);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = truncnorm_raw(2, a, b);
    if c2 != OK { return nan(), c2; }
    let m3, c3 = truncnorm_raw(3, a, b);
    if c3 != OK { return nan(), c3; }
    return fam_skew_raw(m1, m2, m3), OK;
}

fn truncnorm_kurt(a: f64, b: f64) -> f64, i64 {
    let k = truncnorm_check(a, b);
    if k != OK { return nan(), k; }
    let m1, c1 = truncnorm_raw(1, a, b);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = truncnorm_raw(2, a, b);
    if c2 != OK { return nan(), c2; }
    let m3, c3 = truncnorm_raw(3, a, b);
    if c3 != OK { return nan(), c3; }
    let m4, c4 = truncnorm_raw(4, a, b);
    if c4 != OK { return nan(), c4; }
    return fam_kurt_raw(m1, m2, m3, m4), OK;
}
