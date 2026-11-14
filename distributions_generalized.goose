// Three continuous families of the gamma and the beta kind: gengamma,
// genexpon and mielke. The names and the parameters follow SciPy.
//
// Every family has eight calls: pdf, cdf, sf, ppf, mean, var, skew and
// kurt. The kurtosis is the excess kurtosis. Every call returns the result
// first and a fault code last. A value outside the support gives a density
// of 0. The quantile of 0 or 1 gives the bound of the support.
//
// The generalized gamma law takes a shape below 0. Its moments then exist
// only up to a limit that depends on the shape.
namespace statistics;

import std;
import math;

// ---- gengamma: the generalized gamma law --------------------------------

fn gengamma_pdf(x: f64, a: f64, c: f64) -> f64, i64 {
    let k = gengamma_check(a, c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x < 0.0 { return 0.0, OK; }
    if x == 0.0 {
        if c < 0.0 { return 0.0, OK; }
        if c * a - 1.0 > 0.0 { return 0.0, OK; }
        if c * a - 1.0 < 0.0 { return inf(), ERR_OVERFLOW; }
        let g0, cg0 = gamma(a);
        if cg0 != OK { return nan(), cg0; }
        return abs(c) / g0, OK;
    }
    let g, cg = gamma(a);
    if cg != OK { return nan(), cg; }
    let lx = log(x);
    let p = exp((c * a - 1.0) * lx - pow(x, c)) * abs(c) / g;
    if is_inf(p) { return nan(), ERR_OVERFLOW; }
    return p, OK;
}

fn gengamma_cdf(x: f64, a: f64, c: f64) -> f64, i64 {
    let k = gengamma_check(a, c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    let xc = pow(x, c);
    if c > 0.0 {
        let v, cv = gammainc(a, xc);
        if cv != OK { return nan(), cv; }
        return v, OK;
    }
    let v, cv = gammaincc(a, xc);
    if cv != OK { return nan(), cv; }
    return v, OK;
}

fn gengamma_sf(x: f64, a: f64, c: f64) -> f64, i64 {
    let k = gengamma_check(a, c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 1.0, OK; }
    let xc = pow(x, c);
    if c > 0.0 {
        let v, cv = gammaincc(a, xc);
        if cv != OK { return nan(), cv; }
        return v, OK;
    }
    let v, cv = gammainc(a, xc);
    if cv != OK { return nan(), cv; }
    return v, OK;
}

fn gengamma_ppf(q: f64, a: f64, c: f64) -> f64, i64 {
    let k = gengamma_check(a, c);
    if k != OK { return nan(), k; }
    let pc = check_prob(q);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(q, 0.0, inf());
    if e { return ev, ec; }
    if c > 0.0 {
        let v, cv = gammaincinv(a, q);
        if cv != OK { return nan(), cv; }
        return pow(v, 1.0 / c), OK;
    }
    let v, cv = gammainccinv(a, q);
    if cv != OK { return nan(), cv; }
    return pow(v, 1.0 / c), OK;
}

// The check of the shapes of the generalized gamma law. The shape a must be
// above 0 and the shape c must not be 0.
fn gengamma_check(a: f64, c: f64) -> i64 {
    let ka = check_shape(a);
    if ka != OK { return ka; }
    let kc = check_scalar(c);
    if kc != OK { return kc; }
    if c == 0.0 { return ERR_DOMAIN; }
    return OK;
}

// The raw moment of order n of the generalized gamma law. It is the ratio
// of the gamma function of the shape plus n over the shape c, and of the
// shape. The moment exists when that sum is above 0.
fn gengamma_raw(n: i64, a: f64, c: f64) -> f64, i64 {
    let t = a + n as f64 / c;
    if t <= 0.0 { return nan(), ERR_UNDEFINED; }
    let l1, s1, c1 = lgamma(t);
    if c1 != OK { return nan(), c1; }
    let l2, s2, c2 = lgamma(a);
    if c2 != OK { return nan(), c2; }
    let v = exp(l1 - l2);
    if is_inf(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

fn gengamma_mean(a: f64, c: f64) -> f64, i64 {
    let k = gengamma_check(a, c);
    if k != OK { return nan(), k; }
    return gengamma_raw(1, a, c);
}

fn gengamma_var(a: f64, c: f64) -> f64, i64 {
    let k = gengamma_check(a, c);
    if k != OK { return nan(), k; }
    let m1, c1 = gengamma_raw(1, a, c);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = gengamma_raw(2, a, c);
    if c2 != OK { return nan(), c2; }
    return m2 - m1 * m1, OK;
}

fn gengamma_skew(a: f64, c: f64) -> f64, i64 {
    let k = gengamma_check(a, c);
    if k != OK { return nan(), k; }
    let m1, c1 = gengamma_raw(1, a, c);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = gengamma_raw(2, a, c);
    if c2 != OK { return nan(), c2; }
    let m3, c3 = gengamma_raw(3, a, c);
    if c3 != OK { return nan(), c3; }
    return fam_skew_raw(m1, m2, m3), OK;
}

fn gengamma_kurt(a: f64, c: f64) -> f64, i64 {
    let k = gengamma_check(a, c);
    if k != OK { return nan(), k; }
    let m1, c1 = gengamma_raw(1, a, c);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = gengamma_raw(2, a, c);
    if c2 != OK { return nan(), c2; }
    let m3, c3 = gengamma_raw(3, a, c);
    if c3 != OK { return nan(), c3; }
    let m4, c4 = gengamma_raw(4, a, c);
    if c4 != OK { return nan(), c4; }
    return fam_kurt_raw(m1, m2, m3, m4), OK;
}

// ---- genexpon: the generalized exponential law ---------------------------

// The density of the generalized exponential law at the point x. The law is
// a sum of two exponentials. The form below keeps its digits for a small
// point.
fn genexpon_density(x: f64, a: f64, b: f64, c: f64) -> f64, i64 {
    let e, ce = expm1(-c * x);
    if ce != OK { return nan(), ce; }
    // One minus the exponential of minus c times x.
    let u = -e;
    let t = (-a - b) * x + b * u / c;
    return (a + b * u) * exp(t), OK;
}

fn genexpon_pdf(x: f64, a: f64, b: f64, c: f64) -> f64, i64 {
    let k = genexpon_check(a, b, c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x < 0.0 { return 0.0, OK; }
    let f, cf = genexpon_density(x, a, b, c);
    if cf != OK { return nan(), cf; }
    if is_inf(f) { return nan(), ERR_OVERFLOW; }
    return f, OK;
}

fn genexpon_cdf(x: f64, a: f64, b: f64, c: f64) -> f64, i64 {
    let k = genexpon_check(a, b, c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    let e, ce = expm1(-c * x);
    if ce != OK { return nan(), ce; }
    let u = -e;
    let t = (-a - b) * x + b * u / c;
    let v, cv = expm1(t);
    if cv != OK { return nan(), cv; }
    return -v, OK;
}

fn genexpon_sf(x: f64, a: f64, b: f64, c: f64) -> f64, i64 {
    let k = genexpon_check(a, b, c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 1.0, OK; }
    let e, ce = expm1(-c * x);
    if ce != OK { return nan(), ce; }
    let u = -e;
    let t = (-a - b) * x + b * u / c;
    return exp(t), OK;
}

fn genexpon_ppf(p: f64, a: f64, b: f64, c: f64) -> f64, i64 {
    let k = genexpon_check(a, b, c);
    if k != OK { return nan(), k; }
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(p, 0.0, inf());
    if e { return ev, ec; }
    // The quantile needs the Lambert W function. The package holds no such
    // function, so the search finds the root of the cumulative function.
    var lo = 0.0;
    var hi = 1.0;
    // The upper bound grows until the cumulative function passes the
    // probability.
    var count: i64 = 0;
    while count < 200 {
        let v, cv = genexpon_cdf(hi, a, b, c);
        if cv != OK { return nan(), cv; }
        if v >= p { break; }
        hi = hi * 2.0;
        count += 1;
    }
    if count >= 200 { return nan(), ERR_FIT; }
    // The bisection of the root of the cumulative function.
    var mid = (lo + hi) / 2.0;
    count = 0;
    while count < 200 {
        let v, cv = genexpon_cdf(mid, a, b, c);
        if cv != OK { return nan(), cv; }
        if v == p { return mid, OK; }
        if v < p { lo = mid; } else { hi = mid; }
        mid = (lo + hi) / 2.0;
        count += 1;
    }
    return mid, OK;
}

// The check of the shapes of the generalized exponential law. All three
// shapes must be above 0.
fn genexpon_check(a: f64, b: f64, c: f64) -> i64 {
    let ka = check_shape(a);
    if ka != OK { return ka; }
    let kb = check_shape(b);
    if kb != OK { return kb; }
    let kc = check_shape(c);
    if kc != OK { return kc; }
    return OK;
}

// The raw moment of order n of the generalized exponential law. It comes
// from an integration of the density of the law.
fn genexpon_raw(n: i64, a: f64, b: f64, c: f64) -> f64, i64 {
    var gx: f64[>..] = [];
    var gw: f64[>..] = [];
    gauss_legendre(gx, gw, FAM_NODES);
    var total: f64 = 0.0;
    // The range of the integral reaches the point where the density has
    // fallen below a small part of its peak. A law with a slow decay
    // needs a long range.
    let top = genexpon_reach(a, b, c);
    let width = top / FAM_PANELS as f64;
    for k in 0..FAM_PANELS {
        let left = width * k as f64;
        for j in 0..FAM_NODES {
            let x = left + width * (gx[j] + 1.0) / 2.0;
            let d, cd = genexpon_density(x, a, b, c);
            if cd != OK { return nan(), cd; }
            total += pow(x, n as f64) * d * gw[j] * width / 2.0;
        }
    }
    if is_nan(total) || is_inf(total) { return nan(), ERR_OVERFLOW; }
    // The range must hold the whole mass of the law. A law with a slow
    // decay needs a longer range, so the check gives the fit fault.
    if n == 1 {
        let m, cm = genexpon_mass(a, b, c, top);
        if cm != OK { return nan(), cm; }
        if m < 1.0 - 1e-9 { return nan(), ERR_FIT; }
    }
    return total, OK;
}

// The mass of the generalized exponential law up to a point. The form
// with the exponential keeps the digits of a slow decay.
fn genexpon_mass(a: f64, b: f64, c: f64, top: f64) -> f64, i64 {
    let e, ce = expm1(-c * top);
    if ce != OK { return nan(), ce; }
    let u = -e;
    let t = (-a - b) * top + b * u / c;
    return -expm1(t), OK;
}

// The point where the mass of the generalized exponential law is nearly
// one. The search doubles the point until the mass is above that bound.
fn genexpon_reach(a: f64, b: f64, c: f64) -> f64 {
    var top: f64 = 20.0;
    var i: i64 = 0;
    while i < 60 {
        let m, cm = genexpon_mass(a, b, c, top);
        if cm != OK { return top; }
        if m >= 1.0 - 1e-12 { return top; }
        top = top * 2.0;
        i += 1;
    }
    return top;
}

fn genexpon_mean(a: f64, b: f64, c: f64) -> f64, i64 {
    let k = genexpon_check(a, b, c);
    if k != OK { return nan(), k; }
    return genexpon_raw(1, a, b, c);
}

fn genexpon_var(a: f64, b: f64, c: f64) -> f64, i64 {
    let k = genexpon_check(a, b, c);
    if k != OK { return nan(), k; }
    let m1, c1 = genexpon_raw(1, a, b, c);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = genexpon_raw(2, a, b, c);
    if c2 != OK { return nan(), c2; }
    return m2 - m1 * m1, OK;
}

fn genexpon_skew(a: f64, b: f64, c: f64) -> f64, i64 {
    let k = genexpon_check(a, b, c);
    if k != OK { return nan(), k; }
    let m1, c1 = genexpon_raw(1, a, b, c);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = genexpon_raw(2, a, b, c);
    if c2 != OK { return nan(), c2; }
    let m3, c3 = genexpon_raw(3, a, b, c);
    if c3 != OK { return nan(), c3; }
    return fam_skew_raw(m1, m2, m3), OK;
}

fn genexpon_kurt(a: f64, b: f64, c: f64) -> f64, i64 {
    let k = genexpon_check(a, b, c);
    if k != OK { return nan(), k; }
    let m1, c1 = genexpon_raw(1, a, b, c);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = genexpon_raw(2, a, b, c);
    if c2 != OK { return nan(), c2; }
    let m3, c3 = genexpon_raw(3, a, b, c);
    if c3 != OK { return nan(), c3; }
    let m4, c4 = genexpon_raw(4, a, b, c);
    if c4 != OK { return nan(), c4; }
    return fam_kurt_raw(m1, m2, m3, m4), OK;
}

// ---- mielke: the law of Mielke of the beta kind --------------------------

fn mielke_pdf(x: f64, k: f64, s: f64) -> f64, i64 {
    let ck = mielke_check(k, s);
    if ck != OK { return nan(), ck; }
    let v = fam_x(x);
    if v != OK { return nan(), v; }
    if x <= 0.0 { return 0.0, OK; }
    let p = k * pow(x, k - 1.0) / pow(1.0 + pow(x, s), 1.0 + k / s);
    if is_inf(p) { return nan(), ERR_OVERFLOW; }
    return p, OK;
}

fn mielke_cdf(x: f64, k: f64, s: f64) -> f64, i64 {
    let ck = mielke_check(k, s);
    if ck != OK { return nan(), ck; }
    let v = fam_x(x);
    if v != OK { return nan(), v; }
    if x <= 0.0 { return 0.0, OK; }
    return pow(x, k) / pow(1.0 + pow(x, s), k / s), OK;
}

fn mielke_sf(x: f64, k: f64, s: f64) -> f64, i64 {
    let v, k = mielke_cdf(x, k, s);
    if k != OK { return nan(), k; }
    return 1.0 - v, OK;
}

fn mielke_ppf(q: f64, k: f64, s: f64) -> f64, i64 {
    let ck = mielke_check(k, s);
    if ck != OK { return nan(), ck; }
    let pc = check_prob(q);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(q, 0.0, inf());
    if e { return ev, ec; }
    let t = pow(q, s / k);
    return pow(t / (1.0 - t), 1.0 / s), OK;
}

// The check of the shapes of the law of Mielke. Both shapes must be above 0.
fn mielke_check(k: f64, s: f64) -> i64 {
    let ck = check_shape(k);
    if ck != OK { return ck; }
    let cs = check_shape(s);
    if cs != OK { return cs; }
    return OK;
}

// The raw moment of order n of the law of Mielke. It is the beta function
// of (k + n) / s and 1 - n / s, over the beta function of k / s and 1.
// The moment exists when the shape s is above n.
fn mielke_raw(n: f64, k: f64, s: f64) -> f64, i64 {
    if n >= s { return nan(), ERR_UNDEFINED; }
    let top, ct = beta((k + n) / s, 1.0 - n / s);
    if ct != OK { return nan(), ct; }
    let bot, cb = beta(k / s, 1.0);
    if cb != OK { return nan(), cb; }
    return top / bot, OK;
}

fn mielke_mean(k: f64, s: f64) -> f64, i64 {
    let ck = mielke_check(k, s);
    if ck != OK { return nan(), ck; }
    return mielke_raw(1.0, k, s);
}

fn mielke_var(k: f64, s: f64) -> f64, i64 {
    let ck = mielke_check(k, s);
    if ck != OK { return nan(), ck; }
    let m1, c1 = mielke_raw(1.0, k, s);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = mielke_raw(2.0, k, s);
    if c2 != OK { return nan(), c2; }
    return m2 - m1 * m1, OK;
}

fn mielke_skew(k: f64, s: f64) -> f64, i64 {
    let ck = mielke_check(k, s);
    if ck != OK { return nan(), ck; }
    let m1, c1 = mielke_raw(1.0, k, s);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = mielke_raw(2.0, k, s);
    if c2 != OK { return nan(), c2; }
    let m3, c3 = mielke_raw(3.0, k, s);
    if c3 != OK { return nan(), c3; }
    return fam_skew_raw(m1, m2, m3), OK;
}

fn mielke_kurt(k: f64, s: f64) -> f64, i64 {
    let ck = mielke_check(k, s);
    if ck != OK { return nan(), ck; }
    let m1, c1 = mielke_raw(1.0, k, s);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = mielke_raw(2.0, k, s);
    if c2 != OK { return nan(), c2; }
    let m3, c3 = mielke_raw(3.0, k, s);
    if c3 != OK { return nan(), c3; }
    let m4, c4 = mielke_raw(4.0, k, s);
    if c4 != OK { return nan(), c4; }
    return fam_kurt_raw(m1, m2, m3, m4), OK;
}
