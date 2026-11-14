// Three continuous laws of a symmetric shape: gennorm, halfgennorm and
// tukeylambda. The names and the parameters follow SciPy.
//
// Every family has eight calls: pdf, cdf, sf, ppf, mean, var, skew and
// kurt. The kurtosis is the excess kurtosis. Every call returns the result
// first and a fault code last. A value outside the support gives a density
// of 0. The quantile of 0 or 1 gives the bound of the support.
//
// The laws of the generalized normal kind are symmetric around 0. The law
// of Tukey and of lambda has a shape that can bound the law.
namespace statistics;

import std;
import math;


// ---- gennorm: the generalized normal law --------------------------------

fn gennorm_pdf(x: f64, beta: f64) -> f64, i64 {
    let k = check_shape(beta);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    let g, cg = gamma(1.0 / beta);
    if cg != OK { return nan(), cg; }
    let p = beta * exp(-pow(abs(x), beta)) / (2.0 * g);
    if is_inf(p) { return nan(), ERR_OVERFLOW; }
    return p, OK;
}

fn gennorm_cdf(x: f64, beta: f64) -> f64, i64 {
    let k = check_shape(beta);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    let v, cv = gammaincc(1.0 / beta, pow(abs(x), beta));
    if cv != OK { return nan(), cv; }
    // The law is symmetric around 0, so the mass below the point is half of
    // the mass of the absolute value.
    if x >= 0.0 { return 1.0 - v / 2.0, OK; }
    return v / 2.0, OK;
}

fn gennorm_sf(x: f64, beta: f64) -> f64, i64 {
    let v, k = gennorm_cdf(-x, beta);
    if k != OK { return nan(), k; }
    return v, OK;
}

fn gennorm_ppf(q: f64, beta: f64) -> f64, i64 {
    let k = check_shape(beta);
    if k != OK { return nan(), k; }
    let pc = check_prob(q);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(q, -inf(), inf());
    if e { return ev, ec; }
    if q == 0.5 { return 0.0, OK; }
    if q > 0.5 {
        let v, cv = gammainccinv(1.0 / beta, 2.0 * (1.0 - q));
        if cv != OK { return nan(), cv; }
        return pow(v, 1.0 / beta), OK;
    }
    let v, cv = gammainccinv(1.0 / beta, 2.0 * q);
    if cv != OK { return nan(), cv; }
    return -pow(v, 1.0 / beta), OK;
}

fn gennorm_mean(beta: f64) -> f64, i64 {
    let k = check_shape(beta);
    if k != OK { return nan(), k; }
    // The law is symmetric around 0, so its mean is 0.
    return 0.0, OK;
}

fn gennorm_var(beta: f64) -> f64, i64 {
    let k = check_shape(beta);
    if k != OK { return nan(), k; }
    let g1, c1 = gamma(3.0 / beta);
    if c1 != OK { return nan(), c1; }
    let g0, c0 = gamma(1.0 / beta);
    if c0 != OK { return nan(), c0; }
    return g1 / g0, OK;
}

fn gennorm_skew(beta: f64) -> f64, i64 {
    let k = check_shape(beta);
    if k != OK { return nan(), k; }
    // The law is symmetric around 0, so its skewness is 0.
    return 0.0, OK;
}

fn gennorm_kurt(beta: f64) -> f64, i64 {
    let k = check_shape(beta);
    if k != OK { return nan(), k; }
    let g1, c1 = gamma(3.0 / beta);
    if c1 != OK { return nan(), c1; }
    let g0, c0 = gamma(1.0 / beta);
    if c0 != OK { return nan(), c0; }
    let g2, c2 = gamma(5.0 / beta);
    if c2 != OK { return nan(), c2; }
    return g2 * g0 / (g1 * g1) - 3.0, OK;
}

// ---- halfgennorm: the generalized half normal law ------------------------

fn halfgennorm_pdf(x: f64, beta: f64) -> f64, i64 {
    let k = check_shape(beta);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    let g, cg = gamma(1.0 / beta);
    if cg != OK { return nan(), cg; }
    let p = beta * exp(-pow(x, beta)) / g;
    if is_inf(p) { return nan(), ERR_OVERFLOW; }
    return p, OK;
}

fn halfgennorm_cdf(x: f64, beta: f64) -> f64, i64 {
    let k = check_shape(beta);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    let v, cv = gammainc(1.0 / beta, pow(x, beta));
    if cv != OK { return nan(), cv; }
    return v, OK;
}

fn halfgennorm_sf(x: f64, beta: f64) -> f64, i64 {
    let k = check_shape(beta);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 1.0, OK; }
    let v, cv = gammaincc(1.0 / beta, pow(x, beta));
    if cv != OK { return nan(), cv; }
    return v, OK;
}

fn halfgennorm_ppf(q: f64, beta: f64) -> f64, i64 {
    let k = check_shape(beta);
    if k != OK { return nan(), k; }
    let pc = check_prob(q);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(q, 0.0, inf());
    if e { return ev, ec; }
    let v, cv = gammaincinv(1.0 / beta, q);
    if cv != OK { return nan(), cv; }
    return pow(v, 1.0 / beta), OK;
}

// The raw moment of order n of the generalized half normal law. It is the
// ratio of the gamma function of (n + 1) over beta, and of 1 over beta.
fn halfgennorm_raw(n: f64, beta: f64) -> f64, i64 {
    let l1, s1, c1 = lgamma((n + 1.0) / beta);
    if c1 != OK { return nan(), c1; }
    let l2, s2, c2 = lgamma(1.0 / beta);
    if c2 != OK { return nan(), c2; }
    let v = exp(l1 - l2);
    if is_inf(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

fn halfgennorm_mean(beta: f64) -> f64, i64 {
    let k = check_shape(beta);
    if k != OK { return nan(), k; }
    return halfgennorm_raw(1.0, beta);
}

fn halfgennorm_var(beta: f64) -> f64, i64 {
    let k = check_shape(beta);
    if k != OK { return nan(), k; }
    let m1, c1 = halfgennorm_raw(1.0, beta);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = halfgennorm_raw(2.0, beta);
    if c2 != OK { return nan(), c2; }
    return m2 - m1 * m1, OK;
}

fn halfgennorm_skew(beta: f64) -> f64, i64 {
    let k = check_shape(beta);
    if k != OK { return nan(), k; }
    let m1, c1 = halfgennorm_raw(1.0, beta);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = halfgennorm_raw(2.0, beta);
    if c2 != OK { return nan(), c2; }
    let m3, c3 = halfgennorm_raw(3.0, beta);
    if c3 != OK { return nan(), c3; }
    return fam_skew_raw(m1, m2, m3), OK;
}

fn halfgennorm_kurt(beta: f64) -> f64, i64 {
    let k = check_shape(beta);
    if k != OK { return nan(), k; }
    let m1, c1 = halfgennorm_raw(1.0, beta);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = halfgennorm_raw(2.0, beta);
    if c2 != OK { return nan(), c2; }
    let m3, c3 = halfgennorm_raw(3.0, beta);
    if c3 != OK { return nan(), c3; }
    let m4, c4 = halfgennorm_raw(4.0, beta);
    if c4 != OK { return nan(), c4; }
    return fam_kurt_raw(m1, m2, m3, m4), OK;
}
// ---- tukeylambda: the law of Tukey and of lambda ------------------------

// The quantile of the law of Tukey and of lambda. It is the power of the
// probability less the power of one minus the probability, over the shape.
fn tukeylambda_ppf_formula(p: f64, lam: f64) -> f64, i64 {
    if lam == 0.0 {
        // The shape of 0 gives the law of the logistic kind.
        return log(p / (1.0 - p)), OK;
    }
    let a = pow(p, lam);
    let b = pow(1.0 - p, lam);
    return (a - b) / lam, OK;
}

// The bounds of the law of Tukey and of lambda. A shape above 0 bounds the
// law on both sides.
fn tukeylambda_bounds(lam: f64) -> f64, f64 {
    if lam > 0.0 { return -1.0 / lam, 1.0 / lam; }
    return -inf(), inf();
}

fn tukeylambda_pdf(x: f64, lam: f64) -> f64, i64 {
    let k = check_scalar(lam);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    let lo, hi = tukeylambda_bounds(lam);
    if x <= lo || x >= hi { return 0.0, OK; }
    // The density comes from the derivative of the quantile. The law is
    // symmetric, so the sign of the point does not matter.
    let f, cf = tukeylambda_cdf(x, lam);
    if cf != OK { return nan(), cf; }
    if lam == 0.0 {
        // The shape of 0 gives the law of the logistic kind.
        let e = exp(-abs(x));
        return e / ((1.0 + e) * (1.0 + e)), OK;
    }
    let t = pow(f, lam - 1.0) + pow(1.0 - f, lam - 1.0);
    if t <= 0.0 { return 0.0, OK; }
    return 1.0 / t, OK;
}

fn tukeylambda_cdf(x: f64, lam: f64) -> f64, i64 {
    let k = check_scalar(lam);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    let lo, hi = tukeylambda_bounds(lam);
    if x <= lo { return 0.0, OK; }
    if x >= hi { return 1.0, OK; }
    if lam == 0.0 {
        // The shape of 0 gives the law of the logistic kind.
        if x >= 0.0 { return 1.0 / (1.0 + exp(-x)), OK; }
        let e = exp(x);
        return e / (1.0 + e), OK;
    }
    // The cumulative function comes from a bisection of the quantile. The
    // quantile grows with the probability, so the search narrows the
    // probability between two ends.
    var plow = 0.0;
    var phigh = 1.0;
    var pmid = 0.5;
    var count: i64 = 0;
    while count < 60 && abs(pmid - plow) > 1.0e-18 {
        let q, cq = tukeylambda_ppf_formula(pmid, lam);
        if cq != OK { return nan(), cq; }
        if q == x { return pmid, OK; }
        if q > x {
            phigh = pmid;
            pmid = (pmid + plow) / 2.0;
        } else {
            plow = pmid;
            pmid = (pmid + phigh) / 2.0;
        }
        count += 1;
    }
    return pmid, OK;
}

fn tukeylambda_sf(x: f64, lam: f64) -> f64, i64 {
    let v, k = tukeylambda_cdf(x, lam);
    if k != OK { return nan(), k; }
    return 1.0 - v, OK;
}

fn tukeylambda_ppf(p: f64, lam: f64) -> f64, i64 {
    let k = check_scalar(lam);
    if k != OK { return nan(), k; }
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    let lo, hi = tukeylambda_bounds(lam);
    let e, ev, ec = fam_end(p, lo, hi);
    if e { return ev, ec; }
    return tukeylambda_ppf_formula(p, lam);
}

// The variance of the law of Tukey and of lambda. It comes from the beta
// function. The shape of 0 gives the law of the logistic kind.
fn tukeylambda_var(lam: f64) -> f64, i64 {
    let k = check_scalar(lam);
    if k != OK { return nan(), k; }
    if lam == 0.0 { return PI * PI / 3.0, OK; }
    let b, cb = beta(lam + 1.0, lam + 1.0);
    if cb != OK { return nan(), cb; }
    let d = 1.0 / (1.0 + 2.0 * lam) - b;
    if lam < -0.5 { return nan(), ERR_UNDEFINED; }
    if lam == -0.5 { return inf(), ERR_OVERFLOW; }
    return 2.0 / (lam * lam) * d, OK;
}

fn tukeylambda_mean(lam: f64) -> f64, i64 {
    let k = check_scalar(lam);
    if k != OK { return nan(), k; }
    // The law is symmetric around 0, so its mean is 0.
    return 0.0, OK;
}

fn tukeylambda_skew(lam: f64) -> f64, i64 {
    let k = check_scalar(lam);
    if k != OK { return nan(), k; }
    // The law is symmetric around 0, so its skewness is 0.
    return 0.0, OK;
}

fn tukeylambda_kurt(lam: f64) -> f64, i64 {
    let k = check_scalar(lam);
    if k != OK { return nan(), k; }
    if lam == 0.0 { return 1.2, OK; }
    if lam < -0.25 { return nan(), ERR_UNDEFINED; }
    if lam == -0.25 { return inf(), ERR_OVERFLOW; }
    let b1, cb1 = beta(lam + 1.0, lam + 1.0);
    if cb1 != OK { return nan(), cb1; }
    let b2, cb2 = beta(2.0 * lam + 1.0, 2.0 * lam + 1.0);
    if cb2 != OK { return nan(), cb2; }
    let b3, cb3 = beta(3.0 * lam + 1.0, lam + 1.0);
    if cb3 != OK { return nan(), cb3; }
    let num = 1.0 / (4.0 * lam + 1.0) - 4.0 * b3 + 3.0 * b2;
    let den = 2.0 * (1.0 / (2.0 * lam + 1.0) - b1) * (1.0 / (2.0 * lam + 1.0) - b1);
    if den == 0.0 { return nan(), ERR_UNDEFINED; }
    return num / den - 3.0, OK;
}
