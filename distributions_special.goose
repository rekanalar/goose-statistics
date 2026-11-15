// Four continuous families of a special shape: the inverse Gaussian law,
// the Nakagami law, the Kappa law of the third kind and the Pearson law of
// the third kind.
//
// The names and the parameters follow SciPy.
//
// Every family has eight calls: pdf, cdf, sf, ppf, mean, var, skew and
// kurt. The kurtosis is the excess kurtosis. Every call returns the result
// first and a fault code last. A value outside the support gives a density
// of 0. The quantile of 0 or 1 gives the bound of the support.
namespace statistics;

import std;
import math;

// The logarithm of the gamma function for an argument above 0.
fn fam_lngamma(x: f64) -> f64 {
    let v, s, c = lgamma(x);
    if c != OK { return nan(); }
    if s < 0.0 { return nan(); }
    return v;
}

// The logarithm of the cumulative function of the normal law. The far tail
// uses the series of Mills, so the value keeps its digits there.
fn fam_log_ndtr(z: f64) -> f64 {
    if z > -26.0 {
        let e, c = erfc(-z * INV_SQRT2);
        if c != OK { return nan(); }
        if e <= 0.0 { return -inf(); }
        return log(0.5 * e);
    }
    let t = 1.0 / (z * z);
    let s = 1.0 - t * (1.0 - t * (3.0 - t * (15.0 - t * 105.0)));
    return -0.5 * z * z - log(-z) - 0.5 * log(2.0 * PI) + log(s);
}

// ---- invgauss: the inverse Gaussian law -----------------------------------

fn invgauss_pdf(x: f64, mu: f64) -> f64, i64 {
    let c = check_shape(mu);
    if c != OK { return nan(), c; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    let d = x / mu - 1.0;
    let lp = -0.5 * log(2.0 * PI) - 1.5 * log(x) - d * d / (2.0 * x);
    if is_inf(lp) { return 0.0, OK; }
    return exp(lp), OK;
}

// The logarithm of the cumulative function of the inverse Gaussian law.
// The value is a sum of two parts, so the larger part carries it and the
// difference of the parts gives the small correction.
fn invgauss_log_cdf(x: f64, mu: f64) -> f64 {
    let fac = 1.0 / sqrt(x);
    let la = fam_log_ndtr(fac * (x / mu - 1.0));
    let lb = 2.0 / mu + fam_log_ndtr(-fac * (x / mu + 1.0));
    if is_nan(la) || is_nan(lb) { return nan(); }
    let l, cl = log1p(exp(-abs(la - lb)));
    if cl != OK { return nan(); }
    return max(la, lb) + l;
}

// The logarithm of the survival function of the inverse Gaussian law. The
// survival function is the difference of two parts, so the form with expm1
// keeps the digits of a small tail.
fn invgauss_log_sf(x: f64, mu: f64) -> f64 {
    let fac = 1.0 / sqrt(x);
    let la = fam_log_ndtr(-fac * (x / mu - 1.0));
    let lb = 2.0 / mu + fam_log_ndtr(-fac * (x / mu + 1.0));
    if is_nan(la) || is_nan(lb) { return nan(); }
    let d = lb - la;
    if d >= 0.0 { return -inf(); }
    let e, ce = expm1(d);
    if ce != OK { return nan(); }
    let l = log(-e);
    if is_nan(l) { return nan(); }
    return la + l;
}

fn invgauss_cdf(x: f64, mu: f64) -> f64, i64 {
    let c = check_shape(mu);
    if c != OK { return nan(), c; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    let l = invgauss_log_cdf(x, mu);
    if is_nan(l) { return nan(), ERR_NAN; }
    let v = exp(l);
    if v > 1.0 { return 1.0, OK; }
    return v, OK;
}

fn invgauss_sf(x: f64, mu: f64) -> f64, i64 {
    let c = check_shape(mu);
    if c != OK { return nan(), c; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 1.0, OK; }
    let l = invgauss_log_sf(x, mu);
    if is_nan(l) { return nan(), ERR_NAN; }
    let v = exp(l);
    if v > 1.0 { return 1.0, OK; }
    return v, OK;
}

// A root of the logarithm of the cumulative function at the target. The
// function rises with x, so the search needs a point below the target and
// a point above it. The bracket stays within a factor of two of the root.
fn invgauss_root_cdf(mu: f64, lq: f64) -> f64, i64 {
    var hi = mu;
    var iter: i64 = 0;
    while invgauss_log_cdf(hi, mu) < lq && iter < 400 {
        hi = hi * 2.0;
        iter += 1;
    }
    var lo = hi * 0.5;
    iter = 0;
    while lo > 0.0 && invgauss_log_cdf(lo, mu) >= lq && iter < 400 {
        lo = lo * 0.5;
        iter += 1;
    }
    let x, cx = brent_root(lo, hi, 1e-15, 4e-16) { v => invgauss_log_cdf(v, mu) - lq };
    if cx != OK { return nan(), cx; }
    return x, OK;
}

// A root of the logarithm of the survival function at the target. The
// function falls with x, so the search needs a point below the target and
// a point above it. The bracket stays within a factor of two of the root.
fn invgauss_root_sf(mu: f64, lq: f64) -> f64, i64 {
    var hi = mu;
    var iter: i64 = 0;
    while invgauss_log_sf(hi, mu) > lq && iter < 400 {
        hi = hi * 2.0;
        iter += 1;
    }
    var lo = hi * 0.5;
    iter = 0;
    while lo > 0.0 && invgauss_log_sf(lo, mu) <= lq && iter < 400 {
        lo = lo * 0.5;
        iter += 1;
    }
    let x, cx = brent_root(lo, hi, 1e-15, 4e-16) { v => invgauss_log_sf(v, mu) - lq };
    if cx != OK { return nan(), cx; }
    return x, OK;
}

fn invgauss_ppf(q: f64, mu: f64) -> f64, i64 {
    let c = check_shape(mu);
    if c != OK { return nan(), c; }
    let pc = check_prob(q);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(q, 0.0, inf());
    if e { return ev, ec; }
    // The search runs on the logarithm of the cumulative function or of
    // the survival function. A tiny target keeps its digits, and the
    // search never forms one minus the target.
    if q <= 0.5 { return invgauss_root_cdf(mu, log(q)); }
    return invgauss_root_sf(mu, log(1.0 - q));
}

fn invgauss_mean(mu: f64) -> f64, i64 {
    let c = check_shape(mu);
    if c != OK { return nan(), c; }
    return mu, OK;
}

fn invgauss_var(mu: f64) -> f64, i64 {
    let c = check_shape(mu);
    if c != OK { return nan(), c; }
    return mu * mu * mu, OK;
}

fn invgauss_skew(mu: f64) -> f64, i64 {
    let c = check_shape(mu);
    if c != OK { return nan(), c; }
    return 3.0 * sqrt(mu), OK;
}

fn invgauss_kurt(mu: f64) -> f64, i64 {
    let c = check_shape(mu);
    if c != OK { return nan(), c; }
    return 15.0 * mu, OK;
}

// ---- nakagami: the Nakagami law -------------------------------------------

fn nakagami_pdf(x: f64, nu: f64) -> f64, i64 {
    let c = check_shape(nu);
    if c != OK { return nan(), c; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x < 0.0 { return 0.0, OK; }
    if x == 0.0 {
        if nu < 0.5 { return inf(), OK; }
        if nu > 0.5 { return 0.0, OK; }
        // The exponent of x is zero here, so the density has a limit.
        return sqrt(2.0 / PI), OK;
    }
    let lp = log(2.0) + nu * log(nu) - fam_lngamma(nu) + (2.0 * nu - 1.0) * log(x) - nu * x * x;
    if is_inf(lp) { return 0.0, OK; }
    return exp(lp), OK;
}

fn nakagami_cdf(x: f64, nu: f64) -> f64, i64 {
    let c = check_shape(nu);
    if c != OK { return nan(), c; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x < 0.0 { return 0.0, OK; }
    let v, cv = gammainc(nu, nu * x * x);
    if cv != OK { return nan(), cv; }
    return v, OK;
}

fn nakagami_sf(x: f64, nu: f64) -> f64, i64 {
    let c = check_shape(nu);
    if c != OK { return nan(), c; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x < 0.0 { return 1.0, OK; }
    let v, cv = gammaincc(nu, nu * x * x);
    if cv != OK { return nan(), cv; }
    return v, OK;
}

fn nakagami_ppf(q: f64, nu: f64) -> f64, i64 {
    let c = check_shape(nu);
    if c != OK { return nan(), c; }
    let pc = check_prob(q);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(q, 0.0, inf());
    if e { return ev, ec; }
    let v, cv = gammaincinv(nu, q);
    if cv != OK { return nan(), cv; }
    return sqrt(v / nu), OK;
}

fn nakagami_mean(nu: f64) -> f64, i64 {
    let c = check_shape(nu);
    if c != OK { return nan(), c; }
    let m = nakagami_mu(nu);
    if is_nan(m) { return nan(), ERR_UNDEFINED; }
    return m, OK;
}

// The mean of the Nakagami law. It is the ratio of two gamma functions
// over the square root of the shape.
fn nakagami_mu(nu: f64) -> f64 {
    let lg = fam_lngamma(nu + 0.5) - fam_lngamma(nu);
    if is_nan(lg) { return nan(); }
    return exp(lg) / sqrt(nu);
}

fn nakagami_var(nu: f64) -> f64, i64 {
    let c = check_shape(nu);
    if c != OK { return nan(), c; }
    let m = nakagami_mu(nu);
    if is_nan(m) { return nan(), ERR_UNDEFINED; }
    return 1.0 - m * m, OK;
}

fn nakagami_skew(nu: f64) -> f64, i64 {
    let c = check_shape(nu);
    if c != OK { return nan(), c; }
    let m = nakagami_mu(nu);
    if is_nan(m) { return nan(), ERR_UNDEFINED; }
    let v = 1.0 - m * m;
    if v <= 0.0 { return nan(), ERR_UNDEFINED; }
    return m * (1.0 - 4.0 * nu * v) / 2.0 / nu / (v * sqrt(v)), OK;
}

fn nakagami_kurt(nu: f64) -> f64, i64 {
    let c = check_shape(nu);
    if c != OK { return nan(), c; }
    let m = nakagami_mu(nu);
    if is_nan(m) { return nan(), ERR_UNDEFINED; }
    let v = 1.0 - m * m;
    if v <= 0.0 { return nan(), ERR_UNDEFINED; }
    let m2 = m * m;
    let k = (-6.0 * m2 * m2 * nu + (8.0 * nu - 2.0) * m2 - 2.0 * nu + 1.0) / (nu * v * v);
    return k, OK;
}

// ---- kappa3: the Kappa law of the third kind -----------------------------

// The logarithm of a + x^a. The form keeps the digits of a large power.
fn kappa3_log_sum(x: f64, a: f64) -> f64 {
    let u = a * log(x);
    let v = log(a);
    if u > v {
        let l, cl = log1p(exp(v - u));
        if cl != OK { return u; }
        return u + l;
    }
    let l, cl = log1p(exp(u - v));
    if cl != OK { return v; }
    return v + l;
}

fn kappa3_pdf(x: f64, a: f64) -> f64, i64 {
    let c = check_shape(a);
    if c != OK { return nan(), c; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x < 0.0 { return 0.0, OK; }
    if x == 0.0 {
        // The density has a limit at zero, which is a power of the shape.
        return pow(a, -1.0 / a), OK;
    }
    let lp = log(a) - (1.0 / a + 1.0) * kappa3_log_sum(x, a);
    if is_inf(lp) { return 0.0, OK; }
    return exp(lp), OK;
}

fn kappa3_cdf(x: f64, a: f64) -> f64, i64 {
    let c = check_shape(a);
    if c != OK { return nan(), c; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    let lp = log(x) - kappa3_log_sum(x, a) / a;
    if is_inf(lp) { return 1.0, OK; }
    return exp(lp), OK;
}

fn kappa3_sf(x: f64, a: f64) -> f64, i64 {
    let c = check_shape(a);
    if c != OK { return nan(), c; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 1.0, OK; }
    // The form keeps the digits of a small tail.
    let u = a * pow(x, -a);
    if is_inf(u) { return 1.0, OK; }
    let l, cl = log1p(u);
    if cl != OK { return nan(), cl; }
    let e, ce = expm1(-l / a);
    if ce != OK { return nan(), ce; }
    return -e, OK;
}

fn kappa3_ppf(q: f64, a: f64) -> f64, i64 {
    let c = check_shape(a);
    if c != OK { return nan(), c; }
    let pc = check_prob(q);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(q, 0.0, inf());
    if e { return ev, ec; }
    let u = pow(q, -a) - 1.0;
    if u <= 0.0 { return nan(), ERR_OVERFLOW; }
    let r = pow(a / u, 1.0 / a);
    if is_inf(r) { return nan(), ERR_OVERFLOW; }
    return r, OK;
}

// The raw moment of order n of the Kappa law of the third kind. The moment
// exists only when the shape is above n.
fn kappa3_raw(n: i64, a: f64) -> f64, i64 {
    if a <= n { return nan(), ERR_UNDEFINED; }
    let e = n / a;
    let lg = (e - 1.0) * log(a) + fam_lngamma((n + 1.0) / a)
        + fam_lngamma(1.0 - e) - fam_lngamma(1.0 + 1.0 / a);
    if is_nan(lg) { return nan(), ERR_UNDEFINED; }
    return exp(lg), OK;
}

fn kappa3_mean(a: f64) -> f64, i64 {
    let c = check_shape(a);
    if c != OK { return nan(), c; }
    return kappa3_raw(1, a);
}

fn kappa3_var(a: f64) -> f64, i64 {
    let c = check_shape(a);
    if c != OK { return nan(), c; }
    let m1, c1 = kappa3_raw(1, a);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = kappa3_raw(2, a);
    if c2 != OK { return nan(), c2; }
    return m2 - m1 * m1, OK;
}

fn kappa3_skew(a: f64) -> f64, i64 {
    let c = check_shape(a);
    if c != OK { return nan(), c; }
    let m1, c1 = kappa3_raw(1, a);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = kappa3_raw(2, a);
    if c2 != OK { return nan(), c2; }
    let m3, c3 = kappa3_raw(3, a);
    if c3 != OK { return nan(), c3; }
    return fam_skew_raw(m1, m2, m3), OK;
}

fn kappa3_kurt(a: f64) -> f64, i64 {
    let c = check_shape(a);
    if c != OK { return nan(), c; }
    let m1, c1 = kappa3_raw(1, a);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = kappa3_raw(2, a);
    if c2 != OK { return nan(), c2; }
    let m3, c3 = kappa3_raw(3, a);
    if c3 != OK { return nan(), c3; }
    let m4, c4 = kappa3_raw(4, a);
    if c4 != OK { return nan(), c4; }
    return fam_kurt_raw(m1, m2, m3, m4), OK;
}

// ---- pearson3: the Pearson law of the third kind --------------------------

// Below this skew the Pearson law of the third kind is the normal law. The
// shape of the law grows without bound as the skew falls to zero, so the
// law of the normal is the only form that keeps its digits there.
const PEARSON3_NORMAL_SKEW: f64 = 0.000016;

// True when the skew is small enough for the normal law.
fn pearson3_is_normal(skew: f64) -> bool {
    return abs(skew) < PEARSON3_NORMAL_SKEW;
}

// The shape, the scale and the shift of the Pearson law. The skewness of
// the law is the only parameter.
fn pearson3_shape(skew: f64) -> f64, f64, f64 {
    let beta = 2.0 / skew;
    let alpha = beta * beta;
    let zeta = -alpha / beta;
    return beta, alpha, zeta;
}

// The logarithm of the density of the gamma law at a point, for a shape
// above the threshold. The direct form takes the difference of two large
// numbers, so it loses every digit of a small skew. The form with the
// ratio of the point to the shape keeps them.
fn pearson3_log_gamma(tx: f64, alpha: f64) -> f64 {
    let r = tx / alpha;
    return alpha * (log(r) - r + 1.0) - log(r) - 0.5 * log(2.0 * PI * alpha)
        - 1.0 / (12.0 * alpha);
}

fn pearson3_pdf(x: f64, skew: f64) -> f64, i64 {
    let c = check_scalar(skew);
    if c != OK { return nan(), c; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if pearson3_is_normal(skew) { return norm_pdf(x); }
    let beta, alpha, zeta = pearson3_shape(skew);
    // A huge skewness gives a shape without bound. The gamma function of
    // the shape is not a number then.
    if is_inf(alpha) || alpha <= 0.0 { return nan(), ERR_OVERFLOW; }
    let tx = beta * (x - zeta);
    if tx <= 0.0 { return 0.0, OK; }
    var lp: f64;
    if alpha > 1.0e4 {
        lp = log(abs(beta)) + pearson3_log_gamma(tx, alpha);
    } else {
        lp = log(abs(beta)) + (alpha - 1.0) * log(tx) - tx - fam_lngamma(alpha);
    }
    if is_inf(lp) { return 0.0, OK; }
    return exp(lp), OK;
}

fn pearson3_cdf(x: f64, skew: f64) -> f64, i64 {
    let c = check_scalar(skew);
    if c != OK { return nan(), c; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if pearson3_is_normal(skew) { return norm_cdf(x); }
    let beta, alpha, zeta = pearson3_shape(skew);
    let tx = beta * (x - zeta);
    if tx <= 0.0 {
        if skew > 0.0 { return 0.0, OK; }
        return 1.0, OK;
    }
    if skew > 0.0 {
        let v, cv = gammainc(alpha, tx);
        if cv != OK { return nan(), cv; }
        return v, OK;
    }
    let v, cv = gammaincc(alpha, tx);
    if cv != OK { return nan(), cv; }
    return v, OK;
}

fn pearson3_sf(x: f64, skew: f64) -> f64, i64 {
    let c = check_scalar(skew);
    if c != OK { return nan(), c; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if pearson3_is_normal(skew) { return norm_sf(x); }
    let beta, alpha, zeta = pearson3_shape(skew);
    let tx = beta * (x - zeta);
    if tx <= 0.0 {
        if skew > 0.0 { return 1.0, OK; }
        return 0.0, OK;
    }
    if skew > 0.0 {
        let v, cv = gammaincc(alpha, tx);
        if cv != OK { return nan(), cv; }
        return v, OK;
    }
    let v, cv = gammainc(alpha, tx);
    if cv != OK { return nan(), cv; }
    return v, OK;
}

fn pearson3_ppf(q: f64, skew: f64) -> f64, i64 {
    let c = check_scalar(skew);
    if c != OK { return nan(), c; }
    let pc = check_prob(q);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(q, -inf(), inf());
    if e { return ev, ec; }
    if pearson3_is_normal(skew) { return norm_std_ppf(q), OK; }
    let beta, alpha, zeta = pearson3_shape(skew);
    var p = q;
    if skew < 0.0 { p = 1.0 - q; }
    let v, cv = gammaincinv(alpha, p);
    if cv != OK { return nan(), cv; }
    return v / beta + zeta, OK;
}

fn pearson3_mean(skew: f64) -> f64, i64 {
    let c = check_scalar(skew);
    if c != OK { return nan(), c; }
    return 0.0, OK;
}

fn pearson3_var(skew: f64) -> f64, i64 {
    let c = check_scalar(skew);
    if c != OK { return nan(), c; }
    return 1.0, OK;
}

fn pearson3_skew(skew: f64) -> f64, i64 {
    let c = check_scalar(skew);
    if c != OK { return nan(), c; }
    return skew, OK;
}

fn pearson3_kurt(skew: f64) -> f64, i64 {
    let c = check_scalar(skew);
    if c != OK { return nan(), c; }
    let r = 1.5 * skew * skew;
    if is_inf(r) { return nan(), ERR_OVERFLOW; }
    return r, OK;
}
