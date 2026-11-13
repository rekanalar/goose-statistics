// The maximum likelihood fits of the gamma, beta and Weibull families,
// for the goodness-of-fit test. The location is fixed at 0, and the beta
// scale is fixed at 1.
//
// The gamma shape solves log(a) - digamma(a) = log(mean) - mean(log x).
// The beta pair solves the two digamma equations by a nested search: the
// inner search finds b for an a, and the outer search finds a. The
// Weibull shape solves the likelihood equation of the shape. All searches
// use the method of Brent. SciPy stops its beta search at the tolerance
// 1e-4, so the statistic of a beta fit agrees with SciPy to about 1e-4.
namespace statistics;

import std;
import math;

// The sums that the three fits read: the mean, the mean of the logs and
// the mean of the logs of one minus x.
fn mle_means(s: const f64[:], mean: f64&, mlog: f64&, mlog1: f64&) {
    var a: f64 = 0.0;
    var b: f64 = 0.0;
    var c: f64 = 0.0;
    for v in s { a += v; b += log(v); c += log(1.0 - v); }
    let n = s.len as f64;
    mean = a / n;
    mlog = b / n;
    mlog1 = c / n;
}

fn mle_gamma_gap(a: f64, sd: f64) -> f64 {
    let d, c = digamma(a);
    if c != OK { return nan(); }
    return log(a) - d - sd;
}

// The gamma fit: the shape and the scale.
fn mle_gamma(s: const f64[:]) -> f64, f64, i64 {
    var mean: f64 = 0.0;
    var mlog: f64 = 0.0;
    var mlog1: f64 = 0.0;
    mle_means(s, mean, mlog, mlog1);
    let sd = log(mean) - mlog;
    if !(sd > 0.0) || is_inf(sd) { return nan(), nan(), ERR_FIT; }
    let a0 = (3.0 - sd + sqrt((sd - 3.0) * (sd - 3.0) + 24.0 * sd)) / (12.0 * sd);
    var lo = a0 / 2.0;
    var hi = a0 * 2.0;
    for k in 0..60 { if mle_gamma_gap(lo, sd) > 0.0 { break; } lo /= 2.0; }
    for k in 0..60 { if mle_gamma_gap(hi, sd) < 0.0 { break; } hi *= 2.0; }
    let a, c = brent_root(lo, hi, 1e-15, 4e-16) { v => mle_gamma_gap(v, sd) };
    if c != OK { return nan(), nan(), ERR_FIT; }
    let scale = mean / a;
    if !(scale > 0.0) || is_inf(scale) { return nan(), nan(), ERR_FIT; }
    return a, scale, OK;
}

// The b of the inner equation of the beta fit, for one a.
fn mle_beta_b(a: f64, mlog1: f64) -> f64 {
    var lo: f64 = 1.0;
    var hi: f64 = 1.0;
    for k in 0..80 { if mle_beta_inner(a, lo, mlog1) < 0.0 { break; } lo /= 2.0; }
    for k in 0..80 { if mle_beta_inner(a, hi, mlog1) > 0.0 { break; } hi *= 2.0; }
    let b, c = brent_root(lo, hi, 1e-15, 4e-16) { v => mle_beta_inner(a, v, mlog1) };
    if c != OK { return nan(); }
    return b;
}

fn mle_beta_inner(a: f64, b: f64, mlog1: f64) -> f64 {
    let d1, c1 = digamma(b);
    let d2, c2 = digamma(a + b);
    if c1 != OK || c2 != OK { return nan(); }
    return d1 - d2 - mlog1;
}

fn mle_beta_outer(a: f64, mlog: f64, mlog1: f64) -> f64 {
    let b = mle_beta_b(a, mlog1);
    if is_nan(b) { return nan(); }
    let d1, c1 = digamma(a);
    let d2, c2 = digamma(a + b);
    if c1 != OK || c2 != OK { return nan(); }
    return d1 - d2 - mlog;
}

// The beta fit: the two shapes.
fn mle_beta(s: const f64[:]) -> f64, f64, i64 {
    var mean: f64 = 0.0;
    var mlog: f64 = 0.0;
    var mlog1: f64 = 0.0;
    mle_means(s, mean, mlog, mlog1);
    if !(mlog < 0.0) || !(mlog1 < 0.0) { return nan(), nan(), ERR_FIT; }
    var v: f64 = 0.0;
    for x in s { v += (x - mean) * (x - mean); }
    v /= s.len - 1;
    var a0: f64 = 1.0;
    if v > 0.0 { a0 = mean * (mean * (1.0 - mean) / v - 1.0); }
    if !(a0 > 0.0) { a0 = 1.0; }
    var lo = a0;
    var hi = a0;
    for k in 0..60 { if mle_beta_outer(lo, mlog, mlog1) < 0.0 { break; } lo /= 2.0; }
    for k in 0..60 { if mle_beta_outer(hi, mlog, mlog1) > 0.0 { break; } hi *= 2.0; }
    let a, c = brent_root(lo, hi, 1e-15, 4e-16) { t => mle_beta_outer(t, mlog, mlog1) };
    if c != OK { return nan(), nan(), ERR_FIT; }
    let b = mle_beta_b(a, mlog1);
    if is_nan(b) { return nan(), nan(), ERR_FIT; }
    return a, b, OK;
}

// The Weibull fit with the location 0: the shape and the scale.
fn mle_weibull(s: const f64[:]) -> f64, f64, i64 {
    var logs: f64[>..] = [];
    var lmax = -inf();
    var mlog: f64 = 0.0;
    for v in s {
        let l = log(v);
        logs.push(l);
        if l > lmax { lmax = l; }
        mlog += l / s.len;
    }
    if lmax == -inf() || lmax == inf() { return nan(), nan(), ERR_FIT; }
    let m = weibull_shape(logs[..], lmax, mlog);
    if !(m > 0.0) { return nan(), nan(), ERR_FIT; }
    var sp: f64 = 0.0;
    for l in logs { sp += exp(m * (l - lmax)); }
    let scale = exp(lmax + log(sp / s.len) / m);
    if !(scale > 0.0) || is_inf(scale) { return nan(), nan(), ERR_FIT; }
    return m, scale, OK;
}

// The cdf and the survival function of the three families.
fn mle_cdf(family: i64, x: f64, a: f64, b: f64) -> f64, i64 {
    if family == FAMILY_GAMMA { return gammadist_cdf(x, a, b); }
    if family == FAMILY_BETA { return betadist_cdf(x, a, b); }
    return weibull_cdf(x, a, b);
}

fn mle_sf(family: i64, x: f64, a: f64, b: f64) -> f64, i64 {
    if family == FAMILY_GAMMA { return gammadist_sf(x, a, b); }
    if family == FAMILY_BETA { return betadist_sf(x, a, b); }
    return weibull_sf(x, a, b);
}

// One draw from a fitted family by the inverse of the cdf at u.
fn mle_draw(family: i64, a: f64, b: f64, u: f64) -> f64 {
    var v: f64 = nan();
    var c: i64 = OK;
    if u <= 0.0 || u >= 1.0 { return nan(); }
    if family == FAMILY_GAMMA { v, c = gammadist_ppf(u, a, b); }
    else if family == FAMILY_BETA { v, c = betadist_ppf(u, a, b); }
    else { v, c = weibull_ppf(u, a, b); }
    if c != OK { return nan(); }
    return v;
}
