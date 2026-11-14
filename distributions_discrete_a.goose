// Seven discrete families: bernoulli, betabinom, betanbinom, boltzmann,
// dlaplace, logser and nhypergeom. The names and the parameters follow
// SciPy.
//
// Every family has three calls: pmf, cdf and ppf. The mass of a value
// outside the support is 0 and the code is OK. The quantile of 0 or 1
// gives the bound of the support. Every call returns the result first and
// a fault code last.
//
// The laws of this file need the beta law and the Poisson law of the
// package.
namespace statistics;

import std;
import math;

// The check of a count that must be a whole number.
fn disc_whole(x: f64) -> i64 {
    let c = check_scalar(x);
    if c != OK { return c; }
    if x < 0.0 { return ERR_DOMAIN; }
    if x != floor(x) { return ERR_DOMAIN; }
    if x > MAX_COUNT { return ERR_DOMAIN; }
    return OK;
}

// The check of a count of the support of a discrete law.
fn disc_count(k: f64) -> i64 {
    let c = check_scalar(k);
    if c != OK { return c; }
    if k != floor(k) { return ERR_DOMAIN; }
    return OK;
}

// The smallest count whose cumulative function reaches the probability.
// The search walks the counts of the law. The law must have no end at
// the upper side, so the search grows its bound until the cumulative
// function reaches the probability. The block gives the cumulative
// function of a count.
fn disc_find<F>(q: f64, lo: f64, hi: f64) -> f64, i64 {
    var l = lo;
    var r = hi;
    while l < r {
        let mid = floor(0.5 * (l + r));
        if F(mid) >= q { r = mid; } else { l = mid + 1.0; }
    }
    return l, OK;
}

// The smallest count whose cumulative function reaches the probability.
// The law has no end, so the search grows its bound until the cumulative
// function reaches the probability.
fn disc_find_open<F>(q: f64, lo: f64) -> f64, i64 {
    var top = lo + 1.0;
    var steps: i64 = 0;
    var found = false;
    while steps < MAX_ITER_SEARCH * 1000 {
        if F(top) >= q { found = true; break; }
        top = top * 2.0;
        steps += 1;
    }
    if !found { return nan(), ERR_FIT; }
    return disc_find(q, lo, top) { x => F(x) };
}

// ---- bernoulli: the law of 0 and 1 ---------------------------------------

fn bernoulli_pmf(k: f64, p: f64) -> f64, i64 {
    let c = check_prob(p);
    if c != OK { return nan(), c; }
    let kc = disc_count(k);
    if kc != OK { return nan(), kc; }
    if k < 0.0 || k > 1.0 { return 0.0, OK; }
    if k == 0.0 { return 1.0 - p, OK; }
    return p, OK;
}

fn bernoulli_cdf(k: f64, p: f64) -> f64, i64 {
    let c = check_prob(p);
    if c != OK { return nan(), c; }
    let kc = disc_count(k);
    if kc != OK { return nan(), kc; }
    if k < 0.0 { return 0.0, OK; }
    if k < 1.0 { return 1.0 - p, OK; }
    return 1.0, OK;
}

fn bernoulli_ppf(q: f64, p: f64) -> f64, i64 {
    let c = check_prob(q);
    if c != OK { return nan(), c; }
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    if q == 0.0 { return 0.0, OK; }
    if q == 1.0 { return 1.0, OK; }
    if q <= 1.0 - p { return 0.0, OK; }
    return 1.0, OK;
}

// ---- betabinom: the beta binomial law -----------------------------------

// The mass of the beta binomial law in the logarithm. The binomial
// coefficient and the beta function enter in the logarithm.
fn betabinom_logpmf(k: f64, n: f64, a: f64, b: f64) -> f64, i64 {
    let lc, lcc = lchoose(n, k);
    if lcc != OK { return nan(), lcc; }
    let lb1, lcb1 = lbeta(k + a, n - k + b);
    if lcb1 != OK { return nan(), lcb1; }
    let lb2, lcb2 = lbeta(a, b);
    if lcb2 != OK { return nan(), lcb2; }
    return lc + lb1 - lb2, OK;
}

fn betabinom_pmf(k: f64, n: f64, a: f64, b: f64) -> f64, i64 {
    let nc = disc_whole(n);
    if nc != OK { return nan(), nc; }
    let ac = check_shape(a);
    if ac != OK { return nan(), ac; }
    let bc = check_shape(b);
    if bc != OK { return nan(), bc; }
    let kc = disc_count(k);
    if kc != OK { return nan(), kc; }
    if k < 0.0 || k > n { return 0.0, OK; }
    let l, lc = betabinom_logpmf(k, n, a, b);
    if lc != OK { return nan(), lc; }
    return exp(l), OK;
}

fn betabinom_cdf(k: f64, n: f64, a: f64, b: f64) -> f64, i64 {
    let nc = disc_whole(n);
    if nc != OK { return nan(), nc; }
    let ac = check_shape(a);
    if ac != OK { return nan(), ac; }
    let bc = check_shape(b);
    if bc != OK { return nan(), bc; }
    let kc = disc_count(k);
    if kc != OK { return nan(), kc; }
    if k < 0.0 { return 0.0, OK; }
    if k >= n { return 1.0, OK; }
    // The sum of the masses runs from 0 to the count of the law.
    var total: f64 = 0.0;
    var i: f64 = 0.0;
    while i <= k {
        let l, lc = betabinom_logpmf(i, n, a, b);
        if lc != OK { return nan(), lc; }
        total += exp(l);
        i += 1.0;
    }
    return total, OK;
}

fn betabinom_ppf(q: f64, n: f64, a: f64, b: f64) -> f64, i64 {
    let nc = disc_whole(n);
    if nc != OK { return nan(), nc; }
    let ac = check_shape(a);
    if ac != OK { return nan(), ac; }
    let bc = check_shape(b);
    if bc != OK { return nan(), bc; }
    let qc = check_prob(q);
    if qc != OK { return nan(), qc; }
    if q == 0.0 { return 0.0, OK; }
    if q == 1.0 { return n, OK; }
    return disc_find(q, 0.0, n) { x => betabinom_cdf(x, n, a, b) };
}

// ---- betanbinom: the beta negative binomial law -------------------------

fn betanbinom_logpmf(k: f64, n: f64, a: f64, b: f64) -> f64, i64 {
    let lc, lcc = lchoose(n + k - 1.0, k);
    if lcc != OK { return nan(), lcc; }
    let lb1, lcb1 = lbeta(a + n, b + k);
    if lcb1 != OK { return nan(), lcb1; }
    let lb2, lcb2 = lbeta(a, b);
    if lcb2 != OK { return nan(), lcb2; }
    return lc + lb1 - lb2, OK;
}

fn betanbinom_pmf(k: f64, n: f64, a: f64, b: f64) -> f64, i64 {
    let nc = disc_whole(n);
    if nc != OK { return nan(), nc; }
    let ac = check_shape(a);
    if ac != OK { return nan(), ac; }
    let bc = check_shape(b);
    if bc != OK { return nan(), bc; }
    let kc = disc_count(k);
    if kc != OK { return nan(), kc; }
    if k < 0.0 { return 0.0, OK; }
    let l, lc = betanbinom_logpmf(k, n, a, b);
    if lc != OK { return nan(), lc; }
    return exp(l), OK;
}

fn betanbinom_cdf(k: f64, n: f64, a: f64, b: f64) -> f64, i64 {
    let nc = disc_whole(n);
    if nc != OK { return nan(), nc; }
    let ac = check_shape(a);
    if ac != OK { return nan(), ac; }
    let bc = check_shape(b);
    if bc != OK { return nan(), bc; }
    let kc = disc_count(k);
    if kc != OK { return nan(), kc; }
    if k < 0.0 { return 0.0, OK; }
    // The law has no end, so the sum runs from 0 to the count of the law.
    var total: f64 = 0.0;
    var i: f64 = 0.0;
    while i <= k {
        let l, lc = betanbinom_logpmf(i, n, a, b);
        if lc != OK { return nan(), lc; }
        total += exp(l);
        i += 1.0;
    }
    if total > 1.0 { total = 1.0; }
    return total, OK;
}

fn betanbinom_ppf(q: f64, n: f64, a: f64, b: f64) -> f64, i64 {
    let nc = disc_whole(n);
    if nc != OK { return nan(), nc; }
    let ac = check_shape(a);
    if ac != OK { return nan(), ac; }
    let bc = check_shape(b);
    if bc != OK { return nan(), bc; }
    let qc = check_prob(q);
    if qc != OK { return nan(), qc; }
    if q == 0.0 { return 0.0, OK; }
    return disc_find_open(q, 0.0) { x => betanbinom_cdf(x, n, a, b) };
}

// ---- boltzmann: the truncated discrete exponential law -------------------

fn boltzmann_pmf(k: f64, lambda: f64, n: f64) -> f64, i64 {
    let lc = check_shape(lambda);
    if lc != OK { return nan(), lc; }
    let nc = disc_whole(n);
    if nc != OK { return nan(), nc; }
    if n < 1.0 { return nan(), ERR_DOMAIN; }
    let kc = disc_count(k);
    if kc != OK { return nan(), kc; }
    if k < 0.0 || k >= n { return 0.0, OK; }
    let top = 1.0 - exp(-lambda * n);
    if top <= 0.0 { return nan(), ERR_DOMAIN; }
    return (1.0 - exp(-lambda)) * exp(-lambda * k) / top, OK;
}

fn boltzmann_cdf(k: f64, lambda: f64, n: f64) -> f64, i64 {
    let lc = check_shape(lambda);
    if lc != OK { return nan(), lc; }
    let nc = disc_whole(n);
    if nc != OK { return nan(), nc; }
    if n < 1.0 { return nan(), ERR_DOMAIN; }
    let kc = disc_count(k);
    if kc != OK { return nan(), kc; }
    if k < 0.0 { return 0.0, OK; }
    if k >= n - 1.0 { return 1.0, OK; }
    let top = 1.0 - exp(-lambda * n);
    if top <= 0.0 { return nan(), ERR_DOMAIN; }
    return (1.0 - exp(-lambda * (floor(k) + 1.0))) / top, OK;
}

fn boltzmann_ppf(q: f64, lambda: f64, n: f64) -> f64, i64 {
    let lc = check_shape(lambda);
    if lc != OK { return nan(), lc; }
    let nc = disc_whole(n);
    if nc != OK { return nan(), nc; }
    if n < 1.0 { return nan(), ERR_DOMAIN; }
    let qc = check_prob(q);
    if qc != OK { return nan(), qc; }
    if q == 0.0 { return 0.0, OK; }
    if q == 1.0 { return n - 1.0, OK; }
    let top = 1.0 - exp(-lambda * n);
    if top <= 0.0 { return nan(), ERR_DOMAIN; }
    // The count grows while the cumulative function stays under the
    // probability.
    var i: f64 = 0.0;
    while i < n - 1.0 {
        let c = (1.0 - exp(-lambda * (i + 1.0))) / top;
        if c >= q { return i, OK; }
        i += 1.0;
    }
    return n - 1.0, OK;
}

// ---- dlaplace: the discrete Laplacian law --------------------------------

fn dlaplace_pmf(k: f64, a: f64) -> f64, i64 {
    let ac = check_shape(a);
    if ac != OK { return nan(), ac; }
    let kc = disc_count(k);
    if kc != OK { return nan(), kc; }
    return (1.0 - exp(-a)) / (1.0 + exp(-a)) * exp(-a * abs(k)), OK;
}

fn dlaplace_cdf(k: f64, a: f64) -> f64, i64 {
    let ac = check_shape(a);
    if ac != OK { return nan(), ac; }
    let kc = disc_count(k);
    if kc != OK { return nan(), kc; }
    let e = exp(a);
    let f = floor(k);
    if f >= 0.0 { return 1.0 - exp(-a * f) / (e + 1.0), OK; }
    return exp(a * (f + 1.0)) / (e + 1.0), OK;
}

fn dlaplace_ppf(q: f64, a: f64) -> f64, i64 {
    let ac = check_shape(a);
    if ac != OK { return nan(), ac; }
    let qc = check_prob(q);
    if qc != OK { return nan(), qc; }
    if q == 0.0 { return -inf(), OK; }
    if q == 1.0 { return inf(), OK; }
    // The law has no end on the lower side, so the search walks up from
    // the count of the law.
    var l: f64 = -1.0;
    var r: f64 = 1.0;
    while dlaplace_cdf(l, a) >= q { l = l * 2.0; }
    while dlaplace_cdf(r, a) < q { r = r * 2.0; }
    return disc_find(q, l, r) { x => dlaplace_cdf(x, a) };
}

// ---- logser: the logarithmic series law ----------------------------------

fn logser_pmf(k: f64, p: f64) -> f64, i64 {
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    if p <= 0.0 || p >= 1.0 { return nan(), ERR_DOMAIN; }
    let kc = disc_count(k);
    if kc != OK { return nan(), kc; }
    if k < 1.0 { return 0.0, OK; }
    return -pow(p, k) / (k * log1p(-p)), OK;
}

fn logser_cdf(k: f64, p: f64) -> f64, i64 {
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    if p <= 0.0 || p >= 1.0 { return nan(), ERR_DOMAIN; }
    let kc = disc_count(k);
    if kc != OK { return nan(), kc; }
    if k < 1.0 { return 0.0, OK; }
    // The sum of the masses runs from 1 to the count of the law.
    var total: f64 = 0.0;
    var i: f64 = 1.0;
    while i <= k {
        total += -pow(p, i) / (i * log1p(-p));
        i += 1.0;
    }
    if total > 1.0 { total = 1.0; }
    return total, OK;
}

fn logser_ppf(q: f64, p: f64) -> f64, i64 {
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    if p <= 0.0 || p >= 1.0 { return nan(), ERR_DOMAIN; }
    let qc = check_prob(q);
    if qc != OK { return nan(), qc; }
    if q == 0.0 { return 1.0, OK; }
    return disc_find_open(q, 1.0) { x => logser_cdf(x, p) };
}

// ---- nhypergeom: the negative hypergeometric law -------------------------

// The mass of the negative hypergeometric law in the logarithm.
fn nhypergeom_logpmf(k: f64, m: f64, n: f64, r: f64) -> f64, i64 {
    let l1, c1 = lbeta(k + 1.0, r);
    if c1 != OK { return nan(), c1; }
    let l2, c2 = lbeta(k + r, 1.0);
    if c2 != OK { return nan(), c2; }
    let l3, c3 = lbeta(n - k + 1.0, m - r - n + 1.0);
    if c3 != OK { return nan(), c3; }
    let l4, c4 = lbeta(m - r - k + 1.0, 1.0);
    if c4 != OK { return nan(), c4; }
    let l5, c5 = lbeta(n + 1.0, m - n + 1.0);
    if c5 != OK { return nan(), c5; }
    let l6, c6 = lbeta(m + 1.0, 1.0);
    if c6 != OK { return nan(), c6; }
    return -l1 + l2 - l3 + l4 + l5 - l6, OK;
}

fn nhypergeom_pmf(k: f64, m: f64, n: f64, r: f64) -> f64, i64 {
    let mc = disc_whole(m);
    if mc != OK { return nan(), mc; }
    let nc = disc_whole(n);
    if nc != OK { return nan(), nc; }
    let rc = disc_whole(r);
    if rc != OK { return nan(), rc; }
    let kc = disc_count(k);
    if kc != OK { return nan(), kc; }
    if m < 1.0 || n < 0.0 || r < 0.0 { return nan(), ERR_DOMAIN; }
    if n > m { return nan(), ERR_DOMAIN; }
    if r > n { return nan(), ERR_DOMAIN; }
    if k < 0.0 || k > n { return 0.0, OK; }
    if r == 0.0 {
        if k == 0.0 { return 1.0, OK; }
        return 0.0, OK;
    }
    let l, lc = nhypergeom_logpmf(k, m, n, r);
    if lc != OK { return nan(), lc; }
    return exp(l), OK;
}

fn nhypergeom_cdf(k: f64, m: f64, n: f64, r: f64) -> f64, i64 {
    let mc = disc_whole(m);
    if mc != OK { return nan(), mc; }
    let nc = disc_whole(n);
    if nc != OK { return nan(), nc; }
    let rc = disc_whole(r);
    if rc != OK { return nan(), rc; }
    let kc = disc_count(k);
    if kc != OK { return nan(), kc; }
    if m < 1.0 || n < 0.0 || r < 0.0 { return nan(), ERR_DOMAIN; }
    if n > m { return nan(), ERR_DOMAIN; }
    if r > n { return nan(), ERR_DOMAIN; }
    if k < 0.0 { return 0.0, OK; }
    if k >= n { return 1.0, OK; }
    var total: f64 = 0.0;
    var i: f64 = 0.0;
    while i <= k {
        let l, lc = nhypergeom_logpmf(i, m, n, r);
        if lc != OK { return nan(), lc; }
        total += exp(l);
        i += 1.0;
    }
    if total > 1.0 { total = 1.0; }
    return total, OK;
}

fn nhypergeom_ppf(q: f64, m: f64, n: f64, r: f64) -> f64, i64 {
    let mc = disc_whole(m);
    if mc != OK { return nan(), mc; }
    let nc = disc_whole(n);
    if nc != OK { return nan(), nc; }
    let rc = disc_whole(r);
    if rc != OK { return nan(), rc; }
    let qc = check_prob(q);
    if qc != OK { return nan(), qc; }
    if m < 1.0 || n < 0.0 || r < 0.0 { return nan(), ERR_DOMAIN; }
    if n > m { return nan(), ERR_DOMAIN; }
    if r > n { return nan(), ERR_DOMAIN; }
    if q == 0.0 { return 0.0, OK; }
    if q == 1.0 { return n, OK; }
    return disc_find(q, 0.0, n) { x => nhypergeom_cdf(x, m, n, r) };
}
