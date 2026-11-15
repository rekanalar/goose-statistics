// The two non-central hypergeometric laws of Fisher and of Wallenius.
// The names and the parameters follow SciPy 1.18.1: the total m, the
// successes n, the draws d and the odds ratio.
//
// The mass of Fisher is a weighted series over the counts. The mass of
// Wallenius is an integral over the unit interval. Every family has
// three calls: pmf, cdf and ppf. Every call returns the result first and
// a fault code last. The quantile of 0 gives the smallest count of the
// support, as the other discrete families of the package do.
namespace statistics;

import std;
import math;

const NCH_MAX_M: f64 = 100000.0;      // the largest total of the law of Fisher
const NCH_W_MAX_M: f64 = 500.0;       // the largest total of the law of Wallenius

// The check of one count. It rejects a NaN, a fraction, a negative value
// and a value above the bound.
fn nch_count(x: f64, top: f64) -> i64 {
    let c = check_scalar(x);
    if c != OK { return c; }
    if x < 0.0 || x != floor(x) || x > top { return ERR_DOMAIN; }
    return OK;
}

fn nch_check(m: f64, n: f64, d: f64, odds: f64, top: f64) -> i64 {
    let c1 = nch_count(m, top);
    if c1 != OK { return c1; }
    let c2 = nch_count(n, top);
    if c2 != OK { return c2; }
    let c3 = nch_count(d, top);
    if c3 != OK { return c3; }
    let c4 = check_scalar(odds);
    if c4 != OK { return c4; }
    if n > m || d > m || odds <= 0.0 { return ERR_DOMAIN; }
    return OK;
}

// The smallest and the largest count of the support.
fn nch_lo(m: f64, n: f64, d: f64) -> f64 { return max(0.0, d - (m - n)); }
fn nch_hi(n: f64, d: f64) -> f64 { return min(n, d); }

// The log of the number of ways to draw k successes: C(n, k) C(m - n, d - k).
fn nch_lways(m: f64, n: f64, d: f64, k: f64) -> f64, i64 {
    let a, ca = lchoose(n, k);
    if ca != OK { return nan(), ca; }
    let b, cb = lchoose(m - n, d - k);
    if cb != OK { return nan(), cb; }
    return a + b, OK;
}

// ---- nchypergeom_fisher ---------------------------------------------------

// The masses of the law, from the smallest count to the largest, in the
// array. The masses sum to 1. The call gives the smallest count.
fn fisher_table(m: f64, n: f64, d: f64, odds: f64, pm: f64[>..]&) -> f64, i64 {
    pm.resize(0);
    let lo = nch_lo(m, n, d);
    let hi = nch_hi(n, d);
    var top = -inf();
    var k = lo;
    while k <= hi {
        let w, cw = nch_lways(m, n, d, k);
        if cw != OK { return nan(), cw; }
        let l = w + k * log(odds);
        pm.push(l);
        if l > top { top = l; }
        k += 1.0;
    }
    var total: f64 = 0.0;
    for i in 0..pm.len {
        pm[i] = exp(pm[i] - top);
        total += pm[i];
    }
    for i in 0..pm.len { pm[i] = pm[i] / total; }
    return lo, OK;
}

fn nchypergeom_fisher_pmf(k: f64, m: f64, n: f64, d: f64, odds: f64) -> f64, i64 {
    let c = nch_check(m, n, d, odds, NCH_MAX_M);
    if c != OK { return nan(), c; }
    let kc = check_scalar(k);
    if kc != OK { return nan(), kc; }
    if k != floor(k) { return 0.0, OK; }
    let lo = nch_lo(m, n, d);
    if k < lo || k > nch_hi(n, d) { return 0.0, OK; }
    var pm: f64[>..] = [];
    let l0, cl = fisher_table(m, n, d, odds, pm);
    if cl != OK { return nan(), cl; }
    return pm[(k - lo) as i64], OK;
}

fn nchypergeom_fisher_cdf(k: f64, m: f64, n: f64, d: f64, odds: f64) -> f64, i64 {
    let c = nch_check(m, n, d, odds, NCH_MAX_M);
    if c != OK { return nan(), c; }
    let kc = check_scalar(k);
    if kc != OK { return nan(), kc; }
    let lo = nch_lo(m, n, d);
    let hi = nch_hi(n, d);
    if k < lo { return 0.0, OK; }
    if k >= hi { return 1.0, OK; }
    var pm: f64[>..] = [];
    let l0, cl = fisher_table(m, n, d, odds, pm);
    if cl != OK { return nan(), cl; }
    let last = floor(k) - lo;
    var s: f64 = 0.0;
    var i: i64 = 0;
    while i <= last {
        s += pm[i];
        i += 1;
    }
    return min(s, 1.0), OK;
}

fn nchypergeom_fisher_ppf(q: f64, m: f64, n: f64, d: f64, odds: f64) -> f64, i64 {
    let c = nch_check(m, n, d, odds, NCH_MAX_M);
    if c != OK { return nan(), c; }
    let qc = check_prob(q);
    if qc != OK { return nan(), qc; }
    let lo = nch_lo(m, n, d);
    let hi = nch_hi(n, d);
    if q == 0.0 { return lo, OK; }
    if q == 1.0 { return hi, OK; }
    var pm: f64[>..] = [];
    let l0, cl = fisher_table(m, n, d, odds, pm);
    if cl != OK { return nan(), cl; }
    var s: f64 = 0.0;
    for i in 0..pm.len {
        s += pm[i];
        if s >= q { return lo + i, OK; }
    }
    return hi, OK;
}

// ---- nchypergeom_wallenius ------------------------------------------------

// The factor (1 - u^a)^e of the integrand. An exponent of 0 gives 1.
fn wal_factor(u: f64, a: f64, e: f64) -> f64 {
    if e == 0.0 { return 1.0; }
    let p = 1.0 - pow(u, a);
    if p <= 0.0 { return 0.0; }
    return exp(e * log(p));
}

// The mass of the count k, without the checks. The integral over t is
// taken in u = t^(1 / big), which turns the steep end at 0 into a smooth
// one: the integrand is big u^(big - 1) (1 - u^odds)^k (1 - u)^(d - k).
fn wallenius_mass(k: f64, m: f64, n: f64, d: f64, odds: f64) -> f64, i64 {
    let lo = nch_lo(m, n, d);
    let hi = nch_hi(n, d);
    if k < lo || k > hi { return 0.0, OK; }
    // One count only: the mass is 1.
    if lo == hi { return 1.0, OK; }
    let w, cw = nch_lways(m, n, d, k);
    if cw != OK { return nan(), cw; }
    let big = odds * (n - k) + (m - n - (d - k));
    // The integrand lives near u = 1, where the factor u^(big - 1)
    // grows. The panels there halve, so the last one is tiny.
    var total: f64 = 0.0;
    var edge: f64 = 0.0;
    var gap: f64 = 0.5;
    var panels: i64 = 0;
    while edge < 1.0 && panels < 34 {
        let hi_e = min(edge + gap, 1.0);
        total += fam_quad(edge, hi_e) { u => big * exp((big - 1.0) * log(max(u, 1e-300))) * wal_factor(u, odds, k) * wal_factor(u, 1.0, d - k) };
        edge = hi_e;
        // The next panel is half as wide, so the panels meet at 1.
        gap = 0.5 * (1.0 - edge);
        if gap < 1e-12 { break; }
        panels += 1;
    }
    return exp(w) * total, OK;
}

fn nchypergeom_wallenius_pmf(k: f64, m: f64, n: f64, d: f64, odds: f64) -> f64, i64 {
    let c = nch_check(m, n, d, odds, NCH_W_MAX_M);
    if c != OK { return nan(), c; }
    let kc = check_scalar(k);
    if kc != OK { return nan(), kc; }
    if k != floor(k) { return 0.0, OK; }
    let v, cv = wallenius_mass(k, m, n, d, odds);
    if cv != OK { return nan(), cv; }
    return min(max(v, 0.0), 1.0), OK;
}

fn nchypergeom_wallenius_cdf(k: f64, m: f64, n: f64, d: f64, odds: f64) -> f64, i64 {
    let c = nch_check(m, n, d, odds, NCH_W_MAX_M);
    if c != OK { return nan(), c; }
    let kc = check_scalar(k);
    if kc != OK { return nan(), kc; }
    let lo = nch_lo(m, n, d);
    let hi = nch_hi(n, d);
    if k < lo { return 0.0, OK; }
    if k >= hi { return 1.0, OK; }
    var s: f64 = 0.0;
    var j = lo;
    while j <= floor(k) {
        let v, cv = wallenius_mass(j, m, n, d, odds);
        if cv != OK { return nan(), cv; }
        s += v;
        j += 1.0;
    }
    return min(s, 1.0), OK;
}

fn nchypergeom_wallenius_ppf(q: f64, m: f64, n: f64, d: f64, odds: f64) -> f64, i64 {
    let c = nch_check(m, n, d, odds, NCH_W_MAX_M);
    if c != OK { return nan(), c; }
    let qc = check_prob(q);
    if qc != OK { return nan(), qc; }
    let lo = nch_lo(m, n, d);
    let hi = nch_hi(n, d);
    if q == 0.0 { return lo, OK; }
    if q == 1.0 { return hi, OK; }
    var s: f64 = 0.0;
    var j = lo;
    while j <= hi {
        let v, cv = wallenius_mass(j, m, n, d, odds);
        if cv != OK { return nan(), cv; }
        s += v;
        if s >= q { return j, OK; }
        j += 1.0;
    }
    return hi, OK;
}
