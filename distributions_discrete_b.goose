// Seven discrete families: planck, poisson_binom, randint, skellam,
// yulesimon, zipf and zipfian. The names and the parameters follow SciPy.
//
// Every family has three calls: pmf, cdf and ppf. The mass of a value
// outside the support is 0 and the code is OK. The quantile of 0 or 1
// gives the bound of the support. Every call returns the result first and
// a fault code last.
//
// The Poisson binomial law takes a list of probabilities. The other laws
// take a shape or two.
namespace statistics;

import std;
import math;

// The terms of the Euler transform of the alternating series of the
// Riemann zeta function.
const ZIPF_TERMS: i64 = 50;

// ---- planck: the Planck discrete exponential law -------------------------

fn planck_pmf(k: f64, lambda: f64) -> f64, i64 {
    let lc = check_shape(lambda);
    if lc != OK { return nan(), lc; }
    let kc = disc_count(k);
    if kc != OK { return nan(), kc; }
    if k < 0.0 { return 0.0, OK; }
    let m, cm = expm1(-lambda);
    if cm != OK { return nan(), cm; }
    return -m * exp(-lambda * k), OK;
}

fn planck_cdf(k: f64, lambda: f64) -> f64, i64 {
    let lc = check_shape(lambda);
    if lc != OK { return nan(), lc; }
    let kc = disc_count(k);
    if kc != OK { return nan(), kc; }
    if k < 0.0 { return 0.0, OK; }
    let m, cm = expm1(-lambda * (floor(k) + 1.0));
    if cm != OK { return nan(), cm; }
    return -m, OK;
}

fn planck_ppf(q: f64, lambda: f64) -> f64, i64 {
    let lc = check_shape(lambda);
    if lc != OK { return nan(), lc; }
    let qc = check_prob(q);
    if qc != OK { return nan(), qc; }
    if q == 0.0 { return 0.0, OK; }
    return disc_find_open(q, 0.0) { x => planck_cdf(x, lambda) };
}

// ---- poisson_binom: the Poisson binomial law -----------------------------

// The masses of the Poisson binomial law. The list holds the masses of
// the counts 0 to the number of probabilities. The rule adds one
// probability at a time. The mass of a count of the new list is the mass
// of the count without the success plus the mass of the count below with
// the success.
fn poisson_binom_masses(ps: const f64[:]) -> f64[>..] {
    var out: f64[>..] = [];
    out.push(1.0);
    for i in 0..ps.len {
        var next: f64[>..] = [];
        let p = ps[i];
        var k: i64 = 0;
        while k < out.len {
            let old = out[k];
            next.push(old * (1.0 - p));
            k += 1;
        }
        next.push(0.0);
        k = 0;
        while k < out.len {
            let old = out[k];
            next[k + 1] = next[k + 1] + old * p;
            k += 1;
        }
        out = copy(next);
    }
    return out;
}

fn poisson_binom_pmf(k: f64, ps: const f64[:]) -> f64, i64 {
    let c = check_probabilities(ps);
    if c != OK { return nan(), c; }
    let kc = disc_count(k);
    if kc != OK { return nan(), kc; }
    if k < 0.0 || k > ps.len as f64 { return 0.0, OK; }
    let masses = poisson_binom_masses(ps);
    return masses[k as i64], OK;
}

fn poisson_binom_cdf(k: f64, ps: const f64[:]) -> f64, i64 {
    let c = check_probabilities(ps);
    if c != OK { return nan(), c; }
    let kc = disc_count(k);
    if kc != OK { return nan(), kc; }
    if k < 0.0 { return 0.0, OK; }
    if k >= ps.len as f64 { return 1.0, OK; }
    let masses = poisson_binom_masses(ps);
    var total: f64 = 0.0;
    var i: i64 = 0;
    while i <= k as i64 {
        total += masses[i];
        i += 1;
    }
    return total, OK;
}

fn poisson_binom_ppf(q: f64, ps: const f64[:]) -> f64, i64 {
    let c = check_probabilities(ps);
    if c != OK { return nan(), c; }
    let qc = check_prob(q);
    if qc != OK { return nan(), qc; }
    if q == 0.0 { return 0.0, OK; }
    if q == 1.0 { return ps.len as f64, OK; }
    let masses = poisson_binom_masses(ps);
    var total: f64 = 0.0;
    var i: i64 = 0;
    while i < masses.len {
        total += masses[i];
        if total >= q { return i as f64, OK; }
        i += 1;
    }
    return ps.len as f64, OK;
}

// ---- randint: the law of a whole number between two bounds --------------

fn randint_pmf(k: f64, low: f64, high: f64) -> f64, i64 {
    let lc = disc_whole(low);
    if lc != OK { return nan(), lc; }
    let hc = disc_whole(high);
    if hc != OK { return nan(), hc; }
    let kc = disc_count(k);
    if kc != OK { return nan(), kc; }
    if high <= low { return nan(), ERR_DOMAIN; }
    if k < low || k >= high { return 0.0, OK; }
    return 1.0 / (high - low), OK;
}

fn randint_cdf(k: f64, low: f64, high: f64) -> f64, i64 {
    let lc = disc_whole(low);
    if lc != OK { return nan(), lc; }
    let hc = disc_whole(high);
    if hc != OK { return nan(), hc; }
    let kc = disc_count(k);
    if kc != OK { return nan(), kc; }
    if high <= low { return nan(), ERR_DOMAIN; }
    if k < low { return 0.0, OK; }
    if k >= high - 1.0 { return 1.0, OK; }
    return (floor(k) - low + 1.0) / (high - low), OK;
}

fn randint_ppf(q: f64, low: f64, high: f64) -> f64, i64 {
    let lc = disc_whole(low);
    if lc != OK { return nan(), lc; }
    let hc = disc_whole(high);
    if hc != OK { return nan(), hc; }
    let qc = check_prob(q);
    if qc != OK { return nan(), qc; }
    if high <= low { return nan(), ERR_DOMAIN; }
    if q == 0.0 { return low, OK; }
    if q == 1.0 { return high - 1.0, OK; }
    // The count grows while the cumulative function stays under the
    // probability.
    var i: f64 = low;
    while i < high {
        let c = (i - low + 1.0) / (high - low);
        if c >= q { return i, OK; }
        i += 1.0;
    }
    return high - 1.0, OK;
}

// ---- skellam: the difference of two Poisson laws -------------------------

// The mass of the Skellam law. The law needs the modified Bessel
// function of the first kind. The rule computes it from its series.
fn skellam_bessel_i(n: f64, z: f64) -> f64, i64 {
    // The series of the function. Its terms are all of one sign, so the
    // sum keeps its digits.
    let half = z / 2.0;
    var term = pow(half, n) / gamma(n + 1.0);
    var total = term;
    var k: i64 = 1;
    var steps: i64 = 0;
    while steps < MAX_ITER_SERIES {
        term = term * (half * half) / (k as f64 * (k + n));
        total += term;
        if abs(term) < 1e-18 * abs(total) { break; }
        k += 1;
        steps += 1;
    }
    return total, OK;
}

fn skellam_pmf(k: f64, mu1: f64, mu2: f64) -> f64, i64 {
    let c1 = check_shape(mu1);
    if c1 != OK { return nan(), c1; }
    let c2 = check_shape(mu2);
    if c2 != OK { return nan(), c2; }
    let kc = disc_count(k);
    if kc != OK { return nan(), kc; }
    let n = abs(k);
    let b, bc = skellam_bessel_i(n, 2.0 * sqrt(mu1 * mu2));
    if bc != OK { return nan(), bc; }
    let half = k / 2.0;
    let r = exp(-(mu1 + mu2)) * b;
    // The power of the ratio needs the form with the logarithm when the
    // ratio is far from one.
    let l = log(mu1) - log(mu2);
    let p = exp(half * l);
    return r * p, OK;
}

fn skellam_cdf(k: f64, mu1: f64, mu2: f64) -> f64, i64 {
    let c1 = check_shape(mu1);
    if c1 != OK { return nan(), c1; }
    let c2 = check_shape(mu2);
    if c2 != OK { return nan(), c2; }
    let kc = disc_count(k);
    if kc != OK { return nan(), kc; }
    // The law has no end, so the cumulative function comes from the sum
    // of the masses. The sum runs from the lower bound of the law up to
    // the count.
    var total: f64 = 0.0;
    var i: f64 = floor(k);
    // The lower bound of the sum. The mass there falls below the
    // precision of a double.
    var lo: f64 = floor(k) - 1.0;
    var steps: i64 = 0;
    while steps < 200 {
        let m, mc = skellam_pmf(lo, mu1, mu2);
        if mc != OK { return nan(), mc; }
        if m < 1e-18 { break; }
        lo -= 1.0;
        steps += 1;
    }
    lo = floor(lo);
    while lo <= i {
        let m, mc = skellam_pmf(lo, mu1, mu2);
        if mc != OK { return nan(), mc; }
        total += m;
        lo += 1.0;
    }
    if total > 1.0 { total = 1.0; }
    return total, OK;
}

fn skellam_ppf(q: f64, mu1: f64, mu2: f64) -> f64, i64 {
    let c1 = check_shape(mu1);
    if c1 != OK { return nan(), c1; }
    let c2 = check_shape(mu2);
    if c2 != OK { return nan(), c2; }
    let qc = check_prob(q);
    if qc != OK { return nan(), qc; }
    if q == 0.0 { return -inf(), OK; }
    if q == 1.0 { return inf(), OK; }
    // The law has no end on the lower side, so the search walks up from
    // the count of the law.
    var l: f64 = -1.0;
    var r: f64 = 1.0;
    while skellam_cdf(l, mu1, mu2) >= q { l = l * 2.0; }
    while skellam_cdf(r, mu1, mu2) < q { r = r * 2.0; }
    return disc_find(q, l, r) { x => skellam_cdf(x, mu1, mu2) };
}

// ---- yulesimon: the Yule-Simon law ---------------------------------------

fn yulesimon_pmf(k: f64, alpha: f64) -> f64, i64 {
    let ac = check_shape(alpha);
    if ac != OK { return nan(), ac; }
    let kc = disc_count(k);
    if kc != OK { return nan(), kc; }
    if k < 1.0 { return 0.0, OK; }
    let b, bc = beta(k, alpha + 1.0);
    if bc != OK { return nan(), bc; }
    return alpha * b, OK;
}

fn yulesimon_cdf(k: f64, alpha: f64) -> f64, i64 {
    let ac = check_shape(alpha);
    if ac != OK { return nan(), ac; }
    let kc = disc_count(k);
    if kc != OK { return nan(), kc; }
    if k < 1.0 { return 0.0, OK; }
    let b, bc = beta(k, alpha + 1.0);
    if bc != OK { return nan(), bc; }
    let v = 1.0 - k * b;
    if v < 0.0 { return 0.0, OK; }
    return v, OK;
}

fn yulesimon_ppf(q: f64, alpha: f64) -> f64, i64 {
    let ac = check_shape(alpha);
    if ac != OK { return nan(), ac; }
    let qc = check_prob(q);
    if qc != OK { return nan(), qc; }
    if q == 0.0 { return 1.0, OK; }
    return disc_find_open(q, 1.0) { x => yulesimon_cdf(x, alpha) };
}

// ---- zipf: the Zipf law of the infinite kind -----------------------------

// The Riemann zeta function. The shape must be above one.
//
// The function comes from the alternating series of Dirichlet. The
// Euler transform of that series converges much faster than the series
// itself. The transform takes the differences of the terms of the
// series and halves the weight of each step.
fn zipf_zeta(a: f64) -> f64, i64 {
    if a <= 1.0 { return nan(), ERR_DOMAIN; }
    // The terms of the alternating series.
    var d: f64[>..] = [];
    for k in 0..ZIPF_TERMS + 1 {
        d.push(pow(k + 1.0, -a));
    }
    // The Euler transform. Each step takes the differences of the terms
    // and halves the weight.
    var total: f64 = 0.0;
    var weight: f64 = 0.5;
    for n in 0..ZIPF_TERMS + 1 {
        total += d[0] * weight;
        weight = weight * 0.5;
        var next: f64[>..] = [];
        for i in 0..d.len - 1 {
            next.push(d[i] - d[i + 1]);
        }
        d = copy(next);
        if d.len == 0 { break; }
    }
    // The alternating series and the function differ by the factor
    // 1 - 2^(1 - a).
    let den = 1.0 - pow(2.0, 1.0 - a);
    if den == 0.0 { return nan(), ERR_DOMAIN; }
    return total / den, OK;
}

fn zipf_pmf(k: f64, a: f64) -> f64, i64 {
    let ac = check_shape(a);
    if ac != OK { return nan(), ac; }
    if a <= 1.0 { return nan(), ERR_DOMAIN; }
    let kc = disc_count(k);
    if kc != OK { return nan(), kc; }
    if k < 1.0 { return 0.0, OK; }
    let z, zc = zipf_zeta(a);
    if zc != OK { return nan(), zc; }
    return 1.0 / (z * pow(k, a)), OK;
}

fn zipf_cdf(k: f64, a: f64) -> f64, i64 {
    let ac = check_shape(a);
    if ac != OK { return nan(), ac; }
    if a <= 1.0 { return nan(), ERR_DOMAIN; }
    let kc = disc_count(k);
    if kc != OK { return nan(), kc; }
    if k < 1.0 { return 0.0, OK; }
    let z, zc = zipf_zeta(a);
    if zc != OK { return nan(), zc; }
    // The sum of the masses runs from 1 to the count of the law.
    var total: f64 = 0.0;
    var i: f64 = 1.0;
    while i <= k {
        total += 1.0 / (z * pow(i, a));
        i += 1.0;
    }
    if total > 1.0 { total = 1.0; }
    return total, OK;
}

fn zipf_ppf(q: f64, a: f64) -> f64, i64 {
    let ac = check_shape(a);
    if ac != OK { return nan(), ac; }
    if a <= 1.0 { return nan(), ERR_DOMAIN; }
    let qc = check_prob(q);
    if qc != OK { return nan(), qc; }
    if q == 0.0 { return 1.0, OK; }
    return disc_find_open(q, 1.0) { x => zipf_cdf(x, a) };
}

// ---- zipfian: the Zipf law of the finite kind ----------------------------

// The generalized harmonic number of the order a up to the count n.
fn zipfian_harmonic(n: f64, a: f64) -> f64 {
    var total: f64 = 0.0;
    var k: f64 = 1.0;
    while k <= n {
        total += pow(k, -a);
        k += 1.0;
    }
    return total;
}

fn zipfian_pmf(k: f64, a: f64, n: f64) -> f64, i64 {
    let ac = check_shape(a);
    if ac != OK { return nan(), ac; }
    let nc = disc_whole(n);
    if nc != OK { return nan(), nc; }
    if n < 1.0 { return nan(), ERR_DOMAIN; }
    let kc = disc_count(k);
    if kc != OK { return nan(), kc; }
    if k < 1.0 || k > n { return 0.0, OK; }
    let h = zipfian_harmonic(n, a);
    if h <= 0.0 { return nan(), ERR_DOMAIN; }
    return 1.0 / (h * pow(k, a)), OK;
}

fn zipfian_cdf(k: f64, a: f64, n: f64) -> f64, i64 {
    let ac = check_shape(a);
    if ac != OK { return nan(), ac; }
    let nc = disc_whole(n);
    if nc != OK { return nan(), nc; }
    if n < 1.0 { return nan(), ERR_DOMAIN; }
    let kc = disc_count(k);
    if kc != OK { return nan(), kc; }
    if k < 1.0 { return 0.0, OK; }
    if k >= n { return 1.0, OK; }
    let h = zipfian_harmonic(n, a);
    if h <= 0.0 { return nan(), ERR_DOMAIN; }
    return zipfian_harmonic(floor(k), a) / h, OK;
}

fn zipfian_ppf(q: f64, a: f64, n: f64) -> f64, i64 {
    let ac = check_shape(a);
    if ac != OK { return nan(), ac; }
    let nc = disc_whole(n);
    if nc != OK { return nan(), nc; }
    if n < 1.0 { return nan(), ERR_DOMAIN; }
    let qc = check_prob(q);
    if qc != OK { return nan(), qc; }
    if q == 0.0 { return 1.0, OK; }
    if q == 1.0 { return n, OK; }
    let h = zipfian_harmonic(n, a);
    if h <= 0.0 { return nan(), ERR_DOMAIN; }
    // The count grows while the cumulative function stays under the
    // probability.
    var i: f64 = 1.0;
    var total: f64 = 0.0;
    while i <= n {
        total += pow(i, -a) / h;
        if total >= q { return i, OK; }
        i += 1.0;
    }
    return n, OK;
}
