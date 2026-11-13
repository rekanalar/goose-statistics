// The discrete distributions of the package: the binomial, the Poisson,
// the geometric, the negative binomial and the hypergeometric families.
// Each family has a mass, a cumulative function, a survival function, a
// quantile function and the four moments.
//
// Every call returns the result first and a fault code last. A fault gives
// NaN as the result. The checks run in one order: NaN, infinity, domain,
// pole.
namespace statistics;

import std;
import math;

fn chk_count(k: f64) -> i64 {
    if is_nan(k) { return ERR_NAN; }
    if is_inf(k) { return ERR_INFINITY; }
    if k < 0.0 || floor(k) != k { return ERR_DOMAIN; }
    if k > MAX_COUNT { return ERR_DOMAIN; }
    return OK;
}

// ---- the binomial family ----

fn binom_check(n: f64, p: f64) -> i64 {
    if is_nan(n) || is_nan(p) { return ERR_NAN; }
    return first_code(chk_whole(n), chk_unit(p));
}

fn binom_pmf(k: f64, n: f64, p: f64 = 0.5) -> f64, i64 {
    let ck = chk_count(k);
    if ck != OK { return nan(), ck; }
    let c = binom_check(n, p);
    if c != OK { return nan(), c; }
    if k > n { return 0.0, OK; }
    if p == 0.0 { if k == 0.0 { return 1.0, OK; } return 0.0, OK; }
    if p == 1.0 { if k == n { return 1.0, OK; } return 0.0, OK; }
    var l = dd_lchoose(n, k);
    l = dd_add(l, dd_mul_d(dd_log_d(p), k));
    l = dd_add(l, dd_mul_d(dd_log_1m(p), n - k));
    return exp_or_overflow(l);
}

fn binom_cdf(k: f64, n: f64, p: f64 = 0.5) -> f64, i64 {
    let ck = chk_count(k);
    if ck != OK { return nan(), ck; }
    let c = binom_check(n, p);
    if c != OK { return nan(), c; }
    if k >= n { return 1.0, OK; }
    if p == 0.0 { return 1.0, OK; }
    if p == 1.0 { return 0.0, OK; }
    let kf = floor(k);
    return betainc_upper(kf + 1.0, n - kf, p);
}

fn binom_sf(k: f64, n: f64, p: f64 = 0.5) -> f64, i64 {
    let ck = chk_count(k);
    if ck != OK { return nan(), ck; }
    let c = binom_check(n, p);
    if c != OK { return nan(), c; }
    if k >= n { return 0.0, OK; }
    if p == 0.0 { return 0.0, OK; }
    if p == 1.0 { return 1.0, OK; }
    let kf = floor(k);
    return betainc(kf + 1.0, n - kf, p);
}

// ---- the Poisson family ----

fn poisson_pmf(k: f64, rate: f64 = 1.0) -> f64, i64 {
    let ck = chk_count(k);
    if ck != OK { return nan(), ck; }
    let c = check_shape(rate);
    if c != OK { return nan(), c; }
    let l = dd_sub(dd_sub(dd_mul_d(dd_log_d(rate), k), dd_from(rate)), dd_lngamma(dd_from(k + 1.0)));
    return exp_or_overflow(l);
}

fn poisson_cdf(k: f64, rate: f64 = 1.0) -> f64, i64 {
    let ck = chk_count(k);
    if ck != OK { return nan(), ck; }
    let c = check_shape(rate);
    if c != OK { return nan(), c; }
    return gammaincc(floor(k) + 1.0, rate);
}

fn poisson_sf(k: f64, rate: f64 = 1.0) -> f64, i64 {
    let ck = chk_count(k);
    if ck != OK { return nan(), ck; }
    let c = check_shape(rate);
    if c != OK { return nan(), c; }
    return gammainc(floor(k) + 1.0, rate);
}

// ---- the geometric family (the trials up to the first success) ----

fn geom_pmf(k: f64, p: f64 = 0.5) -> f64, i64 {
    let ck = chk_count(k);
    if ck != OK { return nan(), ck; }
    let c = chk_succ(p);
    if c != OK { return nan(), c; }
    if k < 1.0 { return 0.0, OK; }
    if p == 1.0 { if k == 1.0 { return 1.0, OK; } return 0.0, OK; }
    let l = dd_add(dd_mul_d(dd_log_1m(p), k - 1.0), dd_log_d(p));
    return exp_or_overflow(l);
}

fn geom_cdf(k: f64, p: f64 = 0.5) -> f64, i64 {
    let ck = chk_count(k);
    if ck != OK { return nan(), ck; }
    let c = chk_succ(p);
    if c != OK { return nan(), c; }
    if k < 1.0 { return 0.0, OK; }
    if p == 1.0 { return 1.0, OK; }
    let l = dd_mul_d(dd_log_1m(p), floor(k));
    if l.hi < -0.5 { return 1.0 - expf_dd(l), OK; }
    let v, c2 = expm1(dd_to_f64(l));
    return -v, OK;
}

fn geom_sf(k: f64, p: f64 = 0.5) -> f64, i64 {
    let ck = chk_count(k);
    if ck != OK { return nan(), ck; }
    let c = chk_succ(p);
    if c != OK { return nan(), c; }
    if k < 1.0 { return 1.0, OK; }
    if p == 1.0 { return 0.0, OK; }
    return expf_dd(dd_mul_d(dd_log_1m(p), floor(k))), OK;
}

// ---- the negative binomial family (the failures before the r-th success) ----

fn nbinom_check(r: f64, p: f64) -> i64 {
    if is_nan(r) || is_nan(p) { return ERR_NAN; }
    return first_code(check_shape(r), chk_succ(p));
}

fn nbinom_pmf(k: f64, r: f64, p: f64 = 0.5) -> f64, i64 {
    let ck = chk_count(k);
    if ck != OK { return nan(), ck; }
    let c = nbinom_check(r, p);
    if c != OK { return nan(), c; }
    if p == 1.0 { if k == 0.0 { return 1.0, OK; } return 0.0, OK; }
    // Gamma(k + r) / (Gamma(r) k!) = 1 / ((k + r) B(r, k + 1))
    var l = dd_neg(dd_log(dd_add_d(dd_from(r), k)));
    l = dd_sub(l, dd_lnbeta(r, k + 1.0));
    l = dd_add(l, dd_mul_d(dd_log_d(p), r));
    l = dd_add(l, dd_mul_d(dd_log_1m(p), k));
    return exp_or_overflow(l);
}

fn nbinom_cdf(k: f64, r: f64, p: f64 = 0.5) -> f64, i64 {
    let ck = chk_count(k);
    if ck != OK { return nan(), ck; }
    let c = nbinom_check(r, p);
    if c != OK { return nan(), c; }
    if p == 1.0 { return 1.0, OK; }
    return betainc(r, floor(k) + 1.0, p);
}

fn nbinom_sf(k: f64, r: f64, p: f64 = 0.5) -> f64, i64 {
    let ck = chk_count(k);
    if ck != OK { return nan(), ck; }
    let c = nbinom_check(r, p);
    if c != OK { return nan(), c; }
    if p == 1.0 { return 0.0, OK; }
    return betainc_upper(r, floor(k) + 1.0, p);
}

// ---- the hypergeometric family (population m, n successes, d draws) ----

fn hyper_check(m: f64, n: f64, d: f64) -> i64 {
    if is_nan(m) || is_nan(n) || is_nan(d) { return ERR_NAN; }
    let c = first_code(first_code(chk_whole(m), chk_whole(n)), chk_whole(d));
    if c != OK { return c; }
    if m < 1.0 || n > m || d > m || m > MAX_SHAPE { return ERR_DOMAIN; }
    return OK;
}

// The log of the mass at a whole k inside the support, as a pair.
fn hyper_logpmf(k: f64, m: f64, n: f64, d: f64) -> DD {
    let a = dd_lchoose(n, k);
    let b = dd_lchoose(m - n, d - k);
    let t = dd_lchoose(m, d);
    return dd_sub(dd_add(a, b), t);
}

fn hyper_kmin(m: f64, n: f64, d: f64) -> f64 { return max(0.0, d - (m - n)); }
fn hyper_kmax(n: f64, d: f64) -> f64 { return min(n, d); }

fn hypergeom_pmf(k: f64, m: f64, n: f64, d: f64) -> f64, i64 {
    let ck = chk_count(k);
    if ck != OK { return nan(), ck; }
    let c = hyper_check(m, n, d);
    if c != OK { return nan(), c; }
    if k < hyper_kmin(m, n, d) || k > hyper_kmax(n, d) { return 0.0, OK; }
    return expf_dd(hyper_logpmf(k, m, n, d)), OK;
}

// The mode of the hypergeometric family, the count with the largest
// mass.
fn hyper_mode(m: f64, n: f64, d: f64) -> f64 {
    let md = floor((d + 1.0) * (n + 1.0) / (m + 2.0));
    return min(max(md, hyper_kmin(m, n, d)), hyper_kmax(n, d));
}

// The mass at the mode as a pair. It is never below the reciprocal of the
// length of the support, so the pair holds no underflow.
fn hyper_top(m: f64, n: f64, d: f64, mode: f64) -> DD {
    let m2, k2 = dd_exp_parts(hyper_logpmf(mode, m, n, d));
    return dd_scale_dd(m2, k2);
}

// The ratio of the mass at k + 1 to the mass at k as a pair.
fn hyper_up(t: DD, k: f64, m: f64, n: f64, d: f64) -> DD {
    var r = dd_mul_d(dd_mul_d(t, n - k), d - k);
    r = dd_div_d(r, k + 1.0);
    return dd_div_d(r, m - n - d + k + 1.0);
}

// The ratio of the mass at k - 1 to the mass at k as a pair.
fn hyper_down(t: DD, k: f64, m: f64, n: f64, d: f64) -> DD {
    var r = dd_mul_d(dd_mul_d(t, k), m - n - d + k);
    r = dd_div_d(r, n - k + 1.0);
    return dd_div_d(r, d - k + 1.0);
}

// The sums of the masses at and below k, and above k. Both come straight
// from the masses, so a small tail keeps its digits. The support length
// has a bound.
fn hyper_tails(k: f64, m: f64, n: f64, d: f64) -> f64, f64, i64 {
    let kmin = hyper_kmin(m, n, d);
    let kmax = hyper_kmax(n, d);
    if k < kmin { return 0.0, 1.0, OK; }
    if k >= kmax { return 1.0, 0.0, OK; }
    if kmax - kmin + 1.0 > MAX_DIST_SUPPORT { return nan(), nan(), ERR_TOO_LONG; }
    let mode = hyper_mode(m, n, d);
    var t = hyper_top(m, n, d, mode);
    var low = dd_from(0.0);
    var up = dd_from(0.0);
    var j = mode;
    var cur = t;
    while j >= kmin {
        if j <= k { low = dd_add(low, cur); } else { up = dd_add(up, cur); }
        if j == kmin || cur.hi == 0.0 { break; }
        cur = hyper_down(cur, j, m, n, d);
        j -= 1.0;
    }
    j = mode;
    cur = t;
    while j < kmax {
        cur = hyper_up(cur, j, m, n, d);
        j += 1.0;
        if cur.hi == 0.0 { break; }
        if j <= k { low = dd_add(low, cur); } else { up = dd_add(up, cur); }
    }
    return dd_to_f64(low), dd_to_f64(up), OK;
}

fn hypergeom_cdf(k: f64, m: f64, n: f64, d: f64) -> f64, i64 {
    let ck = chk_count(k);
    if ck != OK { return nan(), ck; }
    let c = hyper_check(m, n, d);
    if c != OK { return nan(), c; }
    let lo, up, c2 = hyper_tails(floor(k), m, n, d);
    return lo, c2;
}

fn hypergeom_sf(k: f64, m: f64, n: f64, d: f64) -> f64, i64 {
    let ck = chk_count(k);
    if ck != OK { return nan(), ck; }
    let c = hyper_check(m, n, d);
    if c != OK { return nan(), c; }
    let lo, up, c2 = hyper_tails(floor(k), m, n, d);
    return up, c2;
}

// ---- the quantile of a discrete family ----------------------------------------------------

// The cumulative function of the family number fam at the count k.
// The numbers of the discrete families in the shared helpers.
const D_BINOM: i64 = 0;
const D_POISSON: i64 = 1;
const D_GEOM: i64 = 2;
const D_NBINOM: i64 = 3;
const D_HYPERGEOM: i64 = 4;

fn disc_cdf(fam: i64, k: f64, a: f64, b: f64, c: f64) -> f64, i64 {
    if fam == D_BINOM { return binom_cdf(k, a, b); }
    if fam == D_POISSON { return poisson_cdf(k, a); }
    if fam == D_GEOM { return geom_cdf(k, a); }
    if fam == D_NBINOM { return nbinom_cdf(k, a, b); }
    return hypergeom_cdf(k, a, b, c);
}

// The smallest count whose cumulative function reaches p. The search
// starts at lo and doubles its upper bound until the cumulative function
// reaches p, so a family without an end needs no bound of its own. The
// count of the steps has a bound, and a search that meets its bound gives
// a fault. A probability of 0 gives the smallest count of the support.
fn disc_ppf(fam: i64, p: f64, a: f64, b: f64, c: f64, lo: f64, hi: f64) -> f64, i64 {
    let cp = chk_unit(p);
    if cp != OK { return nan(), cp; }
    if p == 0.0 { return lo, OK; }
    var top = hi;
    if is_inf(hi) {
        top = lo + 1.0;
        var i: i64 = 0;
        var found = false;
        while i < MAX_ITER_SEARCH {
            let v, cc = disc_cdf(fam, top, a, b, c);
            if cc != OK { return nan(), cc; }
            if v >= p { found = true; break; }
            if top > MAX_COUNT * 0.5 { return nan(), ERR_OVERFLOW; }
            top = top * 2.0;
            i += 1;
        }
        if !found { return nan(), ERR_OVERFLOW; }
    } else {
        let v, cc = disc_cdf(fam, top, a, b, c);
        if cc != OK { return nan(), cc; }
        if v < p { return nan(), ERR_OVERFLOW; }
    }
    var l = lo;
    var r = top;
    var steps: i64 = 0;
    while l < r && steps < MAX_ITER_SEARCH {
        let mid = floor(0.5 * (l + r));
        let v, cc = disc_cdf(fam, mid, a, b, c);
        if cc != OK { return nan(), cc; }
        if v >= p { r = mid; } else { l = mid + 1.0; }
        steps += 1;
    }
    if l < r { return nan(), ERR_OVERFLOW; }
    return l, OK;
}

fn binom_ppf(p: f64, n: f64, prob: f64 = 0.5) -> f64, i64 {
    if is_nan(p) { return nan(), ERR_NAN; }
    let c = binom_check(n, prob);
    if c != OK { return nan(), c; }
    return disc_ppf(D_BINOM, p, n, prob, 0.0, 0.0, n);
}

fn poisson_ppf(p: f64, rate: f64 = 1.0) -> f64, i64 {
    if is_nan(p) || is_nan(rate) { return nan(), ERR_NAN; }
    let c = check_shape(rate);
    if c != OK { return nan(), c; }
    return disc_ppf(D_POISSON, p, rate, 0.0, 0.0, 0.0, 1.0 / 0.0);
}

fn geom_ppf(p: f64, prob: f64 = 0.5) -> f64, i64 {
    if is_nan(p) || is_nan(prob) { return nan(), ERR_NAN; }
    let c = chk_succ(prob);
    if c != OK { return nan(), c; }
    return disc_ppf(D_GEOM, p, prob, 0.0, 0.0, 1.0, 1.0 / 0.0);
}

fn nbinom_ppf(p: f64, r: f64, prob: f64 = 0.5) -> f64, i64 {
    if is_nan(p) { return nan(), ERR_NAN; }
    let c = nbinom_check(r, prob);
    if c != OK { return nan(), c; }
    if prob == 1.0 { return 0.0, OK; }
    return disc_ppf(D_NBINOM, p, r, prob, 0.0, 0.0, 1.0 / 0.0);
}

fn hypergeom_ppf(p: f64, m: f64, n: f64, d: f64) -> f64, i64 {
    if is_nan(p) { return nan(), ERR_NAN; }
    let c = hyper_check(m, n, d);
    if c != OK { return nan(), c; }
    return disc_ppf(D_HYPERGEOM, p, m, n, d, hyper_kmin(m, n, d), hyper_kmax(n, d));
}

// ---- the moments of the discrete families --------------------------------------------------

fn binom_mean(n: f64, p: f64 = 0.5) -> f64, i64 {
    let c = binom_check(n, p);
    if c != OK { return nan(), c; }
    return n * p, OK;
}

fn binom_var(n: f64, p: f64 = 0.5) -> f64, i64 {
    let c = binom_check(n, p);
    if c != OK { return nan(), c; }
    return n * p * (1.0 - p), OK;
}

fn binom_skew(n: f64, p: f64 = 0.5) -> f64, i64 {
    let c = binom_check(n, p);
    if c != OK { return nan(), c; }
    let v = n * p * (1.0 - p);
    if v == 0.0 { return nan(), ERR_UNDEFINED; }
    return (1.0 - 2.0 * p) / sqrt(v), OK;
}

fn binom_kurt(n: f64, p: f64 = 0.5) -> f64, i64 {
    let c = binom_check(n, p);
    if c != OK { return nan(), c; }
    let v = n * p * (1.0 - p);
    if v == 0.0 { return nan(), ERR_UNDEFINED; }
    return (1.0 - 6.0 * p * (1.0 - p)) / v, OK;
}

fn poisson_mean(rate: f64 = 1.0) -> f64, i64 {
    let c = check_shape(rate);
    if c != OK { return nan(), c; }
    return rate, OK;
}

fn poisson_var(rate: f64 = 1.0) -> f64, i64 {
    let c = check_shape(rate);
    if c != OK { return nan(), c; }
    return rate, OK;
}

fn poisson_skew(rate: f64 = 1.0) -> f64, i64 {
    let c = check_shape(rate);
    if c != OK { return nan(), c; }
    return 1.0 / sqrt(rate), OK;
}

fn poisson_kurt(rate: f64 = 1.0) -> f64, i64 {
    let c = check_shape(rate);
    if c != OK { return nan(), c; }
    let v = 1.0 / rate;
    if is_inf(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

fn geom_mean(p: f64 = 0.5) -> f64, i64 {
    let c = chk_succ(p);
    if c != OK { return nan(), c; }
    let v = 1.0 / p;
    if is_inf(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

fn geom_var(p: f64 = 0.5) -> f64, i64 {
    let c = chk_succ(p);
    if c != OK { return nan(), c; }
    let v = (1.0 - p) / (p * p);
    if is_inf(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

fn geom_skew(p: f64 = 0.5) -> f64, i64 {
    let c = chk_succ(p);
    if c != OK { return nan(), c; }
    if p == 1.0 { return nan(), ERR_UNDEFINED; }
    return (2.0 - p) / sqrt(1.0 - p), OK;
}

fn geom_kurt(p: f64 = 0.5) -> f64, i64 {
    let c = chk_succ(p);
    if c != OK { return nan(), c; }
    if p == 1.0 { return nan(), ERR_UNDEFINED; }
    return 6.0 + p * p / (1.0 - p), OK;
}

fn nbinom_mean(r: f64, p: f64 = 0.5) -> f64, i64 {
    let c = nbinom_check(r, p);
    if c != OK { return nan(), c; }
    let v = r * (1.0 - p) / p;
    if is_inf(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

fn nbinom_var(r: f64, p: f64 = 0.5) -> f64, i64 {
    let c = nbinom_check(r, p);
    if c != OK { return nan(), c; }
    let v = r * (1.0 - p) / (p * p);
    if is_inf(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

fn nbinom_skew(r: f64, p: f64 = 0.5) -> f64, i64 {
    let c = nbinom_check(r, p);
    if c != OK { return nan(), c; }
    if p == 1.0 { return nan(), ERR_UNDEFINED; }
    return (2.0 - p) / sqrt(r * (1.0 - p)), OK;
}

fn nbinom_kurt(r: f64, p: f64 = 0.5) -> f64, i64 {
    let c = nbinom_check(r, p);
    if c != OK { return nan(), c; }
    if p == 1.0 { return nan(), ERR_UNDEFINED; }
    return 6.0 / r + p * p / (r * (1.0 - p)), OK;
}

fn hypergeom_mean(m: f64, n: f64, d: f64) -> f64, i64 {
    let c = hyper_check(m, n, d);
    if c != OK { return nan(), c; }
    return n * d / m, OK;
}

// The variance d n (m - n) (m - d) / (m^2 (m - 1)) as a pair.
fn hyper_var_dd(m: f64, n: f64, d: f64) -> DD {
    if m == 1.0 { return dd_from(0.0); }
    let num = dd_mul_d(dd_mul_d(dd_mul_d(dd_from(d), n), m - n), m - d);
    let den = dd_mul_d(dd_mul_d(dd_from(m), m), m - 1.0);
    return dd_div(num, den);
}

fn hypergeom_var(m: f64, n: f64, d: f64) -> f64, i64 {
    let c = hyper_check(m, n, d);
    if c != OK { return nan(), c; }
    return dd_to_f64(hyper_var_dd(m, n, d)), OK;
}

// The third and the fourth central moment as sums of the masses of the
// support, in pairs. The sums run outward from the mode, so the large
// terms come first and the loop ends on the bound of the support.
fn hyper_central(m: f64, n: f64, d: f64) -> DD, DD, i64 {
    let kmin = hyper_kmin(m, n, d);
    let kmax = hyper_kmax(n, d);
    if kmax - kmin + 1.0 > MAX_DIST_SUPPORT { return dd_from(nan()), dd_from(nan()), ERR_TOO_LONG; }
    let mean = dd_div_d(dd_mul_d(dd_from(n), d), m);
    let mode = hyper_mode(m, n, d);
    let t = hyper_top(m, n, d, mode);
    var m3 = dd_from(0.0);
    var m4 = dd_from(0.0);
    var j = mode;
    var cur = t;
    while j >= kmin {
        let e = dd_sub(dd_from(j), mean);
        let e2 = dd_mul(e, e);
        m3 = dd_add(m3, dd_mul(cur, dd_mul(e, e2)));
        m4 = dd_add(m4, dd_mul(cur, dd_mul(e2, e2)));
        if j == kmin || cur.hi == 0.0 { break; }
        cur = hyper_down(cur, j, m, n, d);
        j -= 1.0;
    }
    j = mode;
    cur = t;
    while j < kmax {
        cur = hyper_up(cur, j, m, n, d);
        j += 1.0;
        if cur.hi == 0.0 { break; }
        let e = dd_sub(dd_from(j), mean);
        let e2 = dd_mul(e, e);
        m3 = dd_add(m3, dd_mul(cur, dd_mul(e, e2)));
        m4 = dd_add(m4, dd_mul(cur, dd_mul(e2, e2)));
    }
    return m3, m4, OK;
}

fn hypergeom_skew(m: f64, n: f64, d: f64) -> f64, i64 {
    let c = hyper_check(m, n, d);
    if c != OK { return nan(), c; }
    let v = hyper_var_dd(m, n, d);
    if v.hi == 0.0 { return nan(), ERR_UNDEFINED; }
    let m3, m4, c2 = hyper_central(m, n, d);
    if c2 != OK { return nan(), c2; }
    return dd_to_f64(dd_div(m3, dd_mul(v, dd_sqrt(v)))), OK;
}

fn hypergeom_kurt(m: f64, n: f64, d: f64) -> f64, i64 {
    let c = hyper_check(m, n, d);
    if c != OK { return nan(), c; }
    let v = hyper_var_dd(m, n, d);
    if v.hi == 0.0 { return nan(), ERR_UNDEFINED; }
    let m3, m4, c2 = hyper_central(m, n, d);
    if c2 != OK { return nan(), c2; }
    let v2 = dd_mul(v, v);
    return dd_to_f64(dd_div(dd_sub(m4, dd_mul_d(v2, 3.0)), v2)), OK;
}
