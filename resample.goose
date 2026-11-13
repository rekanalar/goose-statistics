// The resampling calls of the package: the bootstrap, the checks and the
// draws.
//
// Every call returns the result first and a fault code last. A fault
// gives NaN as the result.
namespace statistics;

import std;
import math;

// The method of a resampling call. One struct covers the bootstrap, the
// permutation test and the Monte Carlo test. The field interval counts for
// the bootstrap only. The field permutation counts for the permutation
// test only.
struct ResampleMethod {
    kind: i64 = RESAMPLE_BOOTSTRAP,
    resamples: i64 = 9999,
    interval: i64 = INTERVAL_BCA,
    permutation: i64 = PERMUTE_INDEPENDENT
}

// A pair of bounds. A one-sided request gives an infinity for the open end.
struct Interval { low: f64, high: f64 }

// The interval of a bootstrap and the standard error of the statistic.
struct BootstrapResult { interval: Interval, standard_error: f64 }

// The result of a call that faulted.
fn rs_boot_fault(code: i64) -> BootstrapResult, i64 {
    let iv = Interval { low: nan(), high: nan() };
    return BootstrapResult { interval: iv, standard_error: nan() }, code;
}

// ---- checks -----------------------------------------------------------

// The checks of the samples in one list. The sizes array holds the length
// of each sample. Each sample holds two values or more. The call needs at
// least min_groups samples.
fn rs_check_data(data: const f64[:], sizes: const i64[:], min_groups: i64) -> i64 {
    return check_groups_n(data, sizes, min_groups, 2);
}

// The checks of the resample count and of the work. The work is the count
// times the size of one resample. A count below 1 is out of range.
fn rs_check_work(resamples: i64, size: i64) -> i64 {
    if resamples > MAX_RESAMPLES { return ERR_TOO_LONG; }
    if resamples * size > MAX_RESAMPLE_WORK { return ERR_TOO_LONG; }
    if resamples < 1 { return ERR_DOMAIN; }
    return OK;
}

// The value of one statistic. A NaN or an infinity is a fault.
fn rs_stat<F>(xs: const f64[:]) -> f64, i64 {
    let v = F(xs);
    if is_nan(v) { return nan(), ERR_NAN; }
    if is_inf(v) { return nan(), ERR_INFINITY; }
    return v, OK;
}

// True when every value of the list equals the first one.
fn rs_constant(xs: const f64[:]) -> bool {
    for x in xs {
        if x != xs[0] { return false; }
    }
    return true;
}

// ---- draws ------------------------------------------------------------

// A uniform integer in the range from 0 to n, without n. The draw rejects
// the low values of the word, so every result counts the same. The older
// releases of the standard library draw the low values more often, so the
// package makes its own draw. One seed then gives one answer on every
// supported release.
fn rs_rand_int(r: rng&, n: i64) -> i64 {
    let m = n as u64;
    var x = rand_u64(r);
    if x < m {
        let t = (0 - m) % m;
        while x < t { x = rand_u64(r); }
    }
    return (x % m) as i64;
}

// A uniform draw in the open interval from 0 to 1. The draw is a multiple
// of 2^-52 plus half a step, so it is never 0 and never 1.
fn rs_unit(r: rng&) -> f64 {
    let k = rand_u64(r) >> 12;
    return (k + 0.5) * EPS;
}

// A standard normal draw by the Box-Muller form. It uses two uniform draws.
fn rs_normal(r: rng&) -> f64 {
    let u1 = rs_unit(r);
    let u2 = rs_unit(r);
    return sqrt(-2.0 * log(u1)) * cos(6.283185307179586 * u2);
}

// The check of the parameters of a family. Every value of the list must be
// finite, and the list must hold the right count.
fn rs_check_params(family: i64, params: const f64[:]) -> i64 {
    for p in params {
        let cp = check_scalar(p);
        if cp != OK { return cp; }
    }
    if family < 0 || family > FAMILY_F { return ERR_DOMAIN; }
    if params.len != family_arity(family) { return ERR_LENGTH; }
    return OK;
}

// The check of the range of the parameters of a family. The quantile at
// 0.5 finds a parameter that is out of range.
fn rs_check_range(family: i64, params: const f64[:]) -> i64 {
    let v, cq = rs_ppf(family, 0.5, params);
    return cq;
}

// The quantile of a family at one probability.
fn rs_ppf(family: i64, u: f64, params: const f64[:]) -> f64, i64 {
    if family == FAMILY_NORMAL { let v, c = norm_ppf(u, params[0], params[1]); return v, c; }
    if family == FAMILY_UNIFORM { let v, c = unif_ppf(u, params[0], params[1]); return v, c; }
    if family == FAMILY_EXPONENTIAL { let v, c = expon_ppf(u, params[0]); return v, c; }
    if family == FAMILY_LOGNORMAL { let v, c = lognorm_ppf(u, params[0], params[1]); return v, c; }
    if family == FAMILY_WEIBULL { let v, c = weibull_ppf(u, params[0], params[1]); return v, c; }
    if family == FAMILY_GAMMA { let v, c = gammadist_ppf(u, params[0], params[1]); return v, c; }
    if family == FAMILY_BETA { let v, c = betadist_ppf(u, params[0], params[1]); return v, c; }
    if family == FAMILY_CHI2 { let v, c = chi2_ppf(u, params[0]); return v, c; }
    if family == FAMILY_T { let v, c = tdist_ppf(u, params[0]); return v, c; }
    let v, c = fdist_ppf(u, params[0], params[1]);
    return v, c;
}

// One draw from a checked family. The normal family uses the Box-Muller
// form. The other families use the inverse of the cumulative function.
fn rs_draw(family: i64, params: const f64[:], r: rng&) -> f64, i64 {
    if family == FAMILY_NORMAL {
        let v = params[0] + params[1] * rs_normal(r);
        if is_bad(v) { return nan(), ERR_OVERFLOW; }
        return v, OK;
    }
    let v, c = rs_ppf(family, rs_unit(r), params);
    if c != OK { return nan(), c; }
    if is_bad(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

// Fills the buffer with n new draws, one after the other.
fn rs_fill(buf: f64[:], n: i64, family: i64, params: const f64[:], r: rng&) -> i64 {
    for i in 0..n {
        let v, c = rs_draw(family, params, r);
        if c != OK { return c; }
        buf[i] = v;
    }
    return OK;
}

// A Fisher-Yates shuffle of the values from lo to hi, without hi.
fn rs_shuffle(buf: f64[:], lo: i64, hi: i64, r: rng&) {
    var i = hi - 1;
    while i > lo {
        let j = lo + rs_rand_int(r, i - lo + 1);
        let t = buf[i];
        buf[i] = buf[j];
        buf[j] = t;
        i -= 1;
    }
}

// A Fisher-Yates shuffle of the indices from lo to hi, without hi. The
// draw rejects the low values of the word, so the permutation is even on
// the releases where the draw of the standard library is not.
fn rs_shuffle_index(buf: i64[:], lo: i64, hi: i64, r: rng&) {
    var i = hi - 1;
    while i > lo {
        let j = lo + rs_rand_int(r, i - lo + 1);
        let t = buf[i];
        buf[i] = buf[j];
        buf[j] = t;
        i -= 1;
    }
}

// ---- the p-value ------------------------------------------------------

// The p-value of an observed value against a null list. An exact test adds
// nothing to the count. A drawn test adds 1 to the count and to the total.
// The tolerance is 100 eps of the observed value.
fn rs_pvalue(nul: const f64[:], observed: f64, alt: i64, exact: bool) -> f64 {
    let gamma = abs(100.0 * EPS * observed);
    var less: i64 = 0;
    var greater: i64 = 0;
    for v in nul {
        if v <= observed + gamma { less += 1; }
        if v >= observed - gamma { greater += 1; }
    }
    var add: f64 = 0.0;
    if !exact { add = 1.0; }
    let lt: f64 = less;
    let gt: f64 = greater;
    let total: f64 = nul.len;
    let pl = (lt + add) / (total + add);
    let pg = (gt + add) / (total + add);
    if alt == ALT_LESS { return min(pl, 1.0); }
    if alt == ALT_GREATER { return min(pg, 1.0); }
    return min(2.0 * min(pl, pg), 1.0);
}

// ---- the bootstrap ----------------------------------------------------

// The cumulative normal function, exact to 0 and 1 beyond 40 deviations.
fn rs_ndtr(x: f64) -> f64 {
    if x > 40.0 { return 1.0; }
    if x < -40.0 { return 0.0; }
    let v, c = norm_cdf(x);
    return v;
}

// The bias term of the BCa interval: the normal quantile of the share of
// the bootstrap values below the observed value, counting ties as half.
fn rs_bias(boot: const f64[:], observed: f64) -> f64, i64 {
    var lt: i64 = 0;
    var le: i64 = 0;
    for v in boot {
        if v < observed { lt += 1; }
        if v <= observed { le += 1; }
    }
    let p0 = (lt + le) / (2.0 * boot.len);
    if p0 <= 0.0 || p0 >= 1.0 { return nan(), ERR_UNDEFINED; }
    return norm_std_ppf(p0), OK;
}

// The acceleration of the BCa interval from the jackknife of one sample.
// The jackknife leaves out one value at a time. It runs n statistics on
// n - 1 values each.
fn rs_accel<F>(data: const f64[:]) -> f64, i64 {
    let n = data.len;
    var buf: f64[>..] = [];
    for i in 1..n { buf.push(data[i]); }
    var sum: f64 = 0.0;
    var thetas: f64[>..] = [];
    for i in 0..n {
        if i > 0 { buf[i - 1] = data[i - 1]; }
        let v, c = rs_stat(buf[..]) { s => F(s) };
        if c != OK { return nan(), c; }
        thetas.push(v);
        sum += v;
    }
    let dot = sum / n;
    var s2: f64 = 0.0;
    var s3: f64 = 0.0;
    for v in thetas {
        let u = (n - 1) * (dot - v);
        s2 += u * u;
        s3 += u * u * u;
    }
    let num = s3 / (n * n * n);
    let den = s2 / (n * n);
    if den <= 0.0 { return nan(), ERR_UNDEFINED; }
    let a = num / (6.0 * den * sqrt(den));
    if is_bad(a) { return nan(), ERR_UNDEFINED; }
    return a, OK;
}

// The interval of a list of bootstrap values. The kind is one of the
// INTERVAL_ constants. The accel is the acceleration of the BCa form.
fn rs_interval(boot: const f64[:], observed: f64, accel: f64, level: f64, alt: i64, kind: i64) -> f64, f64, i64 {
    var alpha = (1.0 - level) / 2.0;
    if alt != ALT_TWO_SIDED { alpha = 1.0 - level; }
    var q1 = alpha;
    var q2 = 1.0 - alpha;
    if kind == INTERVAL_BCA {
        let z0, cz = rs_bias(boot, observed);
        if cz != OK { return nan(), nan(), cz; }
        let za = norm_std_ppf(alpha);
        let n1 = z0 + za;
        let d1 = 1.0 - accel * n1;
        let n2 = z0 - za;
        let d2 = 1.0 - accel * n2;
        if d1 == 0.0 || d2 == 0.0 { return nan(), nan(), ERR_UNDEFINED; }
        q1 = rs_ndtr(z0 + n1 / d1);
        q2 = rs_ndtr(z0 + n2 / d2);
        if is_nan(q1) || is_nan(q2) { return nan(), nan(), ERR_UNDEFINED; }
    }
    let s = sort_copy(boot);
    var lo = quantile_sorted(s[..], q1, METHOD_LINEAR);
    var hi = quantile_sorted(s[..], q2, METHOD_LINEAR);
    if kind == INTERVAL_BASIC {
        let l2 = 2.0 * observed - hi;
        let h2 = 2.0 * observed - lo;
        lo = l2;
        hi = h2;
    }
    if is_nan(lo) || is_nan(hi) { return nan(), nan(), ERR_UNDEFINED; }
    if is_inf(lo) || is_inf(hi) { return nan(), nan(), ERR_OVERFLOW; }
    if alt == ALT_LESS { lo = -1.0 / 0.0; }
    if alt == ALT_GREATER { hi = 1.0 / 0.0; }
    return lo, hi, OK;
}

// The bootstrap interval of a statistic. The data holds the samples one
// after the other and sizes holds the length of each. A resample draws each
// sample on its own, with replacement. The block gets the resampled list in
// the same layout. The BCa form works for one sample only, because its
// jackknife shortens the sample. The list out holds the resampled values of
// the statistic.
fn bootstrap<F>(out: f64[>..]&, data: const f64[:], sizes: const i64[:], method: ResampleMethod, level: f64, alt: i64, r: rng&) -> BootstrapResult, i64 {
    out.resize(0);
    var c = rs_check_data(data, sizes, 1);
    if c != OK { return rs_boot_fault(c); }
    c = rs_check_work(method.resamples, data.len);
    if c != OK { return rs_boot_fault(c); }
    if method.interval == INTERVAL_BCA && data.len * data.len > MAX_RESAMPLE_WORK {
        return rs_boot_fault(ERR_TOO_LONG);
    }
    if method.kind != RESAMPLE_BOOTSTRAP { return rs_boot_fault(ERR_DOMAIN); }
    if method.interval < INTERVAL_PERCENTILE || method.interval > INTERVAL_BCA { return rs_boot_fault(ERR_DOMAIN); }
    c = check_alt(alt);
    if c != OK { return rs_boot_fault(c); }
    c = check_level(level);
    if c != OK { return rs_boot_fault(c); }
    if method.interval == INTERVAL_BCA && sizes.len > 1 { return rs_boot_fault(ERR_DOMAIN); }
    if method.resamples < 2 { return rs_boot_fault(ERR_TOO_SHORT); }
    if rs_constant(data) { return rs_boot_fault(ERR_UNDEFINED); }
    let observed, co = rs_stat(data) { s => F(s) };
    if co != OK { return rs_boot_fault(co); }
    // The acceleration needs no draw, so a fault shows before the resamples.
    var accel: f64 = 0.0;
    if method.interval == INTERVAL_BCA {
        let a, ca = rs_accel(data) { s => F(s) };
        if ca != OK { return rs_boot_fault(ca); }
        accel = a;
    }
    var buf: f64[>..] = [];
    for x in data { buf.push(x); }
    for b in 0..method.resamples {
        var off: i64 = 0;
        for g in sizes {
            for i in 0..g { buf[off + i] = data[off + rs_rand_int(r, g)]; }
            off += g;
        }
        let v, cv = rs_stat(buf[..]) { s => F(s) };
        if cv != OK { out.resize(0); return rs_boot_fault(cv); }
        out.push(v);
    }
    let se, cs = stdev_core(out[..], 1);
    if cs != OK { out.resize(0); return rs_boot_fault(cs); }
    let lo, hi, ci = rs_interval(out[..], observed, accel, level, alt, method.interval);
    if ci != OK { out.resize(0); return rs_boot_fault(ci); }
    return BootstrapResult { interval: Interval { low: lo, high: hi }, standard_error: se }, OK;
}
