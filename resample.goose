// The resampling calls: the bootstrap interval, the permutation test,
// the Monte Carlo test, the power of a test and the Sobol indices.
//
// The statistic is the trailing block of each call. It takes a list of
// values and returns one f64. The caller owns the random generator and
// passes it in. A call draws one resample at a time, in row-major order,
// sample by sample, so one seed gives one answer. The draws of this package
// differ from the draws of SciPy, so a seeded result agrees with SciPy
// within the Monte Carlo error only.
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

// ---- the permutation test ---------------------------------------------

// The binomial coefficient of n over k, or cap + 1 when it exceeds cap.
fn rs_comb_cap(n: i64, k: i64, cap: i64) -> i64 {
    var kk = k;
    if n - k < kk { kk = n - k; }
    var c: i64 = 1;
    for i in 1..kk + 1 {
        c = c * (n - kk + i) / i;
        if c > cap { return cap + 1; }
    }
    return c;
}

// The factorial of n, or cap + 1 when it exceeds cap.
fn rs_fact_cap(n: i64, cap: i64) -> i64 {
    var f: i64 = 1;
    for i in 2..n + 1 {
        f *= i;
        if f > cap { return cap + 1; }
    }
    return f;
}

// The power b to the power e, or cap + 1 when it exceeds cap.
fn rs_pow_cap(b: i64, e: i64, cap: i64) -> i64 {
    var p: i64 = 1;
    for i in 0..e {
        p *= b;
        if p > cap { return cap + 1; }
    }
    return p;
}

// The count of the distinct permutations of a permutation type, or cap + 1
// when it exceeds cap.
fn rs_nmax(ptype: i64, sizes: const i64[:], cap: i64) -> i64 {
    let m = sizes.len;
    if ptype == PERMUTE_INDEPENDENT {
        var p: i64 = 1;
        var prev: i64 = 0;
        for s in sizes {
            let c = rs_comb_cap(prev + s, s, cap);
            if p > cap / c { return cap + 1; }
            p *= c;
            prev += s;
        }
        return p;
    }
    if ptype == PERMUTE_SAMPLES {
        var k = m;
        if m == 1 { k = 2; }
        return rs_pow_cap(rs_fact_cap(k, cap), sizes[0], cap);
    }
    return rs_pow_cap(rs_fact_cap(sizes[0], cap), m, cap);
}

// The next permutation of a list in the lexicographic order. The call
// returns false and puts the list back in the first order after the last
// permutation.
fn rs_next_perm(p: i64[:]) -> bool {
    var i = p.len - 2;
    while i >= 0 && p[i] >= p[i + 1] { i -= 1; }
    if i < 0 {
        var a: i64 = 0;
        var b = p.len - 1;
        while a < b { let t = p[a]; p[a] = p[b]; p[b] = t; a += 1; b -= 1; }
        return false;
    }
    var j = p.len - 1;
    while p[j] <= p[i] { j -= 1; }
    let t = p[i]; p[i] = p[j]; p[j] = t;
    var a = i + 1;
    var b = p.len - 1;
    while a < b { let t2 = p[a]; p[a] = p[b]; p[b] = t2; a += 1; b -= 1; }
    return true;
}

// The statistic of every distinct split of the samples. The labels give the
// group of each value. The list walks through the permutations of the
// labels in the lexicographic order.
fn rs_exact_independent<F>(out: f64[>..]&, data: const f64[:], sizes: const i64[:], count: i64) -> i64 {
    let n = data.len;
    let m = sizes.len;
    var lab: i64[>..] = [];
    for g in 0..m {
        for i in 0..sizes[g] { lab.push(g); }
    }
    var res: f64[>..] = [];
    for x in data { res.push(x); }
    var cur: i64[>..] = [];
    for g in 0..m { cur.push(0); }
    for t in 0..count {
        var off: i64 = 0;
        for g in 0..m { cur[g] = off; off += sizes[g]; }
        for k in 0..n {
            let g = lab[k];
            res[cur[g]] = data[k];
            cur[g] += 1;
        }
        let v, c = rs_stat(res[..]) { s => F(s) };
        if c != OK { return c; }
        out.push(v);
        rs_next_perm(lab[..]);
    }
    return OK;
}

// The statistic of every combination of the permutations of the pairs or of
// the pairings. The list ix holds the blocks of the permutations. The last
// block advances first.
fn rs_exact_odometer<F>(out: f64[>..]&, data: const f64[:], sizes: const i64[:], ptype: i64, count: i64) -> i64 {
    let m = sizes.len;
    let n = sizes[0];
    // The number of blocks and the size of a block. The pairings permute
    // the observations of each sample. The pairs permute the samples of
    // each observation.
    var l = m;
    var k = n;
    if ptype == PERMUTE_SAMPLES {
        l = n;
        k = m;
        if m == 1 { k = 2; }
    }
    var ix: i64[>..] = [];
    for a in 0..l {
        for b in 0..k { ix.push(b); }
    }
    var res: f64[>..] = [];
    for x in data { res.push(x); }
    for t in 0..count {
        if ptype == PERMUTE_PAIRINGS {
            for j in 0..m {
                for i in 0..n { res[j * n + i] = data[j * n + ix[j * n + i]]; }
            }
        } else if m == 1 {
            for i in 0..n {
                if ix[i * 2] == 0 { res[i] = data[i]; } else { res[i] = -data[i]; }
            }
        } else {
            for j in 0..m {
                for i in 0..n { res[j * n + i] = data[ix[i * m + j] * n + i]; }
            }
        }
        let v, c = rs_stat(res[..]) { s => F(s) };
        if c != OK { return c; }
        out.push(v);
        var q = l - 1;
        while q >= 0 {
            if rs_next_perm(ix[q * k..(q + 1) * k]) { break; }
            q -= 1;
        }
    }
    return OK;
}

// One new random permutation of the working list, of the given type. Each
// call shuffles the list of the last call, so the draw stays uniform.
fn rs_permute_once(res: f64[:], sizes: const i64[:], ptype: i64, r: rng&) {
    let m = sizes.len;
    let n = sizes[0];
    if ptype == PERMUTE_INDEPENDENT {
        rs_shuffle(res, 0, res.len, r);
    } else if ptype == PERMUTE_PAIRINGS {
        for j in 0..m { rs_shuffle(res, j * n, (j + 1) * n, r); }
    } else if m == 1 {
        for i in 0..n {
            if rs_rand_int(r, 2) == 1 { res[i] = -res[i]; }
        }
    } else {
        for i in 0..n {
            var j = m - 1;
            while j > 0 {
                let q = rs_rand_int(r, j + 1);
                let t = res[j * n + i];
                res[j * n + i] = res[q * n + i];
                res[q * n + i] = t;
                j -= 1;
            }
        }
    }
}

// The permutation test of a statistic. The permutation type is a PERMUTE_
// constant of the method. The call lists every permutation when their
// count is at most the resamples of the method, and the p-value is then
// exact. Otherwise it draws the resamples. The list out holds the null
// values of the statistic.
fn permutation_test<F>(out: f64[>..]&, data: const f64[:], sizes: const i64[:], method: ResampleMethod, alt: i64, r: rng&) -> TestResult, i64 {
    out.resize(0);
    let ptype = method.permutation;
    var need: i64 = 1;
    if ptype == PERMUTE_INDEPENDENT { need = 2; }
    var c = rs_check_data(data, sizes, need);
    if c != OK { return fault_result(c); }
    c = rs_check_work(method.resamples, data.len);
    if c != OK { return fault_result(c); }
    if ptype == PERMUTE_SAMPLES || ptype == PERMUTE_PAIRINGS {
        for s in sizes {
            if s != sizes[0] { return fault_result(ERR_LENGTH); }
        }
    }
    if ptype < PERMUTE_INDEPENDENT || ptype > PERMUTE_PAIRINGS { return fault_result(ERR_DOMAIN); }
    if method.kind != RESAMPLE_PERMUTATION { return fault_result(ERR_DOMAIN); }
    c = check_alt(alt);
    if c != OK { return fault_result(c); }
    if rs_constant(data) { return fault_result(ERR_UNDEFINED); }
    let observed, co = rs_stat(data) { s => F(s) };
    if co != OK { return fault_result(co); }
    let nmax = rs_nmax(ptype, sizes, method.resamples);
    let exact = nmax <= method.resamples;
    if exact {
        var ce: i64 = OK;
        if ptype == PERMUTE_INDEPENDENT {
            ce = rs_exact_independent(out, data, sizes, nmax) { s => F(s) };
        } else {
            ce = rs_exact_odometer(out, data, sizes, ptype, nmax) { s => F(s) };
        }
        if ce != OK { out.resize(0); return fault_result(ce); }
    } else {
        var res: f64[>..] = [];
        for x in data { res.push(x); }
        for b in 0..method.resamples {
            rs_permute_once(res[..], sizes, ptype, r);
            let v, cv = rs_stat(res[..]) { s => F(s) };
            if cv != OK { out.resize(0); return fault_result(cv); }
            out.push(v);
        }
    }
    return good_result(observed, rs_pvalue(out[..], observed, alt, exact), nan(), nan());
}

// ---- the Monte Carlo test and the power -------------------------------

// The Monte Carlo test of a statistic. Each resample draws as many values
// as the sample holds from the family. The p-value always adds 1 to the
// count and to the total. The list out holds the null values.
fn monte_carlo_test<F>(out: f64[>..]&, xs: const f64[:], family: i64, params: const f64[:], method: ResampleMethod, alt: i64, r: rng&) -> TestResult, i64 {
    out.resize(0);
    var c = check_sample(xs, 1);
    if c != OK { return fault_result(c); }
    c = rs_check_params(family, params);
    if c != OK { return fault_result(c); }
    c = rs_check_work(method.resamples, xs.len);
    if c != OK { return fault_result(c); }
    c = rs_check_range(family, params);
    if c != OK { return fault_result(c); }
    if method.kind != RESAMPLE_MONTE_CARLO { return fault_result(ERR_DOMAIN); }
    c = check_alt(alt);
    if c != OK { return fault_result(c); }
    let observed, co = rs_stat(xs) { s => F(s) };
    if co != OK { return fault_result(co); }
    var buf: f64[>..] = [];
    for x in xs { buf.push(x); }
    for b in 0..method.resamples {
        let cf = rs_fill(buf[..], xs.len, family, params, r);
        if cf != OK { out.resize(0); return fault_result(cf); }
        let v, cv = rs_stat(buf[..]) { s => F(s) };
        if cv != OK { out.resize(0); return fault_result(cv); }
        out.push(v);
    }
    return good_result(observed, rs_pvalue(out[..], observed, alt, false), nan(), nan());
}

// The power of a test by simulation. Each resample draws n values from the
// family. The block is the test and returns its p-value. The power is the
// share of the p-values below the significance level.
fn power<F>(family: i64, params: const f64[:], n: i64, level: f64, resamples: i64, r: rng&) -> f64, i64 {
    if n > MAX_INPUT { return nan(), ERR_TOO_LONG; }
    if n < 1 { return nan(), ERR_EMPTY; }
    var c = rs_check_params(family, params);
    if c != OK { return nan(), c; }
    c = rs_check_work(resamples, n);
    if c != OK { return nan(), c; }
    c = rs_check_range(family, params);
    if c != OK { return nan(), c; }
    c = check_level(level);
    if c != OK { return nan(), c; }
    var buf: f64[>..] = [];
    buf.resize(n, 0.0);
    var hits: i64 = 0;
    for b in 0..resamples {
        let cf = rs_fill(buf[..], n, family, params, r);
        if cf != OK { return nan(), cf; }
        let p, cp = rs_stat(buf[..]) { s => F(s) };
        if cp != OK { return nan(), cp; }
        if p < 0.0 || p > 1.0 { return nan(), ERR_DOMAIN; }
        if p < level { hits += 1; }
    }
    let h: f64 = hits;
    let rr: f64 = resamples;
    return h / rr, OK;
}

// ---- the Sobol indices ------------------------------------------------

// The mean and the variance of the two lists of evaluations. The lists
// hold the same count of values. A variance at or below 0 gives the
// undefined fault.
fn rs_sobol_center(f_a: const f64[:], f_b: const f64[:]) -> f64, f64, i64 {
    let n = f_a.len;
    var m: f64 = 0.0;
    for i in 0..n { m += f_a[i] + f_b[i]; }
    m = m / (2.0 * n);
    var v: f64 = 0.0;
    for i in 0..n {
        v += (f_a[i] - m) * (f_a[i] - m) + (f_b[i] - m) * (f_b[i] - m);
    }
    v = v / (2.0 * n);
    if is_inf(v) { return nan(), nan(), ERR_OVERFLOW; }
    if v <= 0.0 { return nan(), nan(), ERR_UNDEFINED; }
    return m, v, OK;
}

// The indices from three lists of evaluations: f_a at the points A, f_b at
// the points B and f_ab at the d matrices that take column k from B. The
// list f_ab holds d pages of n values, page by page. The first-order index
// is the Saltelli form and the total index is the Saltelli form of 2010.
// The values are centered by the mean of f_a and f_b. A constant output has
// no variance and gives the undefined fault.
fn sobol_from_evaluations(first: f64[>..]&, total: f64[>..]&, f_a: const f64[:], f_b: const f64[:], f_ab: const f64[:], d: i64) -> i64, i64 {
    first.resize(0);
    total.resize(0);
    if d > MAX_SOBOL_DIMS { return 0, ERR_TOO_LONG; }
    var c = check_pair(f_a, f_b, 2);
    if c != OK { return 0, c; }
    c = check_sample(f_ab, 1);
    if c != OK { return 0, c; }
    if d < 1 { return 0, ERR_DOMAIN; }
    let n = f_a.len;
    if f_ab.len != d * n { return 0, ERR_LENGTH; }
    let m, v, cv = rs_sobol_center(f_a, f_b);
    if cv != OK { return 0, cv; }
    for k in 0..d {
        var s1: f64 = 0.0;
        var s2: f64 = 0.0;
        for i in 0..n {
            let ab = f_ab[k * n + i];
            s1 += (f_b[i] - m) * (ab - f_a[i]);
            s2 += (f_a[i] - ab) * (f_a[i] - ab);
        }
        let a = s1 / n / v;
        let t = 0.5 * s2 / n / v;
        if is_bad(a) || is_bad(t) {
            first.resize(0);
            total.resize(0);
            return 0, ERR_OVERFLOW;
        }
        first.push(a);
        total.push(t);
    }
    return d, OK;
}

// The first-order and total-order Sobol indices of a model. The block is
// the model. It takes a point of d values, each in the open interval from
// 0 to 1, and returns one value. The caller maps the point to its own
// distributions. The call draws n points A and n points B from the
// generator, and runs the model n * (d + 2) times. The draws are
// pseudo-random, not a scrambled Sobol sequence, so the error shrinks more
// slowly than the error of SciPy for the same n.
fn sobol_indices<F>(first: f64[>..]&, total: f64[>..]&, d: i64, n: i64, r: rng&) -> i64, i64 {
    first.resize(0);
    total.resize(0);
    if d > MAX_SOBOL_DIMS { return 0, ERR_TOO_LONG; }
    if n > MAX_SOBOL_POINTS { return 0, ERR_TOO_LONG; }
    if n < 2 { return 0, ERR_TOO_SHORT; }
    if d < 1 { return 0, ERR_DOMAIN; }
    // The state before the first draw. A second pass replays the draws, so
    // the points of the matrices A and B never sit in memory.
    var replay = rng { s: r.s };
    var pa: f64[32];
    var pb: f64[32];
    var pt: f64[32];
    var f_a: f64[>..] = [];
    var f_b: f64[>..] = [];
    for j in 0..n {
        for i in 0..d { pa[i] = rs_unit(r); }
        for i in 0..d { pb[i] = rs_unit(r); }
        let va, ca = rs_stat(pa[0..d]) { s => F(s) };
        if ca != OK { return 0, ca; }
        let vb, cb = rs_stat(pb[0..d]) { s => F(s) };
        if cb != OK { return 0, cb; }
        f_a.push(va);
        f_b.push(vb);
    }
    let m, v, cv = rs_sobol_center(f_a, f_b);
    if cv != OK { return 0, cv; }
    first.resize(d, 0.0);
    total.resize(d, 0.0);
    for j in 0..n {
        for i in 0..d { pa[i] = rs_unit(replay); }
        for i in 0..d { pb[i] = rs_unit(replay); }
        for k in 0..d {
            for i in 0..d { pt[i] = pa[i]; }
            pt[k] = pb[k];
            let ab, cab = rs_stat(pt[0..d]) { s => F(s) };
            if cab != OK { first.resize(0); total.resize(0); return 0, cab; }
            first[k] += (f_b[j] - m) * (ab - f_a[j]);
            total[k] += (f_a[j] - ab) * (f_a[j] - ab);
        }
    }
    for k in 0..d {
        first[k] = first[k] / n / v;
        total[k] = 0.5 * total[k] / n / v;
        if is_bad(first[k]) || is_bad(total[k]) {
            first.resize(0);
            total.resize(0);
            return 0, ERR_OVERFLOW;
        }
    }
    return d, OK;
}
