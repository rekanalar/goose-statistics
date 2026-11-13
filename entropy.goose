// The entropy of a distribution, the differential entropy of a sample
// and the k-statistics.
namespace statistics;

import std;
import math;

// Signature: entropy(pk, base) -> value, code
// The Shannon entropy of a set of probabilities. The values are scaled
// so that they sum to 1. A base of 0 gives the natural logarithm. A
// negative value gives ERR_DOMAIN, as in SciPy, which gives -inf. A sum
// of 0 gives ERR_UNDEFINED.
fn entropy(pk: const f64[:], base: f64) -> f64, i64 {
    let c = check_probabilities_list(pk);
    if c != 0 { return nan(), c; }
    let bc = check_base(base);
    if bc != 0 { return nan(), bc; }
    var acc = acc_new();
    for p in pk { acc_add(acc, p); }
    if acc.fault { return nan(), ERR_UNDEFINED; }
    let total = acc_total(acc);
    if total <= 0.0 { return nan(), ERR_UNDEFINED; }
    var s: f64 = 0.0;
    for p in pk {
        // A share that underflows to 0 adds 0, the limit of q log q.
        let q = p / total;
        if q == 0.0 { continue; }
        s -= q * log(q);
    }
    if base > 0.0 { s = s / log(base); }
    return s, OK;
}

// Signature: relative_entropy(pk, qk, base) -> value, code
// The relative entropy of two sets of probabilities. Each set is scaled
// so that it sums to 1. A zero value of qk with a positive value of pk
// gives ERR_UNDEFINED, as in SciPy, which gives inf. Vectors of another
// length give ERR_LENGTH.
fn relative_entropy(pk: const f64[:], qk: const f64[:], base: f64) -> f64, i64 {
    let c = check_probabilities_list(pk);
    if c != 0 { return nan(), c; }
    let qc = check_probabilities_list(qk);
    if qc != 0 { return nan(), qc; }
    let bc = check_base(base);
    if bc != 0 { return nan(), bc; }
    var pacc = acc_new();
    var qacc = acc_new();
    for p in pk { acc_add(pacc, p); }
    for q in qk { acc_add(qacc, q); }
    if pacc.fault || qacc.fault { return nan(), ERR_UNDEFINED; }
    if pk.len != qk.len { return nan(), ERR_LENGTH; }
    let pt = acc_total(pacc);
    let qt = acc_total(qacc);
    if pt <= 0.0 || qt <= 0.0 { return nan(), ERR_UNDEFINED; }
    var s: f64 = 0.0;
    for i in 0..pk.len {
        let p = pk[i] / pt;
        let q = qk[i] / qt;
        if p == 0.0 { continue; }
        if q == 0.0 { return nan(), ERR_UNDEFINED; }
        s += p * log(p / q);
    }
    if is_bad(s) { return nan(), ERR_OVERFLOW; }
    if base > 0.0 { s = s / log(base); }
    return s, OK;
}

// The check of a list of probabilities. A negative value gives the
// domain fault. An empty list gives the empty fault.
fn check_probabilities_list(ps: const f64[:]) -> i64 {
    let c = check_sample(ps, 1);
    if c != OK { return c; }
    for p in ps {
        if p < 0.0 { return ERR_DOMAIN; }
    }
    return OK;
}

// The check of a base of a logarithm. A base of 0 means the natural
// logarithm. A base below 0 gives the domain fault, and a base of 1
// gives the domain fault, as in SciPy, which gives inf.
fn check_base(base: f64) -> i64 {
    let c = check_scalar(base);
    if c != OK { return c; }
    if base < 0.0 { return ERR_DOMAIN; }
    if base == 1.0 { return ERR_DOMAIN; }
    return OK;
}

// Signature: kstat(xs, order) -> value, code
// The k-statistic of a sample for the order 1 to 4. It is the unbiased
// estimator of the cumulant of that order. The values are centered
// first, so the sum does not lose digits. An order outside 1 to 4 gives
// ERR_DOMAIN. Fewer values than the order gives ERR_TOO_SHORT.
fn kstat(xs: const f64[:], order: i64) -> f64, i64 {
    let c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    if order < 1 || order > 4 { return nan(), ERR_DOMAIN; }
    let n = xs.len;
    if n < order { return nan(), ERR_TOO_SHORT; }
    if order == 1 {
        let m, c2 = exact_mean(xs);
        if c2 != 0 { return nan(), c2; }
        return m, OK;
    }
    let s2, s3, s4, t, c2 = en_centered_sums(xs);
    if c2 != 0 { return nan(), c2; }
    let nf = n as f64;
    if order == 2 {
        // The unbiased variance.
        let v = s2 * t * t / (nf - 1.0);
        if is_bad(v) { return nan(), ERR_OVERFLOW; }
        return v, OK;
    }
    if order == 3 {
        // The unbiased estimator of the third cumulant.
        let v = nf * s3 * t * t * t / ((nf - 1.0) * (nf - 2.0));
        if is_bad(v) { return nan(), ERR_OVERFLOW; }
        return v, OK;
    }
    if n < 4 { return nan(), ERR_TOO_SHORT; }
    // The unbiased estimator of the fourth cumulant.
    let c4 = s4 * t * t * t * t;
    let c2v = s2 * t * t;
    let num = nf * nf * (nf + 1.0) * c4 - 3.0 * nf * (nf - 1.0) * c2v * c2v;
    let v = num / (nf * (nf - 1.0) * (nf - 2.0) * (nf - 3.0));
    if is_bad(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

// The sums of the powers 2, 3 and 4 of the deviations from the mean,
// scaled by the power of two t. A constant sample gives zeros.
fn en_centered_sums(xs: const f64[:]) -> f64, f64, f64, f64, i64 {
    let n = xs.len;
    let m, c = exact_mean(xs);
    if c != 0 { return nan(), nan(), nan(), 1.0, c; }
    var top: f64 = 0.0;
    for x in xs {
        let d = abs(x - m);
        if is_inf(d) { return nan(), nan(), nan(), 1.0, ERR_OVERFLOW; }
        if d > top { top = d; }
    }
    if top == 0.0 { return 0.0, 0.0, 0.0, 1.0, OK; }
    let t = pow2_floor(top);
    var a2 = acc_new();
    var a3 = acc_new();
    var a4 = acc_new();
    for x in xs {
        let d = (x - m) / t;
        let d2 = d * d;
        acc_add(a2, d2);
        acc_add(a3, d2 * d);
        acc_add(a4, d2 * d2);
    }
    if a2.fault || a3.fault || a4.fault { return nan(), nan(), nan(), 1.0, ERR_UNDEFINED; }
    return acc_total(a2), acc_total(a3), acc_total(a4), t, OK;
}

// Signature: kstatvar(xs, order) -> value, code
// The unbiased estimator of the variance of the k-statistic. The order
// is 1 or 2, as in SciPy. Another order gives ERR_DOMAIN.
fn kstatvar(xs: const f64[:], order: i64) -> f64, i64 {
    let c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    if order < 1 || order > 2 { return nan(), ERR_DOMAIN; }
    let n = xs.len;
    if n < order { return nan(), ERR_TOO_SHORT; }
    let k2, c2 = kstat(xs, 2);
    if c2 != 0 { return nan(), c2; }
    let nf = n as f64;
    if order == 1 { return k2 / nf, OK; }
    if n < 4 { return nan(), ERR_TOO_SHORT; }
    let k4, c3 = kstat(xs, 4);
    if c3 != 0 { return nan(), c3; }
    return (2.0 * nf * k2 * k2 + (nf - 1.0) * k4) / (nf * (nf + 1.0)), OK;
}

// Signature: differential_entropy(xs, method, window, base) -> value, code
// The differential entropy of a sample. The method is one of the four
// estimators of SciPy, or ENTROPY_AUTO, which picks the estimator from
// the size of the sample. A window of 0 uses the default window, which
// is floor(sqrt(n) + 0.5). The window must satisfy 2 * window < n, else
// ERR_DOMAIN. A base of 0 gives the natural logarithm.
fn differential_entropy(xs: const f64[:], method: i64, window: i64, base: f64) -> f64, i64 {
    let c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    if method < ENTROPY_AUTO || method > ENTROPY_CORREA { return nan(), ERR_DOMAIN; }
    let bc = check_base(base);
    if bc != 0 { return nan(), bc; }
    let n = xs.len;
    var m = window;
    if m == 0 { m = (sqrt(n as f64) + 0.5) as i64; }
    if m < 1 { return nan(), ERR_DOMAIN; }
    // The check of the window comes before the doubling, so a huge
    // window never wraps the product.
    if m >= n { return nan(), ERR_DOMAIN; }
    if 2 * m >= n { return nan(), ERR_DOMAIN; }
    var use = method;
    if use == ENTROPY_AUTO {
        if n <= ENTROPY_SMALL { use = ENTROPY_VAN_ES; }
        else if n <= ENTROPY_MEDIUM { use = ENTROPY_EBRAHIMI; }
        else { use = ENTROPY_VASICEK; }
    }
    let s = sort_copy(xs);
    var r: f64 = 0.0;
    var c2 = OK;
    if use == ENTROPY_VASICEK {
        r, c2 = en_vasicek(s[..], m);
    } else if use == ENTROPY_VAN_ES {
        r, c2 = en_van_es(s[..], m);
    } else if use == ENTROPY_EBRAHIMI {
        r, c2 = en_ebrahimi(s[..], m);
    } else {
        r, c2 = en_correa(s[..], m);
    }
    if c2 != OK { return nan(), c2; }
    if is_nan(r) { return nan(), ERR_UNDEFINED; }
    if is_inf(r) { return nan(), ERR_UNDEFINED; }
    if base > 0.0 { r = r / log(base); }
    return r, OK;
}

// The padded sample of the estimators. The list holds the window copies
// of the smallest value, the values and the window copies of the
// largest value.
fn en_pad(s: const f64[:], m: i64) -> f64[>..] {
    var p: f64[>..] = [];
    var i: i64 = 0;
    while i < m { p.push(s[0]); i += 1; }
    for x in s { p.push(x); }
    i = 0;
    let n = s.len;
    while i < m { p.push(s[n - 1]); i += 1; }
    return p;
}

// The Vasicek estimator of the differential entropy.
fn en_vasicek(s: const f64[:], m: i64) -> f64, i64 {
    let n = s.len;
    let p = en_pad(s, m);
    var acc = acc_new();
    var i: i64 = 0;
    while i < n {
        let d = p[i + 2 * m] - p[i];
        if d <= 0.0 { return nan(), ERR_UNDEFINED; }
        acc_add(acc, log((n) / (2.0 * m) * d));
        i += 1;
    }
    if acc.fault { return nan(), ERR_UNDEFINED; }
    return acc_total(acc) / (n), OK;
}

// The van Es estimator of the differential entropy. The sample is not
// padded here.
fn en_van_es(s: const f64[:], m: i64) -> f64, i64 {
    let n = s.len;
    var acc = acc_new();
    var i: i64 = 0;
    while i + m < n {
        let d = s[i + m] - s[i];
        if d <= 0.0 { return nan(), ERR_UNDEFINED; }
        acc_add(acc, log(((n + 1) as f64) / (m) * d));
        i += 1;
    }
    if acc.fault { return nan(), ERR_UNDEFINED; }
    let term1 = acc_total(acc) / ((n - m));
    var harmonic: f64 = 0.0;
    var k = m;
    while k <= n {
        harmonic += 1.0 / (k);
        k += 1;
    }
    return term1 + harmonic + log(m as f64) - log((n + 1) as f64), OK;
}

// The Ebrahimi estimator of the differential entropy.
fn en_ebrahimi(s: const f64[:], m: i64) -> f64, i64 {
    let n = s.len;
    let mf = m as f64;
    let p = en_pad(s, m);
    var acc = acc_new();
    var i: i64 = 0;
    while i < n {
        let d = p[i + 2 * m] - p[i];
        if d <= 0.0 { return nan(), ERR_UNDEFINED; }
        // The index of the difference in the padded list runs from 1 to n.
        let idx = (i + 1) as f64;
        var ci: f64 = 2.0;
        if idx <= mf { ci = 1.0 + (idx - 1.0) / mf; }
        if idx >= (n - m + 1) { ci = 1.0 + ((n) - idx) / mf; }
        acc_add(acc, log((n) * d / (ci * mf)));
        i += 1;
    }
    if acc.fault { return nan(), ERR_UNDEFINED; }
    return acc_total(acc) / (n), OK;
}

// The Correa estimator of the differential entropy. The sample is
// padded with the window on each end.
fn en_correa(s: const f64[:], m: i64) -> f64, i64 {
    let n = s.len;
    let p = en_pad(s, m);
    var acc = acc_new();
    var k: i64 = 1;
    while k <= n {
        // The window of 2m + 1 points around the index.
        var sum: f64 = 0.0;
        var num: f64 = 0.0;
        var den: f64 = 0.0;
        var j: i64 = -m;
        while j <= m {
            let q0 = k + j + m - 1;
            let v = p[q0];
            sum += v;
            j += 1;
        }
        let ibar = sum / ((2 * m + 1));
        j = -m;
        while j <= m {
            let q0 = k + j + m - 1;
            let diff = p[q0] - ibar;
            num += diff * (j);
            den += diff * diff;
            j += 1;
        }
        if den == 0.0 { return nan(), ERR_UNDEFINED; }
        let v = num / ((n) * den);
        if v <= 0.0 { return nan(), ERR_UNDEFINED; }
        acc_add(acc, log(v));
        k += 1;
    }
    if acc.fault { return nan(), ERR_UNDEFINED; }
    return -(acc_total(acc) / (n)), OK;
}
