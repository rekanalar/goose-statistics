// The permutation test, the Monte Carlo test and the power of the
// package.
//
// Every call returns the result first and a fault code last. A fault
// gives NaN as the result.
namespace statistics;

import std;
import math;

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
