// The hypothesis tests. Each call returns a TestResult and then a fault
// code. A fault gives NaN in every field of the result.
//
// The file holds the t tests, the analysis of variance, the interval of
// the mean, the count tests and the rank tests.
namespace statistics;

import std;
import math;

// The result of a test. A field that the test does not have holds NaN.
struct TestResult { statistic: f64, pvalue: f64, df: f64, df2: f64 }

// The result of a call that faulted.
fn fault_result(code: i64) -> TestResult, i64 {
    return TestResult { statistic: nan(), pvalue: nan(), df: nan(), df2: nan() }, code;
}

fn good_result(statistic: f64, pvalue: f64, df: f64, df2: f64) -> TestResult, i64 {
    // A NaN p-value never leaves with the code OK.
    if is_nan(pvalue) { return fault_result(ERR_UNDEFINED); }
    return TestResult { statistic: statistic, pvalue: clip01(pvalue), df: df, df2: df2 }, 0;
}

// The check of an alternative.
fn check_alt(alt: i64) -> i64 {
    if alt < 0 || alt > 2 { return ERR_DOMAIN; }
    return 0;
}

// The greatest common divisor of two counts. The older Goose releases
// hold no such call in the standard library.
fn gcd_counts(a: i64, b: i64) -> i64 {
    var x = a;
    var y = b;
    while y != 0 {
        let t = x % y;
        x = y;
        y = t;
    }
    return x;
}

// The p-value of a t statistic. The statistic may be infinite.
fn t_pvalue(t: f64, df: f64, alt: i64) -> f64, i64 {
    if is_inf(t) {
        if alt == ALT_TWO_SIDED { return 0.0, 0; }
        if alt == ALT_LESS { if t < 0.0 { return 0.0, 0; } return 1.0, 0; }
        if t > 0.0 { return 0.0, 0; }
        return 1.0, 0;
    }
    if alt == ALT_TWO_SIDED {
        let p, c = tdist_sf(abs(t), df);
        return 2.0 * p, c;
    }
    if alt == ALT_LESS { return tdist_cdf(t, df); }
    return tdist_sf(t, df);
}

// The t statistic of a difference over a standard error. A zero error
// gives the limit, or the undefined fault when the difference is zero.
fn t_statistic(d: f64, se: f64) -> f64, i64 {
    if is_inf(d) { return nan(), 8; }
    if se == 0.0 || is_nan(se) {
        if d == 0.0 { return nan(), 7; }
        if d > 0.0 { return (1.0 / 0.0), 0; }
        return -(1.0 / 0.0), 0;
    }
    return d / se, 0;
}

fn finish_t(t: f64, df: f64, alt: i64) -> TestResult, i64 {
    let p, c = t_pvalue(t, df, alt);
    if c != 0 { return fault_result(c); }
    return good_result(t, p, df, nan());
}

// ---- the t tests ------------------------------------------------------

// The t test of one sample against a mean.
fn ttest_1samp(xs: const f64[:], popmean: f64 = 0.0, alt: i64 = ALT_TWO_SIDED) -> TestResult, i64 {
    let c = check_sample(xs, 2);
    if c != 0 { return fault_result(c); }
    let ca = check_scalar(popmean);
    if ca != 0 { return fault_result(ca); }
    if check_alt(alt) != 0 { return fault_result(ERR_DOMAIN); }
    return ttest_core(xs, popmean, alt);
}

fn ttest_core(xs: const f64[:], popmean: f64, alt: i64) -> TestResult, i64 {
    let m, cm = mean(xs);
    if cm != 0 { return fault_result(cm); }
    let v, cv = variance_core(xs, 1);
    if cv != 0 { return fault_result(cv); }
    let n = xs.len as f64;
    let t, ct = t_statistic(m - popmean, sqrt(v / n));
    if ct != 0 { return fault_result(ct); }
    return finish_t(t, n - 1.0, alt);
}

// The t test of two independent samples. The Student form pools the
// variances. The Welch form keeps them apart.
fn ttest_ind(xs: const f64[:], ys: const f64[:], alt: i64 = ALT_TWO_SIDED, equal_var: bool = true) -> TestResult, i64 {
    var c = check_sample(xs, 2);
    if c != 0 { return fault_result(c); }
    c = check_sample(ys, 2);
    if c != 0 { return fault_result(c); }
    if check_alt(alt) != 0 { return fault_result(ERR_DOMAIN); }
    let m1, c1 = mean(xs);
    let m2, c2 = mean(ys);
    let v1, c3 = variance_core(xs, 1);
    let v2, c4 = variance_core(ys, 1);
    if c1 != 0 { return fault_result(c1); }
    if c2 != 0 { return fault_result(c2); }
    if c3 != 0 { return fault_result(c3); }
    if c4 != 0 { return fault_result(c4); }
    let n1 = xs.len as f64;
    let n2 = ys.len as f64;
    if equal_var {
        let df = n1 + n2 - 2.0;
        let pooled = v1 * ((n1 - 1.0) / df) + v2 * ((n2 - 1.0) / df);
        let t, ct = t_statistic(m1 - m2, sqrt(pooled * (1.0 / n1 + 1.0 / n2)));
        if ct != 0 { return fault_result(ct); }
        return finish_t(t, df, alt);
    }
    let a = v1 / n1;
    let b = v2 / n2;
    let t, ct = t_statistic(m1 - m2, sqrt(a + b));
    if ct != 0 { return fault_result(ct); }
    var df = n1 + n2 - 2.0;
    let s = a + b;
    if s > 0.0 {
        let ra = a / s;
        let rb = b / s;
        df = 1.0 / (ra * ra / (n1 - 1.0) + rb * rb / (n2 - 1.0));
    }
    let lo = min(n1, n2) - 1.0;
    if df < lo { df = lo; }
    if df > n1 + n2 - 2.0 { df = n1 + n2 - 2.0; }
    return finish_t(t, df, alt);
}

// The t test of two paired samples: the one-sample test of the
// differences against 0.
fn ttest_rel(xs: const f64[:], ys: const f64[:], alt: i64 = ALT_TWO_SIDED) -> TestResult, i64 {
    let c = check_pair(xs, ys, 2);
    if c != 0 { return fault_result(c); }
    if check_alt(alt) != 0 { return fault_result(ERR_DOMAIN); }
    var d: f64[>..] = [];
    for i in 0..xs.len {
        let x = xs[i] - ys[i];
        if is_inf(x) { return fault_result(ERR_OVERFLOW); }
        d.push(x);
    }
    return ttest_core(d[..], 0.0, alt);
}

// ---- the analysis of variance -------------------------------------------

// The check of a set of groups in one buffer. The sizes array gives the
// length of each group. The call returns the code.
fn check_groups(groups: const f64[:], sizes: const i64[:]) -> i64 {
    return check_groups_n(groups, sizes, 2, 1);
}

// The one-way analysis of variance of two or more groups. The groups
// sit one after the other in one buffer.
fn f_oneway(groups: const f64[:], sizes: const i64[:]) -> TestResult, i64 {
    let c = check_groups(groups, sizes);
    if c != 0 { return fault_result(c); }
    let k = sizes.len;
    let total = groups.len;
    if total <= k { return fault_result(ERR_TOO_SHORT); }
    let grand, cg = mean(groups);
    if cg != 0 { return fault_result(cg); }
    var ssb: f64 = 0.0;
    var ssw: f64 = 0.0;
    var at: i64 = 0;
    for g in 0..k {
        let n = sizes[g];
        let part = groups[at..at + n];
        at += n;
        let m, cm = mean(part);
        if cm != 0 { return fault_result(cm); }
        let v, cv = variance_core(part, 0);
        if cv != 0 { return fault_result(cv); }
        let d = m - grand;
        ssb += d * d * n;
        ssw += v * n;
    }
    if is_inf(ssb) || is_inf(ssw) { return fault_result(ERR_OVERFLOW); }
    let df1 = (k - 1) as f64;
    let df2 = (total - k) as f64;
    if ssw == 0.0 {
        if ssb == 0.0 { return fault_result(ERR_UNDEFINED); }
        return good_result((1.0 / 0.0), 0.0, df1, df2);
    }
    let f = (ssb / df1) / (ssw / df2);
    if is_inf(f) { return good_result(f, 0.0, df1, df2); }
    let p, cp = fdist_sf(f, df1, df2);
    if cp != 0 { return fault_result(cp); }
    return good_result(f, p, df1, df2);
}

// ---- the interval of the mean -----------------------------------------

// The interval of the mean of a sample for a confidence from 0 to 1.
// The call gives the lower bound, the upper bound and the code.
fn mean_confidence_interval(xs: const f64[:], confidence: f64 = 0.95) -> f64, f64, i64 {
    var c = check_sample(xs, 2);
    if c != 0 { return nan(), nan(), c; }
    c = check_level(confidence);
    if c != 0 { return nan(), nan(), c; }
    let m, cm = mean(xs);
    if cm != 0 { return nan(), nan(), cm; }
    let se, cs = sem(xs);
    if cs != 0 { return nan(), nan(), cs; }
    let q, cq = tdist_ppf((1.0 - confidence) / 2.0, (xs.len - 1));
    if cq != 0 { return nan(), nan(), cq; }
    let h = -q * se;
    let lo = m - h;
    let hi = m + h;
    if is_inf(h) || is_inf(lo) || is_inf(hi) { return nan(), nan(), 8; }
    return lo, hi, 0;
}

// ---- the count tests --------------------------------------------------

// The p-value of a chi-square statistic: the upper tail.
fn chi_result(stat: f64, df: f64) -> TestResult, i64 {
    let p, c = chi2_sf(stat, df);
    if c != 0 { return fault_result(c); }
    return good_result(stat, p, df, nan());
}

// The chi-square test of fit. An empty expected list means equal
// expected counts. The count ddof lowers the degrees of freedom.
fn chisquare(obs: const f64[:], expected: const f64[:], ddof: i64 = 0) -> TestResult, i64 {
    let c = check_sample(obs, 2);
    if c != 0 { return fault_result(c); }
    if expected.len > 0 {
        let ce = check_sample(expected, 0);
        if ce != 0 { return fault_result(ce); }
        if expected.len != obs.len { return fault_result(ERR_LENGTH); }
    }
    let k = obs.len;
    if ddof < 0 || ddof > k - 2 { return fault_result(ERR_DOMAIN); }
    var total = acc_new();
    for o in obs {
        if o < 0.0 || o > MAX_COUNT { return fault_result(ERR_DOMAIN); }
        acc_add(total, o);
    }
    let sum_obs = acc_total(total);
    var sum_exp = sum_obs;
    if expected.len > 0 {
        var te = acc_new();
        for e in expected {
            if e < 0.0 || e > MAX_COUNT { return fault_result(ERR_DOMAIN); }
            acc_add(te, e);
        }
        sum_exp = acc_total(te);
        if abs(sum_obs - sum_exp) > 1e-8 * max(sum_obs, sum_exp) { return fault_result(ERR_DOMAIN); }
    }
    if sum_obs == 0.0 { return fault_result(ERR_UNDEFINED); }
    var stat = acc_new();
    for i in 0..k {
        var e = sum_obs / k;
        if expected.len > 0 { e = expected[i]; }
        if e == 0.0 { return fault_result(ERR_UNDEFINED); }
        let d = obs[i] - e;
        acc_add(stat, d * d / e);
    }
    return chi_result(acc_total(stat), (k - 1 - ddof));
}

// The chi-square test of a table of counts, stored row by row. A table
// with two rows and two columns gets the continuity correction of Yates
// when correction is true.
fn chi2_contingency(table: const f64[:], rows: i64, cols: i64, correction: bool = true) -> TestResult, i64 {
    let c = check_sample(table, 1);
    if c != 0 { return fault_result(c); }
    if rows < 2 || cols < 2 { return fault_result(ERR_DOMAIN); }
    if rows > table.len || cols > table.len || rows * cols != table.len { return fault_result(ERR_LENGTH); }
    var rs: f64[>..] = [];
    var cs: f64[>..] = [];
    for i in 0..rows { rs.push(0.0); }
    for j in 0..cols { cs.push(0.0); }
    var total = acc_new();
    for i in 0..rows {
        var ra = acc_new();
        for j in 0..cols {
            let o = table[i * cols + j];
            if o < 0.0 || o > MAX_COUNT { return fault_result(ERR_DOMAIN); }
            acc_add(ra, o);
        }
        rs[i] = acc_total(ra);
        acc_add(total, rs[i]);
    }
    for j in 0..cols {
        var ca = acc_new();
        for i in 0..rows { acc_add(ca, table[i * cols + j]); }
        cs[j] = acc_total(ca);
    }
    let n = acc_total(total);
    if n == 0.0 { return fault_result(ERR_UNDEFINED); }
    for i in 0..rows { if rs[i] == 0.0 { return fault_result(ERR_UNDEFINED); } }
    for j in 0..cols { if cs[j] == 0.0 { return fault_result(ERR_UNDEFINED); } }
    let yates = correction && rows == 2 && cols == 2;
    var stat = acc_new();
    for i in 0..rows {
        for j in 0..cols {
            let e = rs[i] * cs[j] / n;
            var o = table[i * cols + j];
            if yates {
                let diff = e - o;
                let mag = min(0.5, abs(diff));
                if diff > 0.0 { o += mag; }
                if diff < 0.0 { o -= mag; }
            }
            let d = o - e;
            acc_add(stat, d * d / e);
        }
    }
    return chi_result(acc_total(stat), ((rows - 1) * (cols - 1)));
}

// The mass at x of the binomial family (kind 0) or of the
// hypergeometric family (kind 1). A count outside the support has mass 0.
fn pmf_at(kind: i64, x: f64, a1: f64, a2: f64, a3: f64) -> f64, i64 {
    if x < 0.0 { return 0.0, 0; }
    if kind == 0 { return binom_pmf(x, a1, a2); }
    return hypergeom_pmf(x, a1, a2, a3);
}

// The binary search of the two-sided exact tests. The mass (or its
// negative) rises with x from lo to hi. The call gives the last x whose
// value is at most d.
fn mass_search(kind: i64, negate: bool, d: f64, lo0: f64, hi0: f64, a1: f64, a2: f64, a3: f64) -> f64, i64 {
    var lo = lo0;
    var hi = hi0;
    var steps: i64 = 0;
    while lo < hi {
        steps += 1;
        if steps > MAX_ITER_SEARCH { return nan(), 8; }
        let mid = lo + floor((hi - lo) / 2.0);
        var v, c = pmf_at(kind, mid, a1, a2, a3);
        if c != 0 { return nan(), c; }
        if negate { v = -v; }
        if v < d { lo = mid + 1.0; }
        else if v > d { hi = mid - 1.0; }
        else { lo = mid; hi = mid; }
    }
    var v, c = pmf_at(kind, lo, a1, a2, a3);
    if c != 0 { return nan(), c; }
    if negate { v = -v; }
    if v <= d { return lo, 0; }
    return lo - 1.0, 0;
}

// The lower and the upper tail of the two families. A count below 0 has
// no mass below it and all the mass above it.
fn tail_cdf(kind: i64, x: f64, a1: f64, a2: f64, a3: f64) -> f64, i64 {
    if x < 0.0 { return 0.0, 0; }
    if kind == 0 { return binom_cdf(x, a1, a2); }
    return hypergeom_cdf(x, a1, a2, a3);
}

fn tail_sf(kind: i64, x: f64, a1: f64, a2: f64, a3: f64) -> f64, i64 {
    if x < 0.0 { return 1.0, 0; }
    if kind == 0 { return binom_sf(x, a1, a2); }
    return hypergeom_sf(x, a1, a2, a3);
}

// The exact binomial test of k successes in n trials against the
// probability p. The statistic is k / n. The two-sided p-value sums the
// masses at most the mass of k, with a relative slack of 1e-7.
fn binomtest(k: f64, n: f64, p: f64 = 0.5, alt: i64 = ALT_TWO_SIDED) -> TestResult, i64 {
    var c = check_scalar(k);
    if c == 0 { c = check_scalar(n); }
    if c == 0 { c = check_scalar(p); }
    if c != 0 { return fault_result(c); }
    if chk_count(k) != 0 || chk_count(n) != 0 || n < 1.0 || k > n { return fault_result(ERR_DOMAIN); }
    if p < 0.0 || p > 1.0 || check_alt(alt) != 0 { return fault_result(ERR_DOMAIN); }
    let stat = k / n;
    if alt == ALT_LESS {
        let v, cv = binom_cdf(k, n, p);
        if cv != 0 { return fault_result(cv); }
        return good_result(stat, v, nan(), nan());
    }
    if alt == ALT_GREATER {
        let v, cv = tail_sf(0, k - 1.0, n, p, 0.0);
        if cv != 0 { return fault_result(cv); }
        return good_result(stat, v, nan(), nan());
    }
    let d, cd = binom_pmf(k, n, p);
    if cd != 0 { return fault_result(cd); }
    let target = d * (1.0 + 1e-7);
    var pval: f64 = 0.0;
    if k < p * n {
        let ix, ci = mass_search(0, true, -target, ceil(p * n), n, n, p, 0.0);
        if ci != 0 { return fault_result(ci); }
        let mx, cm = pmf_at(0, ix, n, p, 0.0);
        if cm != 0 { return fault_result(cm); }
        var cut = ix;
        if target == mx { cut = ix - 1.0; }
        let lo, c1 = binom_cdf(k, n, p);
        let hi, c2 = tail_sf(0, cut, n, p, 0.0);
        if c1 != 0 { return fault_result(c1); }
        if c2 != 0 { return fault_result(c2); }
        pval = lo + hi;
    } else {
        let ix, ci = mass_search(0, false, target, 0.0, floor(p * n), n, p, 0.0);
        if ci != 0 { return fault_result(ci); }
        let lo, c1 = tail_cdf(0, ix, n, p, 0.0);
        let hi, c2 = tail_sf(0, k - 1.0, n, p, 0.0);
        if c1 != 0 { return fault_result(c1); }
        if c2 != 0 { return fault_result(c2); }
        pval = lo + hi;
    }
    return good_result(stat, pval, nan(), nan());
}

// The Fisher exact test of a table with two rows and two columns. The
// four counts come row by row. The statistic is the odds ratio.
fn fisher_exact(table: const f64[:], alt: i64 = ALT_TWO_SIDED) -> TestResult, i64 {
    var c = check_sample(table, 1);
    if c != 0 { return fault_result(c); }
    if table.len != 4 { return fault_result(ERR_LENGTH); }
    for x in table { if chk_count(x) != 0 { return fault_result(ERR_DOMAIN); } }
    if check_alt(alt) != 0 { return fault_result(ERR_DOMAIN); }
    let a = table[0];
    let b = table[1];
    let cc = table[2];
    let d = table[3];
    if a + b == 0.0 || cc + d == 0.0 || a + cc == 0.0 || b + d == 0.0 {
        return good_result(nan(), 1.0, nan(), nan());
    }
    var odds = 1.0 / 0.0;
    if b > 0.0 && cc > 0.0 { odds = a * d / (b * cc); }
    let n1 = a + b;
    let n2 = cc + d;
    let n = a + cc;
    let m = n1 + n2;
    if m > MAX_COUNT { return fault_result(ERR_DOMAIN); }
    var pval: f64 = 0.0;
    if alt == ALT_LESS {
        let v, cv = hypergeom_cdf(a, m, n1, n);
        if cv != 0 { return fault_result(cv); }
        pval = v;
    } else if alt == ALT_GREATER {
        let v, cv = hypergeom_cdf(b, m, n1, b + d);
        if cv != 0 { return fault_result(cv); }
        pval = v;
    } else {
        let r, cr = fisher_two_sided(a, m, n1, n);
        if cr != 0 { return fault_result(cr); }
        pval = r;
    }
    return good_result(odds, pval, nan(), nan());
}

fn fisher_two_sided(a: f64, m: f64, n1: f64, n: f64) -> f64, i64 {
    let mode = floor((n + 1.0) * (n1 + 1.0) / (m + 2.0));
    let pexact, c1 = hypergeom_pmf(a, m, n1, n);
    let pmode, c2 = hypergeom_pmf(mode, m, n1, n);
    if c1 != 0 { return nan(), c1; }
    if c2 != 0 { return nan(), c2; }
    let gamma = 1.0 + 1e-14;
    if abs(pexact - pmode) / max(pexact, pmode) <= 1e-14 { return 1.0, 0; }
    if a < mode {
        let plower, c3 = hypergeom_cdf(a, m, n1, n);
        if c3 != 0 { return nan(), c3; }
        let top, c4 = pmf_at(1, n, m, n1, n);
        if c4 != 0 { return nan(), c4; }
        if top > pexact * gamma { return plower, 0; }
        let guess, c5 = mass_search(1, true, -pexact * gamma, mode, n, m, n1, n);
        if c5 != 0 { return nan(), c5; }
        let up, c6 = tail_sf(1, guess, m, n1, n);
        if c6 != 0 { return nan(), c6; }
        return min(plower + up, 1.0), 0;
    }
    let pupper, c3 = tail_sf(1, a - 1.0, m, n1, n);
    if c3 != 0 { return nan(), c3; }
    let bottom, c4 = pmf_at(1, 0.0, m, n1, n);
    if c4 != 0 { return nan(), c4; }
    if bottom > pexact * gamma { return pupper, 0; }
    let guess, c5 = mass_search(1, false, pexact * gamma, 0.0, mode, m, n1, n);
    if c5 != 0 { return nan(), c5; }
    let lo, c6 = tail_cdf(1, guess, m, n1, n);
    if c6 != 0 { return nan(), c6; }
    return min(pupper + lo, 1.0), 0;
}

// ---- the rank tests ----------------------------------------------------

// The counts of the statistic U for samples of m and n values, m <= n.
// Entry u holds the count of the arrangements with U = u. The table is
// the Gaussian binomial coefficient, built one factor at a time.
fn mw_counts(out: f64[>..]&, m: i64, n: i64) {
    out.resize(0);
    let cells = m * n;
    for u in 0..cells + 1 { out.push(0.0); }
    out[0] = 1.0;
    for i in 1..m + 1 {
        for u in i..cells + 1 { out[u] += out[u - i]; }
        let sh = n + i;
        var u = cells;
        while u >= sh {
            out[u] -= out[u - sh];
            u -= 1;
        }
    }
}

// The exact upper tail P(U >= k) for samples of m and n values.
fn mw_exact_sf(k: i64, m: i64, n: i64) -> f64 {
    var counts: f64[>..] = [];
    mw_counts(counts, m, n);
    var total: f64 = 1.0;
    for i in 1..m + 1 { total = total * (n + i) / i; }
    var tail = acc_new();
    var u = counts.len - 1;
    while u >= k && u >= 0 {
        acc_add(tail, counts[u]);
        u -= 1;
    }
    return acc_total(tail) / total;
}

// The Mann-Whitney U test. The statistic is U of the first sample.
fn mannwhitneyu(xs: const f64[:], ys: const f64[:], continuity: bool = true, alt: i64 = ALT_TWO_SIDED, method: i64 = MW_AUTO) -> TestResult, i64 {
    if xs.len > MAX_INPUT || ys.len > MAX_INPUT { return fault_result(ERR_TOO_LONG); }
    var c = check_sample(xs, 1);
    if c != 0 { return fault_result(c); }
    c = check_sample(ys, 1);
    if c != 0 { return fault_result(c); }
    if check_alt(alt) != 0 || method < 0 || method > MW_ASYMPTOTIC { return fault_result(ERR_DOMAIN); }
    let n1 = xs.len;
    let n2 = ys.len;
    var pooled: f64[>..] = [];
    for x in xs { pooled.push(x); }
    for y in ys { pooled.push(y); }
    var rk: f64[>..] = [];
    let cnt, tt, cr = rankdata_ties(rk, pooled, TIE_AVERAGE);
    if cr != 0 { return fault_result(cr); }
    var r1: f64 = 0.0;
    for i in 0..n1 { r1 += rk[i]; }
    let f1 = n1 as f64;
    let f2 = n2 as f64;
    let u1 = r1 - f1 * (f1 + 1.0) / 2.0;
    let u2 = f1 * f2 - u1;
    let mm = min(n1, n2);
    let nn = max(n1, n2);
    let cells = mm * nn;
    var exact = false;
    if method == MW_EXACT {
        if tt > 0.0 { return fault_result(ERR_DOMAIN); }
        if cells > MAX_MW_CELLS { return fault_result(ERR_OVERFLOW); }
        exact = true;
    }
    if method == MW_AUTO {
        exact = tt == 0.0 && mm <= MAX_EXACT_MW && cells <= MAX_MW_CELLS;
    }
    var u = u1;
    var factor: f64 = 1.0;
    if alt == ALT_LESS { u = u2; }
    if alt == ALT_TWO_SIDED {
        u = max(u1, u2);
        factor = 2.0;
    }
    if exact {
        let p = factor * mw_exact_sf(round(u) as i64, mm, nn);
        return good_result(u1, p, nan(), nan());
    }
    let n = f1 + f2;
    let s = sqrt(f1 * f2 / 12.0 * ((n + 1.0) - tt / (n * (n - 1.0))));
    if s == 0.0 { return fault_result(ERR_UNDEFINED); }
    var num = u - f1 * f2 / 2.0;
    if continuity { num -= 0.5; }
    let sf, cs = norm_sf(num / s);
    if cs != 0 { return fault_result(cs); }
    return good_result(u1, factor * sf, nan(), nan());
}
