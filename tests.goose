// The hypothesis tests. Each call returns a TestResult and then a fault
// code. A fault gives NaN in every field of the result.
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
    let c = check_sample(groups, 1);
    if c != 0 { return c; }
    if sizes.len < 2 { return ERR_TOO_SHORT; }
    var total: i64 = 0;
    for s in sizes {
        if s < 0 { return ERR_DOMAIN; }
        if s == 0 { return ERR_TOO_SHORT; }
        total += s;
        if total > groups.len { return ERR_LENGTH; }
    }
    if total != groups.len { return ERR_LENGTH; }
    return 0;
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
    c = check_scalar(confidence);
    if c != 0 { return nan(), nan(), c; }
    if confidence <= 0.0 || confidence >= 1.0 { return nan(), nan(), 5; }
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

// ---- the Wilcoxon test ---------------------------------------------------

// The distribution of the sum of the chosen ranks. The weights are twice
// the ranks, so they are whole numbers. Entry s holds the count of the
// subsets with weight s. The count of at most 50 weights fits in a double.
fn subset_counts(out: f64[>..]&, w: const i64[:]) {
    out.resize(0);
    var total: i64 = 0;
    for x in w { total += x; }
    for s in 0..total + 1 { out.push(0.0); }
    out[0] = 1.0;
    var top: i64 = 0;
    for x in w {
        top += x;
        var s = top;
        while s >= x {
            out[s] += out[s - x];
            s -= 1;
        }
    }
}

// The Wilcoxon test of a list of differences.
fn wilcoxon(diffs: const f64[:], alt: i64 = ALT_TWO_SIDED, zero_method: i64 = ZERO_WILCOX, correction: bool = false) -> TestResult, i64 {
    let c = check_sample(diffs, 1);
    if c != 0 { return fault_result(c); }
    if check_alt(alt) != 0 || zero_method < 0 || zero_method > ZERO_ZSPLIT { return fault_result(ERR_DOMAIN); }
    var n_zero = 0;
    for v in diffs {
        if v == 0.0 { n_zero += 1; }
    }
    if zero_method == ZERO_WILCOX && n_zero > 0 {
        // The zero differences leave the sample. The copy of the sample
        // is made only then, so the common case holds no copy.
        var kept: f64[>..] = [];
        for v in diffs {
            if v != 0.0 { kept.push(v); }
        }
        return wilcoxon_core(kept[..], diffs.len, n_zero, alt, zero_method, correction);
    }
    return wilcoxon_core(diffs, diffs.len, n_zero, alt, zero_method, correction);
}

// The Wilcoxon test of two paired samples.
fn wilcoxon_paired(xs: const f64[:], ys: const f64[:], alt: i64 = ALT_TWO_SIDED, zero_method: i64 = ZERO_WILCOX, correction: bool = false) -> TestResult, i64 {
    let c = check_pair(xs, ys, 1);
    if c != 0 { return fault_result(c); }
    if check_alt(alt) != 0 || zero_method < 0 || zero_method > ZERO_ZSPLIT { return fault_result(ERR_DOMAIN); }
    var d: f64[>..] = [];
    for i in 0..xs.len {
        let v = xs[i] - ys[i];
        if is_inf(v) { return fault_result(ERR_OVERFLOW); }
        d.push(v);
    }
    return wilcoxon(d[..], alt, zero_method, correction);
}

fn wilcoxon_core(d: const f64[:], n_diff: i64, n_zero: i64, alt: i64, zero_method: i64, correction: bool) -> TestResult, i64 {
    if d.len == 0 || n_zero == n_diff { return fault_result(ERR_UNDEFINED); }
    var absd: f64[>..] = [];
    for v in d { absd.push(abs(v)); }
    var r: f64[>..] = [];
    let cnt, t_all, cr = rankdata_ties(r, absd, TIE_AVERAGE);
    if cr != 0 { return fault_result(cr); }
    var rp: f64 = 0.0;
    var rm: f64 = 0.0;
    var rz: f64 = 0.0;
    for i in 0..d.len {
        if d[i] > 0.0 { rp += r[i]; }
        if d[i] < 0.0 { rm += r[i]; }
        if d[i] == 0.0 { rz += r[i]; }
    }
    var fixed: f64 = 0.0;
    if zero_method == ZERO_ZSPLIT {
        rp += rz / 2.0;
        rm += rz / 2.0;
        fixed = rz / 2.0;
    }
    var stat = rp;
    if alt == ALT_TWO_SIDED { stat = min(rp, rm); }
    let count = d.len;
    // The groups of equal absolute differences give the same groups as the
    // ranks. The zero differences form one group of their own. The sample
    // holds them only when the zero method keeps them.
    let keeps_zero = zero_method != ZERO_WILCOX;
    var t_nz = t_all;
    if keeps_zero && n_zero > 0 {
        let z = n_zero as f64;
        t_nz = t_all - (z * z * z - z);
    }
    var ties = t_nz > 0.0;
    var zeros = n_zero > 0;
    var exact = false;
    if n_diff <= MAX_EXACT_WILCOXON && !ties && !zeros { exact = true; }
    if n_diff <= MAX_ENUM_WILCOXON { exact = true; }
    if exact {
        var weights: i64[>..] = [];
        for i in 0..d.len {
            if d[i] != 0.0 { weights.push(round(2.0 * r[i]) as i64); }
        }
        var cts: f64[>..] = [];
        subset_counts(cts, weights);
        var total: f64 = 1.0;
        for i in 0..weights.len { total *= 2.0; }
        let at = round(2.0 * (rp - fixed)) as i64;
        var le = acc_new();
        var ge = acc_new();
        for s in 0..cts.len {
            if s <= at { acc_add(le, cts[s]); }
            if s >= at { acc_add(ge, cts[s]); }
        }
        let pl = acc_total(le) / total;
        let pg = acc_total(ge) / total;
        var p = pg;
        if alt == ALT_LESS { p = pl; }
        if alt == ALT_TWO_SIDED { p = 2.0 * min(pl, pg); }
        return good_result(stat, p, nan(), nan());
    }
    let cf = count as f64;
    var mn = cf * (cf + 1.0) * 0.25;
    var se = cf * (cf + 1.0) * (2.0 * cf + 1.0);
    if zero_method == ZERO_PRATT {
        let z = n_zero as f64;
        mn -= z * (z + 1.0) * 0.25;
        se -= z * (z + 1.0) * (2.0 * z + 1.0);
    }
    // The Pratt method corrects for the ties of the differences that are
    // not zero. The other methods correct for the ties of all differences.
    var t_corr = t_all;
    if zero_method == ZERO_PRATT && keeps_zero { t_corr = t_nz; }
    se -= 0.5 * t_corr;
    if se <= 0.0 { return fault_result(ERR_UNDEFINED); }
    se = sqrt(se / 24.0);
    var corr: f64 = 0.0;
    if correction {
        if alt == ALT_TWO_SIDED {
            if rp - mn > 0.0 { corr = 0.5; }
            if rp - mn < 0.0 { corr = -0.5; }
        }
        if alt == ALT_LESS { corr = -0.5; }
        if alt == ALT_GREATER { corr = 0.5; }
    }
    let z = (rp - mn - corr) / se;
    var p: f64 = 0.0;
    if alt == ALT_TWO_SIDED {
        let sf, cs = norm_sf(abs(z));
        if cs != 0 { return fault_result(cs); }
        p = 2.0 * sf;
    }
    if alt == ALT_GREATER {
        let sf, cs = norm_sf(z);
        if cs != 0 { return fault_result(cs); }
        p = sf;
    }
    if alt == ALT_LESS {
        let cd, cs = norm_cdf(z);
        if cs != 0 { return fault_result(cs); }
        p = cd;
    }
    return good_result(stat, p, nan(), nan());
}

fn kruskal(groups: const f64[:], sizes: const i64[:]) -> TestResult, i64 {
    let c = check_groups(groups, sizes);
    if c != 0 { return fault_result(c); }
    var rk: f64[>..] = [];
    let cnt, tt, cr = rankdata_ties(rk, groups, TIE_AVERAGE);
    if cr != 0 { return fault_result(cr); }
    let nt = groups.len as f64;
    let corr = 1.0 - tt / (nt * nt * nt - nt);
    if corr <= 0.0 { return fault_result(ERR_UNDEFINED); }
    var h: f64 = 0.0;
    var at: i64 = 0;
    for g in sizes {
        var rs: f64 = 0.0;
        for i in at..at + g { rs += rk[i]; }
        at += g;
        h += rs * rs / g;
    }
    h = 12.0 / (nt * (nt + 1.0)) * h - 3.0 * (nt + 1.0);
    h = h / corr;
    if h < 0.0 { h = 0.0; }
    let df = (sizes.len - 1) as f64;
    let p, cp = chi2_sf(h, df);
    if cp != 0 { return fault_result(cp); }
    return good_result(h, p, df, nan());
}

// ---- the Kolmogorov-Smirnov tests --------------------------------------

// The upper tail of the one-sided statistic: P(D+ >= x) for n values.
// All terms of the sum are positive, so no digits cancel.
fn smirnov_sf(n: i64, x: f64) -> f64 {
    if x <= 0.0 { return 1.0; }
    if x >= 1.0 { return 0.0; }
    let nf = n as f64;
    let jmax = floor(nf * (1.0 - x)) as i64;
    var lc: f64 = 0.0;
    var sum = acc_new();
    for j in 0..jmax + 1 {
        if j > 0 { lc += log((nf - j + 1.0) / j); }
        let a = 1.0 - x - j / nf;
        if a <= 0.0 { continue; }
        let b = x + j / nf;
        acc_add(sum, exp(lc + (nf - j) * log(a) + (j - 1.0) * log(b)));
    }
    return min(1.0, x * acc_total(sum));
}

// The product of two square matrices of side m.
fn mat_mul(out: f64[>..]&, a: const f64[:], b: const f64[:], m: i64) {
    out.resize(0);
    for i in 0..m * m { out.push(0.0); }
    for i in 0..m {
        for k in 0..m {
            let v = a[i * m + k];
            if v == 0.0 { continue; }
            for j in 0..m { out[i * m + j] += v * b[k * m + j]; }
        }
    }
}

// P(D < d) of the two-sided statistic of n values, by the matrix power
// of Marsaglia, Tsang and Wang. The side of the matrix is 2 floor(n d) + 1.
fn mtw_cdf(n: i64, d: f64) -> f64 {
    let nf = n as f64;
    let k = floor(nf * d) as i64 + 1;
    let m = 2 * k - 1;
    let h = k - nf * d;
    var a: f64[>..] = [];
    for i in 0..m {
        for j in 0..m {
            if i - j + 1 < 0 { a.push(0.0); } else { a.push(1.0); }
        }
    }
    for i in 0..m {
        a[i * m] -= pow(h, i + 1.0);
        a[(m - 1) * m + i] -= pow(h, (m - i) * 1.0);
    }
    if 2.0 * h - 1.0 > 0.0 { a[(m - 1) * m] += pow(2.0 * h - 1.0, m * 1.0); }
    for i in 0..m {
        for j in 0..m {
            if i - j + 1 > 0 {
                for g in 1..i - j + 2 { a[i * m + j] /= g; }
            }
        }
    }
    var v: f64[>..] = [];
    for x in a { v.push(x); }
    var ev: i64 = 0;
    var top: i64 = 0;
    while (n >> (top + 1)) > 0 { top += 1; }
    var bit = top - 1;
    var tmp: f64[>..] = [];
    while bit >= 0 {
        mat_mul(tmp, v, v, m);
        ev = 2 * ev;
        v.resize(0);
        for x in tmp { v.push(x); }
        if ((n >> bit) & 1) == 1 {
            mat_mul(tmp, a, v, m);
            v.resize(0);
            for x in tmp { v.push(x); }
        }
        if v[(m / 2) * m + m / 2] > 1e140 {
            for i in 0..m * m { v[i] *= 1e-140; }
            ev += 140;
        }
        bit -= 1;
    }
    var s = v[(k - 1) * m + k - 1];
    for i in 1..n + 1 {
        s = s * i / nf;
        if s < 1e-140 {
            s *= 1e140;
            ev -= 140;
        }
    }
    return s * pow(10.0, ev as f64);
}

// The Pelz-Good series for the two-sided statistic: the upper tail.
fn pelz_good_sf(n: i64, x: f64) -> f64 {
    let nf = n as f64;
    let pi2 = 9.869604401089358;
    let z = sqrt(nf) * x;
    let z2 = z * z;
    let z3 = z2 * z;
    let z4 = z2 * z2;
    let z6 = z4 * z2;
    let z8 = z4 * z4;
    let qlog = -pi2 / 8.0 / z2;
    if qlog < -708.0 { return 1.0; }
    let q = exp(qlog);
    let pi4 = pi2 * pi2;
    let pi6 = pi4 * pi2;
    let k1a = -z2;
    let k1b = pi2 / 4.0;
    let k2a = 6.0 * z6 + 2.0 * z4;
    let k2b = (2.0 * z4 - 5.0 * z2) * pi2 / 4.0;
    let k2c = pi4 * (1.0 - 2.0 * z2) / 16.0;
    let k3d = pi6 * (5.0 - 30.0 * z2) / 64.0;
    let k3c = pi4 * (-60.0 * z2 + 212.0 * z4) / 16.0;
    let k3b = pi2 * (135.0 * z4 - 96.0 * z6) / 4.0;
    let k3a = -30.0 * z6 - 90.0 * z8;
    let maxk = ceil(16.0 * z / 3.141592653589793) as i64;
    var s0: f64 = 0.0;
    var s1: f64 = 0.0;
    var s2: f64 = 0.0;
    var s3: f64 = 0.0;
    var k = maxk;
    while k >= 1 {
        let m = (2 * k - 1) as f64;
        let ms = m * m;
        let m4 = ms * ms;
        let m6 = m4 * ms;
        let qp = pow(q, 8.0 * k);
        s0 = s0 * qp + 1.0;
        s1 = s1 * qp + k1a + k1b * ms;
        s2 = s2 * qp + k2a + k2b * ms + k2c * m4;
        s3 = s3 * qp + k3a + k3b * ms + k3c * m4 + k3d * m6;
        k -= 1;
    }
    let root = sqrt(2.0 * 3.141592653589793);
    let c0 = s0 * q * root / z;
    var c1 = s1 * q * root / (6.0 * z4);
    var c2 = s2 * q * root / (72.0 * z6 * z);
    var c3 = s3 * q * root / (6480.0 * z8 * z2);
    let q2 = exp(-pi2 / 2.0 / z2);
    let sq3 = sqrt(3.0) * z;
    var e2: f64 = 0.0;
    var e3: f64 = 0.0;
    for j in 1..maxk + 1 {
        let jf = j as f64;
        let w = jf * jf * pow(q2, jf * jf);
        e2 += w;
        e3 += (sq3 + 3.141592653589793 * jf) * (sq3 - 3.141592653589793 * jf) * w;
    }
    c2 += e2 * pi2 * root / (-36.0 * z3);
    c3 += e3 * pi2 * root / (216.0 * z6);
    c1 = c1 / sqrt(nf);
    c2 = c2 / nf;
    c3 = c3 / (nf * sqrt(nf));
    return (1.0 - c0) - c1 - c2 - c3;
}

fn clip01(p: f64) -> f64 {
    if p < 0.0 { return 0.0; }
    if p > 1.0 { return 1.0; }
    return p;
}

// The upper tail of the two-sided statistic of n values: P(D >= x).
fn kolmogn_sf(n: i64, x: f64) -> f64 {
    if x >= 1.0 { return 0.0; }
    if x <= 0.0 { return 1.0; }
    let nf = n as f64;
    let t = nf * x;
    if t <= 1.0 {
        if t <= 0.5 { return 1.0; }
        var prob: f64 = 1.0;
        if n <= 140 {
            for i in 1..n + 1 { prob = prob * i / nf * (2.0 * t - 1.0); }
        } else {
            let lf, sg, cl = lgamma(nf + 1.0);
            prob = exp(lf - nf * log(nf) + nf * log(2.0 * t - 1.0));
        }
        return clip01(1.0 - prob);
    }
    if t >= nf - 1.0 { return clip01(2.0 * pow(1.0 - x, nf)); }
    if x >= 0.5 { return clip01(2.0 * smirnov_sf(n, x)); }
    let nx2 = t * x;
    if n <= 140 {
        if nx2 <= 4.0 { return clip01(1.0 - mtw_cdf(n, x)); }
        return clip01(2.0 * smirnov_sf(n, x));
    }
    if nx2 >= 370.0 { return 0.0; }
    if nx2 >= 2.2 { return clip01(2.0 * smirnov_sf(n, x)); }
    if n <= 100000 && nf * pow(x, 1.5) <= 1.4 { return clip01(1.0 - mtw_cdf(n, x)); }
    return clip01(pelz_good_sf(n, x));
}

// The greatest common divisor of two positive values.
// The probability that a lattice path stays inside the band of width h
// around the diagonal from (0, 0) to (m, n). The count of the paths is
// scaled at each column, so the values stay in 0 to 1. This is the
// method of Hodges.
fn ks_inside(m: i64, n: i64, g: i64, h: i64, two_sided: bool) -> f64 {
    let mg = m / g;
    let ng = n / g;
    var minj: i64 = 0;
    var maxj = min(ceil(h as f64 / mg) as i64, n + 1);
    var curlen = maxj - minj;
    var len_a = min(2 * maxj + 2, n + 1);
    var a: f64[>..] = [];
    for i in 0..len_a { a.push(1.0); }
    for i in minj..maxj { a[i] = 0.0; }
    for i in 1..m + 1 {
        let last_minj = minj;
        let last_len = curlen;
        minj = max(floor((ng * i - h) as f64 / mg) as i64 + 1, 0);
        if minj > n { minj = n; }
        if two_sided {
            maxj = min(ceil((ng * i + h) as f64 / mg) as i64, n + 1);
        } else {
            maxj = n + 1;
        }
        if maxj <= minj { return 1.0; }
        var val: f64 = 0.0;
        if minj != 0 { val = 1.0; }
        for jj in minj..maxj {
            let j = jj;
            val = (a[j - last_minj] * i + val * j) / (i + j);
            a[jj - minj] = val;
        }
        curlen = maxj - minj;
        if last_len > curlen {
            for i in curlen..last_len { a[i] = 1.0; }
        }
    }
    return a[maxj - minj - 1];
}

// The two-sided exact p-value of the two-sample test. The two forms
// follow the source of SciPy: the square for equal sizes, the band for
// the others.
fn ks_two_sided_exact(d: f64, n1: i64, n2: i64, g: i64) -> f64 {
    let lcm = (n1 / g) * n2;
    let h = round(d * lcm) as i64;
    if h == 0 { return 1.0; }
    if n1 == n2 { return ks_outside_square(n1, h); }
    return ks_inside(n1, n2, g, h, true);
}

// The probability that a path leaves the band of width h around the
// diagonal of two equal samples.
fn ks_outside_square(n: i64, h: i64) -> f64 {
    var p: f64 = 0.0;
    var k = n / h;
    while k >= 0 {
        var p1: f64 = 1.0;
        for j in 0..h {
            p1 = (n - k * h - j) * p1 / (n + k * h + j + 1);
        }
        p = p1 * (1.0 - p);
        k -= 1;
    }
    return 2.0 * p;
}

// The binomial coefficient of a and b as a float. The loop runs b times.
fn binom_f(a: i64, b: i64) -> f64 {
    if b < 0 || b > a { return 0.0; }
    let k = if b > a - b { a - b } else { b };
    var r: f64 = 1.0;
    for i in 1..k + 1 { r = r * (a - k + i) / i; }
    return r;
}

// The number of lattice paths that touch the one-sided boundary at height
// h, for samples of m >= n values. This is the count of SciPy.
fn ks_paths_outside(m: i64, n: i64, g: i64, h: i64) -> f64 {
    let mg = m / g;
    let ng = n / g;
    // The floor of (mg - h) / mg, as the quotient of a negative number.
    var q = (mg - h) / mg;
    if (mg - h) % mg != 0 && mg - h < 0 { q -= 1; }
    let lxj = n + q;
    if lxj <= 0 { return binom_f(m + n, n); }
    var xj: i64[>..] = [];
    for j in 0..lxj { xj.push((h + mg * j + ng - 1) / ng); }
    var bs: f64[>..] = [];
    for j in 0..lxj { bs.push(0.0); }
    bs[0] = 1.0;
    for j in 1..lxj {
        var bj = binom_f(xj[j] + j, j);
        for i in 0..j { bj -= binom_f(xj[j] - xj[i] + j - i, j - i) * bs[i]; }
        bs[j] = bj;
    }
    var num: f64 = 0.0;
    for j in 0..lxj { num += bs[j] * binom_f(m - xj[j] + n - j, n - j); }
    return num;
}

// The one-sided exact p-value of the two-sample test.
fn ks_one_sided_exact(d: f64, n1: i64, n2: i64, g: i64) -> f64 {
    let lcm = (n1 / g) * n2;
    let h = round(d * lcm) as i64;
    if h == 0 { return 1.0; }
    if n1 == n2 {
        var p: f64 = 1.0;
        for j in 0..h { p = p * (n1 - j) / (n1 + j + 1); }
        return p;
    }
    return ks_paths_outside(n1, n2, g, h) / binom_f(n1 + n2, n2);
}

// The two-sample Kolmogorov-Smirnov test. The statistic is the largest
// gap between the two step functions.
fn ks_2samp(xs: const f64[:], ys: const f64[:], alt: i64 = ALT_TWO_SIDED) -> TestResult, i64 {
    var c = check_sample(xs, 1);
    if c != 0 { return fault_result(c); }
    c = check_sample(ys, 1);
    if c != 0 { return fault_result(c); }
    if check_alt(alt) != 0 { return fault_result(ERR_DOMAIN); }
    let n1 = xs.len;
    let n2 = ys.len;
    let sx = sort_copy(xs);
    let sy = sort_copy(ys);
    var dplus: f64 = 0.0;
    var dminus: f64 = 0.0;
    var i: i64 = 0;
    var j: i64 = 0;
    while i < n1 || j < n2 {
        var v: f64 = 0.0;
        if i < n1 && j < n2 {
            v = if sx[i] <= sy[j] { sx[i] } else { sy[j] };
        } else if i < n1 {
            v = sx[i];
        } else {
            v = sy[j];
        }
        while i < n1 && sx[i] == v { i += 1; }
        while j < n2 && sy[j] == v { j += 1; }
        let gap: f64 = i / (n1 as f64) - j / (n2 as f64);
        if gap > dplus { dplus = gap; }
        if -gap > dminus { dminus = -gap; }
    }
    var d = max(dplus, dminus);
    if alt == ALT_GREATER { d = dplus; }
    if alt == ALT_LESS { d = dminus; }
    if d == 0.0 { return good_result(0.0, 1.0, nan(), nan()); }
    let g = gcd_counts(n1, n2);
    let big = max(n1, n2);
    let small = min(n1, n2);
    let product = n1 * n2;
    let tx, _sx0, _sx1 = tie_sums(sx);
    let ty, _sy0, _sy1 = tie_sums(sy);
    let tied = tx > 0.0 || ty > 0.0;
    var p: f64 = 1.0;
    if product <= MAX_EXACT_KS && !tied {
        if alt == ALT_TWO_SIDED {
            p = ks_two_sided_exact(d, big, small, g);
        } else {
            p = ks_one_sided_exact(d, big, small, g);
        }
        if alt == ALT_TWO_SIDED { p = min(1.0, p); }
        return good_result(d, p, nan(), nan());
    }
    let m = max(n1, n2) as f64;
    let n = min(n1, n2) as f64;
    let en = m * n / (m + n);
    if alt == ALT_TWO_SIDED {
        let q = kolmogn_sf(round(en) as i64, d);
        return good_result(d, q, nan(), nan());
    }
    let z = sqrt(en) * d;
    let expt = -2.0 * z * z - 2.0 * z * (m + 2.0 * n) / (3.0 * sqrt(m * n * (m + n)));
    p = exp(expt);
    return good_result(d, p, nan(), nan());
}

// The Kolmogorov-Smirnov test of one sample against a family. The family
// is one of the FAMILY_ constants. The buffer holds the parameters of
// that family.
fn ks_1samp(xs: const f64[:], family: i64, params: const f64[:], alt: i64 = ALT_TWO_SIDED) -> TestResult, i64 {
    var c = check_sample(xs, 1);
    if c != 0 { return fault_result(c); }
    if check_alt(alt) != 0 { return fault_result(ERR_DOMAIN); }
    if family < 0 || family > FAMILY_F { return fault_result(ERR_DOMAIN); }
    let need = family_arity(family);
    if params.len != need { return fault_result(ERR_LENGTH); }
    for p in params {
        let cp = check_scalar(p);
        if cp != 0 { return fault_result(cp); }
    }
    let s = sort_copy(xs);
    let n = s.len;
    var dplus: f64 = 0.0;
    var dminus: f64 = 0.0;
    var k: i64 = 0;
    while k < n {
        let f, cf = family_cdf(family, s[k], params);
        if cf != 0 { return fault_result(cf); }
        let up = (k + 1) as f64 / n - f;
        let down = f - k as f64 / n;
        if up > dplus { dplus = up; }
        if down > dminus { dminus = down; }
        k += 1;
    }
    var d = max(dplus, dminus);
    if alt == ALT_GREATER { d = dplus; }
    if alt == ALT_LESS { d = dminus; }
    if d == 0.0 { return good_result(0.0, 1.0, nan(), nan()); }
    var p: f64 = 1.0;
    if alt == ALT_TWO_SIDED {
        p = kolmogn_sf(n, d);
    } else {
        p = smirnov_sf(n, d);
    }
    return good_result(d, p, nan(), nan());
}

// The count of the parameters of a family.
fn family_arity(family: i64) -> i64 {
    if family == FAMILY_EXPONENTIAL { return 1; }
    if family == FAMILY_CHI2 { return 1; }
    if family == FAMILY_T { return 1; }
    return 2;
}

// The cumulative function of a family at one value. The value is the
// result of the call and the code is the fault of the call.
fn family_cdf(family: i64, x: f64, params: const f64[:]) -> f64, i64 {
    if family == FAMILY_NORMAL { let v, c = norm_cdf(x, params[0], params[1]); return v, c; }
    if family == FAMILY_UNIFORM { let v, c = unif_cdf(x, params[0], params[1]); return v, c; }
    if family == FAMILY_EXPONENTIAL { let v, c = expon_cdf(x, params[0]); return v, c; }
    if family == FAMILY_LOGNORMAL { let v, c = lognorm_cdf(x, params[0], params[1]); return v, c; }
    if family == FAMILY_WEIBULL { let v, c = weibull_cdf(x, params[0], params[1]); return v, c; }
    if family == FAMILY_GAMMA { let v, c = gammadist_cdf(x, params[0], params[1]); return v, c; }
    if family == FAMILY_BETA { let v, c = betadist_cdf(x, params[0], params[1]); return v, c; }
    if family == FAMILY_CHI2 { let v, c = chi2_cdf(x, params[0]); return v, c; }
    if family == FAMILY_T { let v, c = tdist_cdf(x, params[0]); return v, c; }
    let v, c = fdist_cdf(x, params[0], params[1]);
    return v, c;
}

// ---- the normality test and the correlation tests ---------------------

// The value of a polynomial with the coefficients in rising order.
fn poly_at(c: const f64[:], x: f64) -> f64 {
    var r: f64 = 0.0;
    var k = c.len;
    while k > 0 {
        k -= 1;
        r = r * x + c[k];
    }
    return r;
}

// The Shapiro-Wilk test of Royston (1992). It needs 3 values. A sample
// of equal values gives the undefined fault.
fn shapiro(xs: const f64[:]) -> TestResult, i64 {
    let c = check_sample(xs, 3);
    if c != 0 { return fault_result(c); }
    let s = sort_copy(xs);
    let n = s.len;
    let range = s[n - 1] - s[0];
    if is_inf(range) { return fault_result(ERR_OVERFLOW); }
    if range == 0.0 { return fault_result(ERR_UNDEFINED); }
    var y: f64[>..] = [];
    for v in s { y.push((v - s[0]) / range); }
    let half = n / 2;
    var a: f64[>..] = [];
    if n == 3 {
        a.push(sqrt(0.5));
    } else {
        let an25: f64 = n + 0.25;
        let c1: f64[6] = [0.0, 0.221157, -0.147981, -2.07119, 4.434685, -2.706056];
        let c2: f64[6] = [0.0, 0.042981, -0.293762, -1.752461, 5.682633, -3.582633];
        var m: f64[>..] = [];
        var sum2: f64 = 0.0;
        for i in 1..(half + 1) {
            let q = norm_std_ppf((i - 0.375) / an25);
            m.push(q);
            sum2 += q * q;
        }
        sum2 = 2.0 * sum2;
        let root = sqrt(sum2);
        let nf: f64 = n;
        let rsn = 1.0 / sqrt(nf);
        let a1 = poly_at(c1, rsn) - m[0] / root;
        var first = 1;
        var fac: f64 = 0.0;
        a.push(a1);
        if n > 5 {
            let a2 = poly_at(c2, rsn) - m[1] / root;
            first = 2;
            fac = sqrt((sum2 - 2.0 * m[0] * m[0] - 2.0 * m[1] * m[1]) / (1.0 - 2.0 * a1 * a1 - 2.0 * a2 * a2));
            a.push(a2);
        } else {
            fac = sqrt((sum2 - 2.0 * m[0] * m[0]) / (1.0 - 2.0 * a1 * a1));
        }
        for i in first..half { a.push(-m[i] / fac); }
    }
    var mean_y: f64 = 0.0;
    for v in y { mean_y += v; }
    let nn: f64 = n;
    mean_y = mean_y / nn;
    var ssq: f64 = 0.0;
    var b: f64 = 0.0;
    for v in y { ssq += (v - mean_y) * (v - mean_y); }
    for i in 0..half { b += a[i] * (y[n - 1 - i] - y[i]); }
    var w = b * b / ssq;
    if w > 1.0 { w = 1.0; }
    var p: f64 = 0.0;
    if n == 3 {
        p = 1.90985931710274 * (asin(sqrt(w)) - 1.04719755119660);
        if p < 0.0 { p = 0.0; }
    } else {
        let w1 = 1.0 - w;
        var ly = log(w1);
        let nf: f64 = n;
        let lx = log(nf);
        var mu: f64 = 0.0;
        var sigma: f64 = 0.0;
        if n <= 11 {
            let g: f64[2] = [-2.273, 0.459];
            let gamma = poly_at(g, nf);
            // W stays far below this limit for 4 to 11 values, so the
            // guard only keeps a NaN out of the result.
            if ly >= gamma { return fault_result(ERR_UNDEFINED); }
            ly = -log(gamma - ly);
            let c3: f64[4] = [0.544, -0.39978, 0.025054, -6.714e-4];
            let c4: f64[4] = [1.3822, -0.77857, 0.062767, -0.0020322];
            mu = poly_at(c3, nf);
            sigma = exp(poly_at(c4, nf));
        } else {
            let c5: f64[4] = [-1.5861, -0.31082, -0.083751, 0.0038915];
            let c6: f64[3] = [-0.4803, -0.082676, 0.0030302];
            mu = poly_at(c5, lx);
            sigma = exp(poly_at(c6, lx));
        }
        let q, cq = norm_sf((ly - mu) / sigma);
        if cq != 0 { return fault_result(cq); }
        p = q;
    }
    return good_result(w, p, nan(), nan());
}

// The Pearson test. The p-value uses the beta form of the distribution
// of r, so a p-value near 0 keeps its digits. Two points give p = 1.
fn pearsonr(xs: const f64[:], ys: const f64[:], alt: i64 = ALT_TWO_SIDED) -> TestResult, i64 {
    let ca = check_alt(alt);
    let c = check_pair(xs, ys, 2);
    if c != 0 { return fault_result(c); }
    if ca != 0 { return fault_result(ca); }
    let r, c2 = correlation(xs, ys);
    if c2 != 0 { return fault_result(c2); }
    let n = xs.len;
    let df: f64 = n - 2;
    if n == 2 { return good_result(r, 1.0, df, nan()); }
    let ab = n / 2.0 - 1.0;
    var x: f64 = 0.0;
    if alt == ALT_TWO_SIDED { x = 0.5 * (1.0 - abs(r)); }
    if alt == ALT_LESS { x = 0.5 * (1.0 + r); }
    if alt == ALT_GREATER { x = 0.5 * (1.0 - r); }
    if x < 0.0 { x = 0.0; }
    if x > 1.0 { x = 1.0; }
    let q, cb = betainc(ab, ab, x);
    if cb != 0 { return fault_result(cb); }
    var p = q;
    if alt == ALT_TWO_SIDED { p = 2.0 * q; }
    return good_result(r, p, df, nan());
}

// The Spearman test. The p-value uses the t form with n - 2 degrees of
// freedom. Two points give p = 1.
fn spearmanr(xs: const f64[:], ys: const f64[:], alt: i64 = ALT_TWO_SIDED) -> TestResult, i64 {
    let ca = check_alt(alt);
    let c = check_pair(xs, ys, 2);
    if c != 0 { return fault_result(c); }
    if ca != 0 { return fault_result(ca); }
    let r, c2 = spearman(xs, ys);
    if c2 != 0 { return fault_result(c2); }
    let n = xs.len;
    let df: f64 = n - 2;
    if n == 2 { return good_result(r, 1.0, df, nan()); }
    var t = (1.0 / 0.0);
    let den = (1.0 + r) * (1.0 - r);
    if den > 0.0 { t = r * sqrt(df / den); }
    if den <= 0.0 && r < 0.0 { t = -(1.0 / 0.0); }
    let p, cp = t_pvalue(t, df, alt);
    if cp != 0 { return fault_result(cp); }
    return good_result(r, p, df, nan());
}

// The sums over the runs of equal values of a sorted sample: the pairs
// t(t-1)/2, then t(t-1)(t-2) and t(t-1)(2t+5).
fn tie_sums(s: const f64[:]) -> f64, f64, f64 {
    var pairs: f64 = 0.0;
    var s0: f64 = 0.0;
    var s1: f64 = 0.0;
    var k: i64 = 0;
    while k < s.len {
        var j = k + 1;
        while j < s.len && s[j] == s[k] { j += 1; }
        let t: f64 = j - k;
        pairs += t * (t - 1.0) / 2.0;
        s0 += t * (t - 1.0) * (t - 2.0);
        s1 += t * (t - 1.0) * (2.0 * t + 5.0);
        k = j;
    }
    return pairs, s0, s1;
}

// The exact tails of the number of discordant pairs of n values with no
// ties. The first result is P(D <= c) and the second is P(D = c).
fn kendall_exact_tail(n: i64, c: i64) -> f64, f64 {
    var cur: f64[>..] = [];
    var nxt: f64[>..] = [];
    for k in 0..(c + 1) {
        cur.push(0.0);
        nxt.push(0.0);
    }
    cur[0] = 1.0;
    for j in 2..(n + 1) {
        var window: f64 = 0.0;
        for k in 0..(c + 1) {
            window += cur[k];
            if k >= j { window -= cur[k - j]; }
            nxt[k] = window / j;
        }
        for k in 0..(c + 1) { cur[k] = nxt[k]; }
    }
    var total: f64 = 0.0;
    for k in 0..(c + 1) { total += cur[k]; }
    return total, cur[c];
}

// The Kendall tau-b test. The p-value is exact for a sample with no ties
// of at most MAX_EXACT_KENDALL values, or one pair from a perfect order.
// Any other sample uses the normal form with the variance of the ties.
fn kendalltau(xs: const f64[:], ys: const f64[:], alt: i64 = ALT_TWO_SIDED) -> TestResult, i64 {
    let ca = check_alt(alt);
    let c = check_pair(xs, ys, 2);
    if c != 0 { return fault_result(c); }
    if ca != 0 { return fault_result(ca); }
    let n = xs.len;
    let cd, t = kendall_counts(xs, ys);
    // The tau-b of the pairs, from the same counts.
    let tot = n * (n - 1) / 2;
    let a = tot - t.n1;
    let b = tot - t.n2;
    if a <= 0 || b <= 0 { return fault_result(ERR_UNDEFINED); }
    // Widen before the product. Two i64 values multiply as i64 and wrap
    // above 2 to the power 63, which happens at about 78,000 values.
    let fa: f64 = a;
    let fb: f64 = b;
    let ab = fa * fb;
    let den = sqrt(ab);
    if den == 0.0 { return fault_result(ERR_UNDEFINED); }
    let tau = cd / den;
    if t.n1 == 0 && t.n2 == 0 {
        let con = (tot + cd) / 2;
        let dis = tot - con;
        let low = if dis < con { dis } else { con };
        if n <= MAX_EXACT_KENDALL || low <= 1 {
            if n == 2 { return good_result(tau, 1.0, nan(), nan()); }
            let right = con >= tot - con;
            let cc = if right { tot - con } else { con };
            let cdf, mass = kendall_exact_tail(n, cc);
            var p: f64 = 1.0;
            if alt == ALT_TWO_SIDED {
                if 4 * cc != n * (n - 1) { p = 2.0 * cdf; }
            } else {
                let agree = right == (alt == ALT_GREATER);
                if agree { p = cdf; } else { p = 1.0 - cdf + mass; }
            }
            return good_result(tau, p, nan(), nan());
        }
    }
    let px, x0, x1 = tie_sums(sort_copy(xs));
    let py, y0, y1 = tie_sums(sort_copy(ys));
    let nf: f64 = n;
    let m = nf * (nf - 1.0);
    let v = (m * (2.0 * nf + 5.0) - x1 - y1) / 18.0 + 2.0 * px * py / m + x0 * y0 / (9.0 * m * (nf - 2.0));
    if !(v > 0.0) { return fault_result(ERR_UNDEFINED); }
    let z = cd / sqrt(v);
    var p: f64 = 0.0;
    var cz: i64 = 0;
    if alt == ALT_TWO_SIDED {
        let q, cq = norm_sf(abs(z));
        p = 2.0 * q;
        cz = cq;
    } else if alt == ALT_LESS {
        let q, cq = norm_cdf(z);
        p = q;
        cz = cq;
    } else {
        let q, cq = norm_sf(z);
        p = q;
        cz = cq;
    }
    if cz != 0 { return fault_result(cz); }
    return good_result(tau, p, nan(), nan());
}

// ---- the multiple testing corrections ---------------------------------

// Writes the adjusted p-values into out, in the order of ps. The method
// is one of the ADJUST_ constants. A p-value outside 0 to 1 gives the
// domain fault. Returns the count written and the code.
fn multipletests(out: f64[>..]&, ps: const f64[:], method: i64 = ADJUST_HOLM) -> i64, i64 {
    out.resize(0);
    let c = check_sample(ps, 1);
    if c != 0 { return 0, c; }
    if method < 0 || method > 2 { return 0, 5; }
    for p in ps {
        if p < 0.0 || p > 1.0 { return 0, 5; }
    }
    let n = ps.len;
    var adj: f64[>..] = [];
    for k in 0..n { adj.push(0.0); }
    let nf: f64 = n;
    if method == ADJUST_BONFERRONI {
        // The Bonferroni form keeps the order of the list, so it needs no
        // sort of the indices.
        for k in 0..n { adj[k] = ps[k] * nf; }
    } else {
        var order: i64[>..] = [];
        for k in 0..n { order.push(k); }
        sort(order) { a, b => ps[a] < ps[b] || (ps[a] == ps[b] && a < b) };
        if method == ADJUST_HOLM {
            var run: f64 = 0.0;
            for k in 0..n {
                let v = ps[order[k]] * (nf - k);
                if v > run { run = v; }
                adj[order[k]] = run;
            }
        } else {
            var run: f64 = 1.0;
            var k = n;
            while k > 0 {
                k -= 1;
                let v = ps[order[k]] * nf / (k + 1);
                if v < run { run = v; }
                adj[order[k]] = run;
            }
        }
    }
    for k in 0..n {
        var v = adj[k];
        if v > 1.0 { v = 1.0; }
        out.push(v);
    }
    return n, 0;
}
