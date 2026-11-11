// Timing tests of the hypothesis tests. Each call runs on a crafted input at
// the edge of an exact bound or at a large size. The answer must be sound,
// and the call must finish within ten seconds.
namespace suite;

import std;
import math;
import os;
import statistics from "../..";

const TIME_LIMIT: f64 = 10.0;

fn time_check(bad: i64&, r: statistics::TestResult, c: i64, name: u8[:], t0: f64) {
    let spent = clock() - t0;
    var ok = spent < TIME_LIMIT;
    if c == 0 { ok = ok && !is_nan(r.pvalue) && r.pvalue >= 0.0 && r.pvalue <= 1.0; }
    if !ok {
        bad += 1;
        print("FAIL ", name, " took ", spent, " seconds, code ", c, "\n");
    }
    print("time ", name, " ", spent, "\n");
}

// A seeded sample of n values. With ties the values come from a small set.
fn time_sample(n: i64, ties: bool, st: u64&) -> f64[>..] {
    var out: f64[>..] = [];
    for i in 0..n {
        if ties { out.push(rnd(st) % 7); } else { out.push((rnd(st) % 1000003) / 1000.0 + i * 1e-9); }
    }
    return out;
}

// The count of the calls that take samples. The dispatcher below runs one.
const TIME_OPS: i64 = 14;

fn time_op(op: i64, xs: const f64[:], ys: const f64[:]) -> statistics::TestResult, i64 {
    let half = xs.len / 2;
    let sizes: i64[2] = [half, xs.len - half];
    if op == 0 { return statistics::ttest_ind(xs, ys, statistics::ALT_TWO_SIDED, false); }
    if op == 1 { return statistics::mannwhitneyu(xs, ys); }
    if op == 2 { return statistics::wilcoxon_paired(xs, ys); }
    if op == 3 { return statistics::shapiro(xs); }
    if op == 4 { return statistics::pearsonr(xs, ys); }
    if op == 5 { return statistics::spearmanr(xs, ys); }
    if op == 6 { return statistics::kendalltau(xs, ys); }
    if op == 7 { return statistics::ks_2samp(xs, ys); }
    if op == 8 { return statistics::ks_1samp(xs, statistics::FAMILY_NORMAL, [500.0, 300.0][..]); }
    if op == 9 { return statistics::kruskal(xs, sizes[..]); }
    if op == 10 { return statistics::f_oneway(xs, sizes[..]); }
    if op == 11 { return statistics::ttest_rel(xs, ys); }
    if op == 12 {
        // Counts above 0 whose sums are equal: the second list reverses the first.
        var obs: f64[>..] = [];
        var exp: f64[>..] = [];
        for i in 0..xs.len { obs.push(abs(xs[i]) + 1.0); }
        for i in 0..xs.len { exp.push(obs[xs.len - 1 - i]); }
        return statistics::chisquare(obs[..], exp[..]);
    }
    // Probabilities from 0 to 1.
    var ps: f64[>..] = [];
    for y in ys { ps.push(1.0 / (1.0 + abs(y))); }
    var out: f64[>..] = [];
    let n, c = statistics::multipletests(out, ps[..], statistics::ADJUST_BH);
    return statistics::TestResult { statistic: 0.0, pvalue: 0.5, df: 0.0, df2: 0.0 }, c;
}

fn time_name(op: i64) -> u8[:] {
    if op == 0 { return "ttest_ind"; }
    if op == 1 { return "mannwhitneyu"; }
    if op == 2 { return "wilcoxon_paired"; }
    if op == 3 { return "shapiro"; }
    if op == 4 { return "pearsonr"; }
    if op == 5 { return "spearmanr"; }
    if op == 6 { return "kendalltau"; }
    if op == 7 { return "ks_2samp"; }
    if op == 8 { return "ks_1samp"; }
    if op == 9 { return "kruskal"; }
    if op == 10 { return "f_oneway"; }
    if op == 11 { return "ttest_rel"; }
    if op == 12 { return "chisquare"; }
    return "multipletests";
}

// A crafted sample: 0 sorted, 1 reversed, 2 few distinct values, 3 one
// repeated block after a random start.
fn time_crafted(kind: i64, n: i64, st: u64&) -> f64[>..] {
    var out: f64[>..] = [];
    for i in 0..n {
        if kind == 0 { out.push(i); }
        if kind == 1 { out.push(n - i); }
        if kind == 2 { out.push(1 + rnd(st) % 3); }
        if kind == 3 { out.push(if i < n / 8 { 1 + rnd(st) % 1000 } else { 7.0 }); }
    }
    return out;
}

fn time_median3(a: f64, b: f64, c: f64) -> f64 {
    return max(min(a, b), min(max(a, b), c));
}

// Each call runs on random input three times, then on four crafted inputs of
// the same size. A crafted input takes at most 50 times the median of the
// random runs, with a floor of a tenth of a second for the clock.
fn run_time_ratio(n: i64) -> i64 {
    var bad: i64 = 0;
    var st: u64 = 777;
    for op in 0..TIME_OPS {
        // Chi-square needs positive counts with equal sums, so its random
        // input is a count list.
        var spent: f64[3] = [0.0, 0.0, 0.0];
        for rep in 0..3 {
            let xs = time_sample(n, false, st);
            let ys = time_sample(n, false, st);
            let t0 = clock();
            let r, c = time_op(op, xs[..], ys[..]);
            spent[rep] = clock() - t0;
        }
        let med = time_median3(spent[0], spent[1], spent[2]);
        for kind in 0..4 {
            let xs = time_crafted(kind, n, st);
            let ys = time_crafted((kind + 1) % 4, n, st);
            let t0 = clock();
            let r, c = time_op(op, xs[..], ys[..]);
            let took = clock() - t0;
            if took > 50.0 * med + 0.1 {
                bad += 1;
                print("FAIL ", time_name(op), " crafted ", kind, " took ", took, " against ", med, "\n");
            }
            if c == 0 && (is_nan(r.pvalue) || r.pvalue < 0.0 || r.pvalue > 1.0) {
                bad += 1;
                print("FAIL ", time_name(op), " crafted ", kind, " p-value ", r.pvalue, "\n");
            }
        }
        print("time ratio ", time_name(op), " median ", med, "\n");
    }
    return bad;
}

// Each call runs on the largest input, random and with few distinct values.
// Each run must finish in under 60 seconds.
fn run_time_max() -> i64 {
    var bad: i64 = 0;
    var st: u64 = 31337;
    for op in 0..TIME_OPS {
        // The Mann-Whitney test pools its two samples, so the pool is at
        // the bound when each sample holds half of it.
        let n = if op == 1 { statistics::MAX_INPUT / 2 } else { statistics::MAX_INPUT };
        for kind in 0..2 {
            let xs = if kind == 0 { time_sample(n, false, st) } else { time_crafted(2, n, st) };
            let ys = if kind == 0 { time_sample(n, false, st) } else { time_crafted(2, n, st) };
            let t0 = clock();
            let r, c = time_op(op, xs[..], ys[..]);
            let took = clock() - t0;
            if took >= 60.0 || c != 0 || (c == 0 && (is_nan(r.pvalue) || r.pvalue < 0.0 || r.pvalue > 1.0)) {
                bad += 1;
                print("FAIL ", time_name(op), " at the input bound took ", took, " code ", c, "\n");
            }
            print("time max ", time_name(op), " ", kind, " ", took, "\n");
        }
    }
    return bad;
}

fn run_time_tests() -> i64 {
    var bad: i64 = 0;
    var st: u64 = 4242;
    let big = 200000;
    let xs = time_sample(big, false, st);
    let ys = time_sample(big, false, st);
    let tx = time_sample(big, true, st);
    let ty = time_sample(big, true, st);

    var t0 = clock();
    let r1, c1 = statistics::ttest_ind(xs[..], ys[..], statistics::ALT_TWO_SIDED, false);
    time_check(bad, r1, c1, "ttest_ind at 200000", t0);
    t0 = clock();
    let r2, c2 = statistics::mannwhitneyu(xs[..], ys[..]);
    time_check(bad, r2, c2, "mannwhitneyu at 200000", t0);
    t0 = clock();
    let r3, c3 = statistics::mannwhitneyu(tx[..], ty[..]);
    time_check(bad, r3, c3, "mannwhitneyu ties at 200000", t0);
    t0 = clock();
    let r4, c4 = statistics::wilcoxon_paired(xs[..], ys[..]);
    time_check(bad, r4, c4, "wilcoxon at 200000", t0);
    t0 = clock();
    let r5, c5 = statistics::wilcoxon_paired(tx[..], ty[..]);
    time_check(bad, r5, c5, "wilcoxon ties at 200000", t0);
    t0 = clock();
    let r6, c6 = statistics::shapiro(xs[..]);
    time_check(bad, r6, c6, "shapiro at 200000", t0);
    t0 = clock();
    let r7, c7 = statistics::pearsonr(xs[..], ys[..]);
    time_check(bad, r7, c7, "pearsonr at 200000", t0);
    t0 = clock();
    let r8, c8 = statistics::spearmanr(tx[..], ty[..]);
    time_check(bad, r8, c8, "spearmanr ties at 200000", t0);
    t0 = clock();
    let r9, c9 = statistics::kendalltau(xs[..], ys[..]);
    time_check(bad, r9, c9, "kendalltau at 200000", t0);
    t0 = clock();
    let r10, c10 = statistics::kendalltau(tx[..], ty[..]);
    time_check(bad, r10, c10, "kendalltau ties at 200000", t0);
    t0 = clock();
    let r11, c11 = statistics::ks_2samp(xs[..], ys[..]);
    time_check(bad, r11, c11, "ks_2samp at 200000", t0);
    t0 = clock();
    let r12, c12 = statistics::ks_2samp(tx[..], ty[..]);
    time_check(bad, r12, c12, "ks_2samp ties at 200000", t0);
    t0 = clock();
    let r13, c13 = statistics::ks_1samp(xs[..], statistics::FAMILY_NORMAL, [500.0, 300.0][..]);
    time_check(bad, r13, c13, "ks_1samp at 200000", t0);
    var sizes: i64[>..] = [big / 2, big / 2];
    t0 = clock();
    let r14, c14 = statistics::kruskal(xs[..], sizes[..]);
    time_check(bad, r14, c14, "kruskal at 200000", t0);
    t0 = clock();
    let r15, c15 = statistics::f_oneway(xs[..], sizes[..]);
    time_check(bad, r15, c15, "f_oneway at 200000", t0);
    var out: f64[>..] = [];
    t0 = clock();
    let n16, c16 = statistics::multipletests(out, ys[..], statistics::ADJUST_BH);
    time_check(bad, statistics::TestResult { statistic: 0.0, pvalue: 0.5, df: 0.0, df2: 0.0 }, c16, "multipletests at 200000", t0);

    // The edges of the exact bounds.
    let a8 = time_sample(8, false, st);
    let b12 = time_sample(12500, false, st);
    t0 = clock();
    let r17, c17 = statistics::mannwhitneyu(a8[..], b12[..]);
    time_check(bad, r17, c17, "mannwhitneyu 8 against 12500", t0);
    let w50 = time_sample(50, false, st);
    t0 = clock();
    let r18, c18 = statistics::wilcoxon(w50[..]);
    time_check(bad, r18, c18, "wilcoxon at 50 pairs", t0);
    var w13: f64[>..] = [];
    for i in 0..13 { w13.push((rnd(st) % 5) - 2); }
    t0 = clock();
    let r19, c19 = statistics::wilcoxon(w13[..]);
    time_check(bad, r19, c19, "wilcoxon at 13 pairs with ties", t0);
    let k100 = time_sample(100, false, st);
    let k100b = time_sample(100, false, st);
    t0 = clock();
    let r20, c20 = statistics::ks_2samp(k100[..], k100b[..]);
    time_check(bad, r20, c20, "ks_2samp 100 against 100", t0);
    let k1 = time_sample(1, false, st);
    let k10000 = time_sample(10000, false, st);
    t0 = clock();
    let r21, c21 = statistics::ks_2samp(k1[..], k10000[..], statistics::ALT_GREATER);
    time_check(bad, r21, c21, "ks_2samp 1 against 10000 one-sided", t0);
    let kd = time_sample(33, false, st);
    let ke = time_sample(33, false, st);
    t0 = clock();
    let r22, c22 = statistics::kendalltau(kd[..], ke[..]);
    time_check(bad, r22, c22, "kendalltau at 33", t0);

    // Huge counts.
    t0 = clock();
    let r23, c23 = statistics::binomtest(500000000.0, 1000000000.0, 0.5);
    time_check(bad, r23, c23, "binomtest at 1e9 trials", t0);
    t0 = clock();
    let r24, c24 = statistics::fisher_exact([1e8, 2e8, 3e8, 4e8][..]);
    time_check(bad, r24, c24, "fisher_exact at 1e8 counts", t0);
    t0 = clock();
    let r25, c25 = statistics::chi2_contingency([1e15, 2e15, 3e15, 4e15][..], 2, 2);
    time_check(bad, r25, c25, "chi2_contingency at 1e15 counts", t0);
    t0 = clock();
    let r26, c26 = statistics::normaltest(xs[..]);
    time_check(bad, r26, c26, "normaltest at 200000", t0);
    t0 = clock();
    let r27, c27 = statistics::jarque_bera(tx[..]);
    time_check(bad, r27, c27, "jarque-bera ties at 200000", t0);
    t0 = clock();
    let r28, c28 = statistics::anderson_darling(xs[..]);
    time_check(bad, r28, c28, "anderson-darling at 200000", t0);
    t0 = clock();
    let r29, c29 = statistics::skewtest(xs[..], statistics::ALT_GREATER);
    time_check(bad, r29, c29, "skewtest at 200000", t0);
    t0 = clock();
    let r30, c30 = statistics::kurtosistest(tx[..]);
    time_check(bad, r30, c30, "kurtosistest ties at 200000", t0);
    bad += run_time_ratio(100000);
    // A run at the input bound takes minutes, so the mutation probes skip it.
    var flag: u8[>..] = [];
    if !env("STATISTICS_HOSTILE_N", flag) { bad += run_time_max(); }
    return bad;
}
