// Hostile input tests of the hypothesis tests. Every call gets seeded random
// inputs, often hostile: NaN, infinity, ties, zeros, huge values, wrong
// lengths and wrong options. Each answer must be a value with code 0 and a
// p-value in 0 to 1, or NaN in every field with a known fault code.
namespace suite;

import std;
import math;
import os;
import statistics from "../..";

// The count of seeded inputs for each call: 100000, or the number in the
// variable STATISTICS_HOSTILE_N. The mutation probes use a smaller count.
fn tests_count() -> i64 {
    var buf: u8[>..] = [];
    if !env("STATISTICS_HOSTILE_N", buf) { return 100000; }
    var n: i64 = 0;
    for ch in buf {
        if ch < 48 || ch > 57 { return 100000; }
        n = n * 10 + (ch - 48);
    }
    if n < 1 { return 100000; }
    return n;
}

// One seeded value, often hostile.
fn hostile_value(st: u64&) -> f64 {
    let r = rnd(st) % 16;
    let u = (rnd(st) % 2000001 - 1000000) / 1000000.0;
    if r == 0 { return statistics::nan(); }
    if r == 1 { return (1.0 / 0.0); }
    if r == 2 { return -(1.0 / 0.0); }
    if r == 3 { return 0.0; }
    if r == 4 { return u * 1e300; }
    if r == 5 { return u * 1e-300; }
    if r < 9 { return floor(u * 4.0); }
    return u * 10.0;
}

// A sample of a seeded length from 0 to 14.
fn hostile_sample(st: u64&, hostile: bool) -> f64[>..] {
    var out: f64[>..] = [];
    let n = rnd(st) % 15;
    for k in 0..n {
        if hostile {
            out.push(hostile_value(st));
        } else {
            out.push((rnd(st) % 9 - 4) + (rnd(st) % 5) / 4.0);
        }
    }
    return out;
}

// True when the answer is a value or a clean fault.
fn sound(r: statistics::TestResult, c: i64) -> bool {
    if c == 0 {
        if is_nan(r.statistic) && is_nan(r.pvalue) { return false; }
        if is_nan(r.pvalue) { return false; }
        if r.pvalue < 0.0 || r.pvalue > 1.0 { return false; }
        return true;
    }
    if c < 1 || c > 9 { return false; }
    return is_nan(r.statistic) && is_nan(r.pvalue) && is_nan(r.df) && is_nan(r.df2);
}

fn hostile_report(bad: i64&, ok: bool, name: u8[:], k: i64, c: i64) {
    if ok { return; }
    bad += 1;
    if bad <= 12 { print("FAIL ", name, " input ", k, " code ", c, "\n"); }
}

fn run_hostile_tests() -> i64 {
    var bad: i64 = 0;
    var st: u64 = 987654321;
    let n = tests_count();
    for k in 0..n {
        let hostile = rnd(st) % 3 != 0;
        let xs = hostile_sample(st, hostile);
        var ys = hostile_sample(st, hostile);
        if rnd(st) % 2 == 0 {
            ys.resize(0);
            for v in xs { ys.push(if hostile { hostile_value(st) } else { rnd(st) % 7 - 3 }); }
        }
        var tab: f64[>..] = [];
        let rows = 2 + rnd(st) % 2;
        let cols = 2 + rnd(st) % 2;
        for t in 0..rows * cols {
            tab.push(if rnd(st) % 20 == 0 { hostile_value(st) } else { rnd(st) % 12 });
        }
        var cnt: f64[>..] = [];
        for t in 0..1 + rnd(st) % 5 { cnt.push(rnd(st) % 9); }
        let zs = hostile_sample(st, hostile);
        let alt = rnd(st) % 4;
        let opt = rnd(st) % 4;
        var sizes: i64[>..] = [];
        let groups = rnd(st) % 4;
        for g in 0..groups { sizes.push(rnd(st) % 5 - 1); }
        if rnd(st) % 2 == 0 && xs.len >= 2 {
            sizes.resize(0);
            let cut = 1 + rnd(st) % (xs.len - 1);
            sizes.push(cut);
            sizes.push(xs.len - cut);
        }
        let family = rnd(st) % 11;
        var params: f64[>..] = [];
        if rnd(st) % 3 != 0 {
            params.push(0.5 + (rnd(st) % 5) / 2.0);
            if family != 2 && family != 7 && family != 8 { params.push(1.5 + (rnd(st) % 5)); }
        } else {
            for v in zs { params.push(v); }
        }
        let flag = rnd(st) % 2 == 0;

        let r1, c1 = statistics::ttest_1samp(xs[..], 0.0, alt);
        hostile_report(bad, sound(r1, c1), "ttest_1samp", k, c1);
        let r2, c2 = statistics::ttest_ind(xs[..], ys[..], alt, flag);
        hostile_report(bad, sound(r2, c2), "ttest_ind", k, c2);
        let r3, c3 = statistics::ttest_rel(xs[..], ys[..], alt);
        hostile_report(bad, sound(r3, c3), "ttest_rel", k, c3);
        let r4, c4 = statistics::f_oneway(xs[..], sizes[..]);
        hostile_report(bad, sound(r4, c4), "f_oneway", k, c4);
        let r5, c5 = statistics::kruskal(xs[..], sizes[..]);
        hostile_report(bad, sound(r5, c5), "kruskal", k, c5);
        let r6, c6 = statistics::chisquare(cnt[..], if flag { cnt[..] } else { ys[..] }, opt % 2);
        hostile_report(bad, sound(r6, c6), "chisquare", k, c6);
        let r7, c7 = statistics::chi2_contingency(tab[..], if flag { rows } else { rnd(st) % 4 }, cols, flag);
        hostile_report(bad, sound(r7, c7), "chi2_contingency", k, c7);
        let r8, c8 = statistics::fisher_exact(tab[0..if tab.len >= 4 { 4 } else { tab.len }], alt);
        hostile_report(bad, sound(r8, c8), "fisher_exact", k, c8);
        let r9, c9 = statistics::mannwhitneyu(xs[..], ys[..], flag, alt, opt);
        hostile_report(bad, sound(r9, c9), "mannwhitneyu", k, c9);
        let r10, c10 = statistics::wilcoxon(xs[..], alt, opt, flag);
        hostile_report(bad, sound(r10, c10), "wilcoxon", k, c10);
        let r11, c11 = statistics::wilcoxon_paired(xs[..], ys[..], alt, opt, flag);
        hostile_report(bad, sound(r11, c11), "wilcoxon_paired", k, c11);
        let r12, c12 = statistics::ks_2samp(xs[..], ys[..], alt);
        hostile_report(bad, sound(r12, c12), "ks_2samp", k, c12);
        let r13, c13 = statistics::ks_1samp(xs[..], family, params[..], alt);
        hostile_report(bad, sound(r13, c13), "ks_1samp", k, c13);
        let r14, c14 = statistics::shapiro(xs[..]);
        hostile_report(bad, sound(r14, c14), "shapiro", k, c14);
        let r15, c15 = statistics::pearsonr(xs[..], ys[..], alt);
        hostile_report(bad, sound(r15, c15), "pearsonr", k, c15);
        let r16, c16 = statistics::spearmanr(xs[..], ys[..], alt);
        hostile_report(bad, sound(r16, c16), "spearmanr", k, c16);
        let r17, c17 = statistics::kendalltau(xs[..], ys[..], alt);
        hostile_report(bad, sound(r17, c17), "kendalltau", k, c17);
        let r18, c18 = statistics::binomtest(xs.len, ys.len, hostile_value(st), alt);
        hostile_report(bad, sound(r18, c18), "binomtest", k, c18);
        var out: f64[>..] = [];
        let m, c19 = statistics::multipletests(out, xs[..], opt);
        let ok19 = (c19 == 0 && m == xs.len && out.len == m) || (c19 != 0 && m == 0 && out.len == 0);
        hostile_report(bad, ok19, "multipletests", k, c19);
        let h1, hc1 = statistics::skewtest(xs[..], alt % 3);
        hostile_report(bad, sound(h1, hc1), "skewtest", k, hc1);
        let h2, hc2 = statistics::kurtosistest(xs[..], opt % 3);
        hostile_report(bad, sound(h2, hc2), "kurtosistest", k, hc2);
        let h3, hc3 = statistics::normaltest(xs[..]);
        hostile_report(bad, sound(h3, hc3), "normaltest", k, hc3);
        let h4, hc4 = statistics::jarque_bera(xs[..]);
        hostile_report(bad, sound(h4, hc4), "jarque_bera", k, hc4);
        let h5, hc5 = statistics::anderson_darling(xs[..]);
        hostile_report(bad, sound(h5, hc5), "anderson_darling", k, hc5);
        let lo, hi, c20 = statistics::mean_confidence_interval(xs[..], hostile_value(st));
        let ok20 = (c20 == 0 && lo <= hi) || (c20 != 0 && is_nan(lo) && is_nan(hi));
        hostile_report(bad, ok20, "mean_confidence_interval", k, c20);
    }
    return bad;
}
