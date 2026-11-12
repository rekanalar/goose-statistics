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
    if c < 1 || c > 10 { return false; }
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
        // A clean list, so the fault comes from the method and not from
        // the sample.
        let clean: f64[1] = [0.5];
        let m2, c23 = statistics::multipletests(out, clean[..], 9);
        hostile_report(bad, c23 == 5 && m2 == 0 && out.len == 0, "multipletests, an unknown method", k, c23);
        let m3, c21 = statistics::false_discovery_control(out, xs[..], opt + 2);
        let ok21 = (c21 == 0 && m3 == xs.len && out.len == m3) || (c21 != 0 && m3 == 0 && out.len == 0);
        hostile_report(bad, ok21, "false_discovery_control", k, c21);
        let m4, c22 = statistics::false_discovery_control(out, clean[..], statistics::ADJUST_HOLM);
        hostile_report(bad, c22 == 5 && m4 == 0 && out.len == 0, "false_discovery_control, another method", k, c22);
        let v1, vc1 = statistics::combine_pvalues(xs[..], opt % 3);
        hostile_report(bad, sound(v1, vc1), "combine_pvalues", k, vc1);
        var wts: f64[>..] = [];
        for t in 0..xs.len { wts.push(abs(hostile_value(st))); }
        let v2, vc2 = statistics::combine_pvalues(xs[..], statistics::COMBINE_STOUFFER, wts[..]);
        hostile_report(bad, sound(v2, vc2), "combine_pvalues with weights", k, vc2);
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
        // The variance and rank tests, the transforms and the survival test.
        let w1, wc1 = statistics::levene(xs[..], sizes[..], opt % 3);
        hostile_report(bad, sound(w1, wc1), "levene", k, wc1);
        let w2, wc2 = statistics::fligner(xs[..], sizes[..], opt % 3);
        hostile_report(bad, sound(w2, wc2), "fligner", k, wc2);
        let w3, wc3 = statistics::bartlett(xs[..], sizes[..]);
        hostile_report(bad, sound(w3, wc3), "bartlett", k, wc3);
        let w4, wc4 = statistics::alexandergovern(xs[..], sizes[..]);
        hostile_report(bad, sound(w4, wc4), "alexandergovern", k, wc4);
        let w5, wc5 = statistics::ranksums(xs[..], ys[..], alt);
        hostile_report(bad, sound(w5, wc5), "ranksums", k, wc5);
        let w6, wc6 = statistics::mood(xs[..], ys[..], alt);
        hostile_report(bad, sound(w6, wc6), "mood", k, wc6);
        let w7, wc7 = statistics::brunnermunzel(xs[..], ys[..], alt, flag);
        hostile_report(bad, sound(w7, wc7), "brunnermunzel", k, wc7);
        let w8, wc8 = statistics::ansari(xs[..], ys[..], alt, opt % 3);
        hostile_report(bad, sound(w8, wc8), "ansari", k, wc8);
        let w9, wc9 = statistics::friedmanchisquare(xs[..], opt);
        hostile_report(bad, sound(w9, wc9), "friedmanchisquare", k, wc9);
        var tabout: f64[>..] = [];
        let w10, wc10 = statistics::median_test(tabout, xs[..], sizes[..], opt % 3, flag);
        hostile_report(bad, sound(w10, wc10), "median_test", k, wc10);
        var ci: f64[>..] = [];
        let w11, wc11 = statistics::quantile_test(ci, xs[..], hostile_value(st), hostile_value(st), alt, hostile_value(st));
        hostile_report(bad, sound(w11, wc11), "quantile_test", k, wc11);
        let w12, wc12 = statistics::ttest_ind_from_stats(hostile_value(st), hostile_value(st), hostile_value(st), hostile_value(st), hostile_value(st), hostile_value(st), flag, alt);
        hostile_report(bad, sound(w12, wc12), "ttest_ind_from_stats", k, wc12);
        var tout: f64[>..] = [];
        let tl, tc1 = statistics::boxcox(tout, xs[..], hostile_value(st));
        let okt1 = (tc1 == 0 && tout.len == xs.len) || (tc1 != 0 && is_nan(tl) && tout.len == 0);
        hostile_report(bad, okt1, "boxcox", k, tc1);
        let tl2, tc2 = statistics::yeojohnson(tout, xs[..], hostile_value(st));
        let okt2 = (tc2 == 0 && tout.len == xs.len) || (tc2 != 0 && is_nan(tl2) && tout.len == 0);
        hostile_report(bad, okt2, "yeojohnson", k, tc2);
        let tl3, tc3 = statistics::boxcox_normmax(xs[..], opt);
        hostile_report(bad, (tc3 == 0 && !is_nan(tl3)) || (tc3 != 0 && is_nan(tl3) && tc3 <= 10), "boxcox_normmax", k, tc3);
        let tl4, tc4 = statistics::yeojohnson_normmax(xs[..]);
        hostile_report(bad, (tc4 == 0 && !is_nan(tl4)) || (tc4 != 0 && is_nan(tl4) && tc4 <= 10), "yeojohnson_normmax", k, tc4);
        // The exact tests, the distance tests and the comparisons.
        let x1, xc1 = statistics::barnard_exact(abs(hostile_value(st)), abs(hostile_value(st)), abs(hostile_value(st)), abs(hostile_value(st)), alt, flag, opt);
        hostile_report(bad, sound(x1, xc1), "barnard_exact", k, xc1);
        let x2, xc2 = statistics::boschloo_exact(abs(hostile_value(st)), abs(hostile_value(st)), abs(hostile_value(st)), abs(hostile_value(st)), alt, opt);
        hostile_report(bad, sound(x2, xc2), "boschloo_exact", k, xc2);
        var kout: f64[>..] = [];
        let x3, xc3 = statistics::kstest(kout, xs[..], family, params[..], alt);
        hostile_report(bad, sound(x3, xc3), "kstest", k, xc3);
        let x4, xc4 = statistics::kstest_2samp(kout, xs[..], ys[..], alt);
        hostile_report(bad, sound(x4, xc4), "kstest_2samp", k, xc4);
        let x5, xc5 = statistics::cramervonmises(xs[..], family, params[..]);
        hostile_report(bad, sound(x5, xc5), "cramervonmises", k, xc5);
        let x6, xc6 = statistics::cramervonmises_2samp(xs[..], ys[..], opt % 3);
        hostile_report(bad, sound(x6, xc6), "cramervonmises_2samp", k, xc6);
        let x7, xc7 = statistics::epps_singleton_2samp(xs[..], ys[..]);
        hostile_report(bad, sound(x7, xc7), "epps_singleton_2samp", k, xc7);
        var gen = rng { s: 5 };
        let x8, xc8 = statistics::bws_test(xs[..], ys[..], gen, alt, 20);
        hostile_report(bad, sound(x8, xc8), "bws_test", k, xc8);
        var aout: f64[>..] = [];
        let x9, xc9 = statistics::anderson_ksamp(aout, xs[..], sizes[..], opt % 3);
        hostile_report(bad, sound(x9, xc9), "anderson_ksamp", k, xc9);
        let x10, xc10 = statistics::poisson_means_test(abs(hostile_value(st)), abs(hostile_value(st)), abs(hostile_value(st)), abs(hostile_value(st)), hostile_value(st), alt);
        hostile_report(bad, sound(x10, xc10), "poisson_means_test", k, xc10);
        let x11, xc11 = statistics::page_trend_test(xs[..], opt, flag, params[..], opt % 3);
        hostile_report(bad, sound(x11, xc11), "page_trend_test", k, xc11);
        var hout: f64[>..] = [];
        let x12, xc12 = statistics::tukey_hsd(hout, xs[..], sizes[..], hostile_value(st), flag);
        let okx12 = (xc12 == 0 && hout.len % 4 == 0) || (xc12 != 0 && hout.len == 0);
        hostile_report(bad, okx12, "tukey_hsd", k, xc12);
        var dout: f64[>..] = [];
        let x13, xc13 = statistics::dunnett(dout, xs[..], sizes[..], 0, alt, hostile_value(st));
        let okx13 = (xc13 == 0 && dout.len % 4 == 0) || (xc13 != 0 && dout.len == 0);
        hostile_report(bad, okx13, "dunnett", k, xc13);
        let x14, xc14 = statistics::bayes_mvs(xs[..], hostile_value(st));
        let okx14 = xc14 == 0 || (xc14 != 0 && is_nan(x14.mean) && is_nan(x14.variance) && is_nan(x14.deviation));
        hostile_report(bad, okx14, "bayes_mvs", k, xc14);
        let x15, xc15 = statistics::mvsdist(xs[..]);
        let okx15 = xc15 == 0 || (xc15 != 0 && is_nan(x15.mean) && is_nan(x15.variance) && is_nan(x15.deviation));
        hostile_report(bad, okx15, "mvsdist", k, xc15);
        // The resampling calls. The block returns the mean of its list.
        var gen2 = rng { s: 9 };
        let rm = statistics::ResampleMethod { resamples: 20, interval: statistics::INTERVAL_PERCENTILE, .. };
        var bout: f64[>..] = [];
        let x16, xc16 = statistics::bootstrap(bout, xs[..], sizes[..], rm, hostile_value(st), alt, gen2) { s => mean_of(s) };
        let okx16 = (xc16 == 0 && bout.len == 20) || (xc16 != 0 && bout.len == 0);
        hostile_report(bad, okx16, "bootstrap", k, xc16);
        var gen3 = rng { s: 11 };
        var pout: f64[>..] = [];
        let x17, xc17 = statistics::permutation_test(pout, xs[..], sizes[..], rm, alt, gen3) { s => mean_of(s) };
        let okx17 = (xc17 == 0 && pout.len == 20) || (xc17 != 0 && pout.len == 0);
        hostile_report(bad, okx17, "permutation_test", k, xc17);
        var gen4 = rng { s: 13 };
        var mout: f64[>..] = [];
        let x18, xc18 = statistics::monte_carlo_test(mout, xs[..], family, params[..], rm, alt, gen4) { s => mean_of(s) };
        let okx18 = (xc18 == 0 && mout.len == 20) || (xc18 != 0 && mout.len == 0);
        hostile_report(bad, okx18, "monte_carlo_test", k, xc18);
        var gen5 = rng { s: 17 };
        let x19, xc19 = statistics::power(family, params[..], 10, hostile_value(st), 20, gen5) { s => mean_of(s) };
        let okx19 = xc19 == 0 || (xc19 != 0 && is_nan(x19));
        hostile_report(bad, okx19, "power", k, xc19);
        var gen6 = rng { s: 19 };
        var sfirst: f64[>..] = [];
        var stotal: f64[>..] = [];
        let x20, xc20 = statistics::sobol_indices(sfirst, stotal, 2, 16, gen6) { p => p[0] + p[1] };
        let okx20 = (xc20 == 0 && x20 == 2 && sfirst.len == 2 && stotal.len == 2) || (xc20 != 0 && sfirst.len == 0 && stotal.len == 0);
        hostile_report(bad, okx20, "sobol_indices", k, xc20);
        // The fit and the survival test.
        var gen7 = rng { s: 23 };
        var gout: f64[>..] = [];
        let x21, xc21 = statistics::goodness_of_fit(gout, xs[..], family, opt % 3, 20, gen7);
        let okx21 = (xc21 == 0 && gout.len > 0) || (xc21 != 0 && gout.len == 0);
        hostile_report(bad, okx21, "goodness_of_fit", k, xc21);
        var empty: f64[0] = [];
        let cd1 = statistics::CensoredData { uncensored: xs[..], left: empty[..], right: ys[..], interval: empty[..] };
        let cd2 = statistics::CensoredData { uncensored: ys[..], left: empty[..], right: xs[..], interval: empty[..] };
        let x22, xc22 = statistics::logrank(cd1, cd2, alt);
        hostile_report(bad, sound(x22, xc22), "logrank", k, xc22);
        let lo, hi, c20 = statistics::mean_confidence_interval(xs[..], hostile_value(st));
        let ok20 = (c20 == 0 && lo <= hi) || (c20 != 0 && is_nan(lo) && is_nan(hi));
        hostile_report(bad, ok20, "mean_confidence_interval", k, c20);
    }
    return bad;
}
