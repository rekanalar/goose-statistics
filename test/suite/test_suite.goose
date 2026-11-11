// One test for each documented case of the hypothesis tests. Each case
// states the sample and the answer that SciPy gives.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_accept_tests() -> i64 {
    var bad: i64 = 0;
    let nan_v = statistics::nan();

    // The rank tests.
    let lo: f64[4] = [1.0, 2.0, 3.0, 4.0];
    let hi: f64[4] = [5.0, 6.0, 7.0, 8.5];
    let m1, mc1 = statistics::mannwhitneyu(lo[..], hi[..]);
    check(bad, mc1 == 0 && m1.statistic == 0.0 && near(m1.pvalue, 0.02857142857142857), "mann-whitney, separated samples");
    let m2, mc2 = statistics::mannwhitneyu(hi[..], lo[..], true, statistics::ALT_GREATER);
    check(bad, mc2 == 0 && near(m2.pvalue, 0.014285714285714285), "mann-whitney, greater alternative");
    let t1: f64[4] = [1.0, 2.0, 2.0, 4.0];
    let t2: f64[4] = [2.0, 6.0, 7.0, 8.5];
    let m3, mc3 = statistics::mannwhitneyu(t1[..], t2[..]);
    check(bad, mc3 == 0 && m3.statistic == 2.0 && near(m3.pvalue, 0.10375367752098565), "mann-whitney, a tie uses the normal form");
    let w1: f64[6] = [1.0, -2.0, 0.0, 3.0, 4.0, 5.5];
    let wr1, wc1 = statistics::wilcoxon(w1[..]);
    check(bad, wc1 == 0 && wr1.statistic == 2.0 && near(wr1.pvalue, 0.1875), "wilcoxon, a zero difference");
    let w2: f64[7] = [1.5, -2.1, 0.7, 3.3, 4.4, 5.5, -0.2];
    let wr2, wc2 = statistics::wilcoxon(w2[..]);
    check(bad, wc2 == 0 && wr2.statistic == 5.0 && near(wr2.pvalue, 0.15625), "wilcoxon, exact count");
    let kg: f64[9] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
    let ksz: i64[3] = [3, 3, 3];
    let kr, kc = statistics::kruskal(kg[..], ksz[..]);
    check(bad, kc == 0 && near(kr.statistic, 7.200000000000003) && near(kr.pvalue, 0.02732372244729252), "kruskal, three ordered groups");
    let kx: f64[4] = [0.1, 0.4, 0.7, 0.9];
    let ky: f64[4] = [0.2, 0.5, 0.6, 0.95];
    let ks, ksc = statistics::ks_1samp(kx[..], statistics::FAMILY_UNIFORM, [0.0, 1.0][..]);
    check(bad, ksc == 0 && ks.statistic <= 1.0 && ks.pvalue >= 0.0 && ks.pvalue <= 1.0, "ks one sample, uniform");
    let k2, k2c = statistics::ks_2samp(kx[..], ky[..]);
    check(bad, k2c == 0 && k2.pvalue > 0.0 && k2.pvalue <= 1.0, "ks two samples, small");
    let bad_ks, kbc = statistics::ks_1samp(kx[..], statistics::FAMILY_UNIFORM, [0.0][..]);
    check(bad, kbc == 6 && is_nan(bad_ks.pvalue), "ks one sample, wrong parameter count");

    // Shapiro-Wilk and the correlations.
    let sx: f64[8] = [1.1, 2.3, 3.7, 4.2, 5.9, 6.1, 7.4, 8.8];
    let sy: f64[8] = [2.0, 2.9, 3.1, 5.5, 6.6, 7.2, 9.1, 9.9];
    let sh, shc = statistics::shapiro(sx[..]);
    check(bad, shc == 0 && near(sh.statistic, 0.9800850002460284) && near(sh.pvalue, 0.9633237306489094), "shapiro, eight values");
    let pe, pec = statistics::pearsonr(sx[..], sy[..]);
    check(bad, pec == 0 && near(pe.statistic, 0.9776179086563659) && near(pe.pvalue, 2.7562778385532504e-05) && pe.df == 6.0, "pearson, eight pairs");
    let sp, spc = statistics::spearmanr(sx[..], sy[..]);
    check(bad, spc == 0 && sp.statistic == 1.0 && sp.pvalue == 0.0, "spearman, a perfect order");
    let ke, kec = statistics::kendalltau(sx[..], sy[..]);
    check(bad, kec == 0 && within(ke.statistic, 1.0, 4.0) && near(ke.pvalue, 4.96031746031746e-05), "kendall, exact p-value");
    let two: f64[2] = [1.0, 2.0];
    let s2, s2c = statistics::shapiro(two[..]);
    check(bad, s2c == 2 && is_nan(s2.statistic), "shapiro, two values");
    let p2, p2c = statistics::pearsonr(two[..], [3.0, 1.0][..]);
    check(bad, p2c == 0 && p2.statistic == -1.0 && p2.pvalue == 1.0, "pearson, two points");
    let same: f64[4] = [2.0, 2.0, 2.0, 2.0];
    let ps, psc = statistics::pearsonr(same[..], [1.0, 2.0, 3.0, 4.0][..]);
    check(bad, psc == 7 && is_nan(ps.pvalue), "pearson, a constant sample");
    let ss, ssc = statistics::shapiro(same[..]);
    check(bad, ssc == 7 && is_nan(ss.statistic), "shapiro, a constant sample");

    // The corrections.
    var adj: f64[>..] = [];
    let ps4: f64[4] = [0.01, 0.04, 0.03, 0.005];
    let n1, c1 = statistics::multipletests(adj, ps4[..], statistics::ADJUST_BONFERRONI);
    check(bad, c1 == 0 && n1 == 4 && near(adj[0], 0.04) && near(adj[3], 0.02), "bonferroni");
    let n2, c2 = statistics::multipletests(adj, ps4[..], statistics::ADJUST_HOLM);
    check(bad, c2 == 0 && near(adj[3], 0.02) && near(adj[0], 0.03) && near(adj[2], 0.06) && near(adj[1], 0.06), "holm");
    let n3, c3 = statistics::multipletests(adj, ps4[..], statistics::ADJUST_BH);
    check(bad, c3 == 0 && near(adj[3], 0.02) && near(adj[0], 0.02) && near(adj[2], 0.04) && near(adj[1], 0.04), "benjamini-hochberg");
    let bigp: f64[2] = [0.7, 0.9];
    let n4, c4 = statistics::multipletests(adj, bigp[..], statistics::ADJUST_BONFERRONI);
    check(bad, c4 == 0 && adj[0] == 1.0 && adj[1] == 1.0, "a corrected value clamps to 1");
    let outside: f64[2] = [0.5, 1.5];
    let n5, c5 = statistics::multipletests(adj, outside[..]);
    check(bad, c5 == 5 && n5 == 0 && adj.len == 0, "a p-value above 1");
    let n6, c6 = statistics::multipletests(adj, ps4[..], 9);
    check(bad, c6 == 5, "an unknown correction");

    let alt3, ac3 = statistics::ttest_1samp(sx[..], 0.0, 3);
    check(bad, ac3 == 5 && is_nan(alt3.pvalue), "t test, alternative 3");
    let alt4, ac4 = statistics::ttest_1samp(sx[..], 0.0, -1);
    check(bad, ac4 == 5 && is_nan(alt4.pvalue), "t test, alternative -1");

    // The interval of the mean, a constant sample and the zero margins.
    let ci: f64[5] = [1.1, 2.3, 3.7, 4.2, 5.9];
    let ilo, ihi, ic = statistics::mean_confidence_interval(ci[..], 0.95);
    check(bad, ic == 0 && within(ilo, 1.1612843741312777, 16.0) && within(ihi, 5.7187156258687235, 16.0), "interval of the mean");
    let ib, ib2, ibc = statistics::mean_confidence_interval(ci[..], 1.0);
    check(bad, ibc == 5, "interval, a confidence of 1");
    let same3: f64[3] = [2.0, 2.0, 2.0];
    let c1, c1c = statistics::ttest_1samp(same3[..], 1.0);
    check(bad, c1c == 0 && is_inf(c1.statistic) && c1.pvalue == 0.0, "t test, a constant sample with another mean");
    let c2, c2c = statistics::ttest_1samp(same3[..], 2.0);
    check(bad, c2c == 7 && is_nan(c2.pvalue), "t test, a constant sample with the same mean");
    let zm: f64[4] = [0.0, 0.0, 3.0, 4.0];
    let fz, fzc = statistics::fisher_exact(zm[..]);
    check(bad, fzc == 0 && is_nan(fz.statistic) && fz.pvalue == 1.0, "fisher, a zero margin");
    let cz, czc = statistics::chi2_contingency(zm[..], 2, 2);
    check(bad, czc == 7 && is_nan(cz.pvalue), "chi-square, a zero margin");
    let bf, bfc = statistics::binomtest(2.5, 10.0, 0.5);
    check(bad, bfc == 5 && is_nan(bf.pvalue), "binomial, a fractional count");
    var one: f64[>..] = [];
    let n7, c7 = statistics::multipletests(one, [0.3][..], statistics::ADJUST_HOLM);
    check(bad, c7 == 0 && n7 == 1 && one[0] == 0.3, "a single p-value");

    // A sample above the input bound gives the too long fault in every call
    // and no work starts.
    var long: f64[>..] = [];
    for i in 0..statistics::MAX_INPUT + 1 { long.push(i); }
    let short: f64[2] = [1.0, 2.0];
    let sizes1: i64[1] = [1];
    let l1, lc1 = statistics::ttest_1samp(long[..]);
    let l2, lc2 = statistics::ttest_ind(long[..], short[..]);
    let l3, lc3 = statistics::ttest_rel(long[..], long[..]);
    let l4, lc4 = statistics::f_oneway(long[..], sizes1[..]);
    let l5, lc5 = statistics::kruskal(long[..], sizes1[..]);
    let l6, lc6 = statistics::chisquare(long[..], long[..]);
    let l7, lc7 = statistics::mannwhitneyu(long[..], short[..]);
    let l8, lc8 = statistics::wilcoxon(long[..]);
    let l9, lc9 = statistics::ks_2samp(long[..], short[..]);
    let l10, lc10 = statistics::ks_1samp(long[..], statistics::FAMILY_NORMAL, [0.0, 1.0][..]);
    let l11, lc11 = statistics::shapiro(long[..]);
    let l12, lc12 = statistics::pearsonr(long[..], long[..]);
    let l13, lc13 = statistics::spearmanr(long[..], long[..]);
    let l14, lc14 = statistics::kendalltau(long[..], long[..]);
    let l15, lc15 = statistics::multipletests(one, long[..]);
    let l16, l16b, lc16 = statistics::mean_confidence_interval(long[..]);
    check(bad, lc1 == 9 && lc2 == 9 && lc3 == 9 && lc4 == 9 && lc5 == 9 && lc6 == 9 && lc7 == 9 && lc8 == 9, "too long input, the first eight calls");
    check(bad, lc9 == 9 && lc10 == 9 && lc11 == 9 && lc12 == 9 && lc13 == 9 && lc14 == 9 && lc15 == 9 && lc16 == 9, "too long input, the last eight calls");

    // The faults of every call, with NaN, an empty sample and a bad option.
    let withnan: f64[4] = [1.0, nan_v, 3.0, 4.0];
    let f1, fc1 = statistics::shapiro(withnan[..]);
    check(bad, fc1 == 3 && is_nan(f1.pvalue), "shapiro, NaN");
    let f2, fc2 = statistics::kendalltau(withnan[..], sy[0..4]);
    check(bad, fc2 == 3 && is_nan(f2.pvalue), "kendall, NaN");
    let f3, fc3 = statistics::mannwhitneyu(sx[..], sy[..], true, 7);
    check(bad, fc3 == 5 && is_nan(f3.pvalue), "mann-whitney, bad alternative");
    let f4, fc4 = statistics::spearmanr(sx[0..3], sy[..]);
    check(bad, fc4 == 6 && is_nan(f4.pvalue), "spearman, unequal lengths");
    let empty: f64[0] = [];
    let f5, fc5 = statistics::kruskal(empty[..], [1, 1][..]);
    check(bad, fc5 == 1 && is_nan(f5.pvalue), "kruskal, empty buffer");
    bad += run_value_tests();
    return bad;
}

// The documented samples of the t tests, the count tests and the
// Kolmogorov-Smirnov tests, with the values SciPy 1.18.1 gives.
fn run_value_tests() -> i64 {
    var bad: i64 = 0;
    let tol = 1e-13;
    let a: f64[5] = [2.1, 2.5, 3.2, 2.8, 3.0];
    let r1, c1 = statistics::ttest_1samp(a[..], 2.0);
    check(bad, c1 == 0 && near_tol(r1.statistic, 3.7230313679641367, tol) && r1.df == 4.0 && near_tol(r1.pvalue, 0.020420477739396304, tol), "t test, one sample two-sided");
    let r2, c2 = statistics::ttest_1samp(a[..], 2.0, statistics::ALT_LESS);
    check(bad, c2 == 0 && near_tol(r2.pvalue, 0.9897897611303018, tol), "t test, one sample less");
    let r3, c3 = statistics::ttest_1samp(a[..], 2.0, statistics::ALT_GREATER);
    check(bad, c3 == 0 && near_tol(r3.pvalue, 0.010210238869698152, tol), "t test, one sample greater");
    let x4: f64[4] = [1.0, 2.0, 3.0, 4.0];
    let y5: f64[5] = [2.0, 4.0, 6.0, 8.0, 10.0];
    let r4, c4 = statistics::ttest_ind(x4[..], y5[..], statistics::ALT_TWO_SIDED, true);
    check(bad, c4 == 0 && near_tol(r4.statistic, -2.0578065752724592, tol) && r4.df == 7.0 && near_tol(r4.pvalue, 0.0786192350586938, tol), "t test, Student two samples");
    let r5, c5 = statistics::ttest_ind(x4[..], y5[..], statistics::ALT_TWO_SIDED, false);
    check(bad, c5 == 0 && near_tol(r5.statistic, -2.2514363231593695, tol) && near_tol(r5.df, 5.520787746170677, tol) && near_tol(r5.pvalue, 0.06913359319239237, tol), "t test, Welch");
    let p1: f64[5] = [1.0, 2.0, 3.0, 4.0, 5.0];
    let p2: f64[5] = [1.5, 2.4, 3.9, 4.1, 5.9];
    let r6, c6 = statistics::ttest_rel(p1[..], p2[..]);
    check(bad, c6 == 0 && near_tol(r6.statistic, -3.6452895074630653, tol) && near_tol(r6.pvalue, 0.02186221514033323, tol), "t test, paired");
    let g: f64[9] = [1.0, 2.0, 3.0, 2.0, 4.0, 6.0, 5.0, 6.0, 7.5];
    let gs: i64[3] = [3, 3, 3];
    let r7, c7 = statistics::f_oneway(g[..], gs[..]);
    check(bad, c7 == 0 && near_tol(r7.statistic, 5.936708860759495, tol) && near_tol(r7.pvalue, 0.037829529067205485, tol), "analysis of variance, three groups");
    let one: f64[1] = [2.0];
    let r8, c8 = statistics::ttest_1samp(one[..], 2.0);
    check(bad, c8 == 2 && is_nan(r8.pvalue), "t test, one value");
    let cs: f64[3] = [2.0, 2.0, 2.0];
    let r9, c9 = statistics::ttest_1samp(cs[..], 1.0);
    check(bad, c9 == 0 && is_inf(r9.statistic) && r9.pvalue == 0.0, "t test, constant sample, other mean");
    let r10, c10 = statistics::ttest_1samp(cs[..], 2.0);
    check(bad, c10 == 7 && is_nan(r10.pvalue), "t test, constant sample, same mean");

    let none: f64[0] = [];
    let ob: f64[6] = [16.0, 18.0, 16.0, 14.0, 12.0, 12.0];
    let q1, d1 = statistics::chisquare(ob[..], none[..], 0);
    check(bad, d1 == 0 && near_tol(q1.statistic, 2.0, tol) && near_tol(q1.pvalue, 0.8491450360846096, tol), "chi-square, fit");
    let tb: f64[4] = [10.0, 20.0, 30.0, 25.0];
    let q2, d2 = statistics::chi2_contingency(tb[..], 2, 2, true);
    check(bad, d2 == 0 && near_tol(q2.pvalue, 0.0999616438735349, tol), "chi-square, table with the correction");
    let q3, d3 = statistics::chi2_contingency(tb[..], 2, 2, false);
    check(bad, d3 == 0 && near_tol(q3.pvalue, 0.06115089757606777, tol), "chi-square, table without the correction");
    let b1, e1 = statistics::binomtest(7.0, 10.0, 0.5);
    check(bad, e1 == 0 && near_tol(b1.pvalue, 0.34375, tol), "binomial test, two-sided");
    let b2, e2 = statistics::binomtest(7.0, 10.0, 0.5, statistics::ALT_LESS);
    check(bad, e2 == 0 && near_tol(b2.pvalue, 0.9453125, tol), "binomial test, less");
    let b3, e3 = statistics::binomtest(7.0, 10.0, 0.5, statistics::ALT_GREATER);
    check(bad, e3 == 0 && near_tol(b3.pvalue, 0.171875, tol), "binomial test, greater");
    let ft: f64[4] = [3.0, 1.0, 1.0, 3.0];
    let f1, h1 = statistics::fisher_exact(ft[..]);
    check(bad, h1 == 0 && near_tol(f1.statistic, 9.0, tol) && near_tol(f1.pvalue, 0.48571428571428565, tol), "fisher, two-sided");
    let f2, h2 = statistics::fisher_exact(ft[..], statistics::ALT_LESS);
    check(bad, h2 == 0 && near_tol(f2.pvalue, 0.9857142857142858, tol), "fisher, less");
    let f3, h3 = statistics::fisher_exact(ft[..], statistics::ALT_GREATER);
    check(bad, h3 == 0 && near_tol(f3.pvalue, 0.24285714285714283, tol), "fisher, greater");
    let zt: f64[4] = [0.0, 0.0, 3.0, 4.0];
    let f4, h4 = statistics::fisher_exact(zt[..]);
    check(bad, h4 == 0 && f4.pvalue == 1.0 && is_nan(f4.statistic), "fisher, zero row");
    // The expected counts sum to the observed sum, one of them is 0.
    let ez: f64[3] = [10.0, 0.0, 10.0];
    let ob3: f64[3] = [5.0, 5.0, 10.0];
    let q4, d4 = statistics::chisquare(ob3[..], ez[..], 0);
    check(bad, d4 == 7 && is_nan(q4.pvalue), "chi-square, zero expected count, equal sums");
    let ez2: f64[3] = [10.0, 0.0, 5.0];
    let q5, d5 = statistics::chisquare(ob3[..], ez2[..], 0);
    check(bad, d5 == 5 && is_nan(q5.pvalue), "chi-square, zero expected count, different sums");

    let ux: f64[5] = [0.1, 0.4, 0.7, 0.2, 0.95];
    let k1, kc1 = statistics::ks_1samp(ux[..], statistics::FAMILY_UNIFORM, [0.0, 1.0][..]);
    check(bad, kc1 == 0 && near_tol(k1.statistic, 0.2, tol) && near_tol(k1.pvalue, 0.9616, 1e-4), "ks one sample, documented value");
    let kx: f64[4] = [1.0, 2.0, 3.0, 4.5];
    let ky: f64[4] = [2.5, 6.0, 7.0, 8.0];
    let k2, kc2 = statistics::ks_2samp(kx[..], ky[..]);
    check(bad, kc2 == 0 && near_tol(k2.statistic, 0.75, tol) && near_tol(k2.pvalue, 0.22857142857142862, tol), "ks two samples, documented value");
    // Above about 78,000 values the product of the two pair counts passes
    // 2 to the power 63. A perfect order must still give tau 1.
    var big_x: f64[>..] = [];
    for i in 0..80000 { big_x.push(i); }
    let kb, kc = statistics::kendalltau(big_x[..], big_x[..]);
    check(bad, kc == 0 && near_tol(kb.statistic, 1.0, 1e-12), "Kendall tau of a perfect order of 80000 values");
    let kt, kd = statistics::kendall(big_x[..], big_x[..], 0);
    check(bad, kd == 0 && near_tol(kt, 1.0, 1e-12), "Kendall tau-b of a perfect order of 80000 values");
    return bad;
}
