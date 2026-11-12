// The tests of the variance calls: the measured oracle rows, the tied
// rows, the faults and the no-write checks. The random rows live in
// variance_vectors.goose, which the generator test/make_variance_vectors.py
// writes.
namespace suite;

import std;
import math;
import statistics from "../..";

// One result against the measured row. A want of NaN skips that field.
fn vs_res(bad: i64&, name: u8[:], r: statistics::TestResult, c: i64, stat: f64, p: f64, df: f64, df2: f64, tol: f64) {
    var ok = c == 0 && near_tol(r.statistic, stat, tol) && near_tol(r.pvalue, p, tol);
    if !is_nan(df) { ok = ok && near_tol(r.df, df, 1e-12); }
    if !is_nan(df2) { ok = ok && near_tol(r.df2, df2, 1e-12); }
    check(bad, ok, name);
}

// The calls of several groups, by number, so one loop can test them all.
fn vs_groups(which: i64, groups: const f64[:], sizes: const i64[:]) -> statistics::TestResult, i64 {
    if which == 0 { return statistics::levene(groups, sizes); }
    if which == 1 { return statistics::fligner(groups, sizes); }
    if which == 2 { return statistics::bartlett(groups, sizes); }
    if which == 3 { return statistics::alexandergovern(groups, sizes); }
    if which == 4 { return statistics::levene(groups, sizes, statistics::CENTER_MEAN); }
    if which == 5 { return statistics::fligner(groups, sizes, statistics::CENTER_TRIMMED, 0.2); }
    var out: f64[>..] = [];
    return statistics::median_test(out, groups, sizes);
}

// The calls of two samples, by number.
fn vs_two(which: i64, xs: const f64[:], ys: const f64[:]) -> statistics::TestResult, i64 {
    if which == 0 { return statistics::ranksums(xs, ys); }
    if which == 1 { return statistics::mood(xs, ys); }
    if which == 2 { return statistics::brunnermunzel(xs, ys); }
    return statistics::ansari(xs, ys);
}

fn run_variance_tests() -> i64 {
    var bad: i64 = 0;
    let a: f64[8] = [2.1, 3.4, 1.9, 5.6, 4.4, 3.8, 2.7, 4.9];
    let b: f64[8] = [3.9, 5.2, 4.8, 6.1, 5.5, 4.1, 6.6, 5.0];
    let c3: f64[8] = [1.2, 2.5, 2.2, 3.1, 1.8, 2.9, 2.4, 3.3];
    var g: f64[>..] = [];
    for v in a { g.push(v); }
    for v in b { g.push(v); }
    for v in c3 { g.push(v); }
    let sz: i64[3] = [8, 8, 8];
    let tx: f64[5] = [1.0, 2.0, 2.0, 3.0, 4.0];
    let ty: f64[6] = [2.0, 3.0, 5.0, 6.0, 6.0, 7.0];
    let nan = statistics::nan();
    let inf = 1.0 / 0.0;

    // ---- the measured oracle rows ----
    var r, k = statistics::levene(g[..], sz[..], statistics::CENTER_MEAN);
    vs_res(bad, "levene, mean", r, k, 2.0972332015810293, 0.14777411401301205, 2.0, 21.0, 1e-12);
    r, k = statistics::levene(g[..], sz[..]);
    vs_res(bad, "levene, median is the default", r, k, 2.088976377952756, 0.14879496808938358, 2.0, 21.0, 1e-12);
    r, k = statistics::levene(g[..], sz[..], statistics::CENTER_TRIMMED);
    vs_res(bad, "levene, trimmed cuts nothing of 8 values", r, k, 2.0972332015810293, 0.14777411401301205, 2.0, 21.0, 1e-12);
    r, k = statistics::levene(g[..], sz[..], statistics::CENTER_TRIMMED, 0.25);
    vs_res(bad, "levene, trimmed at 0.25", r, k, 2.0791536050156756, 0.15001950441996406, 2.0, 21.0, 1e-12);
    r, k = statistics::fligner(g[..], sz[..], statistics::CENTER_MEAN);
    vs_res(bad, "fligner, mean", r, k, 4.558521167077154, 0.10235986530988447, 2.0, nan, 1e-12);
    r, k = statistics::fligner(g[..], sz[..]);
    vs_res(bad, "fligner, median", r, k, 4.13604262718082, 0.12643571098066891, 2.0, nan, 1e-12);
    r, k = statistics::fligner(g[..], sz[..], statistics::CENTER_TRIMMED);
    vs_res(bad, "fligner, trimmed cuts nothing of 8 values", r, k, 4.558521167077154, 0.10235986530988447, 2.0, nan, 1e-12);
    // SciPy gives 4.507825712411529 here. Its trimmed mean of the first group
    // is 3.5749999999999997, one rounding step below the exact 3.575. The two
    // centers break a tie of the deviations in different ways.
    r, k = statistics::fligner(g[..], sz[..], statistics::CENTER_TRIMMED, 0.25);
    vs_res(bad, "fligner, trimmed at 0.25 with the exact center", r, k, 4.6973335774878855, 0.09549639425084107, 2.0, nan, 1e-12);
    r, k = statistics::bartlett(g[..], sz[..]);
    vs_res(bad, "bartlett", r, k, 2.744272211811164, 0.2535647395017535, 2.0, nan, 1e-12);
    r, k = statistics::alexandergovern(g[..], sz[..]);
    vs_res(bad, "alexandergovern", r, k, 18.71737507059921, 8.621317624550678e-05, 2.0, nan, 1e-12);
    r, k = statistics::brunnermunzel(a[..], b[..]);
    vs_res(bad, "brunnermunzel, t", r, k, 2.9698484809834995, 0.013032630174299448, 10.761527070103662, nan, 1e-12);
    r, k = statistics::brunnermunzel(a[..], b[..], statistics::ALT_TWO_SIDED, false);
    vs_res(bad, "brunnermunzel, normal", r, k, 2.9698484809834995, 0.0029794666563329845, nan, nan, 1e-12);
    r, k = statistics::mood(a[..], b[..]);
    vs_res(bad, "mood, no ties", r, k, 0.6149400462680908, 0.5385943148264352, nan, nan, 1e-12);
    r, k = statistics::mood(tx[..], ty[..]);
    vs_res(bad, "mood, with ties", r, k, -0.2893493147793367, 0.7723140776921694, nan, nan, 1e-12);
    r, k = statistics::ansari(a[..], b[..], statistics::ALT_TWO_SIDED, statistics::MW_EXACT);
    vs_res(bad, "ansari, exact", r, k, 33.0, 0.6060606060606061, nan, nan, 1e-12);
    r, k = statistics::ansari(a[..], b[..], statistics::ALT_TWO_SIDED, statistics::MW_ASYMPTOTIC);
    vs_res(bad, "ansari, asymptotic", r, k, 33.0, 0.5261684777561433, nan, nan, 1e-12);
    r, k = statistics::ansari(a[..], b[..]);
    vs_res(bad, "ansari, auto takes the exact form", r, k, 33.0, 0.6060606060606061, nan, nan, 1e-12);
    r, k = statistics::ansari(tx[..], ty[..]);
    vs_res(bad, "ansari, auto with ties takes the normal form", r, k, 17.5, 0.6645642667502323, nan, nan, 1e-12);
    r, k = statistics::ranksums(a[..], b[..]);
    vs_res(bad, "ranksums", r, k, -2.2054411323441157, 0.027423154381354665, nan, nan, 1e-12);
    r, k = statistics::ranksums(tx[..], ty[..]);
    vs_res(bad, "ranksums, ties get no correction", r, k, -1.9170289512680814, 0.05523425371806383, nan, nan, 1e-12);
    var d: f64[>..] = [];
    for i in 0..8 { d.push(a[i]); d.push(b[i]); d.push(c3[i]); }
    r, k = statistics::friedmanchisquare(d[..], 3);
    vs_res(bad, "friedmanchisquare", r, k, 14.25, 0.000804733010124613, 2.0, nan, 1e-12);
    let dt: f64[12] = [1.0, 2.0, 3.0, 2.0, 2.0, 1.0, 2.0, 3.0, 1.0, 4.0, 5.0, 2.0];
    r, k = statistics::friedmanchisquare(dt[..], 3);
    vs_res(bad, "friedmanchisquare, ties", r, k, 2.8, 0.24659696394160646, 2.0, nan, 1e-12);
    var out: f64[>..] = [];
    r, k = statistics::median_test(out, g[..], sz[..]);
    vs_res(bad, "median_test", r, k, 16.0, 0.0003354626279025119, 2.0, 3.5999999999999996, 1e-12);
    let tab: f64[6] = [4.0, 8.0, 0.0, 4.0, 0.0, 8.0];
    check(bad, near_list(out[..], tab[..], 0.0), "median_test, the table");
    r, k = statistics::median_test(out, g[..], sz[..], statistics::MEDIAN_ABOVE, false);
    vs_res(bad, "median_test, above and no correction", r, k, 16.0, 0.0003354626279025119, 2.0, 3.5999999999999996, 1e-12);
    r, k = statistics::median_test(out, g[..], sz[..], statistics::MEDIAN_IGNORE);
    check(bad, k == 0 && near_list(out[..], tab[..], 0.0), "median_test, ignore");
    r, k = statistics::quantile_test(out, a[..], 3.5, 0.5);
    vs_res(bad, "quantile_test", r, k, 4.0, 1.0, 1.0, nan, 1e-12);
    let ci95: f64[2] = [1.9, 5.6];
    check(bad, near_list(out[..], ci95[..], 0.0), "quantile_test, interval at 0.95");
    r, k = statistics::quantile_test(out, a[..], 3.5, 0.5, statistics::ALT_TWO_SIDED, 0.9);
    let ci90: f64[2] = [2.1, 4.9];
    check(bad, k == 0 && near_list(out[..], ci90[..], 0.0), "quantile_test, interval at 0.9");
    r, k = statistics::quantile_test(out, a[..], 3.5, 0.5, statistics::ALT_LESS, 0.9);
    check(bad, k == 0 && out.len == 2 && out[0] == -inf && out[1] == 4.9, "quantile_test, less interval");
    r, k = statistics::quantile_test(out, a[..], 3.5, 0.5, statistics::ALT_GREATER, 0.9);
    check(bad, k == 0 && out.len == 2 && out[0] == 2.1 && out[1] == inf, "quantile_test, greater interval");
    r, k = statistics::ttest_ind_from_stats(3.0, 1.2, 8, 4.5, 1.4, 9);
    vs_res(bad, "ttest_ind_from_stats, pooled", r, k, -2.3556227868485973, 0.03252064750046273, 15.0, nan, 1e-12);
    r, k = statistics::ttest_ind_from_stats(3.0, 1.2, 8, 4.5, 1.4, 9, false);
    vs_res(bad, "ttest_ind_from_stats, Welch", r, k, -2.378323901617605, 0.03112605195627196, 14.98793812021785, nan, 1e-12);

    // ---- the summary form equals the raw form ----
    let m1, e1 = statistics::mean(a[..]);
    let m2, e2 = statistics::mean(b[..]);
    let s1, e3 = statistics::stdev(a[..]);
    let s2, e4 = statistics::stdev(b[..]);
    for eq in 0..2 {
        let raw, kr = statistics::ttest_ind(a[..], b[..], statistics::ALT_GREATER, eq == 0);
        let sm, ks = statistics::ttest_ind_from_stats(m1, s1, 8, m2, s2, 8, eq == 0, statistics::ALT_GREATER);
        check(bad, kr == 0 && ks == 0 && near_tol(sm.statistic, raw.statistic, 1e-12)
            && near_tol(sm.pvalue, raw.pvalue, 1e-12) && near_tol(sm.df, raw.df, 1e-12), "ttest_ind_from_stats equals ttest_ind");
    }

    // ---- the one-sided forms ----
    let rl, kl = statistics::ranksums(a[..], b[..], statistics::ALT_LESS);
    let rg, kg = statistics::ranksums(a[..], b[..], statistics::ALT_GREATER);
    check(bad, kl == 0 && kg == 0 && near_tol(rl.pvalue + rg.pvalue, 1.0, 1e-12), "ranksums, the one-sided tails add to 1");
    let bl, kbl = statistics::brunnermunzel(a[..], b[..], statistics::ALT_LESS);
    let bg, kbg = statistics::brunnermunzel(a[..], b[..], statistics::ALT_GREATER);
    check(bad, kbl == 0 && kbg == 0 && near_tol(bl.pvalue + bg.pvalue, 1.0, 1e-12) && bl.pvalue < 0.01 && bg.pvalue > 0.99, "brunnermunzel, the sign of -W");

    // ---- the faults of the group calls ----
    let none: f64[0] = [];
    let nans: f64[4] = [1.0, 2.0, nan, 4.0];
    let infs: f64[4] = [1.0, 2.0, inf, 4.0];
    let two: i64[2] = [2, 2];
    let one_group: i64[1] = [4];
    let bad_len: i64[2] = [2, 3];
    let zero_sz: i64[2] = [4, 0];
    let neg_sz: i64[2] = [-1, 5];
    let vals: f64[4] = [1.0, 2.0, 3.0, 5.0];
    let same: f64[4] = [3.0, 3.0, 3.0, 3.0];
    for w in 0..7 {
        let nm = if w == 0 { "levene" } else if w == 1 { "fligner" } else if w == 2 { "bartlett" } else if w == 3 { "alexandergovern" } else if w == 4 { "levene mean" } else if w == 5 { "fligner trimmed" } else { "median_test" };
        var f, fc = vs_groups(w, none[..], two[..]);
        fault_is(bad, nm, f, fc, 1);
        f, fc = vs_groups(w, nans[..], two[..]);
        fault_is(bad, nm, f, fc, 3);
        f, fc = vs_groups(w, infs[..], two[..]);
        fault_is(bad, nm, f, fc, 4);
        f, fc = vs_groups(w, vals[..], bad_len[..]);
        fault_is(bad, nm, f, fc, 6);
        f, fc = vs_groups(w, vals[..], one_group[..]);
        fault_is(bad, nm, f, fc, 2);
        f, fc = vs_groups(w, vals[..], zero_sz[..]);
        fault_is(bad, nm, f, fc, 2);
        f, fc = vs_groups(w, vals[..], neg_sz[..]);
        fault_is(bad, nm, f, fc, 5);
        f, fc = vs_groups(w, same[..], two[..]);
        fault_is(bad, nm, f, fc, 7);
    }
    var long: f64[>..] = [];
    for i in 0..statistics::MAX_INPUT + 1 { long.push(1.0); }
    let lsz: i64[2] = [2, 2];
    for w in 0..7 {
        let f, fc = vs_groups(w, long[..], lsz[..]);
        fault_is(bad, "a group call refuses a long buffer", f, fc, 9);
    }
    let sizes31: i64[2] = [3, 1];
    let f1, fc1 = statistics::bartlett(vals[..], sizes31[..]);
    fault_is(bad, "bartlett, a group of one value", f1, fc1, 2);
    let f2, fc2 = statistics::alexandergovern(vals[..], sizes31[..]);
    fault_is(bad, "alexandergovern, a group of one value", f2, fc2, 2);
    let f3, fc3 = statistics::levene(vals[..], sizes31[..]);
    check(bad, fc3 == 0, "levene takes a group of one value");
    let eq2: i64[2] = [1, 1];
    let two_vals: f64[2] = [1.0, 2.0];
    let f4, fc4 = statistics::levene(two_vals[..], eq2[..]);
    fault_is(bad, "levene, as many groups as values", f4, fc4, 2);
    let eqs: i64[2] = [2, 2];
    let f5, fc5 = statistics::bartlett(vals[..], eqs[..]);
    check(bad, fc5 == 0, "bartlett, two groups of two values");

    // ---- the center and the trim ----
    let f6, fc6 = statistics::levene(vals[..], two[..], 3);
    fault_is(bad, "levene, an unknown center", f6, fc6, 5);
    let f7, fc7 = statistics::fligner(vals[..], two[..], -1);
    fault_is(bad, "fligner, a negative center", f7, fc7, 5);
    let f8, fc8 = statistics::levene(vals[..], two[..], statistics::CENTER_TRIMMED, 0.5);
    fault_is(bad, "levene, a trim of 0.5", f8, fc8, 5);
    let f9, fc9 = statistics::fligner(vals[..], two[..], statistics::CENTER_TRIMMED, -0.1);
    fault_is(bad, "fligner, a negative trim", f9, fc9, 5);
    let f10, fc10 = statistics::levene(vals[..], two[..], statistics::CENTER_TRIMMED, nan);
    fault_is(bad, "levene, a NaN trim", f10, fc10, 3);
    let f11, fc11 = statistics::levene(vals[..], two[..], statistics::CENTER_TRIMMED, inf);
    fault_is(bad, "levene, an infinite trim", f11, fc11, 4);
    let wide: f64[6] = [1.0, 2.0, 4.0, 7.0, 3.0, 9.0];
    let three2: i64[2] = [3, 3];
    let f12, fc12 = statistics::levene(wide[..], three2[..], statistics::CENTER_MEAN, nan);
    check(bad, fc12 == 0, "levene, the trim does not matter for another center");
    let f13, fc13 = statistics::fligner(nans[..], two[..], 9);
    fault_is(bad, "fligner, the NaN comes before the center", f13, fc13, 3);
    let f14, fc14 = statistics::levene(wide[..], three2[..], statistics::CENTER_TRIMMED, 0.49);
    check(bad, fc14 == 0, "levene, a trim of 0.49");

    // ---- the faults of the two-sample calls ----
    let nan2: f64[3] = [1.0, nan, 3.0];
    let inf2: f64[3] = [1.0, 2.0, -inf];
    let three: f64[3] = [1.0, 4.0, 9.0];
    let single: f64[1] = [2.0];
    for w in 0..4 {
        let nm = if w == 0 { "ranksums" } else if w == 1 { "mood" } else if w == 2 { "brunnermunzel" } else { "ansari" };
        var f, fc = vs_two(w, none[..], three[..]);
        fault_is(bad, nm, f, fc, 1);
        f, fc = vs_two(w, three[..], none[..]);
        fault_is(bad, nm, f, fc, 1);
        f, fc = vs_two(w, nan2[..], three[..]);
        fault_is(bad, nm, f, fc, 3);
        f, fc = vs_two(w, three[..], nan2[..]);
        fault_is(bad, nm, f, fc, 3);
        f, fc = vs_two(w, inf2[..], three[..]);
        fault_is(bad, nm, f, fc, 4);
        f, fc = vs_two(w, three[..], inf2[..]);
        fault_is(bad, nm, f, fc, 4);
        f, fc = vs_two(w, nan2[..], inf2[..]);
        fault_is(bad, nm, f, fc, 3);
        f, fc = vs_two(w, long[..], nan2[..]);
        fault_is(bad, nm, f, fc, 9);
        f, fc = vs_two(w, nan2[..], long[..]);
        fault_is(bad, nm, f, fc, 9);
    }
    let w1, wc1 = statistics::ranksums(three[..], three[..], 3);
    fault_is(bad, "ranksums, an unknown alternative", w1, wc1, 5);
    let w2, wc2 = statistics::mood(three[..], three[..], -1);
    fault_is(bad, "mood, an unknown alternative", w2, wc2, 5);
    let w3, wc3 = statistics::brunnermunzel(three[..], three[..], 7);
    fault_is(bad, "brunnermunzel, an unknown alternative", w3, wc3, 5);
    let w4, wc4 = statistics::ansari(three[..], three[..], 3);
    fault_is(bad, "ansari, an unknown alternative", w4, wc4, 5);
    let w5, wc5 = statistics::ansari(three[..], three[..], statistics::ALT_TWO_SIDED, 3);
    fault_is(bad, "ansari, an unknown method", w5, wc5, 5);
    let w6, wc6 = statistics::ansari(tx[..], ty[..], statistics::ALT_TWO_SIDED, statistics::MW_EXACT);
    fault_is(bad, "ansari, the exact form refuses ties", w6, wc6, 5);
    var big: f64[>..] = [];
    var big2: f64[>..] = [];
    for i in 0..400 { big.push(i * 1.0); big2.push(i + 0.5); }
    let w7, wc7 = statistics::ansari(big[..], big2[..], statistics::ALT_TWO_SIDED, statistics::MW_EXACT);
    fault_is(bad, "ansari, the exact table is too large", w7, wc7, 8);
    let w8, wc8 = statistics::ansari(big[..], big2[..]);
    check(bad, wc8 == 0 && w8.pvalue >= 0.0 && w8.pvalue <= 1.0, "ansari, the auto form takes the normal form for large samples");
    let w9, wc9 = statistics::mood(single[..], single[..]);
    fault_is(bad, "mood, two values", w9, wc9, 2);
    let w10, wc10 = statistics::brunnermunzel(single[..], three[..]);
    fault_is(bad, "brunnermunzel, a sample of one value", w10, wc10, 2);
    let w11, wc11 = statistics::brunnermunzel(same[..], same[..]);
    fault_is(bad, "brunnermunzel, constant samples", w11, wc11, 7);
    let w12, wc12 = statistics::ranksums(single[..], single[..]);
    check(bad, wc12 == 0 && w12.statistic == 0.0 && w12.pvalue == 1.0, "ranksums, one value each");
    let other: f64[1] = [3.0];
    let w13, wc13 = statistics::ansari(single[..], other[..], statistics::ALT_TWO_SIDED, statistics::MW_ASYMPTOTIC);
    fault_is(bad, "ansari, one value each in the normal form", w13, wc13, 7);
    let w14, wc14 = statistics::ansari(single[..], other[..]);
    check(bad, wc14 == 0 && w14.pvalue == 1.0, "ansari, one value each in the exact form");
    let w15, wc15 = statistics::mood(same[..], same[..]);
    check(bad, wc15 == 7 || (wc15 == 0 && !is_nan(w15.pvalue)), "mood, constant samples");

    // ---- the Friedman layout and faults ----
    let f16, fc16 = statistics::friedmanchisquare(none[..], 3);
    fault_is(bad, "friedmanchisquare, no data", f16, fc16, 1);
    let f17, fc17 = statistics::friedmanchisquare(nans[..], 2);
    fault_is(bad, "friedmanchisquare, NaN comes first", f17, fc17, 3);
    let nine: f64[9] = [1.0, 2.0, 3.0, 2.0, 3.0, 1.0, 3.0, 1.0, 2.0];
    let f18, fc18 = statistics::friedmanchisquare(nine[..], 4);
    fault_is(bad, "friedmanchisquare, a length that the columns do not divide", f18, fc18, 6);
    let f19, fc19 = statistics::friedmanchisquare(nine[..], 0);
    fault_is(bad, "friedmanchisquare, zero columns", f19, fc19, 5);
    let f20, fc20 = statistics::friedmanchisquare(nine[..], -3);
    fault_is(bad, "friedmanchisquare, negative columns", f20, fc20, 5);
    let six: f64[6] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    let f21, fc21 = statistics::friedmanchisquare(six[..], 2);
    fault_is(bad, "friedmanchisquare, two samples", f21, fc21, 2);
    let f22, fc22 = statistics::friedmanchisquare(same[..], 2);
    fault_is(bad, "friedmanchisquare, two columns of equal values", f22, fc22, 2);
    let tied: f64[9] = [2.0, 2.0, 2.0, 5.0, 5.0, 5.0, 1.0, 1.0, 1.0];
    let f23, fc23 = statistics::friedmanchisquare(tied[..], 3);
    fault_is(bad, "friedmanchisquare, every row tied", f23, fc23, 7);
    let f24, fc24 = statistics::friedmanchisquare(infs[..], 2);
    fault_is(bad, "friedmanchisquare, infinity comes first", f24, fc24, 4);
    let f25, fc25 = statistics::friedmanchisquare(long[..], 3);
    fault_is(bad, "friedmanchisquare, a long buffer", f25, fc25, 9);
    // One row per block: a perfect order gives the largest statistic.
    let perfect: f64[9] = [1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0];
    let f26, fc26 = statistics::friedmanchisquare(perfect[..], 3);
    vs_res(bad, "friedmanchisquare, a perfect order", f26, fc26, 6.0, 0.049787068367863944, 2.0, nan, 1e-12);

    // ---- the faults of median_test ----
    out.resize(0);
    out.push(1.0);
    let m1b, mc1 = statistics::median_test(out, nans[..], two[..]);
    check(bad, mc1 == 3 && out.len == 0 && is_nan(m1b.pvalue), "median_test, a NaN leaves out empty");
    let m2b, mc2 = statistics::median_test(out, vals[..], two[..], 3);
    fault_is(bad, "median_test, an unknown ties option", m2b, mc2, 5);
    let m3b, mc3 = statistics::median_test(out, vals[..], two[..], -1);
    fault_is(bad, "median_test, a negative ties option", m3b, mc3, 5);
    let m4b, mc4 = statistics::median_test(out, same[..], two[..]);
    check(bad, mc4 == 7 && out.len == 0, "median_test, every value equals the median");
    let m5b, mc5 = statistics::median_test(out, same[..], two[..], statistics::MEDIAN_IGNORE);
    fault_is(bad, "median_test, every value ignored", m5b, mc5, 7);
    let mixed: f64[6] = [1.0, 2.0, 3.0, 3.0, 3.0, 3.0];
    let sz33: i64[2] = [3, 3];
    let m6b, mc6 = statistics::median_test(out, mixed[..], sz33[..], statistics::MEDIAN_IGNORE);
    fault_is(bad, "median_test, a group of ignored values", m6b, mc6, 7);
    let m7b, mc7 = statistics::median_test(out, mixed[..], sz33[..], statistics::MEDIAN_ABOVE);
    let t_above: f64[4] = [1.0, 3.0, 2.0, 0.0];
    check(bad, mc7 == 0 && near_list(out[..], t_above[..], 0.0), "median_test, the table of the above rule");
    let m8b, mc8 = statistics::median_test(out, mixed[..], sz33[..], statistics::MEDIAN_BELOW);
    check(bad, mc8 == 7 && out.len == 0, "median_test, an empty row under the below rule");

    // ---- the faults of quantile_test ----
    out.resize(0);
    out.push(1.0);
    let q1, qc1 = statistics::quantile_test(out, none[..], 0.0, 0.5);
    check(bad, qc1 == 1 && out.len == 0 && is_nan(q1.pvalue), "quantile_test, empty");
    let q2, qc2 = statistics::quantile_test(out, nans[..], 0.0, 0.5);
    fault_is(bad, "quantile_test, a NaN value", q2, qc2, 3);
    let q3, qc3 = statistics::quantile_test(out, infs[..], 0.0, 0.5);
    fault_is(bad, "quantile_test, an infinite value", q3, qc3, 4);
    let q4, qc4 = statistics::quantile_test(out, vals[..], nan, 0.5);
    fault_is(bad, "quantile_test, a NaN quantile", q4, qc4, 3);
    let q5, qc5 = statistics::quantile_test(out, vals[..], 2.0, nan);
    fault_is(bad, "quantile_test, a NaN probability", q5, qc5, 3);
    let q6, qc6 = statistics::quantile_test(out, vals[..], 2.0, 0.0);
    fault_is(bad, "quantile_test, a probability of 0", q6, qc6, 5);
    let q7, qc7 = statistics::quantile_test(out, vals[..], 2.0, 1.0);
    fault_is(bad, "quantile_test, a probability of 1", q7, qc7, 5);
    let q8, qc8 = statistics::quantile_test(out, vals[..], 2.0, 0.5, 3);
    fault_is(bad, "quantile_test, an unknown alternative", q8, qc8, 5);
    let q9, qc9 = statistics::quantile_test(out, vals[..], 2.0, 0.5, statistics::ALT_TWO_SIDED, 1.0);
    fault_is(bad, "quantile_test, a level of 1", q9, qc9, 5);
    let q10, qc10 = statistics::quantile_test(out, vals[..], 2.0, 0.5, statistics::ALT_TWO_SIDED, 0.0);
    fault_is(bad, "quantile_test, a level of 0", q10, qc10, 5);
    let q11, qc11 = statistics::quantile_test(out, vals[..], 2.0, 0.5, statistics::ALT_TWO_SIDED, nan);
    fault_is(bad, "quantile_test, a NaN level", q11, qc11, 3);
    let q12, qc12 = statistics::quantile_test(out, long[..], 2.0, 0.5);
    fault_is(bad, "quantile_test, a long sample", q12, qc12, 9);
    check(bad, out.len == 0, "quantile_test, a fault leaves out empty");
    let q13, qc13 = statistics::quantile_test(out, vals[..], inf, 0.5, statistics::ALT_GREATER);
    check(bad, qc13 == 0 && q13.statistic == 4.0 && q13.pvalue == 1.0, "quantile_test, an infinite quantile above the sample");
    let q14, qc14 = statistics::quantile_test(out, vals[..], -inf, 0.5, statistics::ALT_LESS);
    check(bad, qc14 == 0 && q14.statistic == 0.0 && q14.pvalue == 1.0, "quantile_test, an infinite quantile below the sample");
    // Five values and a level of 0.99 leave no bound on either side.
    let five: f64[5] = [1.0, 2.0, 3.0, 4.0, 5.0];
    let q15, qc15 = statistics::quantile_test(out, five[..], 3.0, 0.5, statistics::ALT_TWO_SIDED, 0.99);
    check(bad, qc15 == 0 && out.len == 2 && out[0] == -inf && out[1] == inf, "quantile_test, a level that gives no bound");
    check(bad, !is_nan(out[0]) && !is_nan(out[1]), "quantile_test, no bound is NaN");

    // ---- the faults of ttest_ind_from_stats ----
    let t1, tc1 = statistics::ttest_ind_from_stats(nan, 1.0, 5, 0.0, 1.0, 5);
    fault_is(bad, "ttest_ind_from_stats, a NaN mean", t1, tc1, 3);
    let t2, tc2 = statistics::ttest_ind_from_stats(0.0, 1.0, 5, inf, 1.0, 5);
    fault_is(bad, "ttest_ind_from_stats, an infinite mean", t2, tc2, 4);
    let t3, tc3 = statistics::ttest_ind_from_stats(0.0, nan, 5, 0.0, 1.0, 5);
    fault_is(bad, "ttest_ind_from_stats, a NaN deviation", t3, tc3, 3);
    let t4, tc4 = statistics::ttest_ind_from_stats(0.0, 1.0, inf, 0.0, 1.0, 5);
    fault_is(bad, "ttest_ind_from_stats, an infinite count", t4, tc4, 4);
    let t5, tc5 = statistics::ttest_ind_from_stats(0.0, 1.0, 1, 0.0, 1.0, 5);
    fault_is(bad, "ttest_ind_from_stats, a count of 1", t5, tc5, 5);
    let t6, tc6 = statistics::ttest_ind_from_stats(0.0, 1.0, 5, 0.0, 1.0, 2.5);
    fault_is(bad, "ttest_ind_from_stats, a count that is no whole number", t6, tc6, 5);
    let t7, tc7 = statistics::ttest_ind_from_stats(0.0, -1.0, 5, 0.0, 1.0, 5);
    fault_is(bad, "ttest_ind_from_stats, a negative deviation", t7, tc7, 5);
    let t8, tc8 = statistics::ttest_ind_from_stats(0.0, 1.0, 5, 0.0, 1.0, 5, true, 3);
    fault_is(bad, "ttest_ind_from_stats, an unknown alternative", t8, tc8, 5);
    let t9, tc9 = statistics::ttest_ind_from_stats(1.0, 0.0, 5, 1.0, 0.0, 5);
    fault_is(bad, "ttest_ind_from_stats, no spread and equal means", t9, tc9, 7);
    let t10, tc10 = statistics::ttest_ind_from_stats(1.0, 0.0, 5, 2.0, 0.0, 5, false);
    fault_is(bad, "ttest_ind_from_stats, no spread and unequal means, Welch", t10, tc10, 7);
    let t11, tc11 = statistics::ttest_ind_from_stats(1.0, 0.0, 5, 2.0, 0.0, 5);
    fault_is(bad, "ttest_ind_from_stats, no spread and unequal means, pooled", t11, tc11, 7);
    let t12, tc12 = statistics::ttest_ind_from_stats(1e308, 1.0, 5, -1e308, 1.0, 5);
    fault_is(bad, "ttest_ind_from_stats, an overflowing difference", t12, tc12, 8);
    let t13, tc13 = statistics::ttest_ind_from_stats(0.0, 1e200, 5, 1.0, 1e200, 5);
    fault_is(bad, "ttest_ind_from_stats, an overflowing variance", t13, tc13, 8);
    let t14, tc14 = statistics::ttest_ind_from_stats(0.0, 1.0, 5, 0.0, 1.0, 5);
    check(bad, tc14 == 0 && t14.statistic == 0.0 && t14.pvalue == 1.0, "ttest_ind_from_stats, equal means give p = 1");
    let t15, tc15 = statistics::ttest_ind_from_stats(0.0, 1.0, 5, 3.0, 1.0, 5, true, statistics::ALT_LESS);
    check(bad, tc15 == 0 && t15.pvalue < 0.01 && t15.statistic < 0.0, "ttest_ind_from_stats, less");

    // ---- no call writes into the lists of the caller ----
    let ka: f64[8] = [2.1, 3.4, 1.9, 5.6, 4.4, 3.8, 2.7, 4.9];
    let kb: f64[8] = [3.9, 5.2, 4.8, 6.1, 5.5, 4.1, 6.6, 5.0];
    var kg: f64[>..] = [];
    for v in ka { kg.push(v); }
    for v in kb { kg.push(v); }
    var kg0: f64[>..] = [];
    for v in kg { kg0.push(v); }
    let ks: i64[2] = [8, 8];
    let ks0: i64[2] = [8, 8];
    for w in 0..7 {
        let nf, nc = vs_groups(w, kg[..], ks[..]);
    }
    let ka0: f64[8] = [2.1, 3.4, 1.9, 5.6, 4.4, 3.8, 2.7, 4.9];
    let kb0: f64[8] = [3.9, 5.2, 4.8, 6.1, 5.5, 4.1, 6.6, 5.0];
    for w in 0..4 {
        let nf, nc = vs_two(w, ka[..], kb[..]);
    }
    let nf1, nc1 = statistics::brunnermunzel(ka[..], kb[..], statistics::ALT_LESS, false);
    let nf2, nc2 = statistics::ansari(ka[..], kb[..], statistics::ALT_GREATER, statistics::MW_EXACT);
    let nf3, nc3 = statistics::quantile_test(out, ka[..], 3.0, 0.4, statistics::ALT_LESS, 0.9);
    let nf4, nc4 = statistics::median_test(out, kg[..], ks[..], statistics::MEDIAN_IGNORE, false);
    let nf5, nc5 = statistics::friedmanchisquare(kg[..], 4);
    let nf6, nc6 = statistics::levene(kg[..], ks[..], statistics::CENTER_TRIMMED, 0.2);
    check(bad, near_list(kg[..], kg0[..], 0.0) && near_list(ka[..], ka0[..], 0.0) && near_list(kb[..], kb0[..], 0.0)
        && ks[0] == ks0[0] && ks[1] == ks0[1], "the calls leave the lists unchanged");

    // ---- no result is NaN or infinite with the code OK ----
    var clean = true;
    for w in 0..7 {
        let rr, cc = vs_groups(w, kg[..], ks[..]);
        if cc == 0 && (is_nan(rr.statistic) || is_inf(rr.statistic) || is_nan(rr.pvalue) || is_inf(rr.pvalue)) { clean = false; }
    }
    for w in 0..4 {
        let rr, cc = vs_two(w, ka[..], kb[..]);
        if cc == 0 && (is_nan(rr.statistic) || is_inf(rr.statistic) || is_nan(rr.pvalue) || is_inf(rr.pvalue)) { clean = false; }
    }
    check(bad, clean, "no result is NaN or infinite with the code OK");

    // ---- the p-value stays in 0 to 1 for samples far apart ----
    var far: f64[>..] = [];
    for i in 0..20 { far.push(i * 1.0); }
    for i in 0..20 { far.push(1000.0 + i); }
    let fsz: i64[2] = [20, 20];
    let fr, fc = statistics::levene(far[..], fsz[..]);
    check(bad, fc == 0 || fc == 7, "levene, samples far apart");
    let fr2, fc2 = statistics::alexandergovern(far[..], fsz[..]);
    check(bad, fc2 == 0 && fr2.pvalue < 1e-6, "alexandergovern, samples far apart");
    return bad;
}
