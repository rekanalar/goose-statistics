// The tests of the multicomp calls: the faults, the edge cases and the
// rules of the calls. The measured rows live in multicomp_vectors.goose.
namespace suite;

import std;
import math;
import statistics from "../..";

fn mc_codes(bad: i64&, got: i64, want: i64, name: u8[:]) {
    if got == want { return; }
    bad += 1;
    print("FAIL ", name, " code ", got, " want ", want, "\n");
}

fn run_multicomp_tests() -> i64 {
    var bad: i64 = 0;
    let g: f64[12] = [2.1, 3.4, 1.9, 5.6, 4.4, 3.8, 3.9, 5.2, 4.8, 6.1, 5.5, 4.1];
    let s2: i64[2] = [6, 6];
    let s3: i64[3] = [4, 4, 4];
    var out: f64[>..] = [];

    // The faults of tukey_hsd, in the order of the checks.
    var n: i64 = 0;
    var c: i64 = 0;
    let w0: f64[0] = [];
    let e0: i64[0] = [];
    n, c = statistics::tukey_hsd(out, w0[..], e0[..], 0.95, true);
    mc_codes(bad, c, statistics::ERR_EMPTY, "tukey empty");
    var nan_g = g;
    nan_g[3] = 0.0 / 0.0;
    n, c = statistics::tukey_hsd(out, nan_g[..], s3[..], 0.95, true);
    mc_codes(bad, c, statistics::ERR_NAN, "tukey nan");
    var inf_g = g;
    inf_g[3] = 1.0 / 0.0;
    n, c = statistics::tukey_hsd(out, inf_g[..], s3[..], 0.95, true);
    mc_codes(bad, c, statistics::ERR_INFINITY, "tukey infinity");
    let one: i64[1] = [12];
    n, c = statistics::tukey_hsd(out, g[..], one[..], 0.95, true);
    mc_codes(bad, c, statistics::ERR_TOO_SHORT, "tukey one group");
    let bad_sum: i64[3] = [4, 4, 5];
    n, c = statistics::tukey_hsd(out, g[..], bad_sum[..], 0.95, true);
    mc_codes(bad, c, statistics::ERR_LENGTH, "tukey sizes do not sum");
    let neg: i64[3] = [-1, 8, 5];
    n, c = statistics::tukey_hsd(out, g[..], neg[..], 0.95, true);
    mc_codes(bad, c, statistics::ERR_DOMAIN, "tukey negative size");
    let single: i64[3] = [1, 5, 6];
    n, c = statistics::tukey_hsd(out, g[..], single[..], 0.95, true);
    mc_codes(bad, c, statistics::ERR_TOO_SHORT, "tukey group of one value");
    n, c = statistics::tukey_hsd(out, g[..], s3[..], 0.0, true);
    mc_codes(bad, c, statistics::ERR_DOMAIN, "tukey level 0");
    n, c = statistics::tukey_hsd(out, g[..], s3[..], 1.0, true);
    mc_codes(bad, c, statistics::ERR_DOMAIN, "tukey level 1");
    n, c = statistics::tukey_hsd(out, g[..], s3[..], 0.0 / 0.0, true);
    mc_codes(bad, c, statistics::ERR_NAN, "tukey level nan");
    check(bad, out.len == 0, "tukey fault leaves out empty");
    let flat: f64[8] = [3.0, 3.0, 3.0, 3.0, 5.0, 5.0, 5.0, 5.0];
    let s22: i64[2] = [4, 4];
    n, c = statistics::tukey_hsd(out, flat[..], s22[..], 0.95, true);
    mc_codes(bad, c, statistics::ERR_UNDEFINED, "tukey constant groups");
    n, c = statistics::tukey_hsd(out, flat[..], s22[..], 0.95, false);
    mc_codes(bad, c, statistics::ERR_UNDEFINED, "tukey constant groups welch");

    // The count of groups and the bound on it.
    var big: f64[>..] = [];
    var bs: i64[>..] = [];
    for i in 0..65 {
        bs.push(2);
        big.push(i * 0.5);
        big.push(i * 0.5 + 1.0 + (i % 3) * 0.1);
    }
    n, c = statistics::tukey_hsd(out, big[..], bs[..], 0.95, true);
    mc_codes(bad, c, statistics::ERR_TOO_LONG, "tukey more than 64 groups");
    n, c = statistics::dunnett(out, big[..], bs[..], 0, statistics::ALT_TWO_SIDED, 0.95);
    mc_codes(bad, c, statistics::ERR_TOO_LONG, "dunnett more than 64 groups");

    // The table of two groups holds one row; three groups hold three.
    n, c = statistics::tukey_hsd(out, g[..], s2[..], 0.95, true);
    check(bad, c == 0 && n == 1 && out.len == 4, "tukey two groups one row");
    check(bad, out[2] < out[0] && out[0] < out[3], "tukey interval holds the statistic");
    n, c = statistics::tukey_hsd(out, g[..], s3[..], 0.95, true);
    check(bad, c == 0 && n == 3 && out.len == 12, "tukey three groups three rows");
    var all_ok = true;
    for i in 0..n {
        let p = out[4 * i + 1];
        if p < 0.0 || p > 1.0 { all_ok = false; }
    }
    check(bad, all_ok, "tukey p-values lie in 0 to 1");

    // A wider level gives a wider interval, and the p-value does not move.
    var wide: f64[>..] = [];
    n, c = statistics::tukey_hsd(wide, g[..], s3[..], 0.99, true);
    check(bad, c == 0 && wide[3] - wide[2] > out[3] - out[2], "tukey wider level wider interval");
    check(bad, wide[1] == out[1], "tukey level leaves the p-value");

    // The sign flips when two groups swap, and the p-value stays.
    let sw: f64[12] = [3.9, 5.2, 4.8, 6.1, 5.5, 4.1, 2.1, 3.4, 1.9, 5.6, 4.4, 3.8];
    var rev: f64[>..] = [];
    var fwd: f64[>..] = [];
    n, c = statistics::tukey_hsd(fwd, g[..], s2[..], 0.95, true);
    n, c = statistics::tukey_hsd(rev, sw[..], s2[..], 0.95, true);
    check(bad, fwd[0] == -rev[0] && near(fwd[1], rev[1]), "tukey swap flips the sign");

    // A tiny p-value keeps its digits: two far groups.
    let far: f64[12] = [0.1, 0.2, 0.15, 0.12, 0.18, 0.11, 9.1, 9.2, 9.15, 9.12, 9.18, 9.11];
    n, c = statistics::tukey_hsd(out, far[..], s2[..], 0.95, true);
    check(bad, c == 0 && out[1] > 0.0 && out[1] < 1e-12, "tukey tiny p-value keeps digits");

    // No call writes into the lists of the caller.
    let keep = g;
    let ks = s3;
    n, c = statistics::tukey_hsd(out, g[..], s3[..], 0.95, false);
    n, c = statistics::dunnett(out, g[..], s3[..], 0, statistics::ALT_TWO_SIDED, 0.95);
    var same = true;
    for i in 0..g.len { if g[i] != keep[i] { same = false; } }
    for i in 0..s3.len { if s3[i] != ks[i] { same = false; } }
    check(bad, same, "no write into the lists of the caller");

    // The faults of dunnett.
    n, c = statistics::dunnett(out, w0[..], e0[..], 0, 0, 0.95);
    mc_codes(bad, c, statistics::ERR_EMPTY, "dunnett empty");
    n, c = statistics::dunnett(out, nan_g[..], s3[..], 0, 0, 0.95);
    mc_codes(bad, c, statistics::ERR_NAN, "dunnett nan");
    n, c = statistics::dunnett(out, inf_g[..], s3[..], 0, 0, 0.95);
    mc_codes(bad, c, statistics::ERR_INFINITY, "dunnett infinity");
    n, c = statistics::dunnett(out, g[..], one[..], 0, 0, 0.95);
    mc_codes(bad, c, statistics::ERR_TOO_SHORT, "dunnett no group beside the control");
    n, c = statistics::dunnett(out, g[..], single[..], 0, 0, 0.95);
    mc_codes(bad, c, statistics::ERR_TOO_SHORT, "dunnett group of one value");
    n, c = statistics::dunnett(out, g[..], s3[..], 3, 0, 0.95);
    mc_codes(bad, c, statistics::ERR_DOMAIN, "dunnett control past the end");
    n, c = statistics::dunnett(out, g[..], s3[..], -1, 0, 0.95);
    mc_codes(bad, c, statistics::ERR_DOMAIN, "dunnett negative control");
    n, c = statistics::dunnett(out, g[..], s3[..], 0, 3, 0.95);
    mc_codes(bad, c, statistics::ERR_DOMAIN, "dunnett bad alternative");
    n, c = statistics::dunnett(out, g[..], s3[..], 0, 0, 1.0);
    mc_codes(bad, c, statistics::ERR_DOMAIN, "dunnett level 1");
    n, c = statistics::dunnett(out, g[..], s3[..], 0, 0, 0.0);
    mc_codes(bad, c, statistics::ERR_DOMAIN, "dunnett level 0");
    n, c = statistics::dunnett(out, flat[..], s22[..], 0, 0, 0.95);
    mc_codes(bad, c, statistics::ERR_UNDEFINED, "dunnett constant groups");
    check(bad, out.len == 0, "dunnett fault leaves out empty");

    // Dunnett is deterministic, and the alternatives fit together.
    var d1: f64[>..] = [];
    var d2: f64[>..] = [];
    n, c = statistics::dunnett(d1, g[..], s3[..], 0, statistics::ALT_TWO_SIDED, 0.95);
    n, c = statistics::dunnett(d2, g[..], s3[..], 0, statistics::ALT_TWO_SIDED, 0.95);
    var rep = c == 0 && d1.len == 8 && d2.len == 8;
    if rep { for i in 0..8 { if d1[i] != d2[i] { rep = false; } } }
    check(bad, rep, "dunnett repeats exactly");
    var gt: f64[>..] = [];
    var lt: f64[>..] = [];
    n, c = statistics::dunnett(gt, g[..], s3[..], 0, statistics::ALT_GREATER, 0.95);
    n, c = statistics::dunnett(lt, g[..], s3[..], 0, statistics::ALT_LESS, 0.95);
    check(bad, gt[3] > 1.0e300 && lt[2] < -1.0e300, "dunnett one-sided intervals reach infinity");
    // The p-values of the two one-sided forms of one group add to about 1
    // when only one group sits beside the control; the single-comparison
    // case has an exact t form.
    let two: i64[2] = [6, 6];
    var o1: f64[>..] = [];
    var o2: f64[>..] = [];
    var o3: f64[>..] = [];
    n, c = statistics::dunnett(o1, g[..], two[..], 0, statistics::ALT_GREATER, 0.95);
    n, c = statistics::dunnett(o2, g[..], two[..], 0, statistics::ALT_LESS, 0.95);
    n, c = statistics::dunnett(o3, g[..], two[..], 0, statistics::ALT_TWO_SIDED, 0.95);
    check(bad, near_tol(o1[1] + o2[1], 1.0, 1e-9), "dunnett one group: the one-sided p-values add to 1");
    let xs1 = [2.1, 3.4, 1.9, 5.6, 4.4, 3.8];
    let xs2 = [3.9, 5.2, 4.8, 6.1, 5.5, 4.1];
    let tr, ct = statistics::ttest_ind(xs2[..], xs1[..], statistics::ALT_TWO_SIDED, true);
    check(bad, ct == 0 && near_tol(o3[1], tr.pvalue, 1e-9) && near_tol(o3[0], tr.statistic, 1e-12), "dunnett one group equals the t test");

    // bayes_mvs and mvsdist: the faults.
    let r, cb = statistics::bayes_mvs(w0[..], 0.9);
    mc_codes(bad, cb, statistics::ERR_EMPTY, "bayes_mvs empty");
    let one_v: f64[1] = [1.0];
    let r1, c1 = statistics::bayes_mvs(one_v[..], 0.9);
    mc_codes(bad, c1, statistics::ERR_TOO_SHORT, "bayes_mvs one value");
    let r2, c2 = statistics::bayes_mvs(nan_g[..], 0.9);
    mc_codes(bad, c2, statistics::ERR_NAN, "bayes_mvs nan");
    let r3, c3 = statistics::bayes_mvs(inf_g[..], 0.9);
    mc_codes(bad, c3, statistics::ERR_INFINITY, "bayes_mvs infinity");
    let r4, c4 = statistics::bayes_mvs(g[..], 0.0);
    mc_codes(bad, c4, statistics::ERR_DOMAIN, "bayes_mvs level 0");
    let r5, c5 = statistics::bayes_mvs(g[..], 1.0);
    mc_codes(bad, c5, statistics::ERR_DOMAIN, "bayes_mvs level 1");
    let r6, c6 = statistics::bayes_mvs(g[..], 0.0 / 0.0);
    mc_codes(bad, c6, statistics::ERR_NAN, "bayes_mvs level nan");
    let same_v: f64[5] = [2.0, 2.0, 2.0, 2.0, 2.0];
    let r7, c7 = statistics::bayes_mvs(same_v[..], 0.9);
    mc_codes(bad, c7, statistics::ERR_UNDEFINED, "bayes_mvs constant sample");
    let three: f64[3] = [1.0, 2.0, 4.0];
    let r8, c8 = statistics::bayes_mvs(three[..], 0.9);
    mc_codes(bad, c8, statistics::ERR_UNDEFINED, "bayes_mvs three values");
    let d0, cd = statistics::mvsdist(one_v[..]);
    mc_codes(bad, cd, statistics::ERR_TOO_SHORT, "mvsdist one value");
    let d9, c9 = statistics::mvsdist(same_v[..]);
    mc_codes(bad, c9, statistics::ERR_UNDEFINED, "mvsdist constant sample");
    check(bad, is_nan(r7.mean) && is_nan(d9.mean), "a fault gives NaN");

    // The intervals nest with the level, and the mean lies inside.
    let ra, ca = statistics::bayes_mvs(g[..], 0.8);
    let rb, cb2 = statistics::bayes_mvs(g[..], 0.95);
    check(bad, ca == 0 && cb2 == 0 && rb.mean_interval.low < ra.mean_interval.low && rb.variance_interval.high > ra.variance_interval.high && rb.deviation_interval.low < ra.deviation_interval.low, "bayes_mvs intervals nest");
    check(bad, ra.mean_interval.low < ra.mean && ra.mean < ra.mean_interval.high, "bayes_mvs mean inside");

    // The default level is 0.9, as in SciPy.
    let rd, cdd = statistics::bayes_mvs(g[..]);
    let re, ce = statistics::bayes_mvs(g[..], 0.9);
    check(bad, rd.mean_interval.low == re.mean_interval.low, "bayes_mvs default level");

    // More than 1000 values: the Gaussian form, and the bound on input.
    var long: f64[>..] = [];
    for i in 0..1500 { long.push((i % 17) * 0.5 + (i % 5) * 0.1); }
    let dl, cl = statistics::mvsdist(long[..]);
    check(bad, cl == 0 && dl.gaussian && dl.mean_df == 0.0, "mvsdist above 1000 values is Gaussian");
    let rl, cl2 = statistics::bayes_mvs(long[..], 0.9);
    check(bad, cl2 == 0 && rl.variance_interval.low < rl.variance && rl.variance < rl.variance_interval.high, "bayes_mvs Gaussian interval");
    return bad;
}
