// One test for each documented case of the multiple-test corrections.
// The values come from measurements against SciPy and statsmodels. The
// random lists live in multipletests_vectors.goose, which the generator
// test/make_multipletests_vectors.py writes.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_multipletests_tests() -> i64 {
    var bad: i64 = 0;
    let acc: f64[25] = [0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205,
                        0.212, 0.216, 0.222, 0.251, 0.269, 0.275, 0.34, 0.341,
                        0.384, 0.569, 0.594, 0.696, 0.762, 0.94, 0.942, 0.975, 0.986];
    var out: f64[>..] = [];

    let n1, c1 = statistics::multipletests(out, acc[..], statistics::ADJUST_HOCHBERG);
    let h5: f64[5] = [0.025, 0.192, 0.882, 0.882, 0.882];
    check(bad, c1 == 0 && n1 == 25 && near_list(out[0..5], h5[..], 1e-14)
        && near_tol(out[24], 0.986, 1e-14), "hochberg, the acceptance list");

    let n2, c2 = statistics::multipletests(out, acc[..], statistics::ADJUST_HOMMEL);
    let m5: f64[5] = [0.025, 0.192, 0.682, 0.6970000000000001, 0.7140000000000001];
    check(bad, c2 == 0 && n2 == 25 && near_list(out[0..5], m5[..], 1e-14)
        && near_tol(out[24], 0.9860000000000001, 1e-14), "hommel, the acceptance list");

    let n3, c3 = statistics::multipletests(out, acc[..], statistics::ADJUST_BY);
    let y3: f64[3] = [0.09539895444383767, 0.3815958177753507, 0.8013512173282364];
    check(bad, c3 == 0 && near_list(out[0..3], y3[..], 1e-14)
        && near_tol(out[24], 1.0, 1e-14), "benjamini-yekutieli, the acceptance list");

    let n4, c4 = statistics::multipletests(out, acc[..], statistics::ADJUST_NONE);
    check(bad, c4 == 0 && n4 == 25 && near_list(out[..], acc[..], 0.0), "none keeps the values");

    let _n5, c5 = statistics::multipletests(out, acc[..], statistics::ADJUST_BH);
    let b3: f64[3] = [0.025, 0.1, 0.21];
    check(bad, c5 == 0 && near_list(out[0..3], b3[..], 1e-14), "benjamini-hochberg, unchanged");

    // The same list in another order gives the same value for each p-value.
    var shuf: f64[25];
    var idx: i64[25];
    for i in 0..25 { idx[i] = (i * 7 + 3) % 25; shuf[i] = acc[idx[i]]; }
    var sorted_out: f64[>..] = [];
    let ns, cs = statistics::multipletests(sorted_out, acc[..], statistics::ADJUST_HOMMEL);
    let nh, ch = statistics::multipletests(out, shuf[..], statistics::ADJUST_HOMMEL);
    var same = cs == 0 && ch == 0;
    for i in 0..25 { if !near_tol(out[i], sorted_out[idx[i]], 1e-15) { same = false; } }
    check(bad, same, "hommel, the shuffled order");

    // The short lists.
    let two: f64[2] = [0.01, 0.04];
    let _n6, c6 = statistics::multipletests(out, two[..], statistics::ADJUST_HOMMEL);
    check(bad, c6 == 0 && near_tol(out[0], 0.02, 1e-14) && near_tol(out[1], 0.04, 1e-14), "hommel, two values");
    let three: f64[3] = [0.03, 0.02, 0.5];
    let _n7, c7 = statistics::multipletests(out, three[..], statistics::ADJUST_HOMMEL);
    let m3: f64[3] = [0.06, 0.045, 0.5];
    check(bad, c7 == 0 && near_list(out[..], m3[..], 1e-14), "hommel, three values");
    let _n8, c8 = statistics::multipletests(out, three[..], statistics::ADJUST_HOCHBERG);
    let h3: f64[3] = [0.06, 0.06, 0.5];
    check(bad, c8 == 0 && near_list(out[..], h3[..], 1e-14), "hochberg, three values");
    let one: f64[1] = [0.3];
    var every = true;
    for m in 0..7 {
        let n9, c9 = statistics::multipletests(out, one[..], m);
        if c9 != 0 || n9 != 1 || out[0] != 0.3 { every = false; }
    }
    check(bad, every, "every method returns one p-value unchanged");

    // The two-step name of SciPy.
    let _f1, fc1 = statistics::false_discovery_control(out, acc[..], statistics::ADJUST_BY);
    check(bad, fc1 == 0 && near_list(out[0..3], y3[..], 1e-14), "false_discovery_control, by");
    let _f2, fc2 = statistics::false_discovery_control(out, acc[..]);
    check(bad, fc2 == 0 && near_list(out[0..3], b3[..], 1e-14), "false_discovery_control, bh");
    let small: f64[3] = [0.3, 0.01, 0.04];
    let _f3, fc3 = statistics::false_discovery_control(out, small[..], statistics::ADJUST_BY);
    let s3: f64[3] = [0.5499999999999999, 0.05499999999999999, 0.10999999999999999];
    check(bad, fc3 == 0 && near_list(out[..], s3[..], 1e-14), "false_discovery_control, the small list");
    let f4, fc4 = statistics::false_discovery_control(out, acc[..], statistics::ADJUST_HOLM);
    check(bad, fc4 == 5 && f4 == 0 && out.len == 0, "false_discovery_control, another method");

    // The combined p-values.
    let r1, k1 = statistics::combine_pvalues(acc[..], statistics::COMBINE_FISHER);
    check(bad, k1 == 0 && near_tol(r1.statistic, 83.82646878826111, 1e-14)
        && near_tol(r1.pvalue, 0.0019244682392491074, 1e-14) && r1.df == 50.0, "fisher, the acceptance list");
    let r2, k2 = statistics::combine_pvalues(acc[..], statistics::COMBINE_STOUFFER);
    check(bad, k2 == 0 && near_tol(r2.statistic, 2.197137995514678, 1e-14)
        && near_tol(r2.pvalue, 0.014005296122329612, 1e-14), "stouffer, the acceptance list");
    var wt: f64[25];
    for i in 0..25 { wt[i] = i + 1; }
    let r3, k3 = statistics::combine_pvalues(acc[..], statistics::COMBINE_STOUFFER, wt[..]);
    check(bad, k3 == 0 && near_tol(r3.statistic, -1.1395380934276564, 1e-13)
        && near_tol(r3.pvalue, 0.872760605365766, 1e-13), "stouffer, weights 1 to 25");
    let r4, k4 = statistics::combine_pvalues(acc[..], statistics::COMBINE_TIPPETT);
    check(bad, k4 == 0 && near_tol(r4.statistic, 0.001, 1e-14)
        && near_tol(r4.pvalue, 0.02470228740295338, 1e-14), "tippett, the acceptance list");
    let p2: f64[2] = [0.2, 0.4];
    let r5, k5 = statistics::combine_pvalues(p2[..], statistics::COMBINE_FISHER);
    check(bad, k5 == 0 && near_tol(r5.statistic, 5.05145728861651, 1e-14)
        && near_tol(r5.pvalue, 0.2820582915446604, 1e-14), "fisher, two values");
    let w2: f64[2] = [1.0, 3.0];
    let r6, k6 = statistics::combine_pvalues(p2[..], statistics::COMBINE_STOUFFER, w2[..]);
    check(bad, k6 == 0 && near_tol(r6.statistic, 0.5064901678795121, 1e-14)
        && near_tol(r6.pvalue, 0.3062562950943614, 1e-14), "stouffer, weights 1 and 3");
    // The weights are scaled by the largest one. A tiny weight keeps its
    // share of the evidence. SciPy loses the digits of such a weight.
    let tiny: f64[2] = [1e-160, 1e-160];
    let r7, k7 = statistics::combine_pvalues(p2[..], statistics::COMBINE_STOUFFER, tiny[..]);
    check(bad, k7 == 0 && near_tol(r7.statistic, 0.7742595360712864, 1e-12)
        && near_tol(r7.pvalue, 0.21938866363752818, 1e-12), "stouffer, tiny weights");
    let huge: f64[2] = [1e300, 1.0];
    let r8, k8 = statistics::combine_pvalues(p2[..], statistics::COMBINE_STOUFFER, huge[..]);
    check(bad, k8 == 0 && near_tol(r8.statistic, 0.8416212335729142, 1e-14)
        && near_tol(r8.pvalue, 0.19999999999999996, 1e-14), "stouffer, a huge weight");

    // The infinite statistics fault, and Tippett takes the edges.
    let z2: f64[2] = [0.0, 0.5];
    let o2: f64[2] = [1.0, 0.5];
    let e1, ec1 = statistics::combine_pvalues(z2[..], statistics::COMBINE_FISHER);
    check(bad, ec1 == 8 && is_nan(e1.pvalue), "fisher, a p-value of 0");
    let e2, ec2 = statistics::combine_pvalues(z2[..], statistics::COMBINE_STOUFFER);
    check(bad, ec2 == 8, "stouffer, a p-value of 0");
    let e3, ec3 = statistics::combine_pvalues(o2[..], statistics::COMBINE_STOUFFER);
    check(bad, ec3 == 8, "stouffer, a p-value of 1");
    let e4, ec4 = statistics::combine_pvalues(z2[..], statistics::COMBINE_TIPPETT);
    check(bad, ec4 == 0 && e4.statistic == 0.0 && e4.pvalue == 0.0, "tippett, a p-value of 0");
    let ones: f64[2] = [1.0, 1.0];
    let e5, ec5 = statistics::combine_pvalues(ones[..], statistics::COMBINE_TIPPETT);
    check(bad, ec5 == 0 && e5.pvalue == 1.0, "tippett, two p-values of 1");
    let e6, ec6 = statistics::combine_pvalues(o2[..], statistics::COMBINE_FISHER);
    check(bad, ec6 == 0 && near_tol(e6.statistic, 1.3862943611198906, 1e-14), "fisher, a p-value of 1");

    // The weights.
    let neg: f64[2] = [1.0, -1.0];
    let zero: f64[2] = [0.0, 0.0];
    let short: f64[1] = [1.0];
    let nanw: f64[2] = [1.0, statistics::nan()];
    let infw: f64[2] = [1.0, 1.0 / 0.0];
    let _g1, gc1 = statistics::combine_pvalues(p2[..], statistics::COMBINE_STOUFFER, neg[..]);
    check(bad, gc1 == 5, "stouffer, a negative weight");
    let _g2, gc2 = statistics::combine_pvalues(p2[..], statistics::COMBINE_STOUFFER, zero[..]);
    check(bad, gc2 == 5, "stouffer, weights of zero");
    let _g3, gc3 = statistics::combine_pvalues(p2[..], statistics::COMBINE_STOUFFER, short[..]);
    check(bad, gc3 == 6, "stouffer, a short list of weights");
    let _g4, gc4 = statistics::combine_pvalues(p2[..], statistics::COMBINE_STOUFFER, nanw[..]);
    check(bad, gc4 == 5, "stouffer, a NaN weight");
    let _g5, gc5 = statistics::combine_pvalues(p2[..], statistics::COMBINE_STOUFFER, infw[..]);
    check(bad, gc5 == 5, "stouffer, an infinite weight");
    let _g6, gc6 = statistics::combine_pvalues(p2[..], statistics::COMBINE_FISHER, w2[..]);
    check(bad, gc6 == 5, "fisher with weights");
    let _g7, gc7 = statistics::combine_pvalues(p2[..], statistics::COMBINE_TIPPETT, w2[..]);
    check(bad, gc7 == 5, "tippett with weights");
    let _g8, gc8 = statistics::combine_pvalues(p2[..], 3);
    check(bad, gc8 == 5, "combine_pvalues, an unknown method");

    // The hostile inputs, for every call.
    let none: f64[0] = [];
    let outside: f64[2] = [0.5, 1.5];
    let below: f64[2] = [-0.1, 0.5];
    let nan_p: f64[2] = [0.5, statistics::nan()];
    let inf_p: f64[2] = [0.5, 1.0 / 0.0];
    let q1, qc1 = statistics::multipletests(out, none[..], statistics::ADJUST_BY);
    check(bad, qc1 == 1 && q1 == 0 && out.len == 0, "multipletests, an empty list");
    let _q2, qc2 = statistics::multipletests(out, outside[..], statistics::ADJUST_HOMMEL);
    check(bad, qc2 == 5 && out.len == 0, "multipletests, a p-value above 1");
    let _q3, qc3 = statistics::multipletests(out, below[..], statistics::ADJUST_HOCHBERG);
    check(bad, qc3 == 5, "multipletests, a negative p-value");
    let _q4, qc4 = statistics::multipletests(out, nan_p[..], statistics::ADJUST_NONE);
    check(bad, qc4 == 3, "multipletests, a NaN");
    let _q5, qc5 = statistics::multipletests(out, inf_p[..], statistics::ADJUST_BY);
    check(bad, qc5 == 4, "multipletests, an infinity");
    let _q6, qc6 = statistics::multipletests(out, acc[..], 7);
    check(bad, qc6 == 5, "multipletests, the method after the last");
    let _q7, qc7 = statistics::multipletests(out, acc[..], -1);
    check(bad, qc7 == 5, "multipletests, a negative method");
    let _q8, qc8 = statistics::false_discovery_control(out, none[..]);
    check(bad, qc8 == 1, "false_discovery_control, an empty list");
    let _q9, qc9 = statistics::false_discovery_control(out, outside[..]);
    check(bad, qc9 == 5, "false_discovery_control, a p-value above 1");
    let x1, xc1 = statistics::combine_pvalues(none[..]);
    check(bad, xc1 == 1 && is_nan(x1.pvalue), "combine_pvalues, an empty list");
    let x2, xc2 = statistics::combine_pvalues(outside[..], statistics::COMBINE_TIPPETT);
    check(bad, xc2 == 5, "combine_pvalues, a p-value above 1");
    let x3, xc3 = statistics::combine_pvalues(nan_p[..], statistics::COMBINE_FISHER);
    check(bad, xc3 == 3, "combine_pvalues, a NaN");
    let x4, xc4 = statistics::combine_pvalues(inf_p[..], statistics::COMBINE_FISHER);
    check(bad, xc4 == 4, "combine_pvalues, an infinity");

    // The caps. The Hommel form takes 10,000 values and refuses 10,001.
    var big: f64[>..] = [];
    for i in 0..10001 { big.push((i * 37) % 10001 / 10001.0); }
    let hb, hbc = statistics::multipletests(out, big[..], statistics::ADJUST_HOMMEL);
    check(bad, hbc == 9 && hb == 0 && out.len == 0, "hommel, 10001 values");
    let hy, hyc = statistics::multipletests(out, big[..], statistics::ADJUST_BY);
    check(bad, hyc == 0 && hy == 10001, "benjamini-yekutieli, 10001 values");
    big.pop();
    let hk, hkc = statistics::multipletests(out, big[..], statistics::ADJUST_HOMMEL);
    var in_range = hkc == 0 && hk == 10000;
    for v in out { if is_nan(v) || v < 0.0 || v > 1.0 { in_range = false; } }
    check(bad, in_range, "hommel, 10000 values");

    var long: f64[>..] = [];
    for i in 0..statistics::MAX_INPUT + 1 { long.push(0.5); }
    let _lg1, lc1 = statistics::multipletests(out, long[..], statistics::ADJUST_NONE);
    let _lg2, lc2 = statistics::false_discovery_control(out, long[..]);
    let _lg3, lc3 = statistics::combine_pvalues(long[..], statistics::COMBINE_TIPPETT);
    check(bad, lc1 == 9 && lc2 == 9 && lc3 == 9, "the calls refuse a list above the input bound");

    // No call writes into the list or the weights of the caller.
    let keep: f64[6] = [0.3, 0.01, 0.5, 0.04, 0.2, 0.9];
    let before: f64[6] = [0.3, 0.01, 0.5, 0.04, 0.2, 0.9];
    let kw: f64[6] = [1.0, 2.0, 3.0, 1.0, 2.0, 3.0];
    let kw0: f64[6] = [1.0, 2.0, 3.0, 1.0, 2.0, 3.0];
    for m in 0..7 {
        let _a, _ac = statistics::multipletests(out, keep[..], m);
    }
    let _b, _bc = statistics::false_discovery_control(out, keep[..], statistics::ADJUST_BY);
    let _c1, _cc1 = statistics::combine_pvalues(keep[..], statistics::COMBINE_FISHER);
    let _c2, _cc2 = statistics::combine_pvalues(keep[..], statistics::COMBINE_STOUFFER, kw[..]);
    let _c3, _cc3 = statistics::combine_pvalues(keep[..], statistics::COMBINE_TIPPETT);
    check(bad, near_list(keep[..], before[..], 0.0) && near_list(kw[..], kw0[..], 0.0),
        "the calls leave the list and the weights unchanged");

    // The old methods keep their values.
    let old: f64[4] = [0.01, 0.04, 0.03, 0.005];
    let _o1, oc1 = statistics::multipletests(out, old[..], statistics::ADJUST_HOLM);
    check(bad, oc1 == 0 && near(out[3], 0.02) && near(out[0], 0.03) && near(out[2], 0.06) && near(out[1], 0.06), "holm, unchanged");
    return bad;
}
