// The fits of the chi-square, the student and the F families, and the
// Filliben statistic of every family.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_fitting_shapes_tests() -> i64 {
    var bad: i64 = 0;
    var out: f64[>..] = [];
    var r3 = rng { s: 3 };
    // a chi-square sample from the quantiles is fitted near its parameters
    var cs: f64[>..] = [];
    for i in 0..150 {
        let v, cv = statistics::chi2_ppf((i + 0.5) / 150.0, 6.0);
        cs.push(v * 0.7);
    }
    let t, c = statistics::goodness_of_fit(out, cs[..], 7, 3, 2, r3);
    chk_code(bad, "the chi-square fit", c, 0);
    check(bad, abs(out[0] - 6.0) < 0.5 && abs(out[2] - 0.7) < 0.05 && out[1] == 0.0, "the chi-square fit recovers its parameters");
    // the Filliben statistic of a sample of quantiles is near one
    check(bad, t.statistic > 0.98, "the Filliben statistic of quantiles is near one");
    // a student sample from the quantiles is fitted near its parameters
    var ts: f64[>..] = [];
    for i in 0..150 {
        let v, cv = statistics::tdist_ppf((i + 0.5) / 150.0, 6.0);
        ts.push(v * 1.2 + 0.5);
    }
    let t2, c2 = statistics::goodness_of_fit(out, ts[..], 8, 3, 2, r3);
    chk_code(bad, "the student fit", c2, 0);
    check(bad, abs(out[0] - 6.0) < 1.0 && abs(out[1] - 0.5) < 0.1 && abs(out[2] - 1.2) < 0.1, "the student fit recovers its parameters");
    check(bad, t2.statistic > 0.98, "the Filliben statistic of student quantiles is near one");
    // an F sample from the quantiles is fitted near its parameters
    var fs: f64[>..] = [];
    for i in 0..150 {
        let v, cv = statistics::fdist_ppf((i + 0.5) / 150.0, 5.0, 9.0);
        fs.push(v * 1.1);
    }
    let t3, c3 = statistics::goodness_of_fit(out, fs[..], 9, 3, 2, r3);
    chk_code(bad, "the F fit", c3, 0);
    check(bad, abs(out[0] - 5.0) < 1.0 && abs(out[1] - 9.0) < 2.0 && abs(out[3] - 1.1) < 0.1 && out[2] == 0.0, "the F fit recovers its parameters");
    check(bad, t3.statistic > 0.98, "the Filliben statistic of F quantiles is near one");
    // a sample from another family gives a small Filliben statistic
    let t4, c4 = statistics::goodness_of_fit(out, cs[..], 8, 3, 2, r3);
    chk_code(bad, "chi-square data against a student fit", c4, 0);
    check(bad, t4.statistic < t.statistic, "the wrong family lowers the Filliben statistic");
    // the Filliben statistic is a correlation, so it stays in the range
    let t5, c5 = statistics::goodness_of_fit(out, cs[..], 7, 3, 2, r3);
    chk_code(bad, "the Filliben statistic", c5, 0);
    check(bad, t5.statistic >= -1.0 && t5.statistic <= 1.0, "the Filliben statistic stays in the range");
    // the same stream gives the same p-value
    var g1 = rng { s: 5 };
    var g2 = rng { s: 5 };
    let a, ca = statistics::goodness_of_fit(out, cs[..], 7, 3, 4, g1);
    let b, cb = statistics::goodness_of_fit(out, cs[..], 7, 3, 4, g2);
    check(bad, ca == 0 && cb == 0 && a.pvalue == b.pvalue && a.statistic == b.statistic, "the stream fixes the p-value");
    // a sample of one value has no fit
    let one: f64[1] = [1.5];
    var oout: f64[>..] = [];
    let or, oc = statistics::goodness_of_fit(oout, one[..], 7, 3, 2, r3);
    chk_code(bad, "a sample of one value", oc, 2);
    // a family outside the list gives the domain fault
    let dr, dc = statistics::goodness_of_fit(oout, cs[..], 10, 3, 2, r3);
    chk_code(bad, "a family outside the list", dc, 5);
    // a statistic outside the list gives the domain fault
    let sr, sc = statistics::goodness_of_fit(oout, cs[..], 7, 4, 2, r3);
    chk_code(bad, "a statistic outside the list", sc, 5);
    // a negative value under the chi-square family gives the domain fault
    let neg: f64[3] = [-1.0, 2.0, 3.0];
    let nr, nc = statistics::goodness_of_fit(oout, neg[..], 7, 3, 2, r3);
    chk_code(bad, "a negative value", nc, 5);
    return bad;
}
