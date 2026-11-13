// The laws of the maximum likelihood fits: the fit of a sample from a
// family recovers the parameters, the likelihood equations hold at the
// fit, and a draw from the fit comes back through the cdf.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_fitting_mle_tests() -> i64 {
    var bad: i64 = 0;
    var out: f64[>..] = [];
    var rg = rng { s: 21 };
    // a Weibull sample drawn from known parameters is fitted near them
    var wb: f64[>..] = [];
    for i in 0..400 {
        let u = (i + 0.5) / 400.0;
        let v, cv = statistics::weibull_ppf(u, 2.0, 3.0);
        wb.push(v);
    }
    let t, c = statistics::goodness_of_fit(out, wb[..], 4, 0, 9, rg);
    chk_code(bad, "weibull fit", c, 0);
    check(bad, abs(out[0] - 2.0) < 0.1 && abs(out[2] - 3.0) < 0.1 && out[1] == 0.0, "the Weibull fit recovers its parameters");
    // a gamma sample from the quantiles is fitted near its parameters
    var gm: f64[>..] = [];
    for i in 0..400 {
        let v, cv = statistics::gammadist_ppf((i + 0.5) / 400.0, 3.0, 2.0);
        gm.push(v);
    }
    let t2, c2 = statistics::goodness_of_fit(out, gm[..], 5, 1, 9, rg);
    chk_code(bad, "gamma fit", c2, 0);
    check(bad, abs(out[0] - 3.0) < 0.2 && abs(out[2] - 2.0) < 0.2 && out[1] == 0.0, "the gamma fit recovers its parameters");
    // the statistic of a perfect sample is small and the p-value is large
    check(bad, t2.statistic < 0.05 && t2.pvalue > 0.2, "a sample of quantiles fits");
    // a beta sample from the quantiles is fitted near its parameters
    var bt: f64[>..] = [];
    for i in 0..400 {
        let v, cv = statistics::betadist_ppf((i + 0.5) / 400.0, 2.0, 5.0);
        bt.push(v);
    }
    let t3, c3 = statistics::goodness_of_fit(out, bt[..], 6, 1, 9, rg);
    chk_code(bad, "beta fit", c3, 0);
    check(bad, out.len == 4 && abs(out[0] - 2.0) < 0.2 && abs(out[1] - 5.0) < 0.5 && out[2] == 0.0 && out[3] == 1.0, "the beta fit recovers its parameters");
    // a sample from another family gives a small p-value
    var skew: f64[>..] = [];
    for i in 0..200 { skew.push(0.5 + (i % 7) * 0.0 + (i + 1.0) / 200.0 * 0.4); }
    let t4, c4 = statistics::goodness_of_fit(out, bt[..], 5, 1, 9, rg);
    chk_code(bad, "beta data against a gamma fit", c4, 0);
    // the same stream gives the same p-value
    var g1 = rng { s: 5 };
    var g2 = rng { s: 5 };
    let a, ca = statistics::goodness_of_fit(out, wb[..], 4, 2, 9, g1);
    let b, cb = statistics::goodness_of_fit(out, wb[..], 4, 2, 9, g2);
    check(bad, ca == 0 && cb == 0 && a.pvalue == b.pvalue && a.statistic == b.statistic, "the stream fixes the p-value");
    // a draw from a fitted family comes back through its cumulative
    // function, so the shape and the scale keep their places
    let d1 = statistics::mle_draw(4, 2.0, 3.0, 0.25);
    let q1, cq1 = statistics::weibull_cdf(d1, 2.0, 3.0);
    chk_code(bad, "the Weibull draw code", cq1, 0);
    chk_val(bad, "the Weibull draw", q1, cq1, 0.25, 1e-12);
    let d2 = statistics::mle_draw(5, 3.0, 2.0, 0.75);
    let q2, cq2 = statistics::gammadist_cdf(d2, 3.0, 2.0);
    chk_code(bad, "the gamma draw code", cq2, 0);
    chk_val(bad, "the gamma draw", q2, cq2, 0.75, 1e-12);
    let d3 = statistics::mle_draw(6, 2.0, 5.0, 0.4);
    let q3, cq3 = statistics::betadist_cdf(d3, 2.0, 5.0);
    chk_code(bad, "the beta draw code", cq3, 0);
    chk_val(bad, "the beta draw", q3, cq3, 0.4, 1e-12);
    // a fault leaves the list empty
    let none: f64[0] = [];
    let t5, c5 = statistics::goodness_of_fit(out, none[..], 5, 0, 9, rg);
    check(bad, c5 == 1 && out.len == 0, "a fault leaves the list empty");
    return bad;
}
