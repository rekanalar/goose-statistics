// Time rows of the families of a special shape: a long run of calls ends
// in a bounded time and the answers stay finite. The quantile of the
// inverse Gaussian law costs a search, so it runs fewer times.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_special_time() -> i64 {
    var bad: i64 = 0;
    var total: f64 = 0.0;
    let n = 50;
    for k in 0..n {
        let x = (k as f64 + 0.5) / n * 6.0 - 1.0;
        let d, cd = statistics::nakagami_cdf(x, 2.5);
        let e, ce = statistics::kappa3_pdf(x, 2.0);
        let f, cf = statistics::pearson3_cdf(x, 1.5);
        if cd != 0 || ce != 0 || cf != 0 { bad += 1; }
        total += d + e + f;
    }
    for k in 0..n / 20 {
        let p = (k as f64 + 0.5) / (n / 20);
        let x, cx = statistics::invgauss_ppf(p, 0.5);
        let d, cd = statistics::invgauss_cdf(x, 0.5);
        if cx != 0 || cd != 0 { bad += 1; }
        total += x + d;
    }
    let m, cm = statistics::nakagami_mean(2.5);
    chk_code(bad, "nakagami mean code", cm, 0);
    let k, ck = statistics::kappa3_mean(3.5);
    chk_code(bad, "kappa3 mean code", ck, 0);
    total += m + k;
    check(bad, !is_nan(total) && !is_inf(total), "the sum of the calls is finite");
    return bad;
}
