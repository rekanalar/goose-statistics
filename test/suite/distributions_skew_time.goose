// Time rows of the families of a skew shape: a long run of calls ends in a
// bounded time and the answers stay finite. The quantile of the skew normal
// law costs a search, so it runs fewer times.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_skew_time() -> i64 {
    var bad: i64 = 0;
    var total: f64 = 0.0;
    let n = 2000;
    for k in 0..n {
        let x = (k as f64 + 0.5) / n as f64 * 6.0 - 1.0;
        let s = ((k as f64 + 0.5) / n as f64 - 0.5) * 4.0;
        let d, cd = statistics::skewnorm_cdf(x, s);
        let e, ce = statistics::skewnorm_pdf(x, s);
        let f, cf = statistics::laplace_asymmetric_cdf(x, 0.5);
        if cd != 0 || ce != 0 || cf != 0 { bad += 1; }
        total += d + e + f;
    }
    for k in 0..n / 20 {
        let p = (k as f64 + 0.5) / (n / 20) as f64;
        let x, cx = statistics::skewnorm_ppf(p, 2.0);
        let d, cd = statistics::skewnorm_cdf(x, 2.0);
        if cx != 0 || cd != 0 { bad += 1; }
        total += x + d;
    }
    let m, cm = statistics::skewnorm_mean(2.0);
    chk_code(bad, "skewnorm mean code", cm, 0);
    let k, ck = statistics::laplace_asymmetric_kurt(0.5);
    chk_code(bad, "laplace_asymmetric kurt code", ck, 0);
    total += m + k;
    check(bad, !is_nan(total) && !is_inf(total), "the sum of the calls is finite");
    return bad;
}
