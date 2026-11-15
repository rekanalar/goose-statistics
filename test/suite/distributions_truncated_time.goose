// Time rows of the families of a truncation: a long run of calls ends in a
// bounded time and the answers stay finite.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_truncated_time() -> i64 {
    var bad: i64 = 0;
    var total: f64 = 0.0;
    let n = 50;
    for k in 0..n {
        let p = (k as f64 + 0.5) / n;
        let x, cx = statistics::truncpareto_ppf(p, 2.0, 5.0);
        let d, cd = statistics::truncpareto_cdf(x, 2.0, 5.0);
        let y, cy = statistics::truncnorm_ppf(p, -1.0, 2.0);
        if cx != 0 || cd != 0 || cy != 0 { bad += 1; }
        total += x + d + y;
    }
    let m, cm = statistics::truncnorm_mean(-1.0, 2.0);
    chk_code(bad, "truncnorm mean code", cm, 0);
    let h, ch = statistics::truncpareto_kurt(2.0, 5.0);
    chk_code(bad, "truncpareto kurt code", ch, 0);
    total += m + h;
    check(bad, !is_nan(total) && !is_inf(total), "the sum of the calls is finite");
    return bad;
}
