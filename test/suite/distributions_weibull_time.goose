// Time rows of the Weibull families: a long run of calls ends in a bounded
// time and the answers stay finite.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_weibull_time() -> i64 {
    var bad: i64 = 0;
    var total: f64 = 0.0;
    let n = 100;
    for k in 0..n {
        let p = (k as f64 + 0.5) / n;
        let q, cq = statistics::invweibull_ppf(p, 1.5);
        let d, cd = statistics::truncweibull_min_cdf(p, 1.5, 1.0, 2.0);
        if cq != 0 || cd != 0 { bad += 1; }
        total += q + d;
    }
    check(bad, !is_nan(total) && !is_inf(total), "the sum of the calls is finite");
    return bad;
}
