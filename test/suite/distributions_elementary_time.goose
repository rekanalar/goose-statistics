// Time rows of the elementary families: a long run of calls ends in a
// bounded time and the answers stay finite.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_elementary_time() -> i64 {
    var bad: i64 = 0;
    var total: f64 = 0.0;
    let n = 20000;
    for k in 0..n {
        let p = (k as f64 + 0.5) / n as f64;
        let q, cq = statistics::cosine_ppf(p);
        let x, cx = statistics::hypsecant_ppf(p);
        let d, cd = statistics::halfnorm_pdf(x);
        if cq != 0 || cx != 0 || cd != 0 { bad += 1; }
        total += q + d;
    }
    check(bad, !is_nan(total) && !is_inf(total), "the sum of the calls is finite");
    return bad;
}
