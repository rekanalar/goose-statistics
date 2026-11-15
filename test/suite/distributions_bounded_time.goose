
// Time rows of the bounded families: a long run of calls ends in a bounded
// time and the answers stay finite.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_bounded_time() -> i64 {
    var bad: i64 = 0;
    var total: f64 = 0.0;
    let n = 200;
    for k in 0..n {
        let p = (k as f64 + 0.5) / n;
        let q, cq = statistics::wald_ppf(p);
        let x, cx = statistics::semicircular_ppf(p);
        let d, cd = statistics::moyal_pdf(x);
        if cq != 0 || cx != 0 || cd != 0 { bad += 1; }
        total += q + d;
    }
    check(bad, !is_nan(total) && !is_inf(total), "the sum of the calls is finite");
    return bad;
}
