// Time rows of the first group of discrete families: a long run of calls
// ends in a bounded time and the answers stay finite. The laws whose
// quantile needs a search cost more than the other calls.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_discrete_a_time() -> i64 {
    var bad: i64 = 0;
    var total: f64 = 0.0;
    let n = 2000;
    for k in 0..n {
        let p = (k as f64 + 0.5) / n as f64;
        let x, cx = statistics::bernoulli_ppf(p, 0.5);
        let d, cd = statistics::betabinom_cdf(x, 10.0, 2.0, 3.0);
        if cx != 0 || cd != 0 { bad += 1; }
        total += x + d;
    }
    // The search of the quantile of the unbounded laws runs a step for
    // each count of the law.
    let m, cm = statistics::betanbinom_ppf(0.5, 10.0, 2.0, 3.0);
    chk_code(bad, "betanbinom ppf code", cm, 0);
    let m2, cm2 = statistics::logser_ppf(0.5, 0.5);
    chk_code(bad, "logser ppf code", cm2, 0);
    let m3, cm3 = statistics::dlaplace_ppf(0.5, 0.5);
    chk_code(bad, "dlaplace ppf code", cm3, 0);
    total += m + m2 + m3;
    check(bad, !is_nan(total) && !is_inf(total), "the sum of the calls is finite");
    return bad;
}
