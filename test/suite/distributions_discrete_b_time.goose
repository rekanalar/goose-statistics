// Time rows of the second group of discrete families: a long run of calls
// ends in a bounded time and the answers stay finite. The Poisson
// binomial law and the laws whose quantile needs a search cost more than
// the other calls.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_discrete_b_time() -> i64 {
    var bad: i64 = 0;
    var total: f64 = 0.0;
    let n = 50;
    for k in 0..n {
        let p = (k as f64 + 0.5) / n;
        let x, cx = statistics::randint_ppf(p, 0.0, 5.0);
        let d, cd = statistics::zipfian_cdf(x, 2.0, 5.0);
        if cx != 0 || cd != 0 { bad += 1; }
        total += x + d;
    }
    // A list of 200 probabilities builds the masses of the law with a
    // quadratic rule.
    var ps: f64[>..] = [];
    for k in 0..200 { ps.push((k % 10) / 10.0); }
    let m, cm = statistics::poisson_binom_ppf(0.5, ps[..]);
    chk_code(bad, "poisson_binom ppf code", cm, 0);
    let m2, cm2 = statistics::yulesimon_ppf(0.5, 2.0);
    chk_code(bad, "yulesimon ppf code", cm2, 0);
    let m3, cm3 = statistics::zipf_ppf(0.5, 2.0);
    chk_code(bad, "zipf ppf code", cm3, 0);
    let m4, cm4 = statistics::skellam_ppf(0.5, 1.0, 2.0);
    chk_code(bad, "skellam ppf code", cm4, 0);
    total += m + m2 + m3 + m4;
    check(bad, !is_nan(total) && !is_inf(total), "the sum of the calls is finite");
    return bad;
}
