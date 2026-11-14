// Time rows of the multivariate families: a long run of calls ends in a
// bounded time and the answers stay finite.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_multivariate_time() -> i64 {
    var bad: i64 = 0;
    var total: f64 = 0.0;
    var g = rng { s: 77 };
    var out: f64[>..] = [];
    for k in 0..10 {
        let c = statistics::multivariate_normal_rvs(out, [1.0, 2.0], [2.0, 0.5, 0.5, 1.0], 2, 5, g);
        if c != 0 { bad += 1; }
        total += out[0];
    }
    let c2 = statistics::multinomial_rvs(out, 20, [0.5, 0.3, 0.2], 10, g);
    if c2 != 0 { bad += 1; }
    let c3 = statistics::wishart_rvs(out, [2.0, 0.5, 0.5, 1.0], 5.0, 2, 5, g);
    if c3 != 0 { bad += 1; }
    total += out[0];
    check(bad, !is_nan(total) && !is_inf(total), "the sum of the calls is finite");
    return bad;
}
