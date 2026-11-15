// Time rows of the two non-central hypergeometric laws.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_nch_time() -> i64 {
    var bad: i64 = 0;
    var total: f64 = 0.0;
    for k in 0..10 {
        let p = (k as f64 + 0.5) / 10.0;
        let a, ka = statistics::nchypergeom_fisher_ppf(p, 500.0, 200.0, 150.0, 2.5);
        let b, kb = statistics::nchypergeom_wallenius_ppf(p, 100.0, 80.0, 60.0, 2.5);
        if ka != 0 || kb != 0 { bad += 1; }
        total += a + b;
    }
    check(bad, !is_nan(total) && !is_inf(total), "the sum of the calls is finite");
    return bad;
}
