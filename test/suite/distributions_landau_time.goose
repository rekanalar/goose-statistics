// Time rows of the Landau law.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_landau_time() -> i64 {
    var bad: i64 = 0;
    var total: f64 = 0.0;
    for k in 0..5 {
        let p = (k as f64 + 0.5) / 5.0;
        let x, c = statistics::landau_ppf(p);
        if c != 0 { bad += 1; }
        total += x;
    }
    check(bad, !is_nan(total) && !is_inf(total), "the sum of the calls is finite");
    return bad;
}
