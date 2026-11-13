// Time tests of the goodness-of-fit test against gamma, beta and Weibull.
namespace suite;

import std;
import math;
import os;
import statistics from "../..";

fn run_time_fitting_mle(n: i64) -> i64 {
    var bad: i64 = 0;
    var out: f64[>..] = [];
    var rg = rng { s: 9 };
    var xs: f64[>..] = [];
    for i in 0..n { xs.push(0.001 + 0.998 * (i + 0.5) / n); }
    for f in 4..7 {
        let t0 = clock();
        let t, c = statistics::goodness_of_fit(out, xs[..], f, 0, 5, rg);
        time_shape(bad, t.statistic, c, "goodness_of_fit mle", t0);
    }
    return bad;
}
