// Time tests of the goodness-of-fit test against the chi-square, student
// and F families.
namespace suite;

import std;
import math;
import os;
import statistics from "../..";

fn run_time_fitting_shapes(n: i64) -> i64 {
    var bad: i64 = 0;
    var out: f64[>..] = [];
    var rg = rng { s: 9 };
    var xs: f64[>..] = [];
    for i in 0..n { xs.push(0.5 + 3.0 * (i + 0.5) / n); }
    for f in 7..10 {
        let t0 = clock();
        let t, c = statistics::goodness_of_fit(out, xs[..], f, 3, 3, rg);
        time_shape(bad, t.statistic, c, "goodness_of_fit shapes", t0);
    }
    return bad;
}
