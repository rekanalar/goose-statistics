// Time tests of the power divergence on a long list of counts.
namespace suite;

import std;
import math;
import os;
import statistics from "../..";

fn run_time_divergence(n: i64) -> i64 {
    var bad: i64 = 0;
    var obs: f64[>..] = [];
    var exp: f64[>..] = [];
    for i in 0..n { obs.push(10.0 + (i % 7)); exp.push(10.0 + (i % 7)); }
    var t0 = clock();
    let r, c = statistics::power_divergence(obs[..], exp[..], 0.0);
    time_shape(bad, r.statistic + 1.0, c, "power_divergence", t0);
    return bad;
}
