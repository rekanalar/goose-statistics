// Time test of the Siegel estimate at the largest accepted sample.
namespace suite;

import std;
import math;
import os;
import statistics from "../..";

fn run_time_line_siegel(n: i64) -> i64 {
    var bad: i64 = 0;
    var st: u64 = 6161;
    let xs = time_sample(n, false, st);
    let ys = time_sample(n, false, st);
    let t0 = clock();
    let s, b, c = statistics::siegelslopes(ys[..], xs[..]);
    time_shape(bad, s, c, "siegelslopes", t0);
    return bad;
}
