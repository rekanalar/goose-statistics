// Time tests of the sigma clip on a large sample.
namespace suite;

import std;
import math;
import os;
import statistics from "../..";

fn run_time_clip(n: i64) -> i64 {
    var bad: i64 = 0;
    var st: u64 = 4141;
    let xs = time_sample(n, false, st);
    var out: f64[>..] = [];
    let t0 = clock();
    let r, c = statistics::sigmaclip(out, xs[..], 3.0, 3.0);
    time_shape(bad, r.high, c, "sigmaclip", t0);
    return bad;
}
