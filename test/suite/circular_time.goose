// Time tests of the circular calls. Each call must run in linear or log-linear
// time on a large sample.
namespace suite;

import std;
import math;
import os;
import statistics from "../..";

fn run_time_circular(n: i64) -> i64 {
    var bad: i64 = 0;
    var st: u64 = 5151;
    let xs = time_sample(n, false, st);
    let ys = time_sample(n, false, st);
    var pos: f64[>..] = [];
    for v in xs { pos.push(abs(v) + 1.0); }
    var out: f64[>..] = [];
    var out2: f64[>..] = [];
    var out3: f64[>..] = [];
    var out4: f64[>..] = [];
    var out5: f64[>..] = [];
    var t0 = clock();
    let v19, c19 = statistics::circmean(xs[..], 0.0, 6.283185307179586);
    time_shape(bad, v19, c19, "circmean", t0);
    return bad;
}
