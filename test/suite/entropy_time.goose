// Time tests of the entropy calls. Each call must run in linear or log-linear
// time on a large sample.
namespace suite;

import std;
import math;
import os;
import statistics from "../..";

fn run_time_entropy(n: i64) -> i64 {
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
    let v9, c9 = statistics::kstat(xs[..], 4);
    time_shape(bad, v9, c9, "kstat", t0);
    t0 = clock();
    let v10, c10 = statistics::differential_entropy(xs[..], statistics::ENTROPY_VASICEK, 0, 0.0);
    time_shape(bad, v10, c10, "differential_entropy", t0);
    return bad;
}
