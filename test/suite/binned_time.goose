// Time tests of the binned calls. Each call must run in linear or log-linear
// time on a large sample.
namespace suite;

import std;
import math;
import os;
import statistics from "../..";

fn run_time_binned(n: i64) -> i64 {
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
    let t20, c20 = statistics::binned_statistic(out, out2, out3, out4, out5, xs[..], ys[..], statistics::STAT_MEDIAN, 1000);
    time_shape(bad, t20.cells, c20, "binned_statistic median", t0);
    t0 = clock();
    let t21, c21 = statistics::binned_statistic_2d(out, out2, out3, out4, out5, xs[..], ys[..], xs[..], statistics::STAT_MEAN, 200, 200);
    time_shape(bad, t21.cells, c21, "binned_statistic_2d", t0);
    return bad;
}
