// Time tests of the descriptive calls. Each call must run in linear or log-linear
// time on a large sample.
namespace suite;

import std;
import math;
import os;
import statistics from "../..";

fn run_time_descriptive(n: i64) -> i64 {
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
    let v1, c1 = statistics::scoreatpercentile(xs[..], 37.0);
    time_shape(bad, v1, c1, "scoreatpercentile", t0);
    t0 = clock();
    let v2, c2 = statistics::percentileofscore(xs[..], 0.5, statistics::PCT_MEAN);
    time_shape(bad, v2, c2, "percentileofscore", t0);
    t0 = clock();
    let v3, c3 = statistics::median_abs_deviation(xs[..], 1.0, false);
    time_shape(bad, v3, c3, "median_abs_deviation", t0);
    t0 = clock();
    let v4, c4 = statistics::pmean(pos[..], 3.0);
    time_shape(bad, v4, c4, "pmean", t0);
    t0 = clock();
    let v5, c5 = statistics::gstd(pos[..], 1);
    time_shape(bad, v5, c5, "gstd", t0);
    t0 = clock();
    let s6, c6 = statistics::describe(xs[..], 1, false);
    time_shape(bad, s6.variance, c6, "describe", t0);
    t0 = clock();
    let _n7, c7 = statistics::lmoment(out, xs[..], 4, true);
    time_shape(bad, out[0], c7, "lmoment", t0);
    t0 = clock();
    let v8, c8 = statistics::expectile(xs[..], 0.3);
    time_shape(bad, v8, c8, "expectile", t0);
    return bad;
}
