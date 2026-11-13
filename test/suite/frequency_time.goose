// Time tests of the frequency calls. Each call must run in linear or log-linear
// time on a large sample.
namespace suite;

import std;
import math;
import os;
import statistics from "../..";

fn run_time_frequency(n: i64) -> i64 {
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
    let v11, c11 = statistics::estimated_cdf(xs[..], 0.1, statistics::CDF_LINEAR);
    time_shape(bad, v11, c11, "estimated_cdf", t0);
    t0 = clock();
    let _n12, c12 = statistics::ecdf(out, out2, xs[..]);
    time_shape(bad, out[0], c12, "ecdf", t0);
    t0 = clock();
    let _n13, c13 = statistics::ecdf_interval(out, out2, xs[..], 0.95);
    time_shape(bad, out[0], c13, "ecdf_interval", t0);
    t0 = clock();
    let _n14, c14 = statistics::ecdf_censored(out, out2, xs[..], ys[..]);
    time_shape(bad, out[0], c14, "ecdf_censored", t0);
    t0 = clock();
    let f15, c15 = statistics::cumfreq(out, xs[..], 1000);
    time_shape(bad, f15.binsize, c15, "cumfreq", t0);
    t0 = clock();
    let v16, c16 = statistics::tiecorrect(xs[..]);
    time_shape(bad, v16, c16, "tiecorrect", t0);
    return bad;
}
