// Time tests of the probability plot calls at large inputs.
namespace suite;

import std;
import math;
import os;
import statistics from "../..";

fn run_time_probplot(n: i64) -> i64 {
    var bad: i64 = 0;
    var st: u64 = 7311;
    let xs = time_sample(n, false, st);
    var osm: f64[>..] = [];
    var osr: f64[>..] = [];
    let t0 = clock();
    let p, c = statistics::probplot(osm, osr, xs[..]);
    time_shape(bad, p.slope, c, "probplot", t0);
    let small = time_sample(2000, false, st);
    var grid: f64[>..] = [];
    let t1 = clock();
    let c1 = statistics::ppcc_plot(grid, small[..], -2.0, 2.0, 200);
    time_shape(bad, 0.0, c1, "ppcc_plot", t1);
    let t2 = clock();
    let s, c2 = statistics::ppcc_max(small[..]);
    time_shape(bad, s, c2, "ppcc_max", t2);
    return bad;
}
