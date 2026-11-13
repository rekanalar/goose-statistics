// Time tests of the distance of two sets of points at the largest size.
namespace suite;

import std;
import math;
import os;
import statistics from "../..";

fn run_time_distance_nd(n: i64) -> i64 {
    var bad: i64 = 0;
    var st: u64 = 7733;
    var u: f64[>..] = [];
    var v: f64[>..] = [];
    for i in 0..n * 2 { u.push((rnd(st) % 1000) as f64 / 10.0); }
    for i in 0..n * 2 { v.push((rnd(st) % 1000) as f64 / 10.0 + 5.0); }
    let t0 = clock();
    let d, c = statistics::wasserstein_distance_nd(u[..], v[..], 2);
    time_shape(bad, d, c, "wasserstein_distance_nd", t0);
    var w: f64[>..] = [];
    for i in 0..n { w.push((rnd(st) % 5) as f64 + 1.0); }
    let t1 = clock();
    let dw, cw = statistics::wasserstein_distance_nd_weighted(u[..], v[..], 2, w[..], w[..]);
    time_shape(bad, dw, cw, "wasserstein_distance_nd weighted", t1);
    return bad;
}
