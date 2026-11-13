// Time tests of the distance calls on a large sample.
namespace suite;

import std;
import math;
import os;
import statistics from "../..";

fn run_time_distance(n: i64) -> i64 {
    var bad: i64 = 0;
    var st: u64 = 5150;
    let xs = time_sample(n, false, st);
    let ys = time_sample(n, true, st);
    var out: f64[>..] = [];
    let t0 = clock();
    let w, c = statistics::wasserstein_distance(xs[..], ys[..]);
    time_shape(bad, w, c, "wasserstein_distance", t0);
    let t1 = clock();
    let e, ce = statistics::energy_distance(xs[..], ys[..]);
    time_shape(bad, e, ce, "energy_distance", t1);
    let sizes: i64[1] = [n];
    let t2 = clock();
    let ct = statistics::obrientransform(out, xs[..], sizes[..]);
    time_shape(bad, 0.0, ct, "obrientransform", t2);
    return bad;
}
