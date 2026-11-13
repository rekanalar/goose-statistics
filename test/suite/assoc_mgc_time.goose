// Time tests of the multiscale graph correlation at the largest sizes.
namespace suite;

import std;
import math;
import os;
import statistics from "../..";

fn run_time_assoc_mgc(n: i64) -> i64 {
    var bad: i64 = 0;
    var st: u64 = 8191;
    var xs: f64[>..] = [];
    var ys: f64[>..] = [];
    for i in 0..n { xs.push((rnd(st) % 100000) as f64 / 100.0); }
    for i in 0..n { ys.push(xs[i] * xs[i] / 1000.0 + (rnd(st) % 1000) as f64 / 100.0); }
    var map: f64[>..] = [];
    var gen = rng { s: 17 };
    let t0 = clock();
    let a, c = statistics::multiscale_graphcorr(map, xs[..], ys[..], 1, 1, 0, gen);
    time_shape(bad, a.statistic, c, "multiscale_graphcorr", t0);
    // the most repetitions that the work bound allows at a size of 100
    var sx: f64[>..] = [];
    var sy: f64[>..] = [];
    for i in 0..100 { sx.push(xs[i]); sy.push(ys[i]); }
    let reps = statistics::MAX_MGC_WORK / 10000 - 1;
    let t1 = clock();
    let b, cb = statistics::multiscale_graphcorr(map, sx[..], sy[..], 1, 1, reps, gen);
    time_shape(bad, b.pvalue, cb, "multiscale_graphcorr with repetitions", t1);
    return bad;
}
