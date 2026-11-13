// Time tests of the rank association calls on a long sample.
namespace suite;

import std;
import math;
import os;
import statistics from "../..";

fn run_time_assoc_rank(n: i64) -> i64 {
    var bad: i64 = 0;
    var st: u64 = 8181;
    let xs = time_sample(n, false, st);
    let ys = time_sample(n, false, st);
    var t0 = clock();
    let a1, c1 = statistics::spearmanrho(xs[..], ys[..]);
    time_shape(bad, a1.statistic, c1, "spearmanrho", t0);
    t0 = clock();
    let a2, c2 = statistics::chatterjeexi(xs[..], ys[..]);
    time_shape(bad, a2.statistic, c2, "chatterjeexi", t0);
    t0 = clock();
    let w1, c3 = statistics::weightedtau(xs[..], ys[..]);
    time_shape(bad, w1, c3, "weightedtau", t0);
    t0 = clock();
    let b1, c4 = statistics::pointbiserialr(xs[..], ys[..]);
    time_shape(bad, b1.statistic, c4, "pointbiserialr", t0);
    return bad;
}
