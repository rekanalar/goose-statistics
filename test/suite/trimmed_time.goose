// Time tests of the trimmed calls. Each call must run in linear or log-linear
// time on a large sample.
namespace suite;

import std;
import math;
import os;
import statistics from "../..";

fn run_time_trimmed(n: i64) -> i64 {
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
    let v17, c17 = statistics::tvar(xs[..], -0.5, 0.5, true, true, true, true, 1);
    time_shape(bad, v17, c17, "tvar", t0);
    t0 = clock();
    let _n18, c18 = statistics::trimboth(out, xs[..], 0.1);
    time_shape(bad, out[0], c18, "trimboth", t0);
    return bad;
}
