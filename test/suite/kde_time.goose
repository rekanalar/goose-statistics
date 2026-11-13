// Time tests of the kernel density estimate on a large sample.
namespace suite;

import std;
import math;
import os;
import statistics from "../..";

fn run_time_kde() -> i64 {
    var bad: i64 = 0;
    var st: u64 = 5521;
    var data: f64[>..] = [];
    for i in 0..5000 { data.push((rnd(st) % 1000) as f64 / 10.0); }
    var points: f64[>..] = [];
    for i in 0..100 { points.push((rnd(st) % 1000) as f64 / 10.0); }
    let t0 = clock();
    var k, ck = statistics::kde_fit(data[..], 1, [], statistics::KDE_SCOTT, 0.0);
    var out: f64[>..] = [];
    let cp = statistics::kde_pdf(k, out, points[..]);
    time_shape(bad, out.len as f64, cp, "kde", t0);
    let t1 = clock();
    var k2, ck2 = statistics::kde_fit(data[..], 2, [], statistics::KDE_SCOTT, 0.0);
    time_shape(bad, k2.factor, ck2, "kde two", t1);
    return bad;
}
