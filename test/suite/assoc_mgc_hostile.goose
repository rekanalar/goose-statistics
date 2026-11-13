// Hostile input tests of the multiscale graph correlation. Each answer is
// a statistic of at most 1 and a p-value in 0 to 1 with code 0, or a fault
// code with NaN.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_assoc_mgc_hostile() -> i64 {
    var bad: i64 = 0;
    var st: u64 = 40213;
    var map: f64[>..] = [];
    var gen = rng { s: 3 };
    let n = tests_count() / 10;
    for k in 0..n {
        let xd = 1 + rnd(st) % 3;
        let yd = 1 + rnd(st) % 2;
        var pts = 3 + rnd(st) % 9;
        let hostile = rnd(st) % 4 == 0;
        var xs: f64[>..] = [];
        var ys: f64[>..] = [];
        for i in 0..pts * xd {
            let r = rnd(st) % 60;
            if hostile && r == 0 { xs.push(0.0 / 0.0); }
            else if hostile && r == 1 { xs.push(1.0 / 0.0); }
            else if hostile && r == 2 { xs.push(1e200); }
            else { xs.push((rnd(st) % 6) as f64); }
        }
        for i in 0..pts * yd { ys.push((rnd(st) % 5) as f64); }
        if rnd(st) % 30 == 0 { ys.push(1.0); }
        let reps = rnd(st) % 4;
        let res, c = statistics::multiscale_graphcorr(map, xs[..], ys[..], xd, yd, reps, gen);
        var ok = false;
        if c == 0 {
            ok = !is_nan(res.statistic) && res.statistic <= 1.0 && res.pvalue > 0.0 && res.pvalue <= 1.0 && map.len == res.rows * res.cols;
        } else {
            ok = c >= 1 && c <= 10 && is_nan(res.statistic) && map.len == 0;
        }
        hostile_report(bad, ok, "multiscale_graphcorr", k, c);
    }
    return bad;
}
