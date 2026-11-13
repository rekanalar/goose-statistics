// Hostile input tests of the circular calls. Each answer must be a value with
// code 0, or NaN with a known fault code.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_circular_hostile() -> i64 {
    var bad: i64 = 0;
    var st: u64 = 246813589;
    let n = tests_count();
    var out: f64[>..] = [];
    var out2: f64[>..] = [];
    for k in 0..n {
        let hostile = rnd(st) % 3 != 0;
        let xs = hostile_sample(st, hostile);
        var ys = hostile_sample(st, hostile);
        if rnd(st) % 2 == 0 {
            ys.resize(0);
            for v in xs { ys.push(if hostile { hostile_value(st) } else { rnd(st) % 7 - 3 }); }
        }
        let per = (rnd(st) % 120) as f64;
        let alpha = (rnd(st) % 120) as f64 / 100.0;
        let level = (rnd(st) % 120) as f64 / 100.0;
        let method = rnd(st) % 12;
        let kind = rnd(st) % 8;
        let stat = rnd(st) % 10;
        let tail = rnd(st) % 4;
        let prop = (rnd(st) % 90) as f64 / 100.0;

        let v30, c30 = statistics::circmean(xs[..], 0.0, alpha);
        hostile_report(bad, sound_value(v30, c30), "circmean", k, c30);
        let v31, c31 = statistics::circvar(xs[..], 0.0, alpha);
        hostile_report(bad, sound_value(v31, c31), "circvar", k, c31);
        let v32, c32 = statistics::circstd(xs[..], 0.0, alpha, rnd(st) % 2 == 0);
        hostile_report(bad, sound_value(v32, c32), "circstd", k, c32);
        let v33, c33 = statistics::directional_stats(out, xs[..], (rnd(st) % 4) as i64, rnd(st) % 2 == 0);
        hostile_report(bad, sound_value(v33, c33), "directional_stats", k, c33);
    }
    return bad;
}
