// Hostile input tests of the entropy calls. Each answer must be a value with
// code 0, or NaN with a known fault code.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_entropy_hostile() -> i64 {
    var bad: i64 = 0;
    var st: u64 = 246813583;
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

        let v12, c12 = statistics::entropy(xs[..], alpha);
        hostile_report(bad, sound_value(v12, c12), "entropy", k, c12);
        let v13, c13 = statistics::relative_entropy(xs[..], ys[..], alpha);
        hostile_report(bad, sound_value(v13, c13), "relative_entropy", k, c13);
        let v14, c14 = statistics::kstat(xs[..], (rnd(st) % 7) as i64);
        hostile_report(bad, sound_value(v14, c14), "kstat", k, c14);
        let v15, c15 = statistics::kstatvar(xs[..], (rnd(st) % 5) as i64);
        hostile_report(bad, sound_value(v15, c15), "kstatvar", k, c15);
        let v16, c16 = statistics::differential_entropy(xs[..], (rnd(st) % 8) as i64, (rnd(st) % 4) as i64, alpha);
        hostile_report(bad, sound_value(v16, c16), "differential_entropy", k, c16);
    }
    return bad;
}
