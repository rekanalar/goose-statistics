// Hostile input tests of the trimmed calls. Each answer must be a value with
// code 0, or NaN with a known fault code.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_trimmed_hostile() -> i64 {
    var bad: i64 = 0;
    var st: u64 = 246813587;
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

        let v25, c25 = statistics::tmin(xs[..], alpha, rnd(st) % 2 == 0);
        hostile_report(bad, sound_value(v25, c25), "tmin", k, c25);
        let v26, c26 = statistics::tmax(xs[..], alpha, rnd(st) % 2 == 0);
        hostile_report(bad, sound_value(v26, c26), "tmax", k, c26);
        let v27, c27 = statistics::tvar(xs[..], alpha, alpha, rnd(st) % 2 == 0, rnd(st) % 2 == 0, true, true, (rnd(st) % 3) as i64);
        hostile_report(bad, sound_value(v27, c27), "tvar", k, c27);
        let n28, c28 = statistics::trimboth(out, xs[..], prop);
        hostile_report(bad, sound_count(n28, c28, out.len), "trimboth", k, c28);
        let n29, c29 = statistics::trim1(out, xs[..], prop, tail);
        hostile_report(bad, sound_count(n29, c29, out.len), "trim1", k, c29);
        let v43, c43 = statistics::tstd(xs[..], alpha, level, rnd(st) % 2 == 0, rnd(st) % 2 == 0, rnd(st) % 2 == 0, rnd(st) % 2 == 0, (rnd(st) % 3) as i64);
        hostile_report(bad, sound_value(v43, c43), "tstd", k, c43);
        let v44, c44 = statistics::tsem(xs[..], alpha, level, rnd(st) % 2 == 0, rnd(st) % 2 == 0, rnd(st) % 2 == 0, rnd(st) % 2 == 0, (rnd(st) % 3) as i64);
        hostile_report(bad, sound_value(v44, c44), "tsem", k, c44);
    }
    return bad;
}
