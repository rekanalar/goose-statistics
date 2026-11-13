// Hostile input tests of the binned calls. Each answer must be a value with
// code 0, or NaN with a known fault code.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_binned_hostile() -> i64 {
    var bad: i64 = 0;
    var st: u64 = 246813591;
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

        let t34, c34 = statistics::binned_statistic(out, out2, out2, out2, out2, xs[..], ys[..], stat, (rnd(st) % 6) as i64);
        hostile_report(bad, sound_table(t34, c34), "binned_statistic", k, c34);
        let t35, c35 = statistics::binned_statistic_2d(out, out2, out2, out2, out2, xs[..], ys[..], ys[..], stat, 2, 2);
        hostile_report(bad, sound_table(t35, c35), "binned_statistic_2d", k, c35);
        let t36, c36 = statistics::binned_statistic_dd(out, out2, out2, out2, out2, xs[..], (rnd(st) % 4) as i64, ys[..], stat, 2);
        hostile_report(bad, sound_table(t36, c36), "binned_statistic_dd", k, c36);
        let ge = edges_of(st);
        let t49, c49 = statistics::binned_statistic_edges(out, out2, out2, out2, out2, xs[..], ys[..], stat, ge[..]);
        hostile_report(bad, sound_table(t49, c49), "binned_statistic_edges", k, c49);
        let t50, c50 = statistics::binned_statistic_with(out, out2, out2, out2, out2, xs[..], ys[..], (rnd(st) % 6) as i64) { v => sum_of(v) };
        hostile_report(bad, sound_table(t50, c50), "binned_statistic_with", k, c50);
        let t51, c51 = statistics::binned_statistic_2d_with(out, out2, out2, out2, out2, xs[..], ys[..], ys[..], 2, 2) { v => sum_of(v) };
        hostile_report(bad, sound_table(t51, c51), "binned_statistic_2d_with", k, c51);
        let t52, c52 = statistics::binned_statistic_dd_with(out, out2, out2, out2, out2, xs[..], (rnd(st) % 4) as i64, ys[..], 2) { v => sum_of(v) };
        hostile_report(bad, sound_table(t52, c52), "binned_statistic_dd_with", k, c52);
    }
    return bad;
}
