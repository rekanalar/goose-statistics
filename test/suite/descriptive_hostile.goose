// Hostile input tests of the descriptive calls. Each answer must be a value with
// code 0, or NaN with a known fault code.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_descriptive_hostile() -> i64 {
    var bad: i64 = 0;
    var st: u64 = 246813581;
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

        let v1, c1 = statistics::scoreatpercentile(xs[..], per);
        hostile_report(bad, sound_value(v1, c1), "scoreatpercentile", k, c1);
        let v2, c2 = statistics::scoreatpercentile_method(xs[..], per, method);
        hostile_report(bad, sound_value(v2, c2), "scoreatpercentile_method", k, c2);
        let v3, c3 = statistics::percentileofscore(xs[..], alpha, kind);
        hostile_report(bad, sound_value(v3, c3), "percentileofscore", k, c3);
        let v4, c4 = statistics::median_abs_deviation(xs[..], 1.0, false);
        hostile_report(bad, sound_value(v4, c4), "median_abs_deviation", k, c4);
        let v5, c5 = statistics::pmean(xs[..], alpha);
        hostile_report(bad, sound_value(v5, c5), "pmean", k, c5);
        let v6, c6 = statistics::gstd(xs[..], (rnd(st) % 3) as i64);
        hostile_report(bad, sound_value(v6, c6), "gstd", k, c6);
        let n7, c7 = statistics::gzscore(out, xs[..], (rnd(st) % 3) as i64);
        hostile_report(bad, sound_count(n7, c7, out.len), "gzscore", k, c7);
        let n8, c8 = statistics::zmap(out, xs[..], ys[..], (rnd(st) % 3) as i64);
        hostile_report(bad, sound_count(n8, c8, out.len), "zmap", k, c8);
        let v9, c9 = statistics::describe(xs[..], (rnd(st) % 3) as i64, rnd(st) % 2 == 0);
        hostile_report(bad, sound_summary(v9, c9), "describe", k, c9);
        let n10, c10 = statistics::lmoment(out, xs[..], (rnd(st) % 8) as i64, rnd(st) % 2 == 0);
        hostile_report(bad, sound_count(n10, c10, out.len), "lmoment", k, c10);
        let v11, c11 = statistics::expectile(xs[..], alpha);
        hostile_report(bad, sound_value(v11, c11), "expectile", k, c11);
        let v37, c37 = statistics::scoreatpercentile_limit(xs[..], per, alpha, level);
        hostile_report(bad, sound_value(v37, c37), "scoreatpercentile_limit", k, c37);
        let v38, c38 = statistics::scoreatpercentile_limit_method(xs[..], per, alpha, level, method);
        hostile_report(bad, sound_value(v38, c38), "scoreatpercentile_limit_method", k, c38);
        let v39, c39 = statistics::median_abs_deviation_normal(xs[..]);
        hostile_report(bad, sound_value(v39, c39), "median_abs_deviation_normal", k, c39);
        let ws = weights_of(st, xs.len);
        let v41, c41 = statistics::pmean_weighted(xs[..], ws[..], alpha);
        hostile_report(bad, sound_value(v41, c41), "pmean_weighted", k, c41);
        let v42, c42 = statistics::expectile_weighted(xs[..], ws[..], alpha);
        hostile_report(bad, sound_value(v42, c42), "expectile_weighted", k, c42);
    }
    return bad;
}
