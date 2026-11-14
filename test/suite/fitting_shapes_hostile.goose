// Hostile input tests of the goodness-of-fit test against the chi-square,
// student and F families and of the Filliben statistic. Each answer is a
// finite statistic with code 0, or a fault code.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_fitting_shapes_hostile() -> i64 {
    var bad: i64 = 0;
    var st: u64 = 61231;
    let n = tests_count() / 400;
    var out: f64[>..] = [];
    var rg = rng { s: 3 };
    for k in 0..n {
        let hostile = rnd(st) % 3 != 0;
        let xs = hostile_sample(st, hostile);
        let fam = rnd(st) % 10;
        let stat = rnd(st) % 4;
        let t, c = statistics::goodness_of_fit(out, xs[..], fam, stat, 3, rg);
        var ok = false;
        if c == 0 { ok = !is_nan(t.statistic) && !is_inf(t.statistic) && t.pvalue > 0.0 && t.pvalue <= 1.0 && out.len >= 2; }
        else { ok = c >= 1 && c <= 10 && is_nan(t.statistic) && out.len == 0; }
        hostile_report(bad, ok, "goodness_of_fit shapes", k, c);
    }
    return bad;
}
