// Hostile input tests of the goodness-of-fit test against gamma, beta and
// Weibull. Each answer is a finite statistic with code 0, or a fault code.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_fitting_mle_hostile() -> i64 {
    var bad: i64 = 0;
    var st: u64 = 51907;
    let n = tests_count() / 8;
    var out: f64[>..] = [];
    var rg = rng { s: 3 };
    for k in 0..n {
        let hostile = rnd(st) % 3 != 0;
        let xs = hostile_sample(st, hostile);
        let fam = 4 + rnd(st) % 3;
        let stat = rnd(st) % 3;
        let t, c = statistics::goodness_of_fit(out, xs[..], fam, stat, 5, rg);
        var ok = false;
        if c == 0 { ok = !is_nan(t.statistic) && !is_inf(t.statistic) && t.pvalue > 0.0 && t.pvalue <= 1.0 && out.len >= 3; }
        else { ok = c >= 1 && c <= 10 && is_nan(t.statistic) && out.len == 0; }
        hostile_report(bad, ok, "goodness_of_fit mle", k, c);
    }
    // a fit above the bound of the repetitions gives the too long fault
    let gx: f64[4] = [1.0, 2.0, 3.0, 4.0];
    var gout: f64[>..] = [];
    var grg = rng { s: 3 };
    let gr, gc = statistics::goodness_of_fit(gout, gx[..], 5, 0, statistics::MAX_GOF_REPS + 1, grg);
    chk_code(bad, "goodness of fit above the repetition bound", gc, 9);
    return bad;
}
