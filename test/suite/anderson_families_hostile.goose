// Hostile input tests of the Anderson-Darling test against four families.
// Each answer is a finite statistic with code 0, or a fault code and empty lists.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_anderson_families_hostile() -> i64 {
    var bad: i64 = 0;
    var st: u64 = 40213;
    let n = tests_count() / 4;
    var crit: f64[>..] = [];
    var lev: f64[>..] = [];
    var par: f64[>..] = [];
    for k in 0..n {
        let hostile = rnd(st) % 3 != 0;
        let xs = hostile_sample(st, hostile);
        let fam = rnd(st) % 9;
        let t, c = statistics::anderson_family(crit, lev, par, xs[..], fam);
        var ok = false;
        if c == 0 { ok = !is_nan(t.statistic) && !is_inf(t.statistic) && t.pvalue >= 0.0 && t.pvalue <= 1.0 && crit.len == lev.len && par.len >= 2; }
        else { ok = c >= 1 && c <= 10 && is_nan(t.statistic) && crit.len == 0 && par.len == 0; }
        hostile_report(bad, ok, "anderson_family", k, c);
    }
    return bad;
}
