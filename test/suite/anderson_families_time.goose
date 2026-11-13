// Time tests of the Anderson-Darling test against four families.
namespace suite;

import std;
import math;
import os;
import statistics from "../..";

fn run_time_anderson_families(n: i64) -> i64 {
    var bad: i64 = 0;
    var st: u64 = 2718;
    let xs = time_sample(n, false, st);
    var crit: f64[>..] = [];
    var lev: f64[>..] = [];
    var par: f64[>..] = [];
    var fam: i64 = 1;
    for f in 1..6 {
        if f == 3 || f == 6 { continue; }
        let t0 = clock();
        let t, c = statistics::anderson_family(crit, lev, par, xs[..], f);
        time_shape(bad, t.statistic, c, "anderson_family", t0);
    }
    return bad;
}
