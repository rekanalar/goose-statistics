// Time rows of the eight laws of the second tail group: a long run of calls ends in a bounded time and
// the answers stay finite.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_noncentral_time() -> i64 {
    var bad: i64 = 0;
    var total: f64 = 0.0;
    let n = 5;
    for k in 0..n {
        let p = (k as f64 + 0.5) / n;
        let x0, c0 = statistics::crystalball_ppf(p, 1.5, 3.0);
        if c0 != 0 { bad += 1; }
        total += x0;
        let x1, c1 = statistics::kstwobign_ppf(p);
        if c1 != 0 { bad += 1; }
        total += x1;
        let x2, c2 = statistics::vonmises_line_ppf(p, 2.0);
        if c2 != 0 { bad += 1; }
        total += x2;
        let x3, c3 = statistics::ncx2_ppf(p, 3.0, 2.0);
        if c3 != 0 { bad += 1; }
        total += x3;
        let x4, c4 = statistics::ncf_ppf(p, 4.0, 7.0, 1.5);
        if c4 != 0 { bad += 1; }
        total += x4;
        let x5, c5 = statistics::nct_ppf(p, 6.0, 1.0);
        if c5 != 0 { bad += 1; }
        total += x5;
        let x6, c6 = statistics::ksone_ppf(p, 50.0);
        if c6 != 0 { bad += 1; }
        total += x6;
        let x7, c7 = statistics::kstwo_ppf(p, 20.0);
        if c7 != 0 { bad += 1; }
        total += x7;
    }
    check(bad, !is_nan(total) && !is_inf(total), "the sum of the calls is finite");
    return bad;
}
