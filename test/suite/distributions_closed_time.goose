// Time rows of the fourteen laws of the first tail group: a long run of calls ends in a bounded time and
// the answers stay finite.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_closed_time() -> i64 {
    var bad: i64 = 0;
    var total: f64 = 0.0;
    let n = 100;
    for k in 0..n {
        let p = (k as f64 + 0.5) / n as f64;
        let x0, c0 = statistics::alpha_ppf(p, 1.5);
        if c0 != 0 { bad += 1; }
        total += x0;
        let x1, c1 = statistics::fatiguelife_ppf(p, 0.7);
        if c1 != 0 { bad += 1; }
        total += x1;
        let x2, c2 = statistics::levy_ppf(p);
        if c2 != 0 { bad += 1; }
        total += x2;
        let x3, c3 = statistics::levy_l_ppf(p);
        if c3 != 0 { bad += 1; }
        total += x3;
        let x4, c4 = statistics::logistic_ppf(p);
        if c4 != 0 { bad += 1; }
        total += x4;
        let x5, c5 = statistics::loglaplace_ppf(p, 2.5);
        if c5 != 0 { bad += 1; }
        total += x5;
        let x6, c6 = statistics::lomax_ppf(p, 2.0);
        if c6 != 0 { bad += 1; }
        total += x6;
        let x7, c7 = statistics::pareto_ppf(p, 2.5);
        if c7 != 0 { bad += 1; }
        total += x7;
        let x8, c8 = statistics::powerlaw_ppf(p, 1.8);
        if c8 != 0 { bad += 1; }
        total += x8;
        let x9, c9 = statistics::reciprocal_ppf(p, 1.0, 4.0);
        if c9 != 0 { bad += 1; }
        total += x9;
        let x10, c10 = statistics::rdist_ppf(p, 3.0);
        if c10 != 0 { bad += 1; }
        total += x10;
        let x11, c11 = statistics::skewcauchy_ppf(p, 0.4);
        if c11 != 0 { bad += 1; }
        total += x11;
        let x12, c12 = statistics::wrapcauchy_ppf(p, 0.6);
        if c12 != 0 { bad += 1; }
        total += x12;
        let x13, c13 = statistics::kappa4_ppf(p, 0.3, 0.2);
        if c13 != 0 { bad += 1; }
        total += x13;
    }
    check(bad, !is_nan(total) && !is_inf(total), "the sum of the calls is finite");
    return bad;
}
