// Time tests of the odds ratio and the relative risk on a large table.
namespace suite;

import std;
import math;
import os;
import statistics from "../..";

fn run_time_contingency_risk() -> i64 {
    var bad: i64 = 0;
    let tab: f64[4] = [2000.0, 1500.0, 2000.0, 4500.0];
    let t0 = clock();
    let s, c = statistics::odds_ratio(tab[..], statistics::ODDS_SAMPLE);
    time_shape(bad, s.statistic, c, "odds_ratio sample", t0);
    let t1 = clock();
    let s2, c2 = statistics::odds_ratio(tab[..], statistics::ODDS_CONDITIONAL);
    time_shape(bad, s2.statistic, c2, "odds_ratio conditional", t1);
    let t2 = clock();
    let r, c3 = statistics::relative_risk(2000.0, 3500.0, 2000.0, 6500.0);
    time_shape(bad, r.statistic, c3, "relative_risk", t2);
    return bad;
}
