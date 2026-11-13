// Time tests of the truncated law, the order statistic and the Dirichlet
// multinomial at their largest accepted input.
namespace suite;

import std;
import math;
import os;
import statistics from "../..";

fn run_time_laws(k: i64) -> i64 {
    var bad: i64 = 0;
    var alpha: f64[>..] = [];
    var counts: f64[>..] = [];
    for i in 0..k { alpha.push(1.0 + (i % 5)); counts.push(1.0); }
    let t0 = clock();
    let v, c = statistics::dirichlet_multinomial_logpmf(counts[..], alpha[..], k);
    time_shape(bad, v, c, "dirichlet_multinomial_logpmf", t0);
    var out: f64[>..] = [];
    let t1 = clock();
    let c1 = statistics::dirichlet_multinomial_var(out, alpha[..], k);
    time_shape(bad, 0.0, c1, "dirichlet_multinomial_var", t1);
    let t2 = clock();
    let m, cm = statistics::order_statistic_mean(0, 0.0, 1.0, 500000, statistics::MAX_ORDER_N);
    time_shape(bad, m, cm, "order_statistic_mean", t2);
    let t3 = clock();
    let tv, ct = statistics::truncate_var(6, 0.0, 1.0, -3.0, 40.0);
    time_shape(bad, tv, ct, "truncate_var", t3);
    return bad;
}
