// Time tests of the Cauchy, Laplace, Rayleigh and Gumbel calls: n calls of
// each form must end in the time limit.
namespace suite;

import std;
import math;
import os;
import statistics from "../..";

fn run_time_distributions_location(n: i64) -> i64 {
    var bad: i64 = 0;
    var s1: f64 = 0.0;
    var s2: f64 = 0.0;
    var s3: f64 = 0.0;
    var c: i64 = 0;
    let t0 = clock();
    for k in 0..n {
        let x = (k % 2001 - 1000) / 100.0;
        let p = (k % 999 + 1) / 1000.0;
        let a, c1 = statistics::cauchy_cdf(x, 0.1, 1.5);
        let b, c2 = statistics::laplace_ppf(p, 0.1, 1.5);
        let d, c3 = statistics::rayleigh_sf(x, -20.0, 1.5);
        let e, c4 = statistics::gumbel_r_pdf(x, 0.1, 1.5);
        let f, c5 = statistics::gumbel_l_ppf(p, 0.1, 1.5);
        s1 += a + b;
        s2 += d + e;
        s3 += f;
        c += c1 + c2 + c3 + c4 + c5;
    }
    time_shape(bad, s1 + s2 + s3, c, "distributions_location", t0);
    return bad;
}
