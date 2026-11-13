// The laws of the distance calls: the distance of a sample with itself is
// zero, the distance grows with the spread of the samples, and the
// transform of a group keeps the mean of the squares of the group.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distance_tests() -> i64 {
    var bad: i64 = 0;
    var out: f64[>..] = [];
    let u: f64[6] = [0.0, 1.0, 3.0, 5.0, 8.0, 13.0];
    let v: f64[5] = [0.5, 2.0, 2.5, 7.0, 20.0];
    let keep: f64[6] = [0.0, 1.0, 3.0, 5.0, 8.0, 13.0];
    // the distance of a sample with itself is zero
    let w, cw = statistics::wasserstein_distance(u[..], u[..]);
    chk_code(bad, "wasserstein of a sample with itself", cw, 0);
    chk_val(bad, "wasserstein of a sample with itself", w, cw, 0.0, 1e-12);
    let e, ce = statistics::energy_distance(u[..], u[..]);
    chk_code(bad, "energy of a sample with itself", ce, 0);
    chk_val(bad, "energy of a sample with itself", e, ce, 0.0, 1e-12);
    // the distance grows with the spread of the second sample
    let near = [0.0, 1.0, 3.0, 5.0, 8.0, 13.0];
    let far = [10.0, 11.0, 13.0, 15.0, 18.0, 23.0];
    let d1, c1 = statistics::wasserstein_distance(u[..], v[..]);
    let d2, c2 = statistics::wasserstein_distance(near[..], far[..]);
    chk_code(bad, "the distance of a shifted sample", c2, 0);
    check(bad, d2 > d1, "the distance grows with the spread");
    // the distance of a sample and a shifted copy of it is the shift
    let shifted = [10.0, 11.0, 13.0, 15.0, 18.0, 23.0];
    let d3, c3 = statistics::wasserstein_distance(u[..], shifted[..]);
    chk_val(bad, "the distance of a shifted sample", d3, c3, 10.0, 1e-12);
    // the energy distance of a shifted sample grows with the shift
    let far2 = [20.0, 21.0, 23.0, 25.0, 28.0, 33.0];
    let d4, c4 = statistics::energy_distance(u[..], shifted[..]);
    let d5, c5 = statistics::energy_distance(u[..], far2[..]);
    chk_code(bad, "the energy of a shifted sample", c5, 0);
    check(bad, d5 > d4, "the energy distance grows with the shift");
    // the transform of a group keeps the spread of the group
    let g: f64[6] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    let sizes: i64[1] = [6];
    let c5 = statistics::obrientransform(out, g[..], sizes[..]);
    chk_code(bad, "the transform of one group", c5, 0);
    check(bad, out.len == 6, "the transform writes one value for each value");
    // the mean of the transformed values is the variance of the group
    let m = mean_of(out[..]);
    let v2, cv = statistics::variance(g[..]);
    chk_code(bad, "the variance of the group", cv, 0);
    chk_val(bad, "the mean of the transformed values", m, c5, v2, 1e-12);
    // the samples of the caller are unchanged
    var same = true;
    for i in 0..6 { if u[i] != keep[i] { same = false; } }
    check(bad, same, "the distance calls leave the samples unchanged");
    // a fault leaves the output empty
    let bad_sizes: i64[1] = [2];
    let c6 = statistics::obrientransform(out, g[..], bad_sizes[..]);
    check(bad, out.len == 0, "a fault leaves the output empty");
    return bad;
}
