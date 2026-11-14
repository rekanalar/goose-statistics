// The rows of the Rice law that the vector file does not hold: the
// support, the end points of the quantile, the law with the shape 0 and
// the round trip of the quantile.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_rice_tests() -> i64 {
    var bad: i64 = 0;
    // Below the support the density is 0 and the cumulative function is 0.
    let v1, c1 = statistics::rice_pdf(-1.0, 2.0);
    chk_val(bad, "rice pdf below the support", v1, c1, 0.0, 1e-12);
    let v2, c2 = statistics::rice_cdf(-1.0, 2.0);
    chk_val(bad, "rice cdf below the support", v2, c2, 0.0, 1e-12);
    let v3, c3 = statistics::rice_sf(-1.0, 2.0);
    chk_val(bad, "rice sf below the support", v3, c3, 1.0, 1e-12);
    // The shape 0 gives the Rayleigh law: the cumulative function at 1 is
    // 1 - exp(-1/2).
    let v4, c4 = statistics::rice_cdf(1.0, 0.0);
    chk_val(bad, "rice cdf with the shape 0", v4, c4, 1.0 - exp(-0.5), 1e-14);
    let v5, c5 = statistics::rice_sf(1.0, 0.0);
    chk_val(bad, "rice sf with the shape 0", v5, c5, exp(-0.5), 1e-14);
    // The quantile of the end points.
    let e1, ce1 = statistics::rice_ppf(0.0, 2.0);
    chk_val(bad, "rice ppf 0", e1, ce1, 0.0, 1e-12);
    let e2, ce2 = statistics::rice_ppf(1.0, 2.0);
    chk_code(bad, "rice ppf 1 is infinite", ce2, 8);
    // The cumulative function and the survival function add up to 1.
    let a1, ka1 = statistics::rice_cdf(2.5, 3.0);
    let a2, ka2 = statistics::rice_sf(2.5, 3.0);
    chk_val(bad, "rice cdf plus sf", a1 + a2, ka1 + ka2, 1.0, 1e-13);
    // The cumulative function of the quantile gives the probability back,
    // in the lower and in the upper half.
    let q1, kq1 = statistics::rice_ppf(0.2, 4.0);
    let r1, kr1 = statistics::rice_cdf(q1, 4.0);
    chk_val(bad, "rice cdf of the quantile", r1, kq1 + kr1, 0.2, 1e-10);
    let q2, kq2 = statistics::rice_ppf(0.95, 4.0);
    let r2, kr2 = statistics::rice_cdf(q2, 4.0);
    chk_val(bad, "rice cdf of the upper quantile", r2, kq2 + kr2, 0.95, 1e-10);
    // A large shape: the mean is near the shape.
    let m, cm = statistics::rice_mean(200.0);
    chk_val(bad, "rice mean of a large shape", m, cm, 200.0025000156, 1e-6);
    // A shape outside its range gives the domain fault.
    let _d1, cd1 = statistics::rice_pdf(1.0, -0.5);
    chk_code(bad, "rice shape below 0", cd1, 5);
    let _d2, cd2 = statistics::rice_pdf(1.0, 4000.0);
    chk_code(bad, "rice shape above the bound", cd2, 5);
    return bad;
}
