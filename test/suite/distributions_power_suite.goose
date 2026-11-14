// The rows of the families of a power that the vector file does not hold:
// the support, the end points of the quantile, the shape bounds and the
// moments that do not exist.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_power_tests() -> i64 {
    var bad: i64 = 0;
    // A value below the support gives a density of 0 and a cumulative of 0.
    let a1, k1 = statistics::fisk_pdf(-1.0, 5.0);
    chk_val(bad, "fisk pdf below the support", a1, k1, 0.0, 1e-12);
    let a2, k2 = statistics::burr_cdf(-1.0, 5.0, 3.0);
    chk_val(bad, "burr cdf below the support", a2, k2, 0.0, 1e-12);
    let a3, k3 = statistics::burr12_sf(-1.0, 2.0, 3.0);
    chk_val(bad, "burr12 sf below the support", a3, k3, 1.0, 1e-12);
    let a4, k4 = statistics::betaprime_pdf(-1.0, 3.0, 5.0);
    chk_val(bad, "betaprime pdf below the support", a4, k4, 0.0, 1e-12);
    let a5, k5 = statistics::genpareto_cdf(-1.0, 0.3);
    chk_val(bad, "genpareto cdf below the support", a5, k5, 0.0, 1e-12);
    // A negative shape bounds the support of the generalized Pareto law.
    let a6, k6 = statistics::genpareto_cdf(5.0, -0.5);
    chk_val(bad, "genpareto cdf above the support", a6, k6, 1.0, 1e-12);
    let a7, k7 = statistics::genpareto_pdf(5.0, -0.5);
    chk_val(bad, "genpareto pdf above the support", a7, k7, 0.0, 1e-12);
    // A shape of 0 gives the exponential law.
    let a8, k8 = statistics::genpareto_cdf(1.0, 0.0);
    chk_val(bad, "genpareto cdf shape 0", a8, k8, 1.0 - exp(-1.0), 1e-12);
    // The quantile of 0 gives the low bound.
    let e1, ke1 = statistics::fisk_ppf(0.0, 5.0);
    chk_val(bad, "fisk ppf 0", e1, ke1, 0.0, 1e-12);
    let e2, ke2 = statistics::betaprime_ppf(0.0, 3.0, 5.0);
    chk_val(bad, "betaprime ppf 0", e2, ke2, 0.0, 1e-12);
    let e3, ke3 = statistics::genpareto_ppf(1.0, -0.5);
    chk_val(bad, "genpareto ppf 1 of a bounded law", e3, ke3, 2.0, 1e-12);
    // The cumulative function of the quantile is the probability.
    let q1, kq1 = statistics::burr_ppf(0.3, 5.0, 3.0);
    chk_code(bad, "burr ppf code", kq1, 0);
    let r1, kr1 = statistics::burr_cdf(q1, 5.0, 3.0);
    chk_val(bad, "burr cdf of the quantile", r1, kr1, 0.3, 1e-9);
    let q2, kq2 = statistics::burr12_ppf(0.7, 2.0, 3.0);
    chk_code(bad, "burr12 ppf code", kq2, 0);
    let r2, kr2 = statistics::burr12_cdf(q2, 2.0, 3.0);
    chk_val(bad, "burr12 cdf of the quantile", r2, kr2, 0.7, 1e-9);
    let q3, kq3 = statistics::betaprime_ppf(0.9, 3.0, 5.0);
    chk_code(bad, "betaprime ppf code", kq3, 0);
    let r3, kr3 = statistics::betaprime_cdf(q3, 3.0, 5.0);
    chk_val(bad, "betaprime cdf of the quantile", r3, kr3, 0.9, 1e-9);
    // A shape of 0 or below gives the domain fault.
    let _z1, cz1 = statistics::fisk_pdf(1.0, 0.0);
    chk_code(bad, "fisk shape 0", cz1, 5);
    let _z2, cz2 = statistics::burr_pdf(1.0, 5.0, -1.0);
    chk_code(bad, "burr second shape below 0", cz2, 5);
    let _z3, cz3 = statistics::betaprime_cdf(1.0, 0.0, 5.0);
    chk_code(bad, "betaprime shape 0", cz3, 5);
    // A moment above the order of the tail does not exist.
    let _u1, cu1 = statistics::fisk_mean(1.0);
    chk_code(bad, "fisk mean of shape 1", cu1, 7);
    let _u2, cu2 = statistics::fisk_var(2.0);
    chk_code(bad, "fisk var of shape 2", cu2, 7);
    let _u3, cu3 = statistics::burr_kurt(4.0, 3.0);
    chk_code(bad, "burr kurt of shape 4", cu3, 7);
    let _u4, cu4 = statistics::betaprime_mean(3.0, 1.0);
    chk_code(bad, "betaprime mean of b 1", cu4, 7);
    let _u5, cu5 = statistics::genpareto_mean(1.0);
    chk_code(bad, "genpareto mean of shape 1", cu5, 7);
    let _u6, cu6 = statistics::genpareto_var(0.5);
    chk_code(bad, "genpareto var of shape 0.5", cu6, 7);
    // The shape of the generalized Pareto law may be negative or 0.
    let m1, km1 = statistics::genpareto_mean(0.0);
    chk_val(bad, "genpareto mean of shape 0", m1, km1, 1.0, 1e-12);
    return bad;
}
