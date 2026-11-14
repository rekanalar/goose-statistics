// The rows of the families that come from the gamma law that the vector
// file does not hold: the support, the end points of the quantile, the
// shape bounds and the moments that do not exist.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_gamma_tests() -> i64 {
    var bad: i64 = 0;
    // A value outside the support gives a density of 0.
    let a1, k1 = statistics::chi_pdf(-1.0, 3.0);
    chk_val(bad, "chi pdf below the support", a1, k1, 0.0, 1e-12);
    let a2, k2 = statistics::chi_cdf(-1.0, 3.0);
    chk_val(bad, "chi cdf below the support", a2, k2, 0.0, 1e-12);
    let a3, k3 = statistics::invgamma_pdf(-1.0, 3.0);
    chk_val(bad, "invgamma pdf below the support", a3, k3, 0.0, 1e-12);
    let a4, k4 = statistics::invgamma_sf(-1.0, 3.0);
    chk_val(bad, "invgamma sf below the support", a4, k4, 1.0, 1e-12);
    // The law of the double gamma law is symmetric around 0.
    let a5, k5 = statistics::dgamma_cdf(0.0, 2.0);
    chk_val(bad, "dgamma cdf at the middle", a5, k5, 0.5, 1e-12);
    let a6, k6 = statistics::dgamma_mean(2.0);
    chk_val(bad, "dgamma mean", a6, k6, 0.0, 1e-12);
    let a7, k7 = statistics::dgamma_skew(2.0);
    chk_val(bad, "dgamma skew", a7, k7, 0.0, 1e-12);
    // The quantile of 0 gives the low bound and of 1 the high bound.
    let e1, ke1 = statistics::chi_ppf(0.0, 3.0);
    chk_val(bad, "chi ppf 0", e1, ke1, 0.0, 1e-12);
    let e2, ke2 = statistics::invgamma_ppf(1.0, 3.0);
    chk_code(bad, "invgamma ppf 1 of an unbounded law", ke2, 8);
    // The cumulative function of the quantile is the probability.
    let q1, kq1 = statistics::chi_ppf(0.3, 3.0);
    chk_code(bad, "chi ppf code", kq1, 0);
    let r1, kr1 = statistics::chi_cdf(q1, 3.0);
    chk_val(bad, "chi cdf of the quantile", r1, kr1, 0.3, 1e-9);
    let q2, kq2 = statistics::invgamma_ppf(0.7, 3.0);
    chk_code(bad, "invgamma ppf code", kq2, 0);
    let r2, kr2 = statistics::invgamma_cdf(q2, 3.0);
    chk_val(bad, "invgamma cdf of the quantile", r2, kr2, 0.7, 1e-9);
    let q3, kq3 = statistics::loggamma_ppf(0.5, 2.0);
    chk_code(bad, "loggamma ppf code", kq3, 0);
    let r3, kr3 = statistics::loggamma_cdf(q3, 2.0);
    chk_val(bad, "loggamma cdf of the quantile", r3, kr3, 0.5, 1e-9);
    // A moment above the order of the tail does not exist.
    let _u1, cu1 = statistics::invgamma_mean(1.0);
    chk_code(bad, "invgamma mean of shape 1", cu1, 7);
    let _u2, cu2 = statistics::invgamma_var(2.0);
    chk_code(bad, "invgamma var of shape 2", cu2, 7);
    let _u3, cu3 = statistics::invgamma_skew(3.0);
    chk_code(bad, "invgamma skew of shape 3", cu3, 7);
    let _u4, cu4 = statistics::invgamma_kurt(4.0);
    chk_code(bad, "invgamma kurt of shape 4", cu4, 7);
    // A shape of 0 gives the domain fault.
    let _z1, cz1 = statistics::chi_pdf(1.0, 0.0);
    chk_code(bad, "chi shape 0", cz1, 5);
    let _z2, cz2 = statistics::dgamma_pdf(1.0, 0.0);
    chk_code(bad, "dgamma shape 0", cz2, 5);
    let _z3, cz3 = statistics::invgamma_pdf(1.0, 0.0);
    chk_code(bad, "invgamma shape 0", cz3, 5);
    let _z4, cz4 = statistics::loggamma_pdf(1.0, 0.0);
    chk_code(bad, "loggamma shape 0", cz4, 5);
    return bad;
}
