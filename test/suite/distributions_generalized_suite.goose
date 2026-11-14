// The rows of the families of a generalized law that the vector file does
// not hold: the support, the end points of the quantile, the shape bounds
// and the moments that do not exist.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_generalized_tests() -> i64 {
    var bad: i64 = 0;
    // A value outside the support gives a density of 0.
    let a1, k1 = statistics::gengamma_pdf(-1.0, 3.0, 2.0);
    chk_val(bad, "gengamma pdf below the support", a1, k1, 0.0, 1e-12);
    let a2, k2 = statistics::gengamma_cdf(-1.0, 3.0, 2.0);
    chk_val(bad, "gengamma cdf below the support", a2, k2, 0.0, 1e-12);
    let a3, k3 = statistics::genexpon_pdf(-1.0, 1.0, 2.0, 3.0);
    chk_val(bad, "genexpon pdf below the support", a3, k3, 0.0, 1e-12);
    let a4, k4 = statistics::genexpon_sf(-1.0, 1.0, 2.0, 3.0);
    chk_val(bad, "genexpon sf below the support", a4, k4, 1.0, 1e-12);
    let a5, k5 = statistics::mielke_pdf(-1.0, 2.0, 3.0);
    chk_val(bad, "mielke pdf below the support", a5, k5, 0.0, 1e-12);
    // A shape below 0 gives the law of a bounded support above 0.
    let a6, k6 = statistics::gengamma_cdf(1.0, 3.0, -2.0);
    chk_code(bad, "gengamma cdf of a negative shape", k6, 0);
    check(bad, a6 >= 0.0 && a6 <= 1.0, "the cumulative of the generalized gamma law is a probability");
    // The quantile of 0 gives the low bound.
    let e1, ke1 = statistics::gengamma_ppf(0.0, 3.0, 2.0);
    chk_val(bad, "gengamma ppf 0", e1, ke1, 0.0, 1e-12);
    let e2, ke2 = statistics::genexpon_ppf(0.0, 1.0, 2.0, 3.0);
    chk_val(bad, "genexpon ppf 0", e2, ke2, 0.0, 1e-12);
    let e3, ke3 = statistics::mielke_ppf(0.0, 2.0, 3.0);
    chk_val(bad, "mielke ppf 0", e3, ke3, 0.0, 1e-12);
    // The cumulative function of the quantile is the probability.
    let q1, kq1 = statistics::gengamma_ppf(0.3, 3.0, 2.0);
    chk_code(bad, "gengamma ppf code", kq1, 0);
    let r1, kr1 = statistics::gengamma_cdf(q1, 3.0, 2.0);
    chk_val(bad, "gengamma cdf of the quantile", r1, kr1, 0.3, 1e-9);
    let q2, kq2 = statistics::genexpon_ppf(0.7, 1.0, 2.0, 3.0);
    chk_code(bad, "genexpon ppf code", kq2, 0);
    let r2, kr2 = statistics::genexpon_cdf(q2, 1.0, 2.0, 3.0);
    chk_val(bad, "genexpon cdf of the quantile", r2, kr2, 0.7, 1e-9);
    let q3, kq3 = statistics::mielke_ppf(0.9, 2.0, 3.0);
    chk_code(bad, "mielke ppf code", kq3, 0);
    let r3, kr3 = statistics::mielke_cdf(q3, 2.0, 3.0);
    chk_val(bad, "mielke cdf of the quantile", r3, kr3, 0.9, 1e-9);
    // A moment above the order of the tail does not exist.
    let _u1, cu1 = statistics::mielke_mean(2.0, 1.0);
    chk_code(bad, "mielke mean of shape s 1", cu1, 7);
    let _u2, cu2 = statistics::mielke_var(2.0, 2.0);
    chk_code(bad, "mielke var of shape s 2", cu2, 7);
    let _u3, cu3 = statistics::gengamma_mean(1.0, -1.0);
    chk_code(bad, "gengamma mean of a negative shape", cu3, 7);
    // A shape of 0 gives the domain fault.
    let _z1, cz1 = statistics::gengamma_pdf(1.0, 0.0, 2.0);
    chk_code(bad, "gengamma shape 0", cz1, 5);
    let _z2, cz2 = statistics::gengamma_pdf(1.0, 3.0, 0.0);
    chk_code(bad, "gengamma shape c 0", cz2, 5);
    let _z3, cz3 = statistics::genexpon_pdf(1.0, 0.0, 2.0, 3.0);
    chk_code(bad, "genexpon shape 0", cz3, 5);
    let _z4, cz4 = statistics::mielke_pdf(1.0, 0.0, 3.0);
    chk_code(bad, "mielke shape 0", cz4, 5);
    return bad;
}
