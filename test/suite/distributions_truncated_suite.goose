// The rows of the families of a truncation that the vector file does not
// hold: the support, the end points of the quantile, the bounds and the
// moments.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_truncated_tests() -> i64 {
    var bad: i64 = 0;
    // A value outside the support gives a density of 0.
    let a1, k1 = statistics::truncexpon_pdf(-1.0, 2.0);
    chk_val(bad, "truncexpon pdf below the support", a1, k1, 0.0, 1e-12);
    let a2, k2 = statistics::truncexpon_pdf(3.0, 2.0);
    chk_val(bad, "truncexpon pdf above the support", a2, k2, 0.0, 1e-12);
    let a3, k3 = statistics::truncexpon_cdf(-1.0, 2.0);
    chk_val(bad, "truncexpon cdf below the support", a3, k3, 0.0, 1e-12);
    let a4, k4 = statistics::truncexpon_sf(3.0, 2.0);
    chk_val(bad, "truncexpon sf above the support", a4, k4, 0.0, 1e-12);
    let a5, k5 = statistics::truncpareto_pdf(0.5, 2.0, 5.0);
    chk_val(bad, "truncpareto pdf below the support", a5, k5, 0.0, 1e-12);
    let a6, k6 = statistics::truncpareto_cdf(7.0, 2.0, 5.0);
    chk_val(bad, "truncpareto cdf above the support", a6, k6, 1.0, 1e-12);
    let a7, k7 = statistics::truncnorm_pdf(-3.0, -1.0, 2.0);
    chk_val(bad, "truncnorm pdf below the support", a7, k7, 0.0, 1e-12);
    let a8, k8 = statistics::truncnorm_sf(3.0, -1.0, 2.0);
    chk_val(bad, "truncnorm sf above the support", a8, k8, 0.0, 1e-12);
    // The quantile of 0 gives the low bound and of 1 the high bound.
    let e1, ke1 = statistics::truncexpon_ppf(0.0, 2.0);
    chk_val(bad, "truncexpon ppf 0", e1, ke1, 0.0, 1e-12);
    let e2, ke2 = statistics::truncexpon_ppf(1.0, 2.0);
    chk_val(bad, "truncexpon ppf 1", e2, ke2, 2.0, 1e-12);
    let e3, ke3 = statistics::truncpareto_ppf(0.0, 2.0, 5.0);
    chk_val(bad, "truncpareto ppf 0", e3, ke3, 1.0, 1e-12);
    let e4, ke4 = statistics::truncpareto_ppf(1.0, 2.0, 5.0);
    chk_val(bad, "truncpareto ppf 1", e4, ke4, 5.0, 1e-12);
    let e5, ke5 = statistics::truncnorm_ppf(0.0, -1.0, 2.0);
    chk_val(bad, "truncnorm ppf 0", e5, ke5, -1.0, 1e-12);
    let e6, ke6 = statistics::truncnorm_ppf(1.0, -1.0, 2.0);
    chk_val(bad, "truncnorm ppf 1", e6, ke6, 2.0, 1e-12);
    // The cumulative function of the quantile is the probability.
    let q1, kq1 = statistics::truncexpon_ppf(0.3, 2.0);
    chk_code(bad, "truncexpon ppf code", kq1, 0);
    let r1, kr1 = statistics::truncexpon_cdf(q1, 2.0);
    chk_val(bad, "truncexpon cdf of the quantile", r1, kr1, 0.3, 1e-9);
    let q2, kq2 = statistics::truncpareto_ppf(0.7, 2.0, 5.0);
    chk_code(bad, "truncpareto ppf code", kq2, 0);
    let r2, kr2 = statistics::truncpareto_cdf(q2, 2.0, 5.0);
    chk_val(bad, "truncpareto cdf of the quantile", r2, kr2, 0.7, 1e-9);
    let q3, kq3 = statistics::truncnorm_ppf(0.9, -1.0, 2.0);
    chk_code(bad, "truncnorm ppf code", kq3, 0);
    let r3, kr3 = statistics::truncnorm_cdf(q3, -1.0, 2.0);
    chk_val(bad, "truncnorm cdf of the quantile", r3, kr3, 0.9, 1e-9);
    // The support of the truncated Pareto law is bounded, so every moment
    // of that law exists.
    let m1, cm1 = statistics::truncpareto_mean(2.0, 5.0);
    chk_code(bad, "truncpareto mean code", cm1, 0);
    check(bad, !is_nan(m1), "the mean of the truncated Pareto law is a number");
    let m2, cm2 = statistics::truncpareto_kurt(2.0, 5.0);
    chk_code(bad, "truncpareto kurt code", cm2, 0);
    check(bad, !is_nan(m2), "the kurtosis of the truncated Pareto law is a number");
    // A bound outside its range gives the domain fault.
    let _z1, cz1 = statistics::truncexpon_pdf(1.0, 0.0);
    chk_code(bad, "truncexpon shape 0", cz1, 5);
    let _z2, cz2 = statistics::truncpareto_pdf(1.5, 0.0, 5.0);
    chk_code(bad, "truncpareto shape 0", cz2, 5);
    let _z3, cz3 = statistics::truncpareto_cdf(1.5, 2.0, 0.5);
    chk_code(bad, "truncpareto end below 1", cz3, 5);
    let _z4, cz4 = statistics::truncnorm_pdf(0.0, 2.0, -1.0);
    chk_code(bad, "truncnorm bounds in the wrong order", cz4, 5);
    return bad;
}
