// The rows of the families of a symmetric law that the vector file does not
// hold: the support, the end points of the quantile and the shape bounds.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_symmetric_tests() -> i64 {
    var bad: i64 = 0;
    // A value outside the support gives a density of 0.
    let a1, k1 = statistics::halfgennorm_pdf(-1.0, 1.5);
    chk_val(bad, "halfgennorm pdf below the support", a1, k1, 0.0, 1e-12);
    let a2, k2 = statistics::halfgennorm_cdf(-1.0, 1.5);
    chk_val(bad, "halfgennorm cdf below the support", a2, k2, 0.0, 1e-12);
    let a3, k3 = statistics::halfgennorm_sf(-1.0, 1.5);
    chk_val(bad, "halfgennorm sf below the support", a3, k3, 1.0, 1e-12);
    // The law of the generalized normal kind is symmetric around 0.
    let a4, k4 = statistics::gennorm_cdf(0.0, 1.5);
    chk_val(bad, "gennorm cdf at the middle", a4, k4, 0.5, 1e-12);
    let a5, k5 = statistics::gennorm_mean(1.5);
    chk_val(bad, "gennorm mean", a5, k5, 0.0, 1e-12);
    let a6, k6 = statistics::gennorm_skew(1.5);
    chk_val(bad, "gennorm skew", a6, k6, 0.0, 1e-12);
    let a7, k7 = statistics::tukeylambda_mean(0.5);
    chk_val(bad, "tukeylambda mean", a7, k7, 0.0, 1e-12);
    let a8, k8 = statistics::tukeylambda_skew(0.5);
    chk_val(bad, "tukeylambda skew", a8, k8, 0.0, 1e-12);
    // The quantile of 0 gives the low bound and of 1 the high bound.
    let e1, ke1 = statistics::halfgennorm_ppf(0.0, 1.5);
    chk_val(bad, "halfgennorm ppf 0", e1, ke1, 0.0, 1e-12);
    let e2, ke2 = statistics::tukeylambda_ppf(0.0, 0.5);
    chk_val(bad, "tukeylambda ppf 0", e2, ke2, -2.0, 1e-12);
    let e3, ke3 = statistics::tukeylambda_ppf(1.0, 0.5);
    chk_val(bad, "tukeylambda ppf 1", e3, ke3, 2.0, 1e-12);
    // A shape above 0 bounds the law of Tukey and of lambda.
    let a9, k9 = statistics::tukeylambda_pdf(3.0, 0.5);
    chk_val(bad, "tukeylambda pdf above the support", a9, k9, 0.0, 1e-12);
    // The cumulative function of the quantile is the probability.
    let q1, kq1 = statistics::gennorm_ppf(0.3, 1.5);
    chk_code(bad, "gennorm ppf code", kq1, 0);
    let r1, kr1 = statistics::gennorm_cdf(q1, 1.5);
    chk_val(bad, "gennorm cdf of the quantile", r1, kr1, 0.3, 1e-9);
    let q2, kq2 = statistics::halfgennorm_ppf(0.7, 1.5);
    chk_code(bad, "halfgennorm ppf code", kq2, 0);
    let r2, kr2 = statistics::halfgennorm_cdf(q2, 1.5);
    chk_val(bad, "halfgennorm cdf of the quantile", r2, kr2, 0.7, 1e-9);
    let q3, kq3 = statistics::tukeylambda_ppf(0.9, 0.5);
    chk_code(bad, "tukeylambda ppf code", kq3, 0);
    let r3, kr3 = statistics::tukeylambda_cdf(q3, 0.5);
    chk_val(bad, "tukeylambda cdf of the quantile", r3, kr3, 0.9, 1e-9);
    // A shape of 0 gives the domain fault.
    let _z1, cz1 = statistics::gennorm_pdf(1.0, 0.0);
    chk_code(bad, "gennorm shape 0", cz1, 5);
    let _z2, cz2 = statistics::halfgennorm_pdf(1.0, 0.0);
    chk_code(bad, "halfgennorm shape 0", cz2, 5);
    let _z3, cz3 = statistics::tukeylambda_pdf(1.0, 0.0 / 0.0);
    chk_code(bad, "tukeylambda shape nan", cz3, 3);
    return bad;
}
