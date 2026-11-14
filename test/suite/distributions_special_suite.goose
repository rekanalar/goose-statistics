// The rows of the families of a special shape that the vector file does not
// hold: the support, the end points of the quantile and the bounds of the
// moments.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_special_tests() -> i64 {
    var bad: i64 = 0;
    // A value outside the support gives a density of 0.
    let a1, k1 = statistics::invgauss_pdf(-1.0, 0.5);
    chk_val(bad, "invgauss pdf below the support", a1, k1, 0.0, 1e-12);
    let a2, k2 = statistics::nakagami_pdf(-1.0, 2.5);
    chk_val(bad, "nakagami pdf below the support", a2, k2, 0.0, 1e-12);
    let a3, k3 = statistics::kappa3_pdf(-1.0, 2.0);
    chk_val(bad, "kappa3 pdf below the support", a3, k3, 0.0, 1e-12);
    // The support of the inverse Gaussian law and of the Kappa law starts
    // at 0.
    let a4, k4 = statistics::invgauss_cdf(-1.0, 0.5);
    chk_val(bad, "invgauss cdf below the support", a4, k4, 0.0, 1e-12);
    let a5, k5 = statistics::invgauss_sf(-1.0, 0.5);
    chk_val(bad, "invgauss sf below the support", a5, k5, 1.0, 1e-12);
    let a6, k6 = statistics::kappa3_cdf(-1.0, 2.0);
    chk_val(bad, "kappa3 cdf below the support", a6, k6, 0.0, 1e-12);
    // The quantile of 0 gives the low bound and of 1 the high bound.
    let e1, ke1 = statistics::invgauss_ppf(0.0, 0.5);
    chk_val(bad, "invgauss ppf 0", e1, ke1, 0.0, 1e-12);
    let e2, ke2 = statistics::nakagami_ppf(0.0, 2.5);
    chk_val(bad, "nakagami ppf 0", e2, ke2, 0.0, 1e-12);
    let e3, ke3 = statistics::kappa3_ppf(0.0, 2.0);
    chk_val(bad, "kappa3 ppf 0", e3, ke3, 0.0, 1e-12);
    // The law of the Pearson kind has no high bound, so its quantile of 1
    // has no value.
    let e4, ke4 = statistics::pearson3_ppf(1.0, 1.5);
    chk_code(bad, "pearson3 ppf 1", ke4, 8);
    // The cumulative function of the quantile is the probability.
    let q1, kq1 = statistics::invgauss_ppf(0.3, 0.5);
    chk_code(bad, "invgauss ppf code", kq1, 0);
    let r1, kr1 = statistics::invgauss_cdf(q1, 0.5);
    chk_val(bad, "invgauss cdf of the quantile", r1, kr1, 0.3, 1e-9);
    let q2, kq2 = statistics::kappa3_ppf(0.3, 2.0);
    chk_code(bad, "kappa3 ppf code", kq2, 0);
    let r2, kr2 = statistics::kappa3_cdf(q2, 2.0);
    chk_val(bad, "kappa3 cdf of the quantile", r2, kr2, 0.3, 1e-9);
    // The moment of order n of the Kappa law of the third kind exists only
    // when the shape is above n.
    let m1, km1 = statistics::kappa3_mean(1.0);
    chk_code(bad, "kappa3 mean below the shape", km1, 7);
    let m2, km2 = statistics::kappa3_var(2.0);
    chk_code(bad, "kappa3 var at the shape", km2, 7);
    let m3, km3 = statistics::kappa3_skew(3.0);
    chk_code(bad, "kappa3 skew at the shape", km3, 7);
    let m4, km4 = statistics::kappa3_kurt(4.0);
    chk_code(bad, "kappa3 kurt at the shape", km4, 7);
    // The law of the Pearson kind is the law of the normal form when the
    // skew falls below the bound. Its moments do not fail there.
    let m9, km9 = statistics::pearson3_mean(0.0);
    chk_code(bad, "pearson3 mean at the skew 0", km9, 0);
    chk_val(bad, "the pearson3 mean at the skew 0", m9, km9, 0.0, 1e-12);
    let m12, km12 = statistics::pearson3_var(0.0);
    chk_val(bad, "the pearson3 var at the skew 0", m12, km12, 1.0, 1e-12);
    let m13, km13 = statistics::pearson3_pdf(1.0, 0.0);
    chk_val(bad, "the pearson3 pdf at the skew 0", m13, km13, 0.24197072451914337, 1e-12);
    let m14, km14 = statistics::pearson3_cdf(1.0, 0.0);
    chk_val(bad, "the pearson3 cdf at the skew 0", m14, km14, 0.84134474606854295, 1e-12);
    let m15, km15 = statistics::pearson3_ppf(0.5, 0.0);
    chk_val(bad, "the pearson3 ppf at the skew 0", m15, km15, 0.0, 1e-12);
    // The density of the Nakagami law at 0 depends on the shape.
    let m10, km10 = statistics::nakagami_pdf(0.0, 0.25);
    check(bad, km10 == 0 && is_inf(m10), "the nakagami pdf at 0 with a small shape");
    let m11, km11 = statistics::nakagami_pdf(0.0, 2.5);
    chk_val(bad, "the nakagami pdf at 0", m11, km11, 0.0, 1e-12);
    // The shape of 0.5 gives a density with a limit at 0.
    let m16, km16 = statistics::nakagami_pdf(0.0, 0.5);
    chk_val(bad, "the nakagami pdf at 0 with the shape 0.5", m16, km16, 0.79788456080286536, 1e-12);
    // The density of the Kappa law of the third kind has a limit at 0.
    let m17, km17 = statistics::kappa3_pdf(0.0, 2.0);
    chk_val(bad, "the kappa3 pdf at 0", m17, km17, 0.70710678118654752, 1e-12);
    // The quantile of the inverse Gaussian law in the deep tail. The answer
    // of mpmath 1.4.1, because the reference library loses the root there.
    let m18, km18 = statistics::invgauss_ppf(1e-300, 0.05);
    chk_val(bad, "the invgauss ppf in the deep tail", m18, km18, 0.00070743290331925498, 1e-9);
    let m19, km19 = statistics::invgauss_cdf(m18, 0.05);
    chk_val(bad, "the invgauss cdf of the deep quantile", m19, km19, 1e-300, 1e-9);
    return bad;
}
