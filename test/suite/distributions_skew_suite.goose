// The rows of the families of a skew shape that the vector file does not
// hold: the support, the end points of the quantile and the bounds of the
// moments.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_skew_tests() -> i64 {
    var bad: i64 = 0;
    // The law of the skew normal form is on the whole line, so its density
    // never falls to 0 outside a support.
    let a1, k1 = statistics::skewnorm_pdf(0.0, 2.0);
    chk_val(bad, "the skewnorm pdf at 0", a1, k1, 0.39894228040143270, 1e-12);
    // The shape of 0 gives the law of the normal form.
    let a2, k2 = statistics::skewnorm_cdf(0.0, 0.0);
    chk_val(bad, "the skewnorm cdf at the shape 0", a2, k2, 0.5, 1e-12);
    let a3, k3 = statistics::skewnorm_mean(0.0);
    chk_val(bad, "the skewnorm mean at the shape 0", a3, k3, 0.0, 1e-12);
    // The quantile of 0 and of 1 have no value, because the law has no
    // bound.
    let e1, ke1 = statistics::skewnorm_ppf(0.0, 2.0);
    chk_code(bad, "skewnorm ppf 0", ke1, 8);
    let e2, ke2 = statistics::laplace_asymmetric_ppf(1.0, 0.5);
    chk_code(bad, "laplace_asymmetric ppf 1", ke2, 8);
    // The cumulative function of the quantile is the probability.
    let q1, kq1 = statistics::skewnorm_ppf(0.3, 2.0);
    chk_code(bad, "skewnorm ppf code", kq1, 0);
    let r1, kr1 = statistics::skewnorm_cdf(q1, 2.0);
    chk_val(bad, "skewnorm cdf of the quantile", r1, kr1, 0.3, 1e-9);
    let q2, kq2 = statistics::laplace_asymmetric_ppf(0.3, 0.5);
    chk_code(bad, "laplace_asymmetric ppf code", kq2, 0);
    let r2, kr2 = statistics::laplace_asymmetric_cdf(q2, 0.5);
    chk_val(bad, "laplace_asymmetric cdf of the quantile", r2, kr2, 0.3, 1e-9);
    // The law of the skew normal form is symmetric in its shape.
    let s1, ks1 = statistics::skewnorm_cdf(1.3, 2.0);
    let s2, ks2 = statistics::skewnorm_sf(-1.3, -2.0);
    chk_val(bad, "the skewnorm sf of the opposite shape", s2, ks2, s1, 1e-12);
    // The shape of the asymmetric Laplace law must be above 0.
    let a4, k4 = statistics::laplace_asymmetric_mean(0.0);
    chk_code(bad, "laplace_asymmetric mean at the shape 0", k4, 5);
    // The deep tail keeps its digits. The answers of mpmath 1.4.1.
    let t1, kt1 = statistics::skewnorm_cdf(-3.0, 2.0);
    chk_val(bad, "the skewnorm cdf in the left tail", t1, kt1, 5.596622287836366e-13, 1e-10);
    let t2, kt2 = statistics::skewnorm_ppf(1e-300, 0.5);
    chk_code(bad, "the skewnorm ppf in the deep tail code", kt2, 0);
    let t2c, kt2c = statistics::skewnorm_cdf(t2, 0.5);
    chk_val(bad, "the skewnorm cdf of the deep quantile", t2c, kt2c, 1e-300, 1e-9);
    let t3, kt3 = statistics::skewnorm_ppf(1e-12, 2.0);
    chk_val(bad, "the skewnorm ppf at 1e-12", t3, kt3, -2.962638342439212, 1e-9);
    let t4, kt4 = statistics::skewnorm_ppf(1.0 - 1e-12, 0.5);
    chk_val(bad, "the skewnorm ppf in the right tail", t4, kt4, 7.130490070377563, 1e-9);
    return bad;
}
