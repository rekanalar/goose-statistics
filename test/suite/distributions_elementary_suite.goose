// The rows of the elementary families that the vector file does not hold:
// the support of each law, the symmetry and the monotonicity of the
// cumulative function.
namespace suite;

import std;
import math;
import statistics from "../..";

// A value against an expected value within an absolute bound. The relative
// bound of chk_val does not work near 0.
fn fam_elem_abs(bad: i64&, name: u8[:], got: f64, c: i64, want: f64, tol: f64) {
    if c == 0 && abs(got - want) <= tol { return; }
    bad += 1;
    print("FAIL ", name, " got ", got, " code ", c, " want ", want, "\n");
}

fn run_distributions_elementary_tests() -> i64 {
    var bad: i64 = 0;
    // The density of a value outside the support is 0 and not a fault.
    let r1, c1 = statistics::anglit_pdf(-2.0);
    chk_val(bad, "anglit pdf below the support", r1, c1, 0.0, 1e-12);
    let r2, c2 = statistics::arcsine_pdf(1.5);
    chk_val(bad, "arcsine pdf above the support", r2, c2, 0.0, 1e-12);
    let r3, c3 = statistics::cosine_pdf(-4.0);
    chk_val(bad, "cosine pdf below the support", r3, c3, 0.0, 1e-12);
    let r4, c4 = statistics::gibrat_pdf(-1.0);
    chk_val(bad, "gibrat pdf below the support", r4, c4, 0.0, 1e-12);
    let r5, c5 = statistics::halfcauchy_pdf(-1.0);
    chk_val(bad, "halfcauchy pdf below the support", r5, c5, 0.0, 1e-12);
    let r6, c6 = statistics::halflogistic_pdf(-1.0);
    chk_val(bad, "halflogistic pdf below the support", r6, c6, 0.0, 1e-12);
    let r7, c7 = statistics::halfnorm_pdf(-1.0);
    chk_val(bad, "halfnorm pdf below the support", r7, c7, 0.0, 1e-12);
    // The cumulative function of a value below the support is 0.
    let r8, c8 = statistics::anglit_cdf(-2.0);
    chk_val(bad, "anglit cdf below the support", r8, c8, 0.0, 1e-12);
    let r9, c9 = statistics::arcsine_cdf(-0.5);
    chk_val(bad, "arcsine cdf below the support", r9, c9, 0.0, 1e-12);
    let r10, c10 = statistics::gibrat_cdf(-0.5);
    chk_val(bad, "gibrat cdf below the support", r10, c10, 0.0, 1e-12);
    let r11, c11 = statistics::halfcauchy_cdf(-0.5);
    chk_val(bad, "halfcauchy cdf below the support", r11, c11, 0.0, 1e-12);
    // The cumulative function of a value above the support is 1.
    let r12, c12 = statistics::anglit_cdf(2.0);
    chk_val(bad, "anglit cdf above the support", r12, c12, 1.0, 1e-12);
    let r13, c13 = statistics::arcsine_cdf(1.5);
    chk_val(bad, "arcsine cdf above the support", r13, c13, 1.0, 1e-12);
    let r14, c14 = statistics::gibrat_cdf(1e10);
    chk_val(bad, "gibrat cdf above the support", r14, c14, 1.0, 1e-12);
    // The survival function and the cumulative function add to 1.
    let c1, k1 = statistics::anglit_cdf(0.3);
    let s1, k2 = statistics::anglit_sf(0.3);
    chk_val(bad, "anglit cdf and sf add to 1", c1 + s1, k1 + k2, 1.0, 1e-12);
    chk_code(bad, "anglit cdf code", k1, 0);
    chk_code(bad, "anglit sf code", k2, 0);
    // The median of a symmetric law is the centre.
    let r15, c15 = statistics::anglit_ppf(0.5);
    fam_elem_abs(bad, "anglit median", r15, c15, 0.0, 1e-12);
    let r16, c16 = statistics::cosine_ppf(0.5);
    fam_elem_abs(bad, "cosine median", r16, c16, 0.0, 1e-9);
    let r17, c17 = statistics::hypsecant_ppf(0.5);
    fam_elem_abs(bad, "hypsecant median", r17, c17, 0.0, 1e-12);
    // The quantile of the end points gives the bounds of the support.
    let r18, c18 = statistics::anglit_ppf(0.0);
    chk_val(bad, "anglit ppf 0", r18, c18, -PI / 4.0, 1e-12);
    let r19, c19 = statistics::anglit_ppf(1.0);
    chk_val(bad, "anglit ppf 1", r19, c19, PI / 4.0, 1e-12);
    let r20, c20 = statistics::arcsine_ppf(0.0);
    chk_val(bad, "arcsine ppf 0", r20, c20, 0.0, 1e-12);
    let r21, c21 = statistics::arcsine_ppf(1.0);
    chk_val(bad, "arcsine ppf 1", r21, c21, 1.0, 1e-12);
    let r22, c22 = statistics::cosine_ppf(0.0);
    chk_val(bad, "cosine ppf 0", r22, c22, -PI, 1e-12);
    let r23, c23 = statistics::cosine_ppf(1.0);
    chk_val(bad, "cosine ppf 1", r23, c23, PI, 1e-12);
    let r24, c24 = statistics::halfcauchy_ppf(0.0);
    chk_val(bad, "halfcauchy ppf 0", r24, c24, 0.0, 1e-12);
    // The cumulative function of the quantile of a probability is that
    // probability.
    let q1, kq = statistics::arcsine_ppf(0.25);
    chk_code(bad, "arcsine ppf code", kq, 0);
    let r25, c25 = statistics::arcsine_cdf(q1);
    chk_val(bad, "arcsine cdf of the quantile", r25, c25, 0.25, 1e-9);
    let q2, kq2 = statistics::cosine_ppf(0.75);
    chk_code(bad, "cosine ppf code", kq2, 0);
    let r26, c26 = statistics::cosine_cdf(q2);
    chk_val(bad, "cosine cdf of the quantile", r26, c26, 0.75, 1e-9);
    let q3, kq3 = statistics::hypsecant_ppf(0.3);
    chk_code(bad, "hypsecant ppf code", kq3, 0);
    let r27, c27 = statistics::hypsecant_cdf(q3);
    chk_val(bad, "hypsecant cdf of the quantile", r27, c27, 0.3, 1e-9);
    // The density of a value rises to the mode and falls after it.
    let m1, km = statistics::halfnorm_pdf(0.0);
    let m2, km2 = statistics::halfnorm_pdf(1.0);
    chk_code(bad, "halfnorm mode code", km, 0);
    chk_code(bad, "halfnorm tail code", km2, 0);
    check(bad, m1 > m2, "the density of the half normal law falls with the value");
    // The mean of the half normal law is below 1 and its variance is below 1.
    let mu, kmu = statistics::halfnorm_mean();
    let va, kva = statistics::halfnorm_var();
    chk_code(bad, "halfnorm mean code", kmu, 0);
    chk_code(bad, "halfnorm var code", kva, 0);
    check(bad, mu > 0.7 && mu < 0.9, "the mean of the half normal law is near 0.8");
    check(bad, va > 0.3 && va < 0.4, "the variance of the half normal law is near 0.36");
    // The half Cauchy law has no moments.
    let _r28, c28 = statistics::halfcauchy_mean();
    chk_code(bad, "halfcauchy mean", c28, 7);
    let _r29, c29 = statistics::halfcauchy_var();
    chk_code(bad, "halfcauchy var", c29, 7);
    let _r30, c30 = statistics::halfcauchy_skew();
    chk_code(bad, "halfcauchy skew", c30, 7);
    let _r31, c31 = statistics::halfcauchy_kurt();
    chk_code(bad, "halfcauchy kurt", c31, 7);
    return bad;
}
