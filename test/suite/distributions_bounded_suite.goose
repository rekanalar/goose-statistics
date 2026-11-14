// The rows of the bounded families that the vector file does not hold: the
// support of each law, the shape bounds and the median of a symmetric law.
namespace suite;

import std;
import math;
import statistics from "../..";

// A value against an expected value within an absolute bound. The relative
// bound of chk_val does not work near 0.
fn fam_bound_abs(bad: i64&, name: u8[:], got: f64, c: i64, want: f64, tol: f64) {
    if c == 0 && abs(got - want) <= tol { return; }
    bad += 1;
    print("FAIL ", name, " got ", got, " code ", c, " want ", want, "\n");
}

fn run_distributions_bounded_tests() -> i64 {
    var bad: i64 = 0;
    // A value outside the support gives a density of 0 and not a fault.
    let v1, c1 = statistics::semicircular_pdf(1.5);
    chk_val(bad, "semicircular pdf above the support", v1, c1, 0.0, 1e-12);
    let v2, c2 = statistics::triang_pdf(1.5, 0.3);
    chk_val(bad, "triang pdf above the support", v2, c2, 0.0, 1e-12);
    let v3, c3 = statistics::trapezoid_pdf(-0.5, 0.2, 0.9);
    chk_val(bad, "trapezoid pdf below the support", v3, c3, 0.0, 1e-12);
    let v4, c4 = statistics::bradford_pdf(2.0, 1.5);
    chk_val(bad, "bradford pdf above the support", v4, c4, 0.0, 1e-12);
    let v5, c5 = statistics::wald_pdf(-1.0);
    chk_val(bad, "wald pdf below the support", v5, c5, 0.0, 1e-12);
    let v6, c6 = statistics::maxwell_pdf(-1.0);
    chk_val(bad, "maxwell pdf below the support", v6, c6, 0.0, 1e-12);
    // The cumulative function of a value below the support is 0 and above
    // the support is 1.
    let v7, c7 = statistics::trapezoid_cdf(-1.0, 0.2, 0.9);
    chk_val(bad, "trapezoid cdf below the support", v7, c7, 0.0, 1e-12);
    let v8, c8 = statistics::bradford_cdf(2.0, 1.5);
    chk_val(bad, "bradford cdf above the support", v8, c8, 1.0, 1e-12);
    let v9, c9 = statistics::wald_cdf(-1.0);
    chk_val(bad, "wald cdf below the support", v9, c9, 0.0, 1e-12);
    // The median of a symmetric law is the centre.
    let m1, cm1 = statistics::semicircular_ppf(0.5);
    fam_bound_abs(bad, "semicircular median", m1, cm1, 0.0, 1e-9);
    let m2, cm2 = statistics::moyal_ppf(0.5);
    check(bad, cm2 == 0 && m2 > 0.0 && m2 < 2.0, "the median of the moyal law is near 1");
    // The cumulative function of the quantile of a probability is that
    // probability.
    let q1, kq1 = statistics::trapezoid_ppf(0.4, 0.2, 0.9);
    chk_code(bad, "trapezoid ppf code", kq1, 0);
    let r1, kr1 = statistics::trapezoid_cdf(q1, 0.2, 0.9);
    chk_val(bad, "trapezoid cdf of the quantile", r1, kr1, 0.4, 1e-9);
    let q2, kq2 = statistics::bradford_ppf(0.3, 1.5);
    chk_code(bad, "bradford ppf code", kq2, 0);
    let r2, kr2 = statistics::bradford_cdf(q2, 1.5);
    chk_val(bad, "bradford cdf of the quantile", r2, kr2, 0.3, 1e-9);
    let q3, kq3 = statistics::wald_ppf(0.25);
    chk_code(bad, "wald ppf code", kq3, 0);
    let r3, kr3 = statistics::wald_cdf(q3);
    chk_val(bad, "wald cdf of the quantile", r3, kr3, 0.25, 1e-8);
    // The quantile of the end points gives the bounds of the support.
    let e1, ce1 = statistics::semicircular_ppf(0.0);
    chk_val(bad, "semicircular ppf 0", e1, ce1, -1.0, 1e-12);
    let e2, ce2 = statistics::triang_ppf(1.0, 0.3);
    chk_val(bad, "triang ppf 1", e2, ce2, 1.0, 1e-12);
    let e3, ce3 = statistics::trapezoid_ppf(0.0, 0.2, 0.9);
    chk_val(bad, "trapezoid ppf 0", e3, ce3, 0.0, 1e-12);
    // A shape of a triangular law outside 0 to 1 gives the domain fault.
    let _z1, cz1 = statistics::triang_pdf(0.5, 1.5);
    chk_code(bad, "triang shape above 1", cz1, 5);
    let _z2, cz2 = statistics::triang_pdf(0.5, -0.5);
    chk_code(bad, "triang shape below 0", cz2, 5);
    let _z3, cz3 = statistics::trapezoid_pdf(0.5, 1.5, 0.9);
    chk_code(bad, "trapezoid shape above 1", cz3, 5);
    let _z4, cz4 = statistics::trapezoid_pdf(0.5, 0.9, 0.2);
    chk_code(bad, "trapezoid shape below d", cz4, 5);
    let _z5, cz5 = statistics::bradford_pdf(0.5, 0.0);
    chk_code(bad, "bradford shape 0", cz5, 5);
    // The moments of the trapezoid law come from the three pieces.
    let t1, kt1 = statistics::trapezoid_mean(0.2, 0.9);
    let t2, kt2 = statistics::trapezoid_var(0.2, 0.9);
    chk_code(bad, "trapezoid mean code", kt1, 0);
    chk_code(bad, "trapezoid var code", kt2, 0);
    check(bad, t1 > 0.52 && t1 < 0.53, "the mean of the trapezoid law is near 0.524");
    check(bad, t2 > 0.062 && t2 < 0.063, "the variance of the trapezoid law is near 0.0623");
    // The moments of the moyal law are the ones of SciPy.
    let y1, ky1 = statistics::moyal_mean();
    let y2, ky2 = statistics::moyal_var();
    chk_code(bad, "moyal mean code", ky1, 0);
    chk_code(bad, "moyal var code", ky2, 0);
    check(bad, y1 > 1.27 && y1 < 1.28, "the mean of the moyal law is near 1.2704");
    check(bad, y2 > 4.93 && y2 < 4.94, "the variance of the moyal law is near 4.9348");
    return bad;
}
