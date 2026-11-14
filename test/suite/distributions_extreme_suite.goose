// The rows of the generalised extreme value law that the vector file does
// not hold: the support, the end points of the quantile and the bounds of
// the moments.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_extreme_tests() -> i64 {
    var bad: i64 = 0;
    // A value outside the support gives a density of 0.
    let a1, k1 = statistics::genextreme_pdf(3.0, 0.5);
    chk_val(bad, "genextreme pdf above the support", a1, k1, 0.0, 1e-12);
    let a2, k2 = statistics::genextreme_pdf(-3.0, -0.5);
    chk_val(bad, "genextreme pdf below the support", a2, k2, 0.0, 1e-12);
    // The cumulative function of the quantile is the probability.
    let q1, kq1 = statistics::genextreme_ppf(0.3, 0.5);
    chk_code(bad, "genextreme ppf code", kq1, 0);
    let r1, kr1 = statistics::genextreme_cdf(q1, 0.5);
    chk_val(bad, "genextreme cdf of the quantile", r1, kr1, 0.3, 1e-9);
    // The moments of the law need a shape above a bound.
    let m5, km5 = statistics::genextreme_mean(-1.5);
    chk_code(bad, "genextreme mean below the bound", km5, 7);
    let m6, km6 = statistics::genextreme_var(-0.6);
    chk_code(bad, "genextreme var below the bound", km6, 7);
    let m7, km7 = statistics::genextreme_skew(-0.4);
    chk_code(bad, "genextreme skew below the bound", km7, 7);
    let m8, km8 = statistics::genextreme_kurt(-0.3);
    chk_code(bad, "genextreme kurt below the bound", km8, 7);
    // The shape of 0 gives the law of Gumbel.
    let m12, km12 = statistics::genextreme_cdf(0.0, 0.0);
    chk_val(bad, "the genextreme cdf at the shape 0", m12, km12, 0.36787944117144233, 1e-12);
    // The density of the Gumbel law at the shape 0. The answer of mpmath 1.4.1.
    let m13, km13 = statistics::genextreme_pdf(1.0, 0.0);
    chk_val(bad, "the genextreme pdf at the shape 0", m13, km13, 0.2546463800435825, 1e-12);
    // The skewness and the kurtosis of a small shape take the Gumbel value.
    let m14, km14 = statistics::genextreme_skew(1e-5);
    chk_val(bad, "the genextreme skew at a small shape", m14, km14, 1.1395470994046486, 1e-12);
    let m15, km15 = statistics::genextreme_kurt(1e-4);
    chk_val(bad, "the genextreme kurt at a small shape", m15, km15, 2.4, 1e-12);
    return bad;
}
