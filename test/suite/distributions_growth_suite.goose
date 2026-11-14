// The rows of the families of a growth or a logistic kind that the vector
// file does not hold: the support, the end points of the quantile and the
// shape bounds.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_growth_tests() -> i64 {
    var bad: i64 = 0;
    // A value below the support gives a density of 0.
    let a1, k1 = statistics::genhalflogistic_pdf(-1.0, 2.0);
    chk_val(bad, "genhalflogistic pdf below the support", a1, k1, 0.0, 1e-12);
    let a2, k2 = statistics::gompertz_pdf(-1.0, 1.5);
    chk_val(bad, "gompertz pdf below the support", a2, k2, 0.0, 1e-12);
    let a3, k3 = statistics::gompertz_cdf(-1.0, 1.5);
    chk_val(bad, "gompertz cdf below the support", a3, k3, 0.0, 1e-12);
    // The support of the generalized half-logistic law ends at 1 / c.
    let a4, k4 = statistics::genhalflogistic_cdf(10.0, 2.0);
    chk_val(bad, "genhalflogistic cdf above the support", a4, k4, 1.0, 1e-12);
    let e1, ke1 = statistics::genhalflogistic_ppf(1.0, 2.0);
    chk_val(bad, "genhalflogistic ppf 1", e1, ke1, 0.5, 1e-12);
    let e2, ke2 = statistics::gompertz_ppf(0.0, 1.5);
    chk_val(bad, "gompertz ppf 0", e2, ke2, 0.0, 1e-12);
    // The cumulative function of the quantile is the probability.
    let q1, kq1 = statistics::genlogistic_ppf(0.3, 2.5);
    chk_code(bad, "genlogistic ppf code", kq1, 0);
    let r1, kr1 = statistics::genlogistic_cdf(q1, 2.5);
    chk_val(bad, "genlogistic cdf of the quantile", r1, kr1, 0.3, 1e-9);
    let q2, kq2 = statistics::gompertz_ppf(0.8, 1.5);
    chk_code(bad, "gompertz ppf code", kq2, 0);
    let r2, kr2 = statistics::gompertz_cdf(q2, 1.5);
    chk_val(bad, "gompertz cdf of the quantile", r2, kr2, 0.8, 1e-9);
    let q3, kq3 = statistics::genhalflogistic_ppf(0.6, 2.0);
    chk_code(bad, "genhalflogistic ppf code", kq3, 0);
    let r3, kr3 = statistics::genhalflogistic_cdf(q3, 2.0);
    chk_val(bad, "genhalflogistic cdf of the quantile", r3, kr3, 0.6, 1e-9);
    // A shape of 0 or below gives the domain fault.
    let _z1, cz1 = statistics::genlogistic_pdf(1.0, 0.0);
    chk_code(bad, "genlogistic shape 0", cz1, 5);
    let _z2, cz2 = statistics::genhalflogistic_cdf(0.1, -1.0);
    chk_code(bad, "genhalflogistic shape below 0", cz2, 5);
    let _z3, cz3 = statistics::gompertz_sf(1.0, 0.0);
    chk_code(bad, "gompertz shape 0", cz3, 5);
    // The law of a shape of 1 of the logistic kind is the logistic law.
    let l1, kl1 = statistics::genlogistic_mean(1.0);
    chk_val(bad, "genlogistic mean of shape 1", l1, kl1, 0.0, 1e-12);
    return bad;
}
