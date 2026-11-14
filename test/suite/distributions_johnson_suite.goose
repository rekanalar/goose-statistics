// The rows of the families of the Johnson kind that the vector file does
// not hold: the support, the end points of the quantile and the shape
// bounds.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_johnson_tests() -> i64 {
    var bad: i64 = 0;
    // A value outside the support gives a density of 0.
    let a1, k1 = statistics::johnsonsb_pdf(-1.0, 1.0, 2.0);
    chk_val(bad, "johnsonsb pdf below the support", a1, k1, 0.0, 1e-12);
    let a2, k2 = statistics::johnsonsb_pdf(2.0, 1.0, 2.0);
    chk_val(bad, "johnsonsb pdf above the support", a2, k2, 0.0, 1e-12);
    let a3, k3 = statistics::johnsonsb_cdf(-1.0, 1.0, 2.0);
    chk_val(bad, "johnsonsb cdf below the support", a3, k3, 0.0, 1e-12);
    let a4, k4 = statistics::johnsonsb_sf(2.0, 1.0, 2.0);
    chk_val(bad, "johnsonsb sf above the support", a4, k4, 0.0, 1e-12);
    // The quantile of 0 gives the low bound and of 1 the high bound.
    let e1, ke1 = statistics::johnsonsb_ppf(0.0, 1.0, 2.0);
    chk_val(bad, "johnsonsb ppf 0", e1, ke1, 0.0, 1e-12);
    let e2, ke2 = statistics::johnsonsb_ppf(1.0, 1.0, 2.0);
    chk_val(bad, "johnsonsb ppf 1", e2, ke2, 1.0, 1e-12);
    // The cumulative function of the quantile is the probability.
    let q1, kq1 = statistics::johnsonsb_ppf(0.3, 1.0, 2.0);
    chk_code(bad, "johnsonsb ppf code", kq1, 0);
    let r1, kr1 = statistics::johnsonsb_cdf(q1, 1.0, 2.0);
    chk_val(bad, "johnsonsb cdf of the quantile", r1, kr1, 0.3, 1e-9);
    let q2, kq2 = statistics::johnsonsu_ppf(0.7, 1.0, 2.0);
    chk_code(bad, "johnsonsu ppf code", kq2, 0);
    let r2, kr2 = statistics::johnsonsu_cdf(q2, 1.0, 2.0);
    chk_val(bad, "johnsonsu cdf of the quantile", r2, kr2, 0.7, 1e-9);
    let q3, kq3 = statistics::jf_skew_t_ppf(0.9, 3.0, 2.0);
    chk_code(bad, "jf_skew_t ppf code", kq3, 0);
    let r3, kr3 = statistics::jf_skew_t_cdf(q3, 3.0, 2.0);
    chk_val(bad, "jf_skew_t cdf of the quantile", r3, kr3, 0.9, 1e-9);
    // A moment above the order of the tail does not exist.
    let _u1, cu1 = statistics::jf_skew_t_mean(0.5, 2.0);
    chk_code(bad, "jf_skew_t mean of shape a 0.5", cu1, 7);
    let _u2, cu2 = statistics::jf_skew_t_var(3.0, 1.0);
    chk_code(bad, "jf_skew_t var of shape b 1", cu2, 7);
    // A shape of 0 gives the domain fault.
    let _z1, cz1 = statistics::jf_skew_t_pdf(1.0, 0.0, 2.0);
    chk_code(bad, "jf_skew_t shape a 0", cz1, 5);
    let _z2, cz2 = statistics::johnsonsb_pdf(0.5, 1.0, 0.0);
    chk_code(bad, "johnsonsb shape b 0", cz2, 5);
    let _z3, cz3 = statistics::johnsonsu_pdf(1.0, 1.0, 0.0);
    chk_code(bad, "johnsonsu shape b 0", cz3, 5);
    return bad;
}
