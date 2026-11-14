// The rows of the Landau law that the vector file does not hold: the
// round trip of the quantile and the left tail.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_landau_tests() -> i64 {
    var bad: i64 = 0;
    let q0, kq0 = statistics::landau_ppf(0.35);
    let r0, kr0 = statistics::landau_cdf(q0);
    chk_val(bad, "landau cdf of the quantile", r0, kq0 + kr0, 0.35, 1e-7);
    let q1, kq1 = statistics::landau_ppf(0.8);
    let r1, kr1 = statistics::landau_cdf(q1);
    chk_val(bad, "landau cdf of the quantile", r1, kq1 + kr1, 0.8, 1e-7);
    let v1, c1 = statistics::landau_pdf(-20.0);
    chk_val(bad, "landau pdf in the far left tail", v1, c1, 0.0, 1e-12);
    let v2, c2 = statistics::landau_cdf(-20.0);
    chk_val(bad, "landau cdf in the far left tail", v2, c2, 0.0, 1e-12);
    return bad;
}
