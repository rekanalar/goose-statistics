// The rows of the histogram law that the vector file does not hold: the
// bounds of the content and the round trip of the quantile.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_histogram_tests() -> i64 {
    var bad: i64 = 0;
    let edges: f64[4] = [0.0, 1.0, 2.0, 3.0];
    let content: f64[3] = [1.0, 2.0, 1.0];
    // The cumulative function of the quantile gives the probability back.
    let q1, kq1 = statistics::histogram_ppf(content, edges, true, 0.35);
    let r1, kr1 = statistics::histogram_cdf(content, edges, true, q1);
    chk_val(bad, "the cdf of the quantile", r1, kq1 + kr1, 0.35, 1e-9);
    // A point under the first edge gives no density.
    let v1, c1 = statistics::histogram_pdf(content, edges, true, -1.0);
    chk_val(bad, "the density under the first edge", v1, c1, 0.0, 1e-12);
    // The last edge gives the whole mass.
    let v2, c2 = statistics::histogram_cdf(content, edges, true, 3.0);
    chk_val(bad, "the cdf at the last edge", v2, c2, 1.0, 1e-12);
    // A probability outside its range gives the domain fault.
    let _d1, cd1 = statistics::histogram_ppf(content, edges, true, 1.5);
    chk_code(bad, "a probability above one", cd1, 5);
    return bad;
}
