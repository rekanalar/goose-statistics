// The rows of the five laws of the third search group that the vector
// file does not hold: the support of each law, the bounds of the shapes
// and the round trip of the quantile.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_search_c_tests() -> i64 {
    var bad: i64 = 0;
    // The cumulative function of the quantile gives the probability back.
    let q1, kq1 = statistics::norminvgauss_ppf(0.4, 1.0, 0.5);
    let r1, kr1 = statistics::norminvgauss_cdf(q1, 1.0, 0.5);
    chk_val(bad, "norminvgauss cdf of the quantile", r1, kq1 + kr1, 0.4, 1e-9);
    let q2, kq2 = statistics::recipinvgauss_ppf(0.6, 2.0);
    let r2, kr2 = statistics::recipinvgauss_cdf(q2, 2.0);
    chk_val(bad, "recipinvgauss cdf of the quantile", r2, kq2 + kr2, 0.6, 1e-9);
    let q3, kq3 = statistics::rel_breitwigner_ppf(0.3, 0.5);
    let r3, kr3 = statistics::rel_breitwigner_cdf(q3, 0.5);
    chk_val(bad, "rel_breitwigner cdf of the quantile", r3, kq3 + kr3, 0.3, 1e-9);
    let q4, kq4 = statistics::studentized_range_ppf(0.7, 3.0, 10.0);
    let r4, kr4 = statistics::studentized_range_cdf(q4, 3.0, 10.0);
    chk_val(bad, "studentized range cdf of the quantile", r4, kq4 + kr4, 0.7, 1e-8);
    let q5, kq5 = statistics::vonmises_ppf(0.8, 2.0);
    let r5, kr5 = statistics::vonmises_cdf(q5, 2.0);
    chk_val(bad, "vonmises cdf of the quantile", r5, kq5 + kr5, 0.8, 1e-9);
    // The quantile of the end points gives the bounds of the law.
    let e1, ce1 = statistics::vonmises_ppf(0.0, 2.0);
    chk_val(bad, "vonmises ppf 0", e1, ce1, (-1.0 / 0.0), 1e-12);
    let e2, ce2 = statistics::vonmises_ppf(1.0, 2.0);
    chk_val(bad, "vonmises ppf 1", e2, ce2, (1.0 / 0.0), 1e-12);
    // A value below the support gives a density of 0.
    let v1, c1 = statistics::recipinvgauss_pdf(-1.0, 2.0);
    chk_val(bad, "recipinvgauss pdf below the support", v1, c1, 0.0, 1e-12);
    let v2, c2 = statistics::recipinvgauss_cdf(-1.0, 2.0);
    chk_val(bad, "recipinvgauss cdf below the support", v2, c2, 0.0, 1e-12);
    // A shape outside its range gives the domain fault.
    let _d1, cd1 = statistics::norminvgauss_pdf(1.0, 1.0, 1.0);
    chk_code(bad, "norminvgauss skewness at the shape", cd1, 5);
    let _d2, cd2 = statistics::recipinvgauss_pdf(1.0, 0.0);
    chk_code(bad, "recipinvgauss shape 0", cd2, 5);
    let _d3, cd3 = statistics::rel_breitwigner_pdf(1.0, 0.0);
    chk_code(bad, "rel_breitwigner shape 0", cd3, 5);
    let _d4, cd4 = statistics::studentized_range_pdf(1.0, 1.0, 10.0);
    chk_code(bad, "studentized range one sample", cd4, 5);
    let _d5, cd5 = statistics::studentized_range_pdf(1.0, 3.0, 0.0);
    chk_code(bad, "studentized range degrees 0", cd5, 5);
    let _d6, cd6 = statistics::vonmises_pdf(1.0, -1.0);
    chk_code(bad, "vonmises shape below 0", cd6, 5);
    return bad;
}
