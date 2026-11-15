// The rows of the five laws of the second search group that the vector
// file does not hold: the support of each law, the bounds of the shapes
// and the round trip of the quantile.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_search_b_tests() -> i64 {
    var bad: i64 = 0;
    // A value below the support gives a density of 0 and a cumulative
    // function of 0.
    let v1, c1 = statistics::loguniform_pdf(0.5, 1.0, 4.0);
    chk_val(bad, "loguniform pdf below the support", v1, c1, 0.0, 1e-12);
    let v2, c2 = statistics::loguniform_cdf(0.5, 1.0, 4.0);
    chk_val(bad, "loguniform cdf below the support", v2, c2, 0.0, 1e-12);
    let v3, c3 = statistics::irwinhall_pdf(-0.5, 3.0);
    chk_val(bad, "irwinhall pdf below the support", v3, c3, 0.0, 1e-12);
    let v4, c4 = statistics::irwinhall_cdf(-0.5, 3.0);
    chk_val(bad, "irwinhall cdf below the support", v4, c4, 0.0, 1e-12);
    let v5, c5 = statistics::gausshyper_pdf(-0.5, 2.0, 3.0, 1.0, 1.0);
    chk_val(bad, "gausshyper pdf below the support", v5, c5, 0.0, 1e-12);
    let v6, c6 = statistics::gausshyper_cdf(-0.5, 2.0, 3.0, 1.0, 1.0);
    chk_val(bad, "gausshyper cdf below the support", v6, c6, 0.0, 1e-12);
    // The bounds of the log-uniform law must rise. Equal bounds give the
    // domain fault.
    let _e0, ce0 = statistics::loguniform_pdf(1.0, 1.0, 1.0);
    chk_code(bad, "loguniform with equal bounds", ce0, 5);
    // The cumulative function of the quantile gives the probability back.
    let q1, kq1 = statistics::loguniform_ppf(0.4, 1.0, 4.0);
    let r1, kr1 = statistics::loguniform_cdf(q1, 1.0, 4.0);
    chk_val(bad, "loguniform cdf of the quantile", r1, kq1 + kr1, 0.4, 1e-9);
    let q2, kq2 = statistics::irwinhall_ppf(0.6, 3.0);
    let r2, kr2 = statistics::irwinhall_cdf(q2, 3.0);
    chk_val(bad, "irwinhall cdf of the quantile", r2, kq2 + kr2, 0.6, 1e-9);
    let q3, kq3 = statistics::gausshyper_ppf(0.3, 2.0, 3.0, 1.0, 1.0);
    let r3, kr3 = statistics::gausshyper_cdf(q3, 2.0, 3.0, 1.0, 1.0);
    chk_val(bad, "gausshyper cdf of the quantile", r3, kq3 + kr3, 0.3, 1e-9);
    let q4, kq4 = statistics::geninvgauss_ppf(0.7, 1.5, 2.0);
    let r4, kr4 = statistics::geninvgauss_cdf(q4, 1.5, 2.0);
    chk_val(bad, "geninvgauss cdf of the quantile", r4, kq4 + kr4, 0.7, 1e-9);
    let q5, kq5 = statistics::genhyperbolic_ppf(0.5, 1.5, 1.0, 0.5);
    let r5, kr5 = statistics::genhyperbolic_cdf(q5, 1.5, 1.0, 0.5);
    chk_val(bad, "genhyperbolic cdf of the quantile", r5, kq5 + kr5, 0.5, 1e-9);
    // The quantile of the end points gives the bounds of the support.
    let e1, ce1 = statistics::loguniform_ppf(0.0, 1.0, 4.0);
    chk_val(bad, "loguniform ppf 0", e1, ce1, 1.0, 1e-12);
    let e2, ce2 = statistics::loguniform_ppf(1.0, 1.0, 4.0);
    chk_val(bad, "loguniform ppf 1", e2, ce2, 4.0, 1e-12);
    let e3, ce3 = statistics::irwinhall_ppf(0.0, 3.0);
    chk_val(bad, "irwinhall ppf 0", e3, ce3, 0.0, 1e-12);
    let e4, ce4 = statistics::irwinhall_ppf(1.0, 3.0);
    chk_val(bad, "irwinhall ppf 1", e4, ce4, 3.0, 1e-12);
    // A shape outside its range gives the domain fault.
    let _d1, cd1 = statistics::loguniform_pdf(2.0, 4.0, 1.0);
    chk_code(bad, "loguniform bounds in the wrong order", cd1, 5);
    let _d2, cd2 = statistics::irwinhall_pdf(2.0, 2.5);
    chk_code(bad, "irwinhall shape with a part", cd2, 5);
    let _d3, cd3 = statistics::gausshyper_pdf(0.5, -1.0, 3.0, 1.0, 1.0);
    chk_code(bad, "gausshyper shape below 0", cd3, 5);
    let _d4, cd4 = statistics::geninvgauss_pdf(1.0, 1.5, 0.0);
    chk_code(bad, "geninvgauss shape 0", cd4, 5);
    let _d5, cd5 = statistics::genhyperbolic_pdf(1.0, 1.5, -1.0, 0.5);
    chk_code(bad, "genhyperbolic shape below 0", cd5, 5);
    return bad;
}
