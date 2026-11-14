// Hostile input for the five laws of the third search group. A NaN gives
// the NaN fault, an infinity gives the infinity fault and a shape
// outside its range gives the domain fault.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_search_c_hostile() -> i64 {
    var bad: i64 = 0;
    let a1, b1 = statistics::norminvgauss_pdf(0.0 / 0.0, 1.0, 0.5);
    chk_code(bad, "a NaN point", b1, 3);
    let a2, b2 = statistics::norminvgauss_cdf(1.0, 0.0 / 0.0, 0.5);
    chk_code(bad, "a NaN shape", b2, 3);
    let a3, b3 = statistics::norminvgauss_cdf(1.0, 1.0, 1.0 / 0.0);
    chk_code(bad, "an infinite skewness", b3, 4);
    let a4, b4 = statistics::recipinvgauss_cdf(1.0, 0.0 / 0.0);
    chk_code(bad, "a NaN shape", b4, 3);
    let a5, b5 = statistics::rel_breitwigner_cdf(1.0, 1.0 / 0.0);
    chk_code(bad, "an infinite shape", b5, 4);
    let a6, b6 = statistics::studentized_range_ppf(1.5, 3.0, 10.0);
    chk_code(bad, "a probability above 1", b6, 5);
    let a7, b7 = statistics::studentized_range_ppf(-0.1, 3.0, 10.0);
    chk_code(bad, "a probability below 0", b7, 5);
    let a8, b8 = statistics::vonmises_ppf(0.5, 0.0 / 0.0);
    chk_code(bad, "a NaN shape", b8, 3);
    // A random sweep: every answer with the code OK is finite. The law
    // of von Mises repeats every turn, so its cumulative function grows
    // above one and falls below zero.
    var st: u64 = 6421;
    for k in 0..tests_count() / 20 {
        let x = (rnd(st) % 801) as f64 / 100.0 - 4.0;
        let a = (rnd(st) % 400 + 1) as f64 / 100.0;
        // The skewness of the normal inverse Gaussian law stays inside
        // the shape.
        let b = (rnd(st) % 99) as f64 / 100.0;
        let v1, cc1 = statistics::norminvgauss_cdf(x, 1.0, b);
        let v2, cc2 = statistics::recipinvgauss_cdf(x, a + 0.5);
        let v3, cc3 = statistics::vonmises_cdf(x, a);
        hostile_report(bad, cc1 == 0 && v1 >= 0.0 && v1 <= 1.0, "norminvgauss cdf", k, cc1);
        hostile_report(bad, cc2 == 0 && v2 >= 0.0 && v2 <= 1.0, "recipinvgauss cdf", k, cc2);
        hostile_report(bad, cc3 == 0 && !is_nan(v3) && !is_inf(v3), "vonmises cdf", k, cc3);
    }
    return bad;
}
