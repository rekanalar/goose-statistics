// Hostile input for the five laws of the second search group. A NaN
// gives the NaN fault, an infinity gives the infinity fault and a shape
// outside its range gives the domain fault.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_search_b_hostile() -> i64 {
    var bad: i64 = 0;
    let a1, b1 = statistics::loguniform_pdf(0.0 / 0.0, 1.0, 4.0);
    chk_code(bad, "a NaN point", b1, 3);
    let a2, b2 = statistics::loguniform_cdf(2.0, 1.0, 0.0 / 0.0);
    chk_code(bad, "a NaN bound", b2, 3);
    let a3, b3 = statistics::loguniform_cdf(2.0, 1.0, 1.0 / 0.0);
    chk_code(bad, "an infinite bound", b3, 4);
    let a4, b4 = statistics::irwinhall_cdf(1.0, 0.0 / 0.0);
    chk_code(bad, "a NaN shape", b4, 3);
    let a5, b5 = statistics::gausshyper_cdf(0.5, 2.0, 3.0, 0.0 / 0.0, 1.0);
    chk_code(bad, "a NaN shape", b5, 3);
    let a6, b6 = statistics::geninvgauss_cdf(1.0, 1.5, 0.0 / 0.0);
    chk_code(bad, "a NaN shape", b6, 3);
    let a7, b7 = statistics::genhyperbolic_ppf(1.5, 1.5, 1.0, 0.5);
    chk_code(bad, "a probability above 1", b7, 5);
    let a8, b8 = statistics::genhyperbolic_ppf(-0.1, 1.5, 1.0, 0.5);
    chk_code(bad, "a probability below 0", b8, 5);
    // A random sweep: every answer with the code OK is finite.
    var st: u64 = 8123;
    for k in 0..tests_count() / 20 {
        let x = (rnd(st) % 401) as f64 / 100.0 - 1.0;
        let a = (rnd(st) % 400 + 1) as f64 / 100.0;
        let b = a + (rnd(st) % 300 + 1) / 100.0;
        let v1, cc1 = statistics::loguniform_pdf(x, a, b);
        let v2, cc2 = statistics::gausshyper_cdf(x, a, b, 1.0, 0.5);
        let v3, cc3 = statistics::geninvgauss_cdf(x, a, b);
        hostile_report(bad, cc1 == 0 && !is_nan(v1) && !is_inf(v1), "loguniform pdf", k, cc1);
        hostile_report(bad, cc2 == 0 && v2 >= 0.0 && v2 <= 1.0, "gausshyper cdf", k, cc2);
        hostile_report(bad, cc3 == 0 && v3 >= 0.0 && v3 <= 1.0, "geninvgauss cdf", k, cc3);
    }
    return bad;
}
