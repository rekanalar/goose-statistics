// Hostile input for the bounded families. A NaN gives the NaN fault, an
// infinity gives the infinity fault and a shape outside its range gives
// the domain fault.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_bounded_hostile() -> i64 {
    var bad: i64 = 0;
    // A NaN value gives the NaN fault.
    let a1, b1 = statistics::semicircular_pdf(0.0 / 0.0);
    chk_code(bad, "a NaN value", b1, 3);
    let a2, b2 = statistics::wald_cdf(0.0 / 0.0);
    chk_code(bad, "a NaN value", b2, 3);
    let a3, b3 = statistics::maxwell_sf(0.0 / 0.0);
    chk_code(bad, "a NaN value", b3, 3);
    // An infinite value gives the infinity fault.
    let a4, b4 = statistics::triang_pdf(1.0 / 0.0, 0.3);
    chk_code(bad, "an infinite value", b4, 4);
    let a5, b5 = statistics::moyal_cdf(-1.0 / 0.0);
    chk_code(bad, "an infinite value", b5, 4);
    // A NaN shape gives the NaN fault.
    let a6, b6 = statistics::triang_pdf(0.5, 0.0 / 0.0);
    chk_code(bad, "a NaN shape", b6, 3);
    let a7, b7 = statistics::bradford_pdf(0.5, 0.0 / 0.0);
    chk_code(bad, "a NaN shape", b7, 3);
    // An infinite shape gives the infinity fault.
    let a8, b8 = statistics::bradford_pdf(0.5, 1.0 / 0.0);
    chk_code(bad, "an infinite shape", b8, 4);
    // A probability outside 0 to 1 gives the domain fault.
    let a9, b9 = statistics::semicircular_ppf(1.5);
    chk_code(bad, "a probability above 1", b9, 5);
    let a10, b10 = statistics::triang_ppf(-0.1, 0.3);
    chk_code(bad, "a probability below 0", b10, 5);
    let a11, b11 = statistics::wald_ppf(0.0 / 0.0);
    chk_code(bad, "a NaN probability", b11, 3);
    // An infinite bound of the support does not fit in a result.
    let a12, b12 = statistics::wald_ppf(1.0);
    chk_code(bad, "the quantile of 1 of an unbounded law", b12, 8);
    let a13, b13 = statistics::maxwell_ppf(1.0);
    chk_code(bad, "the quantile of 1 of an unbounded law", b13, 8);
    // A random sweep: every answer with the code OK is finite.
    var st: u64 = 4471;
    for k in 0..tests_count() {
        let x = ((rnd(st) % 2001) as f64 - 1000.0) / 100.0;
        let p = (rnd(st) % 1001) as f64 / 1000.0;
        let c1 = (rnd(st) % 5) as f64 / 10.0;
        let d1 = c1 + (rnd(st) % 5) / 10.0;
        let v1, cc1 = statistics::trapezoid_pdf(x, c1, d1);
        let v2, cc2 = statistics::bradford_ppf(p, 1.5);
        let v3, cc3 = statistics::moyal_cdf(x);
        hostile_report(bad, cc1 == 0 && !is_nan(v1) && !is_inf(v1), "trapezoid pdf", k, cc1);
        hostile_report(bad, (cc2 == 0 && !is_nan(v2) && !is_inf(v2)) || cc2 == 8, "bradford ppf", k, cc2);
        hostile_report(bad, cc3 == 0 && !is_nan(v3) && !is_inf(v3), "moyal cdf", k, cc3);
    }
    return bad;
}
