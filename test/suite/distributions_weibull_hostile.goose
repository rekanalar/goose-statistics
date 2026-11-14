// Hostile input for the Weibull families. A NaN gives the NaN fault, an
// infinity gives the infinity fault and a shape outside its range gives
// the domain fault.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_weibull_hostile() -> i64 {
    var bad: i64 = 0;
    // A NaN value gives the NaN fault.
    let a1, b1 = statistics::dweibull_pdf(0.0 / 0.0, 1.7);
    chk_code(bad, "a NaN value", b1, 3);
    let a2, b2 = statistics::invweibull_cdf(0.0 / 0.0, 1.5);
    chk_code(bad, "a NaN value", b2, 3);
    // An infinite value gives the infinity fault.
    let a4, b4 = statistics::weibull_max_pdf(-1.0 / 0.0, 1.5);
    chk_code(bad, "an infinite value", b4, 4);
    // A NaN shape gives the NaN fault.
    let a6, b6 = statistics::dweibull_pdf(1.0, 0.0 / 0.0);
    chk_code(bad, "a NaN shape", b6, 3);
    // An infinite shape gives the infinity fault.
    let a8, b8 = statistics::invweibull_pdf(1.0, 1.0 / 0.0);
    chk_code(bad, "an infinite shape", b8, 4);
    // A probability outside 0 to 1 gives the domain fault.
    let a9, b9 = statistics::dweibull_ppf(1.5, 1.7);
    chk_code(bad, "a probability above 1", b9, 5);
    // An infinite bound of the support does not fit in a result.
    let a12, b12 = statistics::dweibull_ppf(0.0, 1.7);
    chk_code(bad, "the quantile of 0 of a law on the whole line", b12, 8);
    let a13, b13 = statistics::invweibull_ppf(1.0, 1.5);
    chk_code(bad, "the quantile of 1 of an unbounded law", b13, 8);
    // A random sweep: every answer with the code OK is finite.
    var st: u64 = 3391;
    for k in 0..tests_count() {
        let x = ((rnd(st) % 2001) as f64 - 1000.0) / 100.0;
        let c1 = 0.5 + (rnd(st) % 30) as f64 / 10.0;
        let v1, cc1 = statistics::dweibull_pdf(x, c1);
        let v2, cc2 = statistics::weibull_max_cdf(x, c1);
        hostile_report(bad, cc1 == 0 && !is_nan(v1) && !is_inf(v1), "dweibull pdf", k, cc1);
        hostile_report(bad, cc2 == 0 && !is_nan(v2) && !is_inf(v2), "weibull_max cdf", k, cc2);
    }
    return bad;
}
