// Hostile input for the families of a power or an exponent. A NaN gives
// the NaN fault, an infinity gives the infinity fault and a shape outside
// its range gives the domain fault.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_exponent_hostile() -> i64 {
    var bad: i64 = 0;
    // A NaN value gives the NaN fault.
    let a3, b3 = statistics::exponpow_sf(0.0 / 0.0, 1.5);
    chk_code(bad, "a NaN value", b3, 3);
    // An infinite value gives the infinity fault.
    let a5, b5 = statistics::powernorm_cdf(1.0 / 0.0, 1.6);
    chk_code(bad, "an infinite value", b5, 4);
    // A NaN shape gives the NaN fault.
    let a7, b7 = statistics::powerlognorm_pdf(1.0, 0.0 / 0.0, 2.0);
    chk_code(bad, "a NaN shape", b7, 3);
    // A probability outside 0 to 1 gives the domain fault.
    let a10, b10 = statistics::exponweib_ppf(-0.1, 1.5, 2.0);
    chk_code(bad, "a probability below 0", b10, 5);
    let a11, b11 = statistics::powernorm_ppf(0.0 / 0.0, 1.6);
    chk_code(bad, "a NaN probability", b11, 3);
    // The quantile of 1 of an unbounded law does not fit in a result.
    let a12, b12 = statistics::exponpow_ppf(1.0, 1.5);
    chk_code(bad, "the quantile of 1 of an unbounded law", b12, 8);
    // A random sweep: every answer with the code OK is finite.
    var st: u64 = 8117;
    for k in 0..tests_count() {
        let x = ((rnd(st) % 2001) as f64 - 1000.0) / 100.0;
        let p = (rnd(st) % 1001) as f64 / 1000.0;
        let c1 = 0.5 + (rnd(st) % 30) as f64 / 10.0;
        let v1, cc1 = statistics::powernorm_pdf(x, c1);
        let v2, cc2 = statistics::exponweib_pdf(x, c1, 2.0);
        let v3, cc3 = statistics::powerlognorm_ppf(p, c1, 2.0);
        hostile_report(bad, cc1 == 0 && !is_nan(v1) && !is_inf(v1), "powernorm pdf", k, cc1);
        hostile_report(bad, cc2 == 0 && !is_nan(v2) && !is_inf(v2), "exponweib pdf", k, cc2);
        hostile_report(bad, (cc3 == 0 && !is_nan(v3) && !is_inf(v3)) || cc3 == 8, "powerlognorm ppf", k, cc3);
    }
    return bad;
}
