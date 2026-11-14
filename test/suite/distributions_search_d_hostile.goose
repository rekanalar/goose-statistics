// Hostile input for the stable law. A NaN gives the NaN fault, an
// infinity gives the infinity fault and a shape outside its range gives
// the domain fault.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_search_d_hostile() -> i64 {
    var bad: i64 = 0;
    let a1, b1 = statistics::levy_stable_pdf(0.0 / 0.0, 1.5, 0.0);
    chk_code(bad, "a NaN point", b1, 3);
    let a2, b2 = statistics::levy_stable_cdf(1.0, 0.0 / 0.0, 0.0);
    chk_code(bad, "a NaN shape", b2, 3);
    let a3, b3 = statistics::levy_stable_cdf(1.0, 1.5, 1.0 / 0.0);
    chk_code(bad, "an infinite skewness", b3, 4);
    let a4, b4 = statistics::levy_stable_ppf(1.5, 1.5, 0.0);
    chk_code(bad, "a probability above 1", b4, 5);
    let a5, b5 = statistics::levy_stable_ppf(-0.1, 1.5, 0.0);
    chk_code(bad, "a probability below 0", b5, 5);
    // An infinite point gives the limit of the law.
    let v1, c1 = statistics::levy_stable_cdf(1.0 / 0.0, 1.5, 0.0);
    chk_val(bad, "an infinite point above", v1, c1, 1.0, 1e-12);
    let v2, c2 = statistics::levy_stable_cdf(-1.0 / 0.0, 1.5, 0.0);
    chk_val(bad, "an infinite point below", v2, c2, 0.0, 1e-12);
    let v3, c3 = statistics::levy_stable_pdf(1.0 / 0.0, 1.5, 0.0);
    chk_val(bad, "an infinite point of the density", v3, c3, 0.0, 1e-12);
    // A random sweep: every answer with the code OK is finite.
    var st: u64 = 3187;
    for k in 0..tests_count() / 20 {
        let x = (rnd(st) % 801) as f64 / 100.0 - 4.0;
        let a = (rnd(st) % 195 + 5) as f64 / 100.0;
        let b = (rnd(st) % 201 - 100) as f64 / 100.0;
        let v4, cc4 = statistics::levy_stable_cdf(x, a, b);
        let v5, cc5 = statistics::levy_stable_pdf(x, a, b);
        hostile_report(bad, cc4 == 0 && v4 >= 0.0 && v4 <= 1.0, "levy cdf", k, cc4);
        hostile_report(bad, cc5 == 0 && !is_nan(v5) && !is_inf(v5), "levy pdf", k, cc5);
    }
    return bad;
}
