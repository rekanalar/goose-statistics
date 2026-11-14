// Hostile input for the first group of discrete families. A NaN gives the
// NaN fault, an infinity gives the infinity fault and a shape outside its
// range gives the domain fault.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_discrete_a_hostile() -> i64 {
    var bad: i64 = 0;
    // A NaN count gives the NaN fault.
    let a1, b1 = statistics::bernoulli_pmf(0.0 / 0.0, 0.5);
    chk_code(bad, "a NaN count", b1, 3);
    let a2, b2 = statistics::betabinom_cdf(0.0 / 0.0, 10.0, 2.0, 3.0);
    chk_code(bad, "a NaN count", b2, 3);
    // An infinite count gives the infinity fault.
    let a3, b3 = statistics::boltzmann_pmf(1.0 / 0.0, 2.0, 5.0);
    chk_code(bad, "an infinite count", b3, 4);
    // A NaN shape gives the NaN fault.
    let a4, b4 = statistics::betabinom_pmf(1.0, 10.0, 0.0 / 0.0, 3.0);
    chk_code(bad, "a NaN shape", b4, 3);
    let a5, b5 = statistics::logser_pmf(1.0, 0.0 / 0.0);
    chk_code(bad, "a NaN shape", b5, 3);
    // An infinite shape gives the infinity fault.
    let a6, b6 = statistics::dlaplace_pmf(1.0, 1.0 / 0.0);
    chk_code(bad, "an infinite shape", b6, 4);
    // A probability outside 0 to 1 gives the domain fault.
    let a7, b7 = statistics::bernoulli_ppf(1.5, 0.5);
    chk_code(bad, "a probability above 1", b7, 5);
    let a8, b8 = statistics::betabinom_ppf(-0.1, 10.0, 2.0, 3.0);
    chk_code(bad, "a probability below 0", b8, 5);
    let a9, b9 = statistics::logser_ppf(0.0 / 0.0, 0.5);
    chk_code(bad, "a NaN probability", b9, 3);
    // A count that is not a whole number gives the domain fault.
    let a10, b10 = statistics::betabinom_pmf(1.5, 10.0, 2.0, 3.0);
    chk_code(bad, "a count with a part", b10, 5);
    let a11, b11 = statistics::boltzmann_ppf(0.5, 2.0, 5.5);
    chk_code(bad, "a count with a part", b11, 5);
    // A random sweep: every answer with the code OK is finite.
    var st: u64 = 8821;
    for k in 0..tests_count() {
        let x = (rnd(st) % 21) as f64 - 10.0;
        let p = (rnd(st) % 999 + 1) as f64 / 1000.0;
        let a = (rnd(st) % 100) as f64 / 10.0 + 0.1;
        let b = (rnd(st) % 100) as f64 / 10.0 + 0.1;
        let n = (rnd(st) % 20) as f64 + 1.0;
        let v1, cc1 = statistics::betabinom_pmf(x, n, a, b);
        let v2, cc2 = statistics::betanbinom_pmf(x, n, a, b);
        let v3, cc3 = statistics::logser_cdf(x, p);
        hostile_report(bad, cc1 == 0 && !is_nan(v1) && !is_inf(v1), "betabinom pmf", k, cc1);
        hostile_report(bad, cc2 == 0 && !is_nan(v2) && !is_inf(v2), "betanbinom pmf", k, cc2);
        hostile_report(bad, cc3 == 0 && !is_nan(v3) && !is_inf(v3), "logser cdf", k, cc3);
    }
    return bad;
}
