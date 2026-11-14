// Hostile input for the second group of discrete families. A NaN gives
// the NaN fault, an infinity gives the infinity fault and a shape outside
// its range gives the domain fault.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_discrete_b_hostile() -> i64 {
    var bad: i64 = 0;
    // A NaN count gives the NaN fault.
    let a1, b1 = statistics::planck_pmf(0.0 / 0.0, 2.0);
    chk_code(bad, "a NaN count", b1, 3);
    let a2, b2 = statistics::randint_cdf(0.0 / 0.0, 0.0, 5.0);
    chk_code(bad, "a NaN count", b2, 3);
    // An infinite count gives the infinity fault.
    let a3, b3 = statistics::yulesimon_pmf(1.0 / 0.0, 2.0);
    chk_code(bad, "an infinite count", b3, 4);
    // A NaN shape gives the NaN fault.
    let a4, b4 = statistics::planck_pmf(1.0, 0.0 / 0.0);
    chk_code(bad, "a NaN shape", b4, 3);
    let a5, b5 = statistics::zipf_pmf(1.0, 0.0 / 0.0);
    chk_code(bad, "a NaN shape", b5, 3);
    // An infinite shape gives the infinity fault.
    let a6, b6 = statistics::zipf_pmf(1.0, 1.0 / 0.0);
    chk_code(bad, "an infinite shape", b6, 4);
    // A probability outside 0 to 1 gives the domain fault.
    let a7, b7 = statistics::planck_ppf(1.5, 2.0);
    chk_code(bad, "a probability above 1", b7, 5);
    let a8, b8 = statistics::zipf_ppf(-0.1, 2.0);
    chk_code(bad, "a probability below 0", b8, 5);
    // A count that is not a whole number gives the domain fault.
    let a9, b9 = statistics::planck_pmf(1.5, 2.0);
    chk_code(bad, "a count with a part", b9, 5);
    let a10, b10 = statistics::randint_pmf(1.0, 0.5, 5.0);
    chk_code(bad, "a bound with a part", b10, 5);
    // A probability of the list of the Poisson binomial outside 0 to 1
    // gives the domain fault.
    let bad_ps: f64[2] = [0.5, 1.5];
    let a11, b11 = statistics::poisson_binom_pmf(1.0, bad_ps);
    chk_code(bad, "a probability above 1 in the list", b11, 5);
    // A random sweep: every answer with the code OK is finite.
    var st: u64 = 9241;
    for k in 0..tests_count() {
        let x = (rnd(st) % 21) as f64 - 10.0;
        let p = (rnd(st) % 1001) as f64 / 1000.0;
        let a = (rnd(st) % 100) as f64 / 10.0 + 1.1;
        let n = (rnd(st) % 20) as f64;
        let ps: f64[2] = [p, 0.5];
        let v1, cc1 = statistics::planck_pmf(x, a);
        let v2, cc2 = statistics::zipf_cdf(x, a);
        let v3, cc3 = statistics::poisson_binom_cdf(x, ps);
        hostile_report(bad, cc1 == 0 && !is_nan(v1) && !is_inf(v1), "planck pmf", k, cc1);
        hostile_report(bad, cc2 == 0 && !is_nan(v2) && !is_inf(v2), "zipf cdf", k, cc2);
        hostile_report(bad, cc3 == 0 && !is_nan(v3) && !is_inf(v3), "poisson_binom cdf", k, cc3);
    }
    return bad;
}
