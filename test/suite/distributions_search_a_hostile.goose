// Hostile input for the six laws of the first search group. A NaN gives
// the NaN fault, an infinity gives the infinity fault and a shape
// outside its range gives the domain fault.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_search_a_hostile() -> i64 {
    var bad: i64 = 0;
    let a1, b1 = statistics::argus_pdf(0.0 / 0.0, 1.0);
    chk_code(bad, "a NaN point", b1, 3);
    let a2, b2 = statistics::argus_cdf(0.5, 0.0 / 0.0);
    chk_code(bad, "a NaN shape", b2, 3);
    let a3, b3 = statistics::argus_cdf(0.5, 1.0 / 0.0);
    chk_code(bad, "an infinite shape", b3, 4);
    let a4, b4 = statistics::dpareto_lognorm_cdf(1.0, 0.0 / 0.0, 1.0, 1.0, 1.0);
    chk_code(bad, "a NaN location", b4, 3);
    let a5, b5 = statistics::dpareto_lognorm_pdf(1.0, 0.0, 1.0, 1.0, 1.0 / 0.0);
    chk_code(bad, "an infinite shape", b5, 4);
    let a6, b6 = statistics::erlang_cdf(1.0, 0.0 / 0.0);
    chk_code(bad, "a NaN shape", b6, 3);
    let a7, b7 = statistics::exponnorm_cdf(1.0, 1.0 / 0.0);
    chk_code(bad, "an infinite shape", b7, 4);
    let a8, b8 = statistics::foldcauchy_ppf(1.5, 1.0);
    chk_code(bad, "a probability above 1", b8, 5);
    let a9, b9 = statistics::foldnorm_ppf(-0.1, 0.5);
    chk_code(bad, "a probability below 0", b9, 5);
    // A random sweep: every answer with the code OK is finite.
    var st: u64 = 4711;
    for k in 0..tests_count() / 20 {
        let x = (rnd(st) % 401) as f64 / 100.0 - 1.0;
        let c = (rnd(st) % 400 + 1) as f64 / 100.0;
        // The Erlang law wants a whole number of terms.
        let a = (rnd(st) % 9 + 1) as f64;
        let v1, cc1 = statistics::argus_pdf(x, c);
        let v2, cc2 = statistics::foldnorm_cdf(x, c);
        let v3, cc3 = statistics::erlang_cdf(x, a);
        hostile_report(bad, cc1 == 0 && !is_nan(v1) && !is_inf(v1), "argus pdf", k, cc1);
        hostile_report(bad, cc2 == 0 && v2 >= 0.0 && v2 <= 1.0, "foldnorm cdf", k, cc2);
        hostile_report(bad, cc3 == 0 && v3 >= 0.0 && v3 <= 1.0, "erlang cdf", k, cc3);
    }
    return bad;
}
