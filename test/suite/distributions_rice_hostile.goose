// Hostile input for the Rice law. A NaN gives the NaN fault, an infinity
// gives the infinity fault and a shape outside its range gives the domain
// fault.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_rice_hostile() -> i64 {
    var bad: i64 = 0;
    let a1, b1 = statistics::rice_pdf(0.0 / 0.0, 2.0);
    chk_code(bad, "a NaN point", b1, 3);
    let a2, b2 = statistics::rice_cdf(1.0 / 0.0, 2.0);
    chk_code(bad, "an infinite point", b2, 4);
    let a3, b3 = statistics::rice_sf(1.0, 0.0 / 0.0);
    chk_code(bad, "a NaN shape", b3, 3);
    let a4, b4 = statistics::rice_cdf(1.0, 1.0 / 0.0);
    chk_code(bad, "an infinite shape", b4, 4);
    let a5, b5 = statistics::rice_ppf(1.5, 2.0);
    chk_code(bad, "a probability above 1", b5, 5);
    let a6, b6 = statistics::rice_ppf(-0.1, 2.0);
    chk_code(bad, "a probability below 0", b6, 5);
    let a7, b7 = statistics::rice_ppf(0.0 / 0.0, 2.0);
    chk_code(bad, "a NaN probability", b7, 3);
    let a8, b8 = statistics::rice_mean(-1.0);
    chk_code(bad, "a negative shape of the mean", b8, 5);
    // A random sweep: every answer with the code OK is finite.
    var st: u64 = 5531;
    for k in 0..tests_count() / 50 {
        let x = (rnd(st) % 4001) as f64 / 100.0 - 5.0;
        let b = (rnd(st) % 401) as f64 / 10.0;
        let v1, cc1 = statistics::rice_pdf(x, b);
        let v2, cc2 = statistics::rice_cdf(x, b);
        let v3, cc3 = statistics::rice_sf(x, b);
        hostile_report(bad, cc1 == 0 && !is_nan(v1) && !is_inf(v1), "rice pdf", k, cc1);
        hostile_report(bad, cc2 == 0 && v2 >= 0.0 && v2 <= 1.0, "rice cdf", k, cc2);
        hostile_report(bad, cc3 == 0 && v3 >= 0.0 && v3 <= 1.0, "rice sf", k, cc3);
    }
    return bad;
}
