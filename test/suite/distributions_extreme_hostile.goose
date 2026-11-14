// Hostile input for the generalised extreme value law. A NaN gives the NaN
// fault, an infinity gives the infinity fault and a shape outside its range
// gives the domain fault.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_extreme_hostile() -> i64 {
    var bad: i64 = 0;
    let a1, b1 = statistics::genextreme_cdf(0.0 / 0.0, 0.5);
    chk_code(bad, "a NaN value", b1, 3);
    let a2, b2 = statistics::genextreme_sf(1.0 / 0.0, 0.5);
    chk_code(bad, "an infinite value", b2, 4);
    let a3, b3 = statistics::genextreme_pdf(1.0, 0.0 / 0.0);
    chk_code(bad, "a NaN shape", b3, 3);
    let a4, b4 = statistics::genextreme_ppf(1.5, 0.5);
    chk_code(bad, "a probability above 1", b4, 5);
    let a5, b5 = statistics::genextreme_ppf(0.0 / 0.0, 0.5);
    chk_code(bad, "a NaN probability", b5, 3);
    // A random sweep: every answer with the code OK is finite.
    var st: u64 = 20261115;
    for k in 0..tests_count() {
        let x = ((rnd(st) % 2001) as f64 - 1000.0) / 100.0;
        let s = ((rnd(st) % 801) as f64 - 400.0) / 100.0;
        let v1, cc1 = statistics::genextreme_cdf(x, s);
        let v2, cc2 = statistics::genextreme_pdf(x, s);
        let v3, cc3 = statistics::genextreme_sf(x, s);
        hostile_report(bad, cc1 == 0 && !is_nan(v1) && !is_inf(v1), "genextreme cdf", k, cc1);
        hostile_report(bad, cc2 == 0 && !is_nan(v2) && !is_inf(v2), "genextreme pdf", k, cc2);
        hostile_report(bad, cc3 == 0 && !is_nan(v3) && !is_inf(v3), "genextreme sf", k, cc3);
    }
    return bad;
}
