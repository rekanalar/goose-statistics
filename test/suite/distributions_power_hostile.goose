// Hostile input for the families of a power. A NaN gives the NaN fault, an
// infinity gives the infinity fault and a shape outside its range gives the
// domain fault.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_power_hostile() -> i64 {
    var bad: i64 = 0;
    let a1, b1 = statistics::fisk_sf(0.0 / 0.0, 5.0);
    chk_code(bad, "a NaN value", b1, 3);
    let a2, b2 = statistics::burr_cdf(1.0 / 0.0, 5.0, 3.0);
    chk_code(bad, "an infinite value", b2, 4);
    let a3, b3 = statistics::burr12_pdf(1.0, 0.0 / 0.0, 3.0);
    chk_code(bad, "a NaN shape", b3, 3);
    let a4, b4 = statistics::genpareto_pdf(1.0, 0.0 / 0.0);
    chk_code(bad, "a NaN shape of genpareto", b4, 3);
    let a5, b5 = statistics::betaprime_ppf(-0.1, 3.0, 5.0);
    chk_code(bad, "a probability below 0", b5, 5);
    let a6, b6 = statistics::fisk_ppf(0.0 / 0.0, 5.0);
    chk_code(bad, "a NaN probability", b6, 3);
    let a7, b7 = statistics::fisk_ppf(1.0, 5.0);
    chk_code(bad, "the quantile of 1 of an unbounded law", b7, 8);
    // A random sweep: every answer with the code OK is finite.
    var st: u64 = 9127;
    for k in 0..tests_count() {
        let x = ((rnd(st) % 2001) as f64 - 1000.0) / 100.0;
        let p = (rnd(st) % 1001) as f64 / 1000.0;
        let c1 = 0.5 + (rnd(st) % 40) as f64 / 10.0;
        let v1, cc1 = statistics::fisk_pdf(x, c1);
        let v2, cc2 = statistics::burr12_cdf(x, c1, 2.0);
        let v3, cc3 = statistics::genpareto_ppf(p, c1 - 2.0);
        let v4, cc4 = statistics::betaprime_ppf(p, c1, 3.0);
        hostile_report(bad, cc1 == 0 && !is_nan(v1) && !is_inf(v1), "fisk pdf", k, cc1);
        hostile_report(bad, cc2 == 0 && !is_nan(v2) && !is_inf(v2), "burr12 cdf", k, cc2);
        hostile_report(bad, (cc3 == 0 && !is_nan(v3) && !is_inf(v3)) || cc3 == 8, "genpareto ppf", k, cc3);
        hostile_report(bad, (cc4 == 0 && !is_nan(v4) && !is_inf(v4)) || cc4 == 8, "betaprime ppf", k, cc4);
    }
    return bad;
}
