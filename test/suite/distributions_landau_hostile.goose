// Hostile input for the Landau law.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_landau_hostile() -> i64 {
    var bad: i64 = 0;
    let a, ka = statistics::landau_pdf(0.0 / 0.0);
    chk_code(bad, "landau a NaN point", ka, 3);
    let b, kb = statistics::landau_cdf(1.0 / 0.0);
    chk_code(bad, "landau an infinite point", kb, 4);
    let c, kc = statistics::landau_ppf(-0.1);
    chk_code(bad, "landau a probability below 0", kc, 5);
    var st: u64 = 99;
    for k in 0..tests_count() / 200 {
        let x = (rnd(st) % 1201) as f64 / 100.0 - 8.0;
        let v, cc = statistics::landau_cdf(x);
        hostile_report(bad, cc == 0 && v >= 0.0 && v <= 1.0, "landau cdf", k, cc);
    }
    return bad;
}
