// Hostile input for the families that come from the gamma law. A NaN gives
// the NaN fault, an infinity gives the infinity fault and a shape outside
// its range gives the domain fault.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_gamma_hostile() -> i64 {
    var bad: i64 = 0;
    let a1, b1 = statistics::chi_sf(0.0 / 0.0, 3.0);
    chk_code(bad, "a NaN value", b1, 3);
    let a2, b2 = statistics::dgamma_cdf(1.0 / 0.0, 2.0);
    chk_code(bad, "an infinite value", b2, 4);
    let a3, b3 = statistics::invgamma_pdf(1.0, 0.0 / 0.0);
    chk_code(bad, "a NaN shape", b3, 3);
    let a4, b4 = statistics::loggamma_ppf(1.5, 2.0);
    chk_code(bad, "a probability above 1", b4, 5);
    let a5, b5 = statistics::chi_ppf(0.0 / 0.0, 3.0);
    chk_code(bad, "a NaN probability", b5, 3);
    let a6, b6 = statistics::loggamma_pdf(1.0, 1.0 / 0.0);
    chk_code(bad, "an infinite shape", b6, 4);
    // A random sweep: every answer with the code OK is finite.
    var st: u64 = 4271;
    for k in 0..tests_count() {
        let x = ((rnd(st) % 2001) as f64 - 1000.0) / 100.0;
        let c1 = 0.5 + (rnd(st) % 40) as f64 / 10.0;
        let v1, cc1 = statistics::chi_pdf(x, c1);
        let v2, cc2 = statistics::dgamma_cdf(x, c1);
        // The density of the chi law does not fit for a large point.
        hostile_report(bad, (cc1 == 0 && !is_nan(v1) && !is_inf(v1)) || cc1 == 8, "chi pdf", k, cc1);
        hostile_report(bad, cc2 == 0 && !is_nan(v2) && !is_inf(v2), "dgamma cdf", k, cc2);
    }
    // The quantiles of this group cost a search each. The sweep uses a
    // fiftieth of the count of the other calls.
    for k in 0..tests_count() / 50 {
        let p = (rnd(st) % 1001) as f64 / 1000.0;
        let c1 = 0.5 + (rnd(st) % 40) as f64 / 10.0;
        let v3, cc3 = statistics::invgamma_ppf(p, c1);
        let v4, cc4 = statistics::loggamma_ppf(p, c1);
        hostile_report(bad, (cc3 == 0 && !is_nan(v3) && !is_inf(v3)) || cc3 == 8, "invgamma ppf", k, cc3);
        hostile_report(bad, (cc4 == 0 && !is_nan(v4) && !is_inf(v4)) || cc4 == 8, "loggamma ppf", k, cc4);
    }
    return bad;
}
