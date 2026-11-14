// Hostile input for the families of a special shape. A NaN gives the NaN
// fault, an infinity gives the infinity fault and a shape outside its range
// gives the domain fault.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_special_hostile() -> i64 {
    var bad: i64 = 0;
    let a1, b1 = statistics::invgauss_cdf(0.0 / 0.0, 0.5);
    chk_code(bad, "a NaN value", b1, 3);
    let a2, b2 = statistics::nakagami_sf(1.0 / 0.0, 2.5);
    chk_code(bad, "an infinite value", b2, 4);
    let a3, b3 = statistics::kappa3_pdf(1.0, 0.0 / 0.0);
    chk_code(bad, "a NaN shape", b3, 3);
    let a4, b4 = statistics::invgauss_ppf(1.5, 0.5);
    chk_code(bad, "a probability above 1", b4, 5);
    let a5, b5 = statistics::kappa3_ppf(0.0 / 0.0, 2.0);
    chk_code(bad, "a NaN probability", b5, 3);
    let a6, b6 = statistics::pearson3_pdf(1.0, 1.0 / 0.0);
    chk_code(bad, "an infinite shape", b6, 4);
    // A random sweep: every answer with the code OK is finite.
    var st: u64 = 20261113;
    for k in 0..tests_count() {
        let x = ((rnd(st) % 2001) as f64 - 1000.0) / 100.0;
        let s = (rnd(st) % 400) as f64 / 100.0 + 0.05;
        let v1, cc1 = statistics::invgauss_pdf(x, s);
        let v2, cc2 = statistics::nakagami_cdf(x, s);
        let v3, cc3 = statistics::kappa3_pdf(x, s);
        let v4, cc4 = statistics::pearson3_cdf(x, s - 1.987);
        hostile_report(bad, cc1 == 0 && !is_nan(v1) && !is_inf(v1), "invgauss pdf", k, cc1);
        hostile_report(bad, cc2 == 0 && !is_nan(v2) && !is_inf(v2), "nakagami cdf", k, cc2);
        hostile_report(bad, cc3 == 0 && !is_nan(v3) && !is_inf(v3), "kappa3 pdf", k, cc3);
        hostile_report(bad, cc4 == 0 && !is_nan(v4) && !is_inf(v4), "pearson3 cdf", k, cc4);
    }
    // The quantile of the inverse Gaussian law costs a search each. The
    // sweep uses a fiftieth of the count of the other calls.
    for k in 0..tests_count() / 50 {
        let p = (rnd(st) % 1001) as f64 / 1000.0;
        let s = (rnd(st) % 400) as f64 / 100.0 + 0.05;
        let v6, cc6 = statistics::invgauss_ppf(p, s);
        hostile_report(bad, (cc6 == 0 && !is_nan(v6) && !is_inf(v6)) || cc6 == 8, "invgauss ppf", k, cc6);
    }
    return bad;
}
