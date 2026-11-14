// Hostile input for the families of a truncation. A NaN gives the NaN
// fault, an infinity gives the infinity fault and a bound outside its range
// gives the domain fault.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_truncated_hostile() -> i64 {
    var bad: i64 = 0;
    let a1, b1 = statistics::truncexpon_sf(0.0 / 0.0, 2.0);
    chk_code(bad, "a NaN value", b1, 3);
    let a2, b2 = statistics::truncpareto_cdf(1.0 / 0.0, 2.0, 5.0);
    chk_code(bad, "an infinite value", b2, 4);
    let a3, b3 = statistics::truncnorm_pdf(0.0, 0.0 / 0.0, 2.0);
    chk_code(bad, "a NaN bound", b3, 3);
    let a4, b4 = statistics::truncexpon_ppf(1.5, 2.0);
    chk_code(bad, "a probability above 1", b4, 5);
    let a5, b5 = statistics::truncpareto_ppf(0.0 / 0.0, 2.0, 5.0);
    chk_code(bad, "a NaN probability", b5, 3);
    let a6, b6 = statistics::truncnorm_pdf(0.0, 0.0 / 0.0, 2.0);
    chk_code(bad, "a NaN bound of truncnorm", b6, 3);
    // A random sweep: every answer with the code OK is finite.
    var st: u64 = 6413;
    for k in 0..tests_count() {
        let x = ((rnd(st) % 2001) as f64 - 1000.0) / 100.0;
        let b1 = 0.5 + (rnd(st) % 40) as f64 / 10.0;
        let b2 = 1.5 + (rnd(st) % 90) as f64 / 10.0;
        let v1, cc1 = statistics::truncexpon_pdf(x, b1);
        let v2, cc2 = statistics::truncpareto_cdf(x, b1, b2);
        hostile_report(bad, cc1 == 0 && !is_nan(v1) && !is_inf(v1), "truncexpon pdf", k, cc1);
        hostile_report(bad, cc2 == 0 && !is_nan(v2) && !is_inf(v2), "truncpareto cdf", k, cc2);
    }
    // The quantile of this group costs a search each. The sweep uses a
    // fiftieth of the count of the other calls.
    for k in 0..tests_count() / 50 {
        let p = (rnd(st) % 1001) as f64 / 1000.0;
        let b1 = 0.5 + (rnd(st) % 40) as f64 / 10.0;
        let v3, cc3 = statistics::truncnorm_ppf(p, -b1, b1);
        hostile_report(bad, cc3 == 0 && !is_nan(v3) && !is_inf(v3), "truncnorm ppf", k, cc3);
    }
    return bad;
}
