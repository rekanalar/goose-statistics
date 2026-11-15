// The rows of the two non-central hypergeometric laws that the vector
// file does not hold: the sum of the masses, the round trip of the
// quantile and the bounds.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_nch_tests() -> i64 {
    var bad: i64 = 0;
    var total_nchypergeom_fisher: f64 = 0.0;
    for k in 0..16 {
        let v, c = statistics::nchypergeom_fisher_pmf(k, 30.0, 10.0, 12.0, 0.5);
        if c == 0 { total_nchypergeom_fisher += v; } else { bad += 1; }
    }
    chk_val(bad, "nchypergeom_fisher masses sum to 1", total_nchypergeom_fisher, 0, 1.0, 1e-9);
    let q, kq = statistics::nchypergeom_fisher_ppf(0.5, 30.0, 10.0, 12.0, 0.5);
    let r, kr = statistics::nchypergeom_fisher_cdf(q, 30.0, 10.0, 12.0, 0.5);
    check(bad, kq == 0 && kr == 0 && r >= 0.5, "nchypergeom_fisher cdf of the quantile reaches the probability");
    let e1, ce1 = statistics::nchypergeom_fisher_pmf(0.5, 30.0, 10.0, 12.0, 0.5);
    chk_val(bad, "nchypergeom_fisher pmf of a fraction", e1, ce1, 0.0, 1e-12);
    let _d1, cd1 = statistics::nchypergeom_fisher_pmf(1.0, 30.0, 40.0, 12.0, 0.5);
    chk_code(bad, "nchypergeom_fisher successes above the total", cd1, 5);
    let _d2, cd2 = statistics::nchypergeom_fisher_pmf(1.0, 30.0, 10.0, 12.0, 0.0);
    chk_code(bad, "nchypergeom_fisher odds 0", cd2, 5);
    var total_nchypergeom_wallenius: f64 = 0.0;
    for k in 0..16 {
        let v, c = statistics::nchypergeom_wallenius_pmf(k, 30.0, 10.0, 12.0, 0.5);
        if c == 0 { total_nchypergeom_wallenius += v; } else { bad += 1; }
    }
    chk_val(bad, "nchypergeom_wallenius masses sum to 1", total_nchypergeom_wallenius, 0, 1.0, 1e-9);
    let q, kq = statistics::nchypergeom_wallenius_ppf(0.5, 30.0, 10.0, 12.0, 0.5);
    let r, kr = statistics::nchypergeom_wallenius_cdf(q, 30.0, 10.0, 12.0, 0.5);
    check(bad, kq == 0 && kr == 0 && r >= 0.5, "nchypergeom_wallenius cdf of the quantile reaches the probability");
    let e1, ce1 = statistics::nchypergeom_wallenius_pmf(0.5, 30.0, 10.0, 12.0, 0.5);
    chk_val(bad, "nchypergeom_wallenius pmf of a fraction", e1, ce1, 0.0, 1e-12);
    let _d1, cd1 = statistics::nchypergeom_wallenius_pmf(1.0, 30.0, 40.0, 12.0, 0.5);
    chk_code(bad, "nchypergeom_wallenius successes above the total", cd1, 5);
    let _d2, cd2 = statistics::nchypergeom_wallenius_pmf(1.0, 30.0, 10.0, 12.0, 0.0);
    chk_code(bad, "nchypergeom_wallenius odds 0", cd2, 5);
    return bad;
}
