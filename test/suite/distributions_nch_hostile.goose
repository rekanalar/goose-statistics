// Hostile input for the two non-central hypergeometric laws.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_nch_hostile() -> i64 {
    var bad: i64 = 0;
    let a, ka = statistics::nchypergeom_fisher_pmf(0.0 / 0.0, 30.0, 10.0, 12.0, 0.5);
    chk_code(bad, "nchypergeom_fisher a NaN count", ka, 3);
    let b, kb = statistics::nchypergeom_fisher_cdf(3.0, 30.0, 10.0, 12.0, 1.0 / 0.0);
    chk_code(bad, "nchypergeom_fisher an infinite odds ratio", kb, 4);
    let c, kc = statistics::nchypergeom_fisher_ppf(1.5, 30.0, 10.0, 12.0, 0.5);
    chk_code(bad, "nchypergeom_fisher a probability above 1", kc, 5);
    let d, kd = statistics::nchypergeom_fisher_pmf(1.0, 1e9, 10.0, 12.0, 0.5);
    chk_code(bad, "nchypergeom_fisher a total above the bound", kd, 5);
    let a, ka = statistics::nchypergeom_wallenius_pmf(0.0 / 0.0, 30.0, 10.0, 12.0, 0.5);
    chk_code(bad, "nchypergeom_wallenius a NaN count", ka, 3);
    let b, kb = statistics::nchypergeom_wallenius_cdf(3.0, 30.0, 10.0, 12.0, 1.0 / 0.0);
    chk_code(bad, "nchypergeom_wallenius an infinite odds ratio", kb, 4);
    let c, kc = statistics::nchypergeom_wallenius_ppf(1.5, 30.0, 10.0, 12.0, 0.5);
    chk_code(bad, "nchypergeom_wallenius a probability above 1", kc, 5);
    let d, kd = statistics::nchypergeom_wallenius_pmf(1.0, 1e9, 10.0, 12.0, 0.5);
    chk_code(bad, "nchypergeom_wallenius a total above the bound", kd, 5);
    return bad;
}
