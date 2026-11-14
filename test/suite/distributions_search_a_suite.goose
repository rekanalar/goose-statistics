// The rows of the six laws of the first search group that the vector
// file does not hold: the support of each law, the bounds of the shapes
// and the round trip of the quantile.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_search_a_tests() -> i64 {
    var bad: i64 = 0;
    // A value below the support gives a density of 0 and a cumulative
    // function of 0.
    let v1, c1 = statistics::argus_pdf(-1.0, 1.0);
    chk_val(bad, "argus pdf below the support", v1, c1, 0.0, 1e-12);
    let v2, c2 = statistics::argus_cdf(-1.0, 1.0);
    chk_val(bad, "argus cdf below the support", v2, c2, 0.0, 1e-12);
    let v3, c3 = statistics::dpareto_lognorm_pdf(-1.0, 0.0, 1.0, 1.0, 1.0);
    chk_val(bad, "dpl pdf below the support", v3, c3, 0.0, 1e-12);
    let v4, c4 = statistics::foldnorm_pdf(-1.0, 0.5);
    chk_val(bad, "foldnorm pdf below the support", v4, c4, 0.0, 1e-12);
    let v5, c5 = statistics::foldcauchy_cdf(-1.0, 1.0);
    chk_val(bad, "foldcauchy cdf below the support", v5, c5, 0.0, 1e-12);
    let v6, c6 = statistics::erlang_cdf(-1.0, 3.0);
    chk_val(bad, "erlang cdf below the support", v6, c6, 0.0, 1e-12);
    // The cumulative function of the quantile gives the probability back.
    let q1, kq1 = statistics::argus_ppf(0.4, 1.0);
    let r1, kr1 = statistics::argus_cdf(q1, 1.0);
    chk_val(bad, "argus cdf of the quantile", r1, kq1 + kr1, 0.4, 1e-9);
    let q2, kq2 = statistics::dpareto_lognorm_ppf(0.6, 0.0, 1.0, 1.0, 1.0);
    let r2, kr2 = statistics::dpareto_lognorm_cdf(q2, 0.0, 1.0, 1.0, 1.0);
    chk_val(bad, "dpl cdf of the quantile", r2, kq2 + kr2, 0.6, 1e-9);
    let q3, kq3 = statistics::erlang_ppf(0.3, 3.0);
    let r3, kr3 = statistics::erlang_cdf(q3, 3.0);
    chk_val(bad, "erlang cdf of the quantile", r3, kq3 + kr3, 0.3, 1e-9);
    let q4, kq4 = statistics::exponnorm_ppf(0.7, 1.0);
    let r4, kr4 = statistics::exponnorm_cdf(q4, 1.0);
    chk_val(bad, "exponnorm cdf of the quantile", r4, kq4 + kr4, 0.7, 1e-9);
    let q5, kq5 = statistics::foldcauchy_ppf(0.5, 1.0);
    let r5, kr5 = statistics::foldcauchy_cdf(q5, 1.0);
    chk_val(bad, "foldcauchy cdf of the quantile", r5, kq5 + kr5, 0.5, 1e-9);
    let q6, kq6 = statistics::foldnorm_ppf(0.8, 0.5);
    let r6, kr6 = statistics::foldnorm_cdf(q6, 0.5);
    chk_val(bad, "foldnorm cdf of the quantile", r6, kq6 + kr6, 0.8, 1e-9);
    // The quantile of the end points gives the bounds of the support.
    let e1, ce1 = statistics::argus_ppf(0.0, 1.0);
    chk_val(bad, "argus ppf 0", e1, ce1, 0.0, 1e-12);
    let e2, ce2 = statistics::argus_ppf(1.0, 1.0);
    chk_val(bad, "argus ppf 1", e2, ce2, 1.0, 1e-12);
    let e3, ce3 = statistics::erlang_ppf(0.0, 3.0);
    chk_val(bad, "erlang ppf 0", e3, ce3, 0.0, 1e-12);
    let e4, ce4 = statistics::foldnorm_ppf(0.0, 0.5);
    chk_val(bad, "foldnorm ppf 0", e4, ce4, 0.0, 1e-12);
    // A shape outside its range gives the domain fault.
    let _d1, cd1 = statistics::argus_pdf(0.5, 0.0);
    chk_code(bad, "argus shape 0", cd1, 5);
    let _d2, cd2 = statistics::dpareto_lognorm_pdf(2.0, 0.0, 0.0, 1.0, 1.0);
    chk_code(bad, "dpl scale 0", cd2, 5);
    let _d3, cd3 = statistics::erlang_pdf(2.0, 2.5);
    chk_code(bad, "erlang shape with a part", cd3, 5);
    let _d4, cd4 = statistics::exponnorm_pdf(1.0, -1.0);
    chk_code(bad, "exponnorm shape below 0", cd4, 5);
    let _d5, cd5 = statistics::foldcauchy_pdf(1.0, -1.0);
    chk_code(bad, "foldcauchy shape below 0", cd5, 5);
    let _d6, cd6 = statistics::foldnorm_pdf(1.0, -1.0);
    chk_code(bad, "foldnorm shape below 0", cd6, 5);
    return bad;
}
