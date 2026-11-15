// The rows of the first group of discrete families that the vector file
// does not hold: the support of each law, the bounds of the shapes and
// the round trip of the quantile.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_discrete_a_tests() -> i64 {
    var bad: i64 = 0;
    // A count outside the support gives a mass of 0 and not a fault.
    let v1, c1 = statistics::bernoulli_pmf(2.0, 0.5);
    chk_val(bad, "bernoulli pmf above the support", v1, c1, 0.0, 1e-12);
    let v2, c2 = statistics::betabinom_pmf(11.0, 10.0, 2.0, 3.0);
    chk_val(bad, "betabinom pmf above the support", v2, c2, 0.0, 1e-12);
    let v3, c3 = statistics::betanbinom_pmf(-1.0, 10.0, 2.0, 3.0);
    chk_val(bad, "betanbinom pmf below the support", v3, c3, 0.0, 1e-12);
    let v4, c4 = statistics::boltzmann_pmf(5.0, 2.0, 5.0);
    chk_val(bad, "boltzmann pmf above the support", v4, c4, 0.0, 1e-12);
    let v5, c5 = statistics::logser_pmf(0.0, 0.5);
    chk_val(bad, "logser pmf below the support", v5, c5, 0.0, 1e-12);
    let v6, c6 = statistics::nhypergeom_pmf(6.0, 10.0, 5.0, 3.0);
    chk_val(bad, "nhypergeom pmf above the support", v6, c6, 0.0, 1e-12);
    // The cumulative function of a count below the support is 0 and above
    // the support is 1.
    let v7, c7 = statistics::betabinom_cdf(-1.0, 10.0, 2.0, 3.0);
    chk_val(bad, "betabinom cdf below the support", v7, c7, 0.0, 1e-12);
    let v8, c8 = statistics::boltzmann_cdf(10.0, 2.0, 5.0);
    chk_val(bad, "boltzmann cdf above the support", v8, c8, 1.0, 1e-12);
    let v9, c9 = statistics::nhypergeom_cdf(-1.0, 10.0, 5.0, 3.0);
    chk_val(bad, "nhypergeom cdf below the support", v9, c9, 0.0, 1e-12);
    // The cumulative function of the quantile of a probability is that
    // probability.
    let q1, kq1 = statistics::betabinom_ppf(0.4, 10.0, 2.0, 3.0);
    chk_code(bad, "betabinom ppf code", kq1, 0);
    let r1, kr1 = statistics::betabinom_cdf(q1, 10.0, 2.0, 3.0);
    check(bad, kr1 == 0 && r1 >= 0.4 - 1e-9, "betabinom cdf of the quantile");
    let q2, kq2 = statistics::betanbinom_ppf(0.3, 10.0, 2.0, 3.0);
    chk_code(bad, "betanbinom ppf code", kq2, 0);
    let r2, kr2 = statistics::betanbinom_cdf(q2, 10.0, 2.0, 3.0);
    check(bad, kr2 == 0 && r2 >= 0.3 - 1e-9, "betanbinom cdf of the quantile");
    let q3, kq3 = statistics::logser_ppf(0.5, 0.5);
    chk_code(bad, "logser ppf code", kq3, 0);
    let r3, kr3 = statistics::logser_cdf(q3, 0.5);
    check(bad, kr3 == 0 && r3 >= 0.5 - 1e-9, "logser cdf of the quantile");
    let q4, kq4 = statistics::dlaplace_ppf(0.25, 0.5);
    chk_code(bad, "dlaplace ppf code", kq4, 0);
    let r4, kr4 = statistics::dlaplace_cdf(q4, 0.5);
    check(bad, kr4 == 0 && r4 >= 0.25 - 1e-9, "dlaplace cdf of the quantile");
    let q5, kq5 = statistics::nhypergeom_ppf(0.7, 10.0, 5.0, 3.0);
    chk_code(bad, "nhypergeom ppf code", kq5, 0);
    let r5, kr5 = statistics::nhypergeom_cdf(q5, 10.0, 5.0, 3.0);
    check(bad, kr5 == 0 && r5 >= 0.7 - 1e-9, "nhypergeom cdf of the quantile");
    // The quantile of the cumulative function of a count gives that count
    // back. The search stops at the first count that reaches the
    // probability, so a count with a mass of zero gives the next one.
    let t1, kt1 = statistics::nhypergeom_cdf(0.0, 10.0, 5.0, 3.0);
    let u1, ku1 = statistics::nhypergeom_ppf(t1, 10.0, 5.0, 3.0);
    chk_val(bad, "nhypergeom ppf of its own cdf", u1, ku1, 0.0, 1e-12);
    let t2, kt2 = statistics::dlaplace_cdf(0.0, 0.5);
    let u2, ku2 = statistics::dlaplace_ppf(t2, 0.5);
    chk_val(bad, "dlaplace ppf of its own cdf", u2, ku2, 0.0, 1e-12);
    // The quantile of the end points gives the bounds of the support.
    let e1, ce1 = statistics::bernoulli_ppf(0.0, 0.5);
    chk_val(bad, "bernoulli ppf 0", e1, ce1, 0.0, 1e-12);
    let e2, ce2 = statistics::bernoulli_ppf(1.0, 0.5);
    chk_val(bad, "bernoulli ppf 1", e2, ce2, 1.0, 1e-12);
    let e3, ce3 = statistics::betabinom_ppf(0.0, 10.0, 2.0, 3.0);
    chk_val(bad, "betabinom ppf 0", e3, ce3, 0.0, 1e-12);
    let e4, ce4 = statistics::betabinom_ppf(1.0, 10.0, 2.0, 3.0);
    chk_val(bad, "betabinom ppf 1", e4, ce4, 10.0, 1e-12);
    // The law of the negative hypergeometric with no success holds all its
    // mass at 0.
    let z1, cz1 = statistics::nhypergeom_pmf(0.0, 10.0, 5.0, 0.0);
    chk_val(bad, "nhypergeom pmf with no success", z1, cz1, 1.0, 1e-12);
    let z2, cz2 = statistics::nhypergeom_pmf(1.0, 10.0, 5.0, 0.0);
    chk_val(bad, "nhypergeom pmf with no success", z2, cz2, 0.0, 1e-12);
    // A shape outside its range gives the domain fault.
    let _d1, cd1 = statistics::betabinom_pmf(1.0, 10.0, 0.0, 3.0);
    chk_code(bad, "betabinom shape 0", cd1, 5);
    let _d2, cd2 = statistics::betabinom_pmf(1.0, 10.0, 2.0, -1.0);
    chk_code(bad, "betabinom shape below 0", cd2, 5);
    let _d3, cd3 = statistics::betabinom_pmf(1.0, 10.5, 2.0, 3.0);
    chk_code(bad, "betabinom count not whole", cd3, 5);
    let _d4, cd4 = statistics::logser_pmf(1.0, 1.0);
    chk_code(bad, "logser shape 1", cd4, 5);
    let _d5, cd5 = statistics::logser_pmf(1.0, 0.0);
    chk_code(bad, "logser shape 0", cd5, 5);
    let _d6, cd6 = statistics::dlaplace_pmf(1.0, 0.0);
    chk_code(bad, "dlaplace shape 0", cd6, 5);
    let _d7, cd7 = statistics::boltzmann_pmf(1.0, 2.0, 0.0);
    chk_code(bad, "boltzmann count 0", cd7, 5);
    let _d8, cd8 = statistics::nhypergeom_pmf(1.0, 10.0, 20.0, 3.0);
    chk_code(bad, "nhypergeom count above the population", cd8, 5);
    let _d9, cd9 = statistics::nhypergeom_pmf(1.0, 10.0, 5.0, 8.0);
    chk_code(bad, "nhypergeom success above the count", cd9, 5);
    return bad;
}
