// The rows of the second group of discrete families that the vector file
// does not hold: the support of each law, the bounds of the shapes and
// the round trip of the quantile.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_discrete_b_tests() -> i64 {
    var bad: i64 = 0;
    // A count outside the support gives a mass of 0 and not a fault.
    let v1, c1 = statistics::planck_pmf(-1.0, 2.0);
    chk_val(bad, "planck pmf below the support", v1, c1, 0.0, 1e-12);
    let v2, c2 = statistics::randint_pmf(5.0, 0.0, 5.0);
    chk_val(bad, "randint pmf above the support", v2, c2, 0.0, 1e-12);
    let v3, c3 = statistics::yulesimon_pmf(0.0, 2.0);
    chk_val(bad, "yulesimon pmf below the support", v3, c3, 0.0, 1e-12);
    let v4, c4 = statistics::zipf_pmf(0.0, 2.0);
    chk_val(bad, "zipf pmf below the support", v4, c4, 0.0, 1e-12);
    let v5, c5 = statistics::zipfian_pmf(6.0, 2.0, 5.0);
    chk_val(bad, "zipfian pmf above the support", v5, c5, 0.0, 1e-12);
    // The cumulative function of a count below the support is 0 and above
    // the support is 1.
    let v6, c6 = statistics::planck_cdf(-1.0, 2.0);
    chk_val(bad, "planck cdf below the support", v6, c6, 0.0, 1e-12);
    let v7, c7 = statistics::randint_cdf(-1.0, 0.0, 5.0);
    chk_val(bad, "randint cdf below the support", v7, c7, 0.0, 1e-12);
    let v8, c8 = statistics::zipfian_cdf(10.0, 2.0, 5.0);
    chk_val(bad, "zipfian cdf above the support", v8, c8, 1.0, 1e-12);
    // The cumulative function of the quantile of a probability is that
    // probability.
    let q1, kq1 = statistics::planck_ppf(0.4, 2.0);
    chk_code(bad, "planck ppf code", kq1, 0);
    let r1, kr1 = statistics::planck_cdf(q1, 2.0);
    check(bad, kr1 == 0 && r1 >= 0.4 - 1e-9, "planck cdf of the quantile");
    let q2, kq2 = statistics::randint_ppf(0.3, 0.0, 5.0);
    chk_code(bad, "randint ppf code", kq2, 0);
    let r2, kr2 = statistics::randint_cdf(q2, 0.0, 5.0);
    check(bad, kr2 == 0 && r2 >= 0.3 - 1e-9, "randint cdf of the quantile");
    let q3, kq3 = statistics::yulesimon_ppf(0.6, 2.0);
    chk_code(bad, "yulesimon ppf code", kq3, 0);
    let r3, kr3 = statistics::yulesimon_cdf(q3, 2.0);
    check(bad, kr3 == 0 && r3 >= 0.6 - 1e-9, "yulesimon cdf of the quantile");
    let q4, kq4 = statistics::zipf_ppf(0.8, 2.0);
    chk_code(bad, "zipf ppf code", kq4, 0);
    let r4, kr4 = statistics::zipf_cdf(q4, 2.0);
    check(bad, kr4 == 0 && r4 >= 0.8 - 1e-9, "zipf cdf of the quantile");
    let q5, kq5 = statistics::zipfian_ppf(0.4, 2.0, 5.0);
    chk_code(bad, "zipfian ppf code", kq5, 0);
    let r5, kr5 = statistics::zipfian_cdf(q5, 2.0, 5.0);
    check(bad, kr5 == 0 && r5 >= 0.4 - 1e-9, "zipfian cdf of the quantile");
    let q6, kq6 = statistics::skellam_ppf(0.3, 1.0, 2.0);
    chk_code(bad, "skellam ppf code", kq6, 0);
    let r6, kr6 = statistics::skellam_cdf(q6, 1.0, 2.0);
    check(bad, kr6 == 0 && r6 >= 0.3 - 1e-9, "skellam cdf of the quantile");
    // The quantile of the end points gives the bounds of the support.
    let e1, ce1 = statistics::planck_ppf(0.0, 2.0);
    chk_val(bad, "planck ppf 0", e1, ce1, 0.0, 1e-12);
    let e2, ce2 = statistics::randint_ppf(0.0, 0.0, 5.0);
    chk_val(bad, "randint ppf 0", e2, ce2, 0.0, 1e-12);
    let e3, ce3 = statistics::randint_ppf(1.0, 0.0, 5.0);
    chk_val(bad, "randint ppf 1", e3, ce3, 4.0, 1e-12);
    let e4, ce4 = statistics::zipfian_ppf(0.0, 2.0, 5.0);
    chk_val(bad, "zipfian ppf 0", e4, ce4, 1.0, 1e-12);
    let e5, ce5 = statistics::zipfian_ppf(1.0, 2.0, 5.0);
    chk_val(bad, "zipfian ppf 1", e5, ce5, 5.0, 1e-12);
    // The law of the Poisson binomial with one probability is the law of
    // Bernoulli.
    let ps: f64[1] = [0.4];
    let v7, c7 = statistics::poisson_binom_pmf(0.0, ps);
    chk_val(bad, "poisson_binom pmf of 0", v7, c7, 0.6, 1e-12);
    let v8, c8 = statistics::poisson_binom_pmf(1.0, ps);
    chk_val(bad, "poisson_binom pmf of 1", v8, c8, 0.4, 1e-12);
    // A shape outside its range gives the domain fault.
    let _d1, cd1 = statistics::planck_pmf(1.0, 0.0);
    chk_code(bad, "planck shape 0", cd1, 5);
    let _d2, cd2 = statistics::yulesimon_pmf(1.0, 0.0);
    chk_code(bad, "yulesimon shape 0", cd2, 5);
    let _d3, cd3 = statistics::zipf_pmf(1.0, 1.0);
    chk_code(bad, "zipf shape 1", cd3, 5);
    let _d4, cd4 = statistics::zipf_pmf(1.0, 0.5);
    chk_code(bad, "zipf shape below 1", cd4, 5);
    let _d5, cd5 = statistics::randint_pmf(1.0, 5.0, 5.0);
    chk_code(bad, "randint bound above the low bound", cd5, 5);
    let _d6, cd6 = statistics::zipfian_pmf(1.0, 2.0, 0.0);
    chk_code(bad, "zipfian count 0", cd6, 5);
    let _d7, cd7 = statistics::skellam_pmf(1.0, 0.0, 2.0);
    chk_code(bad, "skellam shape 0", cd7, 5);
    return bad;
}
