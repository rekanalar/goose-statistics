// The rows of the fourteen laws of the first tail group that the vector file does not hold: the round
// trip of the quantile and the bounds of the shapes.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_closed_tests() -> i64 {
    var bad: i64 = 0;
    let q0, kq0 = statistics::alpha_ppf(0.35, 1.5);
    let r0, kr0 = statistics::alpha_cdf(q0, 1.5);
    chk_val(bad, "alpha cdf of the quantile", r0, kq0 + kr0, 0.35, 1e-9);
    let q1, kq1 = statistics::fatiguelife_ppf(0.35, 0.7);
    let r1, kr1 = statistics::fatiguelife_cdf(q1, 0.7);
    chk_val(bad, "fatiguelife cdf of the quantile", r1, kq1 + kr1, 0.35, 1e-9);
    let q2, kq2 = statistics::levy_ppf(0.35);
    let r2, kr2 = statistics::levy_cdf(q2);
    chk_val(bad, "levy cdf of the quantile", r2, kq2 + kr2, 0.35, 1e-9);
    let q3, kq3 = statistics::levy_l_ppf(0.35);
    let r3, kr3 = statistics::levy_l_cdf(q3);
    chk_val(bad, "levy_l cdf of the quantile", r3, kq3 + kr3, 0.35, 1e-9);
    let q4, kq4 = statistics::logistic_ppf(0.35);
    let r4, kr4 = statistics::logistic_cdf(q4);
    chk_val(bad, "logistic cdf of the quantile", r4, kq4 + kr4, 0.35, 1e-9);
    let q5, kq5 = statistics::loglaplace_ppf(0.35, 2.5);
    let r5, kr5 = statistics::loglaplace_cdf(q5, 2.5);
    chk_val(bad, "loglaplace cdf of the quantile", r5, kq5 + kr5, 0.35, 1e-9);
    let q6, kq6 = statistics::lomax_ppf(0.35, 2.0);
    let r6, kr6 = statistics::lomax_cdf(q6, 2.0);
    chk_val(bad, "lomax cdf of the quantile", r6, kq6 + kr6, 0.35, 1e-9);
    let q7, kq7 = statistics::pareto_ppf(0.35, 2.5);
    let r7, kr7 = statistics::pareto_cdf(q7, 2.5);
    chk_val(bad, "pareto cdf of the quantile", r7, kq7 + kr7, 0.35, 1e-9);
    let q8, kq8 = statistics::powerlaw_ppf(0.35, 1.8);
    let r8, kr8 = statistics::powerlaw_cdf(q8, 1.8);
    chk_val(bad, "powerlaw cdf of the quantile", r8, kq8 + kr8, 0.35, 1e-9);
    let q9, kq9 = statistics::reciprocal_ppf(0.35, 1.0, 4.0);
    let r9, kr9 = statistics::reciprocal_cdf(q9, 1.0, 4.0);
    chk_val(bad, "reciprocal cdf of the quantile", r9, kq9 + kr9, 0.35, 1e-9);
    let q10, kq10 = statistics::rdist_ppf(0.35, 3.0);
    let r10, kr10 = statistics::rdist_cdf(q10, 3.0);
    chk_val(bad, "rdist cdf of the quantile", r10, kq10 + kr10, 0.35, 1e-9);
    let q11, kq11 = statistics::skewcauchy_ppf(0.35, 0.4);
    let r11, kr11 = statistics::skewcauchy_cdf(q11, 0.4);
    chk_val(bad, "skewcauchy cdf of the quantile", r11, kq11 + kr11, 0.35, 1e-9);
    let q12, kq12 = statistics::wrapcauchy_ppf(0.35, 0.6);
    let r12, kr12 = statistics::wrapcauchy_cdf(q12, 0.6);
    chk_val(bad, "wrapcauchy cdf of the quantile", r12, kq12 + kr12, 0.35, 1e-9);
    let q13, kq13 = statistics::kappa4_ppf(0.35, 0.3, 0.2);
    let r13, kr13 = statistics::kappa4_cdf(q13, 0.3, 0.2);
    chk_val(bad, "kappa4 cdf of the quantile", r13, kq13 + kr13, 0.35, 1e-9);
    // A shape outside its range gives the domain fault.
    let _d1, cd1 = statistics::alpha_pdf(1.0, 0.0);
    chk_code(bad, "alpha shape 0", cd1, 5);
    let _d2, cd2 = statistics::skewcauchy_pdf(1.0, 1.0);
    chk_code(bad, "skewcauchy shape 1", cd2, 5);
    let _d3, cd3 = statistics::wrapcauchy_pdf(1.0, 1.0);
    chk_code(bad, "wrapcauchy shape 1", cd3, 5);
    // The support of each bounded law.
    let s1, k1 = statistics::pareto_pdf(0.5, 2.5);
    chk_val(bad, "pareto pdf below 1", s1, k1, 0.0, 1e-12);
    let s2, k2 = statistics::powerlaw_cdf(2.0, 1.8);
    chk_val(bad, "powerlaw cdf above 1", s2, k2, 1.0, 1e-12);
    let s3, k3 = statistics::levy_cdf(-1.0);
    chk_val(bad, "levy cdf below 0", s3, k3, 0.0, 1e-12);
    let s4, k4 = statistics::levy_l_cdf(1.0);
    chk_val(bad, "levy_l cdf above 0", s4, k4, 1.0, 1e-12);
    let s5, k5 = statistics::wrapcauchy_ppf(0.5, 0.6);
    chk_val(bad, "wrapcauchy median", s5, k5, 3.141592653589793, 1e-12);
    return bad;
}
