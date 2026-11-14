// The rows of the Weibull families that the vector file does not hold: the
// support of each law, the end points of quantile and the shape bounds.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_weibull_tests() -> i64 {
    var bad: i64 = 0;
    // A value outside the support gives a density of 0 and not a fault.
    let v3, c3 = statistics::weibull_max_pdf(1.0, 1.5);
    chk_val(bad, "weibull_max pdf above the support", v3, c3, 0.0, 1e-12);
    let v4, c4 = statistics::truncweibull_min_pdf(3.0, 1.5, 1.0, 2.0);
    chk_val(bad, "truncweibull_min pdf above the support", v4, c4, 0.0, 1e-12);
    // The cumulative function of a value below the support is 0.
    let v6, c6 = statistics::invweibull_cdf(-1.0, 1.5);
    chk_val(bad, "invweibull cdf below the support", v6, c6, 0.0, 1e-12);
    let v7, c7 = statistics::truncweibull_min_cdf(0.5, 1.5, 1.0, 2.0);
    chk_val(bad, "truncweibull_min cdf below the support", v7, c7, 0.0, 1e-12);
    // The cumulative function of a value above the support is 1.
    let v8, c8 = statistics::truncweibull_min_cdf(3.0, 1.5, 1.0, 2.0);
    chk_val(bad, "truncweibull_min cdf above the support", v8, c8, 1.0, 1e-12);
    let v9, c9 = statistics::weibull_max_cdf(1.0, 1.5);
    chk_val(bad, "weibull_max cdf above the support", v9, c9, 1.0, 1e-12);
    // The quantile of the end points gives the bounds of the support.
    let e1, ce1 = statistics::truncweibull_min_ppf(0.0, 1.5, 1.0, 2.0);
    chk_val(bad, "truncweibull_min ppf 0", e1, ce1, 1.0, 1e-12);
    let e2, ce2 = statistics::truncweibull_min_ppf(1.0, 1.5, 1.0, 2.0);
    chk_val(bad, "truncweibull_min ppf 1", e2, ce2, 2.0, 1e-12);
    let e3, ce3 = statistics::weibull_max_ppf(1.0, 1.5);
    chk_val(bad, "weibull_max ppf 1", e3, ce3, 0.0, 1e-12);
    // The cumulative function of the quantile of a probability is that
    // probability.
    let q1, kq1 = statistics::invweibull_ppf(0.3, 1.5);
    chk_code(bad, "invweibull ppf code", kq1, 0);
    let r1, kr1 = statistics::invweibull_cdf(q1, 1.5);
    chk_val(bad, "invweibull cdf of the quantile", r1, kr1, 0.3, 1e-9);
    let q3, kq3 = statistics::truncweibull_min_ppf(0.6, 1.5, 1.0, 2.0);
    chk_code(bad, "truncweibull_min ppf code", kq3, 0);
    let r3, kr3 = statistics::truncweibull_min_cdf(q3, 1.5, 1.0, 2.0);
    chk_val(bad, "truncweibull_min cdf of the quantile", r3, kr3, 0.6, 1e-9);
    // A shape of 0 or below gives the domain fault.
    let _z1, cz1 = statistics::dweibull_pdf(1.0, 0.0);
    chk_code(bad, "dweibull shape 0", cz1, 5);
    let _z2, cz2 = statistics::invweibull_pdf(1.0, -1.0);
    chk_code(bad, "invweibull shape below 0", cz2, 5);
    // A bound that is not above the lower bound gives the domain fault.
    let _z5, cz5 = statistics::truncweibull_min_pdf(1.0, 1.5, 1.0, 1.0);
    chk_code(bad, "truncweibull_min bounds equal", cz5, 5);
    let _z6, cz6 = statistics::truncweibull_min_pdf(1.0, 1.5, 2.0, 1.0);
    chk_code(bad, "truncweibull_min bounds below", cz6, 5);
    let _z7, cz7 = statistics::truncweibull_min_pdf(1.0, 1.5, -1.0, 2.0);
    chk_code(bad, "truncweibull_min lower bound below 0", cz7, 5);
    // The law of the double Weibull is symmetric around 0.
    let s1, cs1 = statistics::dweibull_cdf(0.5, 1.7);
    let s2, cs2 = statistics::dweibull_cdf(-0.5, 1.7);
    chk_code(bad, "dweibull cdf code", cs1, 0);
    chk_code(bad, "dweibull cdf code", cs2, 0);
    check(bad, abs(s1 + s2 - 1.0) < 1e-12, "the double Weibull law is symmetric");
    // The Weibull law of the negative line is the mirror of the positive
    // half of the double Weibull law.
    let w1, cw1 = statistics::weibull_max_cdf(-0.7, 1.5);
    let w2, cw2 = statistics::dweibull_cdf(0.7, 1.5);
    chk_code(bad, "weibull_max cdf code", cw1, 0);
    chk_code(bad, "dweibull cdf code", cw2, 0);
    check(bad, abs(w2 - (1.0 - 0.5 * w1)) < 1e-12, "the Weibull law of the negative line is the mirror");
    // The mean of the double Weibull law is 0 and its law is symmetric.
    let m1, cm1 = statistics::dweibull_mean(1.7);
    chk_code(bad, "dweibull mean code", cm1, 0);
    chk_val(bad, "dweibull mean", m1, cm1, 0.0, 1e-12);
    let k1, ck1 = statistics::dweibull_skew(1.7);
    chk_code(bad, "dweibull skew code", ck1, 0);
    chk_val(bad, "dweibull skew", k1, ck1, 0.0, 1e-12);
    return bad;
}
