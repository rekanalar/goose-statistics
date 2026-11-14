// The rows of the mixture that the vector file does not hold: the bounds
// of the weights and the round trip of the quantile.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_mixture_tests() -> i64 {
    var bad: i64 = 0;
    let comps: statistics::Dist[2] = [statistics::Dist { kind: 1, a: 0.0, b: 1.0 },
                                      statistics::Dist { kind: 1, a: 3.0, b: 1.0 }];
    let w: f64[2] = [0.3, 0.7];
    // The cumulative function of the quantile gives the probability back.
    let q1, kq1 = statistics::mixture_icdf(comps, w, 0.4);
    let r1, kr1 = statistics::mixture_cdf(comps, w, q1);
    chk_val(bad, "the cdf of the quantile", r1, kq1 + kr1, 0.4, 1e-9);
    // The support of the mixture is the span of the components.
    let lo, hi, c1 = statistics::mixture_support(comps, w);
    check(bad, c1 == 0 && lo < 0.0 && hi > 3.0, "the support of the mixture");
    // The mean of the mixture is the weighted mean of the components.
    let m1, cm1 = statistics::mixture_mean(comps, w);
    chk_val(bad, "the mean of the mixture", m1, cm1, 2.0999999999999996, 1e-9);
    // A discrete component gives the domain fault.
    let mixed: statistics::Dist[2] = [statistics::Dist { kind: 4, a: 10.0, b: 0.5 },
                                    statistics::Dist { kind: 1, a: 0.0, b: 1.0 }];
    let _a1, b1 = statistics::mixture_pdf(mixed, w, 0.0);
    chk_code(bad, "a discrete component", b1, 5);
    // An empty mixture gives the empty fault.
    let empty: statistics::Dist[0] = [];
    let we: f64[0] = [];
    let _a2, b2 = statistics::mixture_pdf(empty, we, 0.0);
    chk_code(bad, "an empty mixture", b2, 1);
    return bad;
}
