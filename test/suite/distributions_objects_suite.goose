// The rows of the distribution objects that the vector file does not
// hold: the bounds of the shapes, the support and the fault codes.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_objects_tests() -> i64 {
    var bad: i64 = 0;
    // The support of each kind.
    let normal = statistics::Dist { kind: 1, a: 1.0, b: 2.0 };
    let uniform = statistics::Dist { kind: 2, a: 0.0, b: 2.0 };
    let logistic = statistics::Dist { kind: 3, a: 0.0, b: 1.0 };
    let binomial = statistics::Dist { kind: 4, a: 10.0, b: 0.5 };
    let lo1, hi1, c1 = statistics::dist_support(normal);
    check(bad, c1 == 0 && lo1 < 0.0 && hi1 > 0.0, "the support of the normal law");
    let lo2, hi2, c2 = statistics::dist_support(uniform);
    check(bad, c2 == 0 && lo2 == 0.0 && hi2 == 2.0, "the support of the uniform law");
    let lo3, hi3, c3 = statistics::dist_support(binomial);
    check(bad, c3 == 0 && lo3 == 0.0 && hi3 == 10.0, "the support of the binomial law");
    // The cumulative function of the quantile gives the probability back.
    let q1, kq1 = statistics::dist_icdf(normal, 0.3);
    let r1, kr1 = statistics::dist_cdf(normal, q1);
    chk_val(bad, "the cdf of the quantile", r1, kq1 + kr1, 0.3, 1e-9);
    let q2, kq2 = statistics::dist_icdf(uniform, 0.75);
    let r2, kr2 = statistics::dist_cdf(uniform, q2);
    chk_val(bad, "the cdf of the quantile", r2, kq2 + kr2, 0.75, 1e-9);
    let q3, kq3 = statistics::dist_icdf(binomial, 0.4);
    let r3, kr3 = statistics::dist_cdf(binomial, q3);
    // A discrete law gives the first point where the function reaches the
    // probability, so the function there is at least the probability.
    check(bad, kq3 + kr3 == 0 && r3 >= 0.4, "the cdf of the discrete quantile");
    // The complement of the cumulative function is the survival function.
    let s1, ks1 = statistics::dist_ccdf(normal, 1.0);
    let t1, kt1 = statistics::dist_cdf(normal, 1.0);
    chk_val(bad, "the survival function", s1, ks1 + kt1, 1.0 - t1, 1e-12);
    // A probability outside its range gives the domain fault.
    let _d1, cd1 = statistics::dist_icdf(normal, 1.5);
    chk_code(bad, "a probability above one", cd1, 5);
    let _d2, cd2 = statistics::dist_logpdf(normal, 0.0 / 0.0);
    chk_code(bad, "a NaN point", cd2, 3);
    // A shape outside its range gives the domain fault.
    let _n1, cn1 = statistics::Normal(0.0, -1.0);
    chk_code(bad, "a negative scale", cn1, 5);
    let _n2, cn2 = statistics::Uniform(2.0, 1.0);
    chk_code(bad, "a range that ends before it starts", cn2, 5);
    let _n3, cn3 = statistics::Binomial(-1.0, 0.5);
    chk_code(bad, "a negative count", cn3, 5);
    let _n4, cn4 = statistics::Binomial(10.0, 1.5);
    chk_code(bad, "a chance above one", cn4, 5);
    // A kind that the package does not know gives the domain fault.
    let _d3, cd3 = statistics::dist_pdf(statistics::Dist { kind: 9, a: 0.0, b: 1.0 }, 0.0);
    chk_code(bad, "an unknown kind", cd3, 5);
    // The logistic law has a mean of zero and a median of zero.
    let m1, cm1 = statistics::dist_mean(logistic);
    chk_val(bad, "the mean of the logistic law", m1, cm1, 0.0, 1e-12);
    return bad;
}
