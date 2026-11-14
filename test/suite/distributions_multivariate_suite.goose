// The rows of the multivariate families that the vector file does not
// hold: the moments, the bounds of the shapes and the sum of the draws.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_multivariate_tests() -> i64 {
    var bad: i64 = 0;
    // The density of the normal law of a dimension.
    let v1, c1 = statistics::multivariate_normal_pdf([0.5, 1.5], [1.0, 2.0], [2.0, 0.5, 0.5, 1.0], 2);
    chk_val(bad, "multivariate normal pdf", v1, c1, 0.10429393134976603, 1e-9);
    // The mean and the covariance of the law.
    var m1: f64[>..] = [];
    let k1 = statistics::multivariate_normal_mean(m1, [1.0, 2.0], [2.0, 0.5, 0.5, 1.0], 2);
    chk_val(bad, "multivariate normal mean", m1[0], k1, 1.0, 1e-12);
    // The mass of the multinomial law.
    let v2, c2 = statistics::multinomial_pmf([2.0, 1.0, 1.0], 4, [0.5, 0.3, 0.2]);
    chk_val(bad, "multinomial pmf", v2, c2, 0.17999999999999997, 1e-9);
    // The mass of the multivariate hypergeometric law.
    let v3, c3 = statistics::multivariate_hypergeom_pmf([2.0, 1.0, 1.0], [5.0, 4.0, 3.0], 4.0);
    chk_val(bad, "multivariate hypergeom pmf", v3, c3, 0.24242424242424238, 1e-9);
    // The density of the Dirichlet law.
    let v4, c4 = statistics::dirichlet_pdf([0.5, 0.3, 0.2], [2.0, 3.0, 1.5]);
    chk_val(bad, "dirichlet pdf", v4, c4, 3.2686771817967584, 1e-9);
    var v5: f64[>..] = [];
    let k5 = statistics::dirichlet_mean(v5, [2.0, 3.0, 1.5]);
    chk_val(bad, "dirichlet mean", v5[0], k5, 0.3076923076923077, 1e-9);
    // The log density of the von Mises-Fisher law.
    let v6, c6 = statistics::vonmises_fisher_logpdf([0.0, 0.6, 0.8], [0.0, 0.0, 1.0], 2.0);
    chk_val(bad, "vonmises fisher logpdf", v6, c6, -1.5262444390235141, 1e-9);
    let v7, c7 = statistics::vonmises_fisher_logpdf([0.6, 0.8], [1.0, 0.0], 3.0);
    chk_val(bad, "vonmises fisher logpdf of a dimension 2", v7, c7, -1.6231846882227665, 1e-9);
    // The draws of each law give a finite value.
    var g = rng { s: 21 };
    var out: f64[>..] = [];
    let k8 = statistics::multivariate_normal_rvs(out, [1.0, 2.0], [2.0, 0.5, 0.5, 1.0], 2, 5, g);
    check(bad, k8 == 0 && out.len == 10, "the normal draw gives its values");
    let k9 = statistics::dirichlet_rvs(out, [2.0, 3.0, 1.5], 5, g);
    check(bad, k9 == 0 && out.len == 15, "the dirichlet draw gives its values");
    let k10 = statistics::multinomial_rvs(out, 10, [0.5, 0.3, 0.2], 5, g);
    check(bad, k10 == 0 && out.len == 15, "the multinomial draw gives its values");
    let k11 = statistics::wishart_rvs(out, [2.0, 0.5, 0.5, 1.0], 5.0, 2, 3, g);
    check(bad, k11 == 0 && out.len == 12, "the wishart draw gives its values");
    // A shape outside its range gives the domain fault.
    let _d1, cd1 = statistics::multivariate_normal_pdf([0.5, 1.5], [1.0, 2.0], [2.0, 0.5, 0.5, -1.0], 2);
    chk_code(bad, "a covariance that is not positive definite", cd1, 5);
    let _d2, cd2 = statistics::multinomial_pmf([2.0, 1.0, 2.0], 4, [0.5, 0.3, 0.2]);
    chk_code(bad, "a count that does not sum to the total", cd2, 5);
    let _d3, cd3 = statistics::dirichlet_pdf([0.5, 0.3, 0.3], [2.0, 3.0, 1.5]);
    chk_code(bad, "a point off the simplex", cd3, 5);
    // A list of another length than the call needs gives the length fault.
    let _d5, cd5 = statistics::multinomial_pmf([2.0, 1.0], 4, [0.5, 0.3, 0.2]);
    chk_code(bad, "a count list of another length", cd5, 6);
    let _d6, cd6 = statistics::dirichlet_pdf([0.5, 0.6], [2.0, 3.0, 1.5]);
    chk_code(bad, "a point list of another length", cd6, 6);
    let _d4, cd4 = statistics::vonmises_fisher_logpdf([0.0, 0.6, 0.8], [0.0, 0.0, 2.0], 2.0);
    chk_code(bad, "a mean that is not a unit vector", cd4, 5);
    return bad;
}
