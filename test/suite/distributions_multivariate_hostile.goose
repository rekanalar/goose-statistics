// Hostile input for the multivariate families. A NaN gives the NaN
// fault, an infinity gives the infinity fault and a shape outside its
// range gives the domain fault.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_multivariate_hostile() -> i64 {
    var bad: i64 = 0;
    let a1, b1 = statistics::multivariate_normal_pdf([0.0 / 0.0, 1.0], [1.0, 2.0], [2.0, 0.5, 0.5, 1.0], 2);
    chk_code(bad, "a NaN point", b1, 3);
    let a2, b2 = statistics::multivariate_normal_pdf([0.5, 1.5], [1.0, 2.0], [2.0, 0.5, 0.5, 1.0 / 0.0], 2);
    chk_code(bad, "an infinite covariance", b2, 4);
    let a3, b3 = statistics::multivariate_t_df(0.0 / 0.0);
    chk_code(bad, "a NaN degrees of freedom", b3, 3);
    let a4, b4 = statistics::multinomial_pmf([2.0, 1.0, 1.0], 4, [0.5, 0.3, 0.0 / 0.0]);
    chk_code(bad, "a NaN chance", b4, 3);
    let a5, b5 = statistics::dirichlet_pdf([0.5, 0.3, 0.2], [2.0, 0.0 / 0.0, 1.5]);
    chk_code(bad, "a NaN shape", b5, 3);
    let a6, b6 = statistics::multivariate_hypergeom_pmf([2.0, 1.0, 1.0], [5.0, 4.0, 3.0], -1.0);
    chk_code(bad, "a negative total", b6, 5);
    return bad;
}
