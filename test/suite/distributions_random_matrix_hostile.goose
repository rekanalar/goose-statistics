// Hostile input for the families of the random matrices.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_random_matrix_hostile() -> i64 {
    var bad: i64 = 0;
    var g = rng { s: 55 };
    var out: f64[>..] = [];
    let a1, b1 = statistics::vonmises_fisher_logpdf([0.0 / 0.0, 0.6, 0.8], [0.0, 0.0, 1.0], 2.0);
    chk_code(bad, "a NaN point", b1, 3);
    let a2, b2 = statistics::vonmises_fisher_logpdf([0.0, 0.6, 0.8], [0.0, 0.0, 1.0 / 0.0], 2.0);
    chk_code(bad, "an infinite mean", b2, 4);
    let k3 = statistics::random_correlation_rvs(out, [0.5, 0.5, 1.0 / 0.0, 1.0], g);
    chk_code(bad, "an infinite eigenvalue", k3, 4);
    let k4 = statistics::random_correlation_rvs(out, [-0.5, 0.5, 1.0, 1.0], g);
    chk_code(bad, "a negative eigenvalue", k4, 5);
    let k5 = statistics::unitary_group_rvs(out, 100, 1, g);
    chk_code(bad, "a dimension above the bound", k5, 5);
    let k6 = statistics::random_table_rvs(out, [3.0, 2.0], [2.0, 3.0], -1, g);
    chk_code(bad, "a negative count of draws", k6, 5);
    return bad;
}
