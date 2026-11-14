// Time rows of the families of the random matrices.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_random_matrix_time() -> i64 {
    var bad: i64 = 0;
    var total: f64 = 0.0;
    var g = rng { s: 99 };
    var out: f64[>..] = [];
    for k in 0..5 {
        let c = statistics::ortho_group_rvs(out, 4, 2, g);
        if c != 0 { bad += 1; }
        total += out[0];
    }
    let c2 = statistics::unitary_group_rvs(out, 4, 2, g);
    if c2 != 0 { bad += 1; }
    let c3 = statistics::vonmises_fisher_rvs(out, [0.0, 0.0, 1.0], 2.0, 5, g);
    if c3 != 0 { bad += 1; }
    total += out[0];
    check(bad, !is_nan(total) && !is_inf(total), "the sum of the calls is finite");
    return bad;
}
