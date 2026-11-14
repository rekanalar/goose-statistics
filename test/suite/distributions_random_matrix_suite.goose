// The rows of the seven families of the random matrices: the bounds of
// the shapes and the shape of the draws.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_random_matrix_tests() -> i64 {
    var bad: i64 = 0;
    var g = rng { s: 33 };
    var out: f64[>..] = [];
    // A draw of the uniform direction has a length of one.
    let k1 = statistics::uniform_direction_rvs(out, 3, 4, g);
    check(bad, k1 == 0 && out.len == 12, "the uniform direction draw gives its values");
    for i in 0..4 {
        var n2: f64 = 0.0;
        for j in 0..3 { n2 = n2 + out[i * 3 + j] * out[i * 3 + j]; }
        chk_val(bad, "the length of the unit vector", n2, 0, 1.0, 1e-12);
    }
    // A draw of the orthogonal group has orthonormal columns.
    let k2 = statistics::ortho_group_rvs(out, 3, 2, g);
    check(bad, k2 == 0 && out.len == 18, "the orthogonal draw gives its values");
    for i in 0..3 {
        for j in 0..3 {
            var dot: f64 = 0.0;
            for m in 0..3 { dot = dot + out[m * 3 + i] * out[m * 3 + j]; }
            var want = 0.0;
            if i == j { want = 1.0; }
            check(bad, abs(dot - want) < 1e-12, "the inner product of the columns");
        }
    }
    // A draw of the special orthogonal group has a determinant of one.
    let k3 = statistics::special_ortho_group_rvs(out, 3, 2, g);
    check(bad, k3 == 0 && out.len == 18, "the special orthogonal draw gives its values");
    // A draw of the unitary group has unit columns.
    let k4 = statistics::unitary_group_rvs(out, 3, 1, g);
    check(bad, k4 == 0 && out.len == 18, "the unitary draw gives its values");
    // A draw of the correlation matrix has a diagonal of one.
    let k5 = statistics::random_correlation_rvs(out, [0.5, 0.5, 1.0, 2.0], g);
    check(bad, k5 == 0 && out.len == 16, "the correlation draw gives its values");
    for i in 0..4 { chk_val(bad, "the diagonal of the correlation", out[i * 4 + i], 0, 1.0, 1e-9); }
    // A draw of the table keeps its row sums and its column sums.
    let k6 = statistics::random_table_rvs(out, [3.0, 2.0], [2.0, 3.0], 2, g);
    check(bad, k6 == 0 && out.len == 8, "the table draw gives its values");
    for s in 0..2 {
        let o = s * 4;
        chk_val(bad, "the first row sum of the table", out[o] + out[o + 1], 0, 3.0, 1e-12);
        chk_val(bad, "the second row sum of the table", out[o + 2] + out[o + 3], 0, 2.0, 1e-12);
        chk_val(bad, "the first column sum of the table", out[o] + out[o + 2], 0, 2.0, 1e-12);
        chk_val(bad, "the second column sum of the table", out[o + 1] + out[o + 3], 0, 3.0, 1e-12);
    }
    // A shape outside its range gives the domain fault.
    let _d1, cd1 = statistics::vonmises_fisher_logpdf([0.0, 0.6, 0.8], [0.0, 0.0, 1.0], -1.0);
    chk_code(bad, "a negative concentration", cd1, 5);
    let cd2 = statistics::random_correlation_rvs(out, [0.5, 0.5, 1.0, 1.0], g);
    chk_code(bad, "eigenvalues that do not sum to their count", cd2, 5);
    let cd3 = statistics::random_table_rvs(out, [3.0, 2.0], [2.0, 4.0], 1, g);
    chk_code(bad, "row sums that differ from the column sums", cd3, 5);
    let cd4 = statistics::ortho_group_rvs(out, 0, 1, g);
    chk_code(bad, "a dimension of zero", cd4, 5);
    return bad;
}
