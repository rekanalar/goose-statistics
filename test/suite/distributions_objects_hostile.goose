// Hostile input for the distribution objects, the mixture, the
// covariance forms and the histogram law.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_objects_hostile() -> i64 {
    var bad: i64 = 0;
    var out: f64[>..] = [];
    let n = statistics::Dist { kind: 1, a: 0.0, b: 1.0 };
    let a1, b1 = statistics::dist_pdf(n, 0.0 / 0.0);
    chk_code(bad, "a NaN point", b1, 3);
    let a2, b2 = statistics::dist_ilogcdf(n, 1.0);
    chk_code(bad, "a log probability above zero", b2, 5);
    let a3, b3 = statistics::dist_moment(n, -1);
    chk_code(bad, "a negative order", b3, 5);
    let a4, b4 = statistics::dist_moment(n, 100000);
    chk_code(bad, "an order above the bound", b4, 5);
    let _n1, c1 = statistics::Normal(0.0 / 0.0, 1.0);
    chk_code(bad, "a NaN location", c1, 3);
    let _n2, c2 = statistics::Normal(0.0, 1.0 / 0.0);
    chk_code(bad, "an infinite scale", c2, 4);
    // The mixture checks its weights.
    let comps: statistics::Dist[2] = [n, n];
    let w1: f64[2] = [0.5, 0.6];
    let a5, b5 = statistics::mixture_pdf(comps, w1, 0.0);
    chk_code(bad, "weights that do not sum to one", b5, 5);
    let w2: f64[1] = [1.0];
    let a6, b6 = statistics::mixture_pdf(comps, w2, 0.0);
    chk_code(bad, "weights of another length", b6, 6);
    let w3: f64[2] = [-0.5, 1.5];
    let a7, b7 = statistics::mixture_pdf(comps, w3, 0.0);
    chk_code(bad, "a negative weight", b7, 5);
    // The covariance forms check their buffers.
    let dg: f64[3] = [1.0, 2.0, 3.0];
    let x3: f64[3] = [1.0, 2.0, 3.0];
    let k1 = statistics::covariance_matrix(out, 1, dg, 4);
    chk_code(bad, "a buffer of the wrong length", k1, 6);
    let k2 = statistics::covariance_matrix(out, 9, dg, 3);
    chk_code(bad, "an unknown kind", k2, 5);
    let k3 = statistics::covariance_whiten(out, 1, dg, 3, [1.0, 2.0]);
    chk_code(bad, "a point of the wrong length", k3, 6);
    let bad_chol: f64[4] = [1.0, 0.5, 0.0, 1.0];
    let k4 = statistics::covariance_matrix(out, 3, bad_chol, 2);
    chk_code(bad, "a factor that is not lower triangular", k4, 5);
    let a8, b8 = statistics::covariance_log_pdet(1, dg, 0);
    chk_code(bad, "a dimension of zero", b8, 5);
    // The histogram law checks its bins.
    let edges: f64[3] = [0.0, 1.0, 1.0];
    let content: f64[2] = [1.0, 1.0];
    let a9, b9 = statistics::histogram_pdf(content, edges, true, 0.5);
    chk_code(bad, "edges that do not rise", b9, 5);
    let e2: f64[3] = [0.0, 1.0, 2.0];
    let c3: f64[2] = [0.0, 0.0];
    let a10, b10 = statistics::histogram_pdf(c3, e2, true, 0.5);
    chk_code(bad, "a histogram with no mass", b10, 5);
    let c4: f64[2] = [1.0, -1.0];
    let a11, b11 = statistics::histogram_pdf(c4, e2, true, 0.5);
    chk_code(bad, "a negative bin", b11, 5);
    return bad;
}
