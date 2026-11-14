// The rows of the covariance forms that the vector file does not hold:
// the fault codes and the round trip of the whitening.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_covariance_tests() -> i64 {
    var bad: i64 = 0;
    var out: f64[>..] = [];
    let x: f64[3] = [1.0, -2.0, 0.5];
    let dg: f64[3] = [1.0, 2.0, 3.0];
    // The whitening and the colouring cancel each other.
    let k1 = statistics::covariance_whiten(out, 1, dg, 3, x);
    var w: f64[>..] = [];
    for v in out { w.push(v); }
    let k2 = statistics::covariance_colorize(out, 1, dg, 3, w[..]);
    check(bad, k1 == 0 && k2 == 0, "the whitening and the colouring give their values");
    for i in 0..3 {
        check(bad, abs(out[i] - x[i]) < 1e-10, "the round trip of the whitening");
    }
    // A zero on the diagonal drops out of the log determinant and of the
    // rank, and the whitening gives zero for it.
    let sing: f64[3] = [1.0, 0.0, 3.0];
    let v1, c1 = statistics::covariance_log_pdet(1, sing, 3);
    chk_val(bad, "the log pdet of a singular form", v1, c1, log(3.0), 1e-12);
    let r1, cr1 = statistics::covariance_rank(1, sing, 3);
    chk_val(bad, "the rank of a singular form", r1 as f64, cr1, 2.0, 1e-12);
    // A dimension of zero gives the domain fault.
    let _v3, k3 = statistics::covariance_log_pdet(1, dg, 0);
    chk_code(bad, "a dimension of zero", k3, 5);
    return bad;
}
