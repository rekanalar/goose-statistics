// Time rows for the covariance forms: a form of a large dimension must
// end in a bounded time.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_covariance_time() -> i64 {
    var bad: i64 = 0;
    var out: f64[>..] = [];
    var diag: f64[>..] = [];
    var x: f64[>..] = [];
    for i in 0..100 {
        diag.push(1.0 + (i % 5));
        x.push((i - 50.0) * 0.01);
    }
    let k1 = statistics::covariance_whiten(out, 1, diag[..], 100, x[..]);
    check(bad, k1 == 0, "the whitening of a large form");
    var w: f64[>..] = [];
    for v in out { w.push(v); }
    let k2 = statistics::covariance_colorize(out, 1, diag[..], 100, w[..]);
    check(bad, k2 == 0, "the colouring of a large form");
    let v1, c1 = statistics::covariance_log_pdet(1, diag[..], 100);
    check(bad, c1 == 0 && v1 > 0.0, "the log pdet of a large form");
    let r1, cr1 = statistics::covariance_rank(1, diag[..], 100);
    check(bad, cr1 == 0 && r1 == 100, "the rank of a large form");
    return bad;
}
