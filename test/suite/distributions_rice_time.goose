// Time rows of the Rice law: a long run of calls ends in a bounded time
// and the answers stay finite. A large shape needs the most terms of the
// sum of the cumulative function.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_rice_time() -> i64 {
    var bad: i64 = 0;
    var total: f64 = 0.0;
    let n = 20;
    for k in 0..n {
        let p = (k as f64 + 0.5) / n;
        let x, cx = statistics::rice_ppf(p, 5.0);
        let d, cd = statistics::rice_pdf(x, 5.0);
        if cx != 0 || cd != 0 { bad += 1; }
        total += x + d;
    }
    // A shape of 3000 is the largest. The sum runs about 40000 terms.
    let c, cc = statistics::rice_cdf(3000.0, 3000.0);
    chk_code(bad, "rice cdf of the largest shape", cc, 0);
    let q, cq = statistics::rice_ppf(0.5, 3000.0);
    chk_code(bad, "rice ppf of the largest shape", cq, 0);
    let m, cm = statistics::rice_mean(100.0);
    chk_code(bad, "rice mean code", cm, 0);
    total += c + q + m;
    check(bad, !is_nan(total) && !is_inf(total), "the sum of the calls is finite");
    return bad;
}
