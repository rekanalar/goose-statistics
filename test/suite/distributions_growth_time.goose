// Time rows of the families of a growth or a logistic kind: a long run of
// calls ends in a bounded time and the answers stay finite. The moments
// integrate, so a few calls of them run too.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_growth_time() -> i64 {
    var bad: i64 = 0;
    var total: f64 = 0.0;
    let n = 50;
    for k in 0..n {
        let p = (k as f64 + 0.5) / n;
        let x, cx = statistics::gompertz_ppf(p, 1.5);
        let d, cd = statistics::gompertz_cdf(x, 1.5);
        if cx != 0 || cd != 0 { bad += 1; }
        total += x + d;
    }
    let m, cm = statistics::gompertz_kurt(1.5);
    chk_code(bad, "gompertz kurt code", cm, 0);
    let h, ch = statistics::genhalflogistic_kurt(2.0);
    chk_code(bad, "genhalflogistic kurt code", ch, 0);
    total += m + h;
    check(bad, !is_nan(total) && !is_inf(total), "the sum of the calls is finite");
    return bad;
}
