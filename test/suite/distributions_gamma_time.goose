// Time rows of the families that come from the gamma law: a long run of
// calls ends in a bounded time and the answers stay finite. The moments of
// the law of the logarithm of the gamma law need the zeta function, so a
// few calls of them run too.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_gamma_time() -> i64 {
    var bad: i64 = 0;
    var total: f64 = 0.0;
    let n = 2000;
    for k in 0..n {
        let p = (k as f64 + 0.5) / n as f64;
        let x, cx = statistics::chi_ppf(p, 3.0);
        let d, cd = statistics::chi_cdf(x, 3.0);
        let y, cy = statistics::loggamma_ppf(p, 2.0);
        if cx != 0 || cd != 0 || cy != 0 { bad += 1; }
        total += x + d + y;
    }
    let m, cm = statistics::loggamma_mean(2.0);
    chk_code(bad, "loggamma mean code", cm, 0);
    let h, ch = statistics::loggamma_kurt(2.0);
    chk_code(bad, "loggamma kurt code", ch, 0);
    total += m + h;
    check(bad, !is_nan(total) && !is_inf(total), "the sum of the calls is finite");
    return bad;
}
