// Time rows of the families of a power or an exponent: a long run of calls
// ends in a bounded time and the answers stay finite. The moments use the
// rule of the file distributions_numeric.goose, so they cost more than the
// other calls.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_exponent_time() -> i64 {
    var bad: i64 = 0;
    var total: f64 = 0.0;
    let n = 50;
    for k in 0..n {
        let p = (k as f64 + 0.5) / n;
        let x, cx = statistics::exponweib_ppf(p, 1.5, 2.0);
        let d, cd = statistics::powernorm_cdf(x, 1.6);
        if cx != 0 || cd != 0 { bad += 1; }
        total += x + d;
    }
    // The rule of the moments runs on each call of the moments.
    let m, cm = statistics::powerlognorm_mean(1.4, 2.0);
    chk_code(bad, "powerlognorm mean code", cm, 0);
    let k2, ck2 = statistics::powerlognorm_kurt(1.4, 2.0);
    chk_code(bad, "powerlognorm kurt code", ck2, 0);
    total += m + k2;
    check(bad, !is_nan(total) && !is_inf(total), "the sum of the calls is finite");
    return bad;
}
