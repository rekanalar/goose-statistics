// Time rows of the families of a symmetric law: a long run of calls ends in
// a bounded time and the answers stay finite. The cumulative function of
// the law of Tukey and of lambda searches, so it runs too.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_symmetric_time() -> i64 {
    var bad: i64 = 0;
    var total: f64 = 0.0;
    let n = 2000;
    for k in 0..n {
        let p = (k as f64 + 0.5) / n as f64;
        let x, cx = statistics::gennorm_ppf(p, 1.5);
        let d, cd = statistics::gennorm_cdf(x, 1.5);
        let y, cy = statistics::tukeylambda_ppf(p, 0.5);
        if cx != 0 || cd != 0 || cy != 0 { bad += 1; }
        total += x + d + y;
    }
    let m, cm = statistics::halfgennorm_kurt(1.5);
    chk_code(bad, "halfgennorm kurt code", cm, 0);
    let h, ch = statistics::gennorm_var(1.5);
    chk_code(bad, "gennorm var code", ch, 0);
    total += m + h;
    check(bad, !is_nan(total) && !is_inf(total), "the sum of the calls is finite");
    return bad;
}
