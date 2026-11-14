// Time rows of the families of the Johnson kind: a long run of calls ends
// in a bounded time and the answers stay finite. The moments of the
// Johnson law of the B kind integrate, so a few calls of them run too.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_johnson_time() -> i64 {
    var bad: i64 = 0;
    var total: f64 = 0.0;
    let n = 2000;
    for k in 0..n {
        let p = (k as f64 + 0.5) / n as f64;
        let x, cx = statistics::johnsonsb_ppf(p, 1.0, 2.0);
        let d, cd = statistics::johnsonsb_cdf(x, 1.0, 2.0);
        let y, cy = statistics::johnsonsu_ppf(p, 1.0, 2.0);
        if cx != 0 || cd != 0 || cy != 0 { bad += 1; }
        total += x + d + y;
    }
    let m, cm = statistics::johnsonsb_mean(1.0, 2.0);
    chk_code(bad, "johnsonsb mean code", cm, 0);
    let h, ch = statistics::johnsonsb_kurt(1.0, 2.0);
    chk_code(bad, "johnsonsb kurt code", ch, 0);
    total += m + h;
    check(bad, !is_nan(total) && !is_inf(total), "the sum of the calls is finite");
    return bad;
}
