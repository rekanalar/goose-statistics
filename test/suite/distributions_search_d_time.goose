// Time rows of the stable law: a long run of calls ends in a bounded
// time and the answers stay finite. The quantile of the law costs a
// search, and the density costs an integral over an angle.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_search_d_time() -> i64 {
    var bad: i64 = 0;
    var total: f64 = 0.0;
    let n = 60;
    for k in 0..n {
        let p = (k as f64 + 0.5) / n as f64;
        let x, cx = statistics::levy_stable_ppf(p, 1.5, 0.0);
        if cx != 0 { bad += 1; }
        total += x;
    }
    let y, cy = statistics::levy_stable_pdf(1.0, 1.05, 0.4);
    chk_code(bad, "levy pdf code", cy, 0);
    let z, cz = statistics::levy_stable_cdf(1.0, 0.8, 0.9);
    chk_code(bad, "levy cdf code", cz, 0);
    total += y + z;
    check(bad, !is_nan(total) && !is_inf(total), "the sum of the calls is finite");
    return bad;
}
