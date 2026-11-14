// Time rows of the generalised extreme value law: a long run of calls ends
// in a bounded time and the answers stay finite.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_extreme_time() -> i64 {
    var bad: i64 = 0;
    var total: f64 = 0.0;
    let n = 2000;
    for k in 0..n {
        let x = (k as f64 + 0.5) / n as f64 * 6.0 - 1.0;
        let s = ((k as f64 + 0.5) / n as f64 - 0.5) * 4.0;
        let d, cd = statistics::genextreme_cdf(x, s);
        let e, ce = statistics::genextreme_pdf(x, s);
        let f, cf = statistics::genextreme_sf(x, s);
        if cd != 0 || ce != 0 || cf != 0 { bad += 1; }
        total += d + e + f;
    }
    let m, cm = statistics::genextreme_mean(0.5);
    chk_code(bad, "genextreme mean code", cm, 0);
    let k, ck = statistics::genextreme_kurt(0.5);
    chk_code(bad, "genextreme kurt code", ck, 0);
    total += m + k;
    check(bad, !is_nan(total) && !is_inf(total), "the sum of the calls is finite");
    return bad;
}
