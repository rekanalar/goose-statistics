// Time rows of the six laws of the first search group: a long run of
// calls ends in a bounded time and the answers stay finite. The
// quantile of each law costs a search.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_search_a_time() -> i64 {
    var bad: i64 = 0;
    var total: f64 = 0.0;
    let n = 50;
    for k in 0..n {
        let p = (k as f64 + 0.5) / n;
        let x, cx = statistics::argus_ppf(p, 1.0);
        let y, cy = statistics::foldnorm_ppf(p, 0.5);
        if cx != 0 || cy != 0 { bad += 1; }
        total += x + y;
    }
    let z, cz = statistics::dpareto_lognorm_ppf(0.5, 0.0, 1.0, 1.0, 1.0);
    chk_code(bad, "dpl ppf code", cz, 0);
    let w, cw = statistics::exponnorm_ppf(0.5, 1.0);
    chk_code(bad, "exponnorm ppf code", cw, 0);
    total += z + w;
    check(bad, !is_nan(total) && !is_inf(total), "the sum of the calls is finite");
    return bad;
}
