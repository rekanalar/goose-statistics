// Time rows of the five laws of the second search group: a long run of
// calls ends in a bounded time and the answers stay finite. The
// quantile of each law costs a search.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_search_b_time() -> i64 {
    var bad: i64 = 0;
    var total: f64 = 0.0;
    let n = 50;
    for k in 0..n {
        let p = (k as f64 + 0.5) / n;
        let x, cx = statistics::loguniform_ppf(p, 1.0, 4.0);
        let y, cy = statistics::irwinhall_ppf(p, 3.0);
        if cx != 0 || cy != 0 { bad += 1; }
        total += x + y;
    }
    let z, cz = statistics::gausshyper_ppf(0.5, 2.0, 3.0, 1.0, 1.0);
    chk_code(bad, "gausshyper ppf code", cz, 0);
    let w, cw = statistics::geninvgauss_ppf(0.5, 1.5, 2.0);
    chk_code(bad, "geninvgauss ppf code", cw, 0);
    let v, cv = statistics::genhyperbolic_ppf(0.5, 1.5, 1.0, 0.5);
    chk_code(bad, "genhyperbolic ppf code", cv, 0);
    total += z + w + v;
    check(bad, !is_nan(total) && !is_inf(total), "the sum of the calls is finite");
    return bad;
}
