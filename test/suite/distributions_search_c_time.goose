// Time rows of the five laws of the third search group: a long run of
// calls ends in a bounded time and the answers stay finite. The
// quantile of each law costs a search.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_search_c_time() -> i64 {
    var bad: i64 = 0;
    var total: f64 = 0.0;
    let n = 10;
    for k in 0..n {
        let p = (k as f64 + 0.5) / n;
        let x, cx = statistics::norminvgauss_ppf(p, 1.0, 0.5);
        let y, cy = statistics::recipinvgauss_ppf(p, 2.0);
        if cx != 0 || cy != 0 { bad += 1; }
        total += x + y;
    }
    let z, cz = statistics::rel_breitwigner_ppf(0.5, 0.5);
    chk_code(bad, "rel_breitwigner ppf code", cz, 0);
    let w, cw = statistics::studentized_range_ppf(0.5, 3.0, 10.0);
    chk_code(bad, "studentized range ppf code", cw, 0);
    let v, cv = statistics::vonmises_ppf(0.5, 2.0);
    chk_code(bad, "vonmises ppf code", cv, 0);
    total += z + w + v;
    check(bad, !is_nan(total) && !is_inf(total), "the sum of the calls is finite");
    return bad;
}
