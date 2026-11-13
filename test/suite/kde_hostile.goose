// Hostile input tests of the kernel density estimate. A good answer
// holds a density that is finite and above zero. A fault gives an empty
// list.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_kde_hostile() -> i64 {
    var bad: i64 = 0;
    var st: u64 = 9173;
    let n = tests_count();
    for k in 0..n {
        let dims = 1 + rnd(st) % 3;
        let points = 2 + rnd(st) % 8;
        var data: f64[>..] = [];
        for i in 0..points * dims {
            if rnd(st) % 5 == 0 { data.push((rnd(st) % 21) as f64 - 10.0); }
            else { data.push((rnd(st) % 21) as f64); }
        }
        let bw = rnd(st) % 3;
        let kd, ck = statistics::kde_fit(data[..], dims, [], bw, 0.5);
        if ck == 0 {
            var out: f64[>..] = [];
            let cp = statistics::kde_pdf(kd, out, data[..]);
            var ok = cp == 0 && out.len == points;
            for i in 0..out.len {
                if is_nan(out[i]) || out[i] < 0.0 { ok = false; }
            }
            hostile_report(bad, ok, "kde", k, cp);
        } else {
            hostile_report(bad, ck >= 1 && ck <= 10, "kde", k, ck);
        }
    }
    // a sample whose covariance overflows gives the overflow fault
    let huge_d: f64[4] = [1e200, -1e200, 5e199, -3e199];
    let hk, hc = statistics::kde_fit(huge_d[..], 1, [], statistics::KDE_SCOTT, 0.0);
    chk_code(bad, "kde with an overflowing covariance", hc, 8);
    return bad;
}
