// Hostile input tests of the sigma clip. Each answer is a clip with
// finite bounds and code 0, or a fault code with an empty output.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_clip_hostile() -> i64 {
    var bad: i64 = 0;
    var st: u64 = 777001;
    let n = tests_count();
    var out: f64[>..] = [];
    for k in 0..n {
        let hostile = rnd(st) % 3 != 0;
        let xs = hostile_sample(st, hostile);
        let low = hostile_value(st);
        let high = hostile_value(st);
        let r, c = statistics::sigmaclip(out, xs[..], low, high);
        var ok = false;
        if c == 0 { ok = !is_nan(r.low) && !is_nan(r.high) && r.kept == out.len && r.kept > 0; }
        else { ok = c >= 1 && c <= 10 && out.len == 0 && r.kept == 0; }
        hostile_report(bad, ok, "sigmaclip", k, c);
    }
    return bad;
}
