// Hostile input tests of the contingency calls. Each answer is a number
// with code 0, or NaN or an empty list with a known fault code.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_contingency_hostile() -> i64 {
    var bad: i64 = 0;
    var st: u64 = 616263;
    let n = tests_count();
    var out: f64[>..] = [];
    for k in 0..n {
        let hostile = rnd(st) % 3 != 0;
        let xs = hostile_sample(st, hostile);
        let rows = rnd(st) % 5 - 1;
        let cols = rnd(st) % 5 - 1;
        let method = rnd(st) % 5 - 1;
        let lam = (rnd(st) % 5 - 2) / 2.0;
        let axis = rnd(st) % 5 - 1;
        var shape: i64[>..] = [];
        let dims = rnd(st) % 4;
        for d in 0..dims { shape.push(rnd(st) % 5 - 1); }
        let v1, c1 = statistics::association(xs[..], rows, cols, method, rnd(st) % 2 == 0, lam);
        hostile_report(bad, sound_finite(v1, c1), "association", k, c1);
        let c2 = statistics::expected_freq(out, xs[..], shape[..]);
        hostile_report(bad, (c2 == 0 && out.len == xs.len) || (c2 >= 1 && c2 <= 10 && out.len == 0), "expected_freq", k, c2);
        let c3 = statistics::margins(out, xs[..], shape[..], axis);
        hostile_report(bad, c3 == 0 || (c3 >= 1 && c3 <= 10 && out.len == 0), "margins", k, c3);
    }
    // a count above the exact bound gives the domain fault
    let hg: f64[4] = [1e160, 1e160, 1e160, 1e160];
    let hv, hc = statistics::association(hg[..], 2, 2);
    chk_code(bad, "association with counts above the bound", hc, 5);
    return bad;
}
