// Hostile input tests of the Siegel estimate. Each answer is a finite
// slope and intercept with code 0, or a fault code with NaN in both.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_line_siegel_hostile() -> i64 {
    var bad: i64 = 0;
    var st: u64 = 771203;
    let n = tests_count();
    for k in 0..n {
        let hostile = rnd(st) % 3 != 0;
        let xs = hostile_sample(st, hostile);
        let ys = hostile_sample(st, hostile);
        let method = rnd(st) % 2;
        let s, b, c = statistics::siegelslopes(ys[..], xs[..], method);
        var ok = false;
        if c == 0 { ok = !is_nan(s) && !is_inf(s) && !is_nan(b) && !is_inf(b); }
        else { ok = c >= 1 && c <= 10 && is_nan(s) && is_nan(b); }
        hostile_report(bad, ok, "siegelslopes", k, c);
    }
    // a sample above the point bound gives the too long fault
    var big_x: f64[>..] = [];
    var big_y: f64[>..] = [];
    for i in 0..statistics::MAX_SIEGEL_POINTS + 1 {
        big_x.push(i as f64);
        big_y.push(2.0 * i as f64);
    }
    let _bs, _bb, bc = statistics::siegelslopes(big_y[..], big_x[..]);
    chk_code(bad, "siegelslopes too long", bc, 9);
    return bad;
}
