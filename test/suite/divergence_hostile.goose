// Hostile input tests of the power divergence. Each answer is a result
// with a p-value in 0 to 1 and code 0, or NaN with a known fault code.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_divergence_hostile() -> i64 {
    var bad: i64 = 0;
    var st: u64 = 515253;
    let n = tests_count();
    for k in 0..n {
        let hostile = rnd(st) % 3 != 0;
        let xs = hostile_sample(st, hostile);
        var ys = hostile_sample(st, hostile);
        if rnd(st) % 2 == 0 {
            ys.resize(0);
            for v in xs { ys.push(if hostile { hostile_value(st) } else { v }); }
        }
        let lam = (rnd(st) % 9 - 4) / 2.0;
        let lam2 = hostile_value(st);
        let ddof = rnd(st) % 5 - 1;
        let r1, c1 = statistics::power_divergence(xs[..], ys[..], lam, ddof);
        hostile_report(bad, sound(r1, c1), "power_divergence", k, c1);
        let r2, c2 = statistics::power_divergence(xs[..], ys[..], lam2, ddof);
        hostile_report(bad, sound(r2, c2), "power_divergence odd power", k, c2);
        let r3, c3 = statistics::power_divergence_uniform(xs[..], lam, ddof);
        hostile_report(bad, sound(r3, c3), "power_divergence_uniform", k, c3);
    }
    return bad;
}
