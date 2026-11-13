// Hostile input tests of the probability plot calls. Each answer is a
// finite result with code 0, or a fault code with an empty output.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_probplot_hostile() -> i64 {
    var bad: i64 = 0;
    var st: u64 = 90417;
    let n = tests_count();
    var osm: f64[>..] = [];
    var osr: f64[>..] = [];
    for k in 0..n {
        let hostile = rnd(st) % 3 != 0;
        let xs = hostile_sample(st, hostile);
        let fam = rnd(st) % 12;
        let shape = (rnd(st) % 9) / 3.0 - 0.5;
        let p, c = statistics::probplot(osm, osr, xs[..], fam, shape);
        var ok = false;
        if c == 0 { ok = !is_nan(p.slope) && !is_inf(p.slope) && !is_nan(p.r) && osm.len == xs.len && osr.len == xs.len; }
        else { ok = c >= 1 && c <= 10 && is_nan(p.slope) && osm.len == 0 && osr.len == 0; }
        hostile_report(bad, ok, "probplot", k, c);
        let a = (rnd(st) % 9) / 3.0 - 1.0;
        let b = (rnd(st) % 9) / 3.0 - 1.0;
        let c2 = statistics::ppcc_plot(osm, xs[..], a, b, rnd(st) % 12, 3 + (rnd(st) % 9) % 8);
        var ok2 = false;
        if c2 == 0 { ok2 = true; for r in osm { if is_nan(r) || is_inf(r) { ok2 = false; } } }
        else { ok2 = c2 >= 1 && c2 <= 10 && osm.len == 0; }
        hostile_report(bad, ok2, "ppcc_plot", k, c2);
        if k % 20 == 0 {
            let s, c3 = statistics::ppcc_max(xs[..], 3, a, b);
            var ok3 = false;
            if c3 == 0 { ok3 = !is_nan(s) && !is_inf(s); }
            else { ok3 = c3 >= 1 && c3 <= 10 && is_nan(s); }
            hostile_report(bad, ok3, "ppcc_max", k, c3);
        }
    }
    // a sample above the bound gives the too long fault in the searches
    var big: f64[>..] = [];
    for i in 0..statistics::MAX_PPCC_POINTS + 1 { big.push(i as f64); }
    let cb = statistics::ppcc_plot(osm, big[..], 0.0, 1.0, 3);
    chk_code(bad, "ppcc_plot too many values", cb, 9);
    let sb, cs = statistics::ppcc_max(big[..]);
    chk_code(bad, "ppcc_max too many values", cs, 9);
    return bad;
}
