// Hostile input tests of the distance calls. Each answer is a finite
// distance with code 0, or a fault code with no distance.
namespace suite;

import std;
import math;
import statistics from "../..";

// A list of positive weights with one weight for each value.
fn distance_weights(st: u64&, n: i64) -> f64[>..] {
    var out: f64[>..] = [];
    for k in 0..n { out.push((rnd(st) % 5) + 1.0); }
    return out;
}

fn run_distance_hostile() -> i64 {
    var bad: i64 = 0;
    var st: u64 = 883301;
    let n = tests_count();
    var out: f64[>..] = [];
    for k in 0..n {
        let hostile = rnd(st) % 3 != 0;
        let xs = hostile_sample(st, hostile);
        let ys = hostile_sample(st, hostile);
        let w, c = statistics::wasserstein_distance(xs[..], ys[..]);
        var ok = false;
        if c == 0 { ok = !is_nan(w) && !is_inf(w) && w >= 0.0; }
        else { ok = c >= 1 && c <= 10; }
        hostile_report(bad, ok, "wasserstein_distance", k, c);
        let e, ce = statistics::energy_distance(xs[..], ys[..]);
        var oke = false;
        if ce == 0 { oke = !is_nan(e) && !is_inf(e) && e >= 0.0; }
        else { oke = ce >= 1 && ce <= 10; }
        hostile_report(bad, oke, "energy_distance", k, ce);
        // the weighted calls with a weight for each value
        let ws = distance_weights(st, xs.len);
        let wv = distance_weights(st, ys.len);
        let ww, cw = statistics::wasserstein_distance_weighted(xs[..], ys[..], ws[..], wv[..]);
        var okw = false;
        if cw == 0 { okw = !is_nan(ww) && !is_inf(ww) && ww >= 0.0; }
        else { okw = cw >= 1 && cw <= 10; }
        hostile_report(bad, okw, "wasserstein_distance_weighted", k, cw);
        // the transform of a hostile group
        let sizes: i64[1] = [xs.len];
        let ct = statistics::obrientransform(out, xs[..], sizes[..]);
        var okt = false;
        if ct == 0 { okt = out.len == xs.len; }
        else { okt = ct >= 1 && ct <= 10 && out.len == 0; }
        hostile_report(bad, okt, "obrientransform", k, ct);
    }
    // weights whose sum overflows give the overflow fault
    let ov_u: f64[2] = [0.0, 1.0];
    let ov_big: f64[2] = [1e308, 1e308];
    let ov_one: f64[2] = [1.0, 1.0];
    let ov1, ovc1 = statistics::wasserstein_distance_weighted(ov_u[..], ov_u[..], ov_big[..], ov_one[..]);
    chk_code(bad, "wasserstein weights that overflow", ovc1, 8);
    let ov2, ovc2 = statistics::energy_distance_weighted(ov_u[..], ov_u[..], ov_one[..], ov_big[..]);
    chk_code(bad, "energy weights that overflow", ovc2, 8);
    return bad;
}
