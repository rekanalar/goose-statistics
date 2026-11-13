// Hostile input tests of the distance of two sets of points. Each answer
// is a finite distance with code 0, or a fault code with no distance.
namespace suite;

import std;
import math;
import statistics from "../..";

// A set of points. A point may hold a value that has no defined distance.
fn nd_points(st: u64&, n: i64, dims: i64, hostile: bool) -> f64[>..] {
    var out: f64[>..] = [];
    for k in 0..n * dims {
        let r = rnd(st) % 40;
        if hostile && r == 0 { out.push(0.0 / 0.0); }
        else if hostile && r == 1 { out.push(1.0 / 0.0); }
        else if hostile && r == 2 { out.push(1e200); }
        else { out.push((rnd(st) % 11) as f64 - 5.0); }
    }
    return out;
}

fn run_distance_nd_hostile() -> i64 {
    var bad: i64 = 0;
    var st: u64 = 61207;
    let n = tests_count();
    for k in 0..n {
        let dims = 1 + rnd(st) % 3;
        let hostile = rnd(st) % 3 != 0;
        var nu = rnd(st) % 7;
        var nv = rnd(st) % 7;
        if rnd(st) % 20 == 0 { nu = nu + 1; }
        let u = nd_points(st, nu, dims, hostile);
        let v = nd_points(st, nv, dims, hostile);
        let d, c = statistics::wasserstein_distance_nd(u[..], v[..], dims + rnd(st) % 2);
        var ok = false;
        if c == 0 { ok = !is_nan(d) && !is_inf(d) && d >= 0.0; }
        else { ok = c >= 1 && c <= 10 && is_nan(d); }
        hostile_report(bad, ok, "wasserstein_distance_nd", k, c);
        var uw: f64[>..] = [];
        var vw: f64[>..] = [];
        for j in 0..nu { uw.push((rnd(st) % 4) as f64); }
        for j in 0..nv { vw.push((rnd(st) % 4) as f64); }
        let dw, cw = statistics::wasserstein_distance_nd_weighted(u[..], v[..], dims, uw[..], vw[..]);
        var okw = false;
        if cw == 0 { okw = !is_nan(dw) && !is_inf(dw) && dw >= 0.0; }
        else { okw = cw >= 1 && cw <= 10 && is_nan(dw); }
        hostile_report(bad, okw, "wasserstein_distance_nd_weighted", k, cw);
    }
    return bad;
}
