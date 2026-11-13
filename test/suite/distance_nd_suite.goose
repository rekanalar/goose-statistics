// The laws of the distance of two sets of points: it is 0 for equal
// sets, it is symmetric, it keeps the value of the one-dimensional call,
// it moves with a shift of both sets and it grows with the scale.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distance_nd_tests() -> i64 {
    var bad: i64 = 0;
    let a: f64[8] = [0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 3.0, 3.0];
    let b: f64[6] = [1.0, 1.0, 2.0, 0.0, 4.0, 4.0];
    // equal sets
    let z, cz = statistics::wasserstein_distance_nd(a[..], a[..], 2);
    chk_val(bad, "equal sets", z, cz, 0.0, 1e-12);
    // symmetry
    let ab, c1 = statistics::wasserstein_distance_nd(a[..], b[..], 2);
    let ba, c2 = statistics::wasserstein_distance_nd(b[..], a[..], 2);
    chk_val(bad, "symmetry", ab, c1, ba, 1e-12);
    check(bad, c2 == 0, "symmetry code");
    // one dimension equals the call of one dimension
    let u: f64[5] = [0.0, 1.0, 3.0, 5.0, 8.0];
    let v: f64[4] = [0.5, 2.0, 2.5, 7.0];
    let d1, c3 = statistics::wasserstein_distance_nd(u[..], v[..], 1);
    let d0, c4 = statistics::wasserstein_distance(u[..], v[..]);
    chk_val(bad, "one dimension", d1, c3, d0, 1e-12);
    check(bad, c4 == 0, "one dimension code");
    // a shift of both sets keeps the distance
    var sa: f64[8] = [0.0; 8];
    var sb: f64[6] = [0.0; 6];
    for k in 0..8 { sa[k] = a[k] + 10.0; }
    for k in 0..6 { sb[k] = b[k] + 10.0; }
    let sh, c5 = statistics::wasserstein_distance_nd(sa[..], sb[..], 2);
    chk_val(bad, "shift", sh, c5, ab, 1e-12);
    // a scale of both sets scales the distance
    var ta: f64[8] = [0.0; 8];
    var tb: f64[6] = [0.0; 6];
    for k in 0..8 { ta[k] = a[k] * 3.0; }
    for k in 0..6 { tb[k] = b[k] * 3.0; }
    let sc, c6 = statistics::wasserstein_distance_nd(ta[..], tb[..], 2);
    chk_val(bad, "scale", sc, c6, 3.0 * ab, 1e-12);
    // the weights are scaled to a sum of 1
    let w1: f64[4] = [1.0, 2.0, 1.0, 1.0];
    let w2: f64[3] = [2.0, 1.0, 1.0];
    var x1: f64[4] = [0.0; 4];
    var x2: f64[3] = [0.0; 3];
    for k in 0..4 { x1[k] = w1[k] * 7.0; }
    for k in 0..3 { x2[k] = w2[k] * 0.25; }
    let wa, c7 = statistics::wasserstein_distance_nd_weighted(a[..], b[..], 2, w1[..], w2[..]);
    let wb, c8 = statistics::wasserstein_distance_nd_weighted(a[..], b[..], 2, x1[..], x2[..]);
    chk_val(bad, "weights scale", wb, c8, wa, 1e-12);
    check(bad, c7 == 0, "weights code");
    // equal weights equal the call with no weights
    let e1: f64[4] = [1.0, 1.0, 1.0, 1.0];
    let e2: f64[3] = [1.0, 1.0, 1.0];
    let we, c9 = statistics::wasserstein_distance_nd_weighted(a[..], b[..], 2, e1[..], e2[..]);
    chk_val(bad, "equal weights", we, c9, ab, 1e-12);
    // a point of weight 0 does not count
    let p: f64[10] = [0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 3.0, 3.0, 50.0, 50.0];
    let pw: f64[5] = [1.0, 1.0, 1.0, 1.0, 0.0];
    let pe: f64[3] = [1.0, 1.0, 1.0];
    let zq, c10 = statistics::wasserstein_distance_nd_weighted(p[..], b[..], 2, pw[..], pe[..]);
    chk_val(bad, "a zero weight", zq, c10, ab, 1e-12);
    // the triangle inequality on three sets
    let c: f64[4] = [2.0, 2.0, 0.0, 5.0];
    let ac, c11 = statistics::wasserstein_distance_nd(a[..], c[..], 2);
    let bc, c12 = statistics::wasserstein_distance_nd(b[..], c[..], 2);
    check(bad, c11 == 0 && c12 == 0 && ac <= ab + bc + 1e-12, "triangle inequality");
    return bad;
}
