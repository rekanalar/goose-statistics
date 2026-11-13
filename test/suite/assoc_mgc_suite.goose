// The laws of the multiscale graph correlation: the map lies in the range
// of a correlation, a shift or a scale of a sample keeps the statistic,
// the order of the points does not matter, the two samples swap their
// scales, and the p-value follows the permutations.
namespace suite;

import std;
import math;
import statistics from "../..";

fn mgc_perm_table(out: f64[>..]&, d: const f64[:], order: const i64[:], n: i64) {
    out.resize(0);
    for i in 0..n {
        for j in 0..n { out.push(d[order[i] * n + order[j]]); }
    }
}

fn run_assoc_mgc_tests() -> i64 {
    var bad: i64 = 0;
    var map: f64[>..] = [];
    var gen = rng { s: 11 };
    let x: f64[12] = [2.1, 3.4, 1.9, 5.6, 4.4, 3.3, 6.1, 2.8, 4.9, 3.7, 5.1, 1.2];
    let y: f64[12] = [1.2, 3.9, 2.2, 5.1, 4.8, 2.9, 6.6, 3.1, 4.1, 3.5, 5.9, 1.0];
    let a, ca = statistics::multiscale_graphcorr(map, x[..], y[..], 1, 1, 0, gen);
    chk_code(bad, "base code", ca, 0);
    check(bad, a.statistic <= 1.0 && a.statistic > 0.0, "the statistic of a related pair is positive");
    check(bad, a.rows * a.cols == map.len, "the map holds rows times columns values");
    var inside = true;
    for v in map {
        if is_nan(v) || v > 1.0 || v < -1.0 { inside = false; }
    }
    check(bad, inside, "every local correlation lies in the range of a correlation");
    var edge = true;
    for l in 0..a.cols { if map[l] != 0.0 { edge = false; } }
    for k in 0..a.rows { if map[k * a.cols] != 0.0 { edge = false; } }
    check(bad, edge, "scale 1 gives 0, because a point is its own nearest neighbor");
    check(bad, a.pvalue == 1.0, "no repetition gives the p-value 1");
    // a shift and a scale of a sample keep the statistic
    var sx: f64[12] = [0.0; 12];
    var sy: f64[12] = [0.0; 12];
    for k in 0..12 { sx[k] = x[k] * 3.0 + 100.0; sy[k] = y[k] * 0.5 - 7.0; }
    let s, cs = statistics::multiscale_graphcorr(map, sx[..], sy[..], 1, 1, 0, gen);
    chk_val(bad, "shift and scale", s.statistic, cs, a.statistic, 1e-9);
    check(bad, s.kx == a.kx && s.ky == a.ky, "shift and scale keep the scale");
    // a very small and a very large scale keep the statistic: the product
    // of two variances must not underflow or overflow
    for k in 0..12 { sx[k] = x[k] * 1e-100; sy[k] = y[k] * 1e-100; }
    let t1, ct1 = statistics::multiscale_graphcorr(map, sx[..], sy[..], 1, 1, 0, gen);
    chk_val(bad, "tiny scale", t1.statistic, ct1, a.statistic, 1e-9);
    for k in 0..12 { sx[k] = x[k] * 1e80; sy[k] = y[k] * 1e80; }
    let t2, ct2 = statistics::multiscale_graphcorr(map, sx[..], sy[..], 1, 1, 0, gen);
    chk_val(bad, "huge scale", t2.statistic, ct2, a.statistic, 1e-9);
    // a joint permutation of the points keeps the statistic
    let order: i64[12] = [5, 3, 11, 0, 8, 1, 7, 2, 10, 4, 9, 6];
    var px: f64[12] = [0.0; 12];
    var py: f64[12] = [0.0; 12];
    for k in 0..12 { px[k] = x[order[k]]; py[k] = y[order[k]]; }
    let p, cp = statistics::multiscale_graphcorr(map, px[..], py[..], 1, 1, 0, gen);
    chk_val(bad, "joint permutation", p.statistic, cp, a.statistic, 1e-9);
    // the two samples swap
    let w, cw = statistics::multiscale_graphcorr(map, y[..], x[..], 1, 1, 0, gen);
    chk_val(bad, "swap", w.statistic, cw, a.statistic, 1e-9);
    check(bad, w.kx == a.ky && w.ky == a.kx, "swap swaps the scale");
    // the p-value takes the values of the permutation count
    var g1 = rng { s: 5 };
    let q, cq = statistics::multiscale_graphcorr(map, x[..], y[..], 1, 1, 19, g1);
    let steps = q.pvalue * 20.0;
    check(bad, cq == 0 && abs(steps - round(steps)) < 1e-9 && q.pvalue >= 0.05 && q.pvalue <= 1.0, "the p-value counts permutations");
    var g2 = rng { s: 5 };
    let q2, cq2 = statistics::multiscale_graphcorr(map, x[..], y[..], 1, 1, 19, g2);
    chk_val(bad, "one seed gives one p-value", q2.pvalue, cq2, q.pvalue, 0.0);
    // the shortcut of the permutations equals the full statistic of a permuted table
    let n = 12;
    var dx: f64[>..] = [];
    var dy: f64[>..] = [];
    for i in 0..n {
        for j in 0..n { dx.push(abs(x[i] - x[j])); dy.push(abs(y[i] - y[j])); }
    }
    var ident: i64[>..] = [];
    for k in 0..n { ident.push(k); }
    var g3 = rng { s: 9 };
    var count: i64 = 0;
    for rep in 0..30 {
        statistics::rs_shuffle_index(ident[..], 0, n, g3);
        var tab: f64[>..] = [];
        mgc_perm_table(tab, dy[..], ident[..], n);
        let r, cr = statistics::multiscale_graphcorr_dist(map, dx[..], tab[..], 0, gen);
        if r.statistic >= a.statistic { count += 1; }
    }
    var g4 = rng { s: 9 };
    let full, cf = statistics::multiscale_graphcorr_dist(map, dx[..], dy[..], 30, g4);
    chk_val(bad, "the permutation shortcut", full.pvalue, cf, (1.0 + count as f64) / 31.0, 1e-12);
    // a strong curve gives a small p-value, and a shuffled sample gives a large one
    var cx: f64[>..] = [];
    var cy: f64[>..] = [];
    for k in 0..30 {
        let t = (k as f64 - 15.0) / 5.0;
        cx.push(t);
        cy.push(t * t + 0.01 * (k % 3) as f64);
    }
    var g5 = rng { s: 21 };
    let c1, cc1 = statistics::multiscale_graphcorr(map, cx[..], cy[..], 1, 1, 99, g5);
    check(bad, cc1 == 0 && c1.pvalue <= 0.05, "a curve gives a small p-value");
    // a constant sample gives equal statistics, so every permutation counts
    let flat: f64[6] = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0];
    let ramp: f64[6] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    var g6 = rng { s: 2 };
    let fl, cfl = statistics::multiscale_graphcorr(map, flat[..], ramp[..], 1, 1, 9, g6);
    chk_val(bad, "equal statistics count as extreme", fl.pvalue, cfl, 1.0, 1e-12);
    // the work bound refuses many repetitions on a large sample
    var big: f64[>..] = [];
    for k in 0..200 { big.push(k as f64); }
    let wr, cwr = statistics::multiscale_graphcorr(map, big[..], big[..], 1, 1, 300, gen);
    chk_code(bad, "more work than the bound allows", cwr, 9);
    check(bad, is_nan(wr.statistic) && map.len == 0, "a refused call gives an empty result");
    // the two-sample form pools the samples
    let u: f64[10] = [0.1, 0.4, 0.2, 0.9, 0.3, 0.7, 0.5, 0.6, 0.8, 0.0];
    let v: f64[8] = [1.1, 1.4, 1.2, 1.9, 1.3, 1.7, 1.5, 1.6];
    let t, ct = statistics::multiscale_graphcorr_twosample(map, u[..], v[..], 1, 0, gen);
    chk_code(bad, "two-sample code", ct, 0);
    check(bad, t.statistic > 0.5, "two separate samples give a large statistic");
    return bad;
}
