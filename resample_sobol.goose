// The Sobol indices of the package.
//
// Every call returns the result first and a fault code last. A fault
// gives NaN as the result.
namespace statistics;

import std;
import math;

fn rs_sobol_center(f_a: const f64[:], f_b: const f64[:]) -> f64, f64, i64 {
    let n = f_a.len;
    var m: f64 = 0.0;
    for i in 0..n { m += f_a[i] + f_b[i]; }
    m = m / (2.0 * n);
    var v: f64 = 0.0;
    for i in 0..n {
        v += (f_a[i] - m) * (f_a[i] - m) + (f_b[i] - m) * (f_b[i] - m);
    }
    v = v / (2.0 * n);
    if is_inf(v) { return nan(), nan(), ERR_OVERFLOW; }
    if v <= 0.0 { return nan(), nan(), ERR_UNDEFINED; }
    return m, v, OK;
}

// The indices from three lists of evaluations: f_a at the points A, f_b at
// the points B and f_ab at the d matrices that take column k from B. The
// list f_ab holds d pages of n values, page by page. The first-order index
// is the Saltelli form and the total index is the Saltelli form of 2010.
// The values are centered by the mean of f_a and f_b. A constant output has
// no variance and gives the undefined fault.
fn sobol_from_evaluations(first: f64[>..]&, total: f64[>..]&, f_a: const f64[:], f_b: const f64[:], f_ab: const f64[:], d: i64) -> i64, i64 {
    first.resize(0);
    total.resize(0);
    if d > MAX_SOBOL_DIMS { return 0, ERR_TOO_LONG; }
    var c = check_pair(f_a, f_b, 2);
    if c != OK { return 0, c; }
    c = check_sample(f_ab, 1);
    if c != OK { return 0, c; }
    if d < 1 { return 0, ERR_DOMAIN; }
    let n = f_a.len;
    if f_ab.len != d * n { return 0, ERR_LENGTH; }
    let m, v, cv = rs_sobol_center(f_a, f_b);
    if cv != OK { return 0, cv; }
    for k in 0..d {
        var s1: f64 = 0.0;
        var s2: f64 = 0.0;
        for i in 0..n {
            let ab = f_ab[k * n + i];
            s1 += (f_b[i] - m) * (ab - f_a[i]);
            s2 += (f_a[i] - ab) * (f_a[i] - ab);
        }
        let a = s1 / n / v;
        let t = 0.5 * s2 / n / v;
        if is_bad(a) || is_bad(t) {
            first.resize(0);
            total.resize(0);
            return 0, ERR_OVERFLOW;
        }
        first.push(a);
        total.push(t);
    }
    return d, OK;
}

// The first-order and total-order Sobol indices of a model. The block is
// the model. It takes a point of d values, each in the open interval from
// 0 to 1, and returns one value. The caller maps the point to its own
// distributions. The call draws n points A and n points B from the
// generator, and runs the model n * (d + 2) times. The draws are
// pseudo-random, not a scrambled Sobol sequence, so the error shrinks more
// slowly than the error of SciPy for the same n.
fn sobol_indices<F>(first: f64[>..]&, total: f64[>..]&, d: i64, n: i64, r: rng&) -> i64, i64 {
    first.resize(0);
    total.resize(0);
    if d > MAX_SOBOL_DIMS { return 0, ERR_TOO_LONG; }
    if n > MAX_SOBOL_POINTS { return 0, ERR_TOO_LONG; }
    if n < 2 { return 0, ERR_TOO_SHORT; }
    if d < 1 { return 0, ERR_DOMAIN; }
    // The state before the first draw. A second pass replays the draws, so
    // the points of the matrices A and B never sit in memory.
    var replay = rng { s: r.s };
    var pa: f64[32];
    var pb: f64[32];
    var pt: f64[32];
    var f_a: f64[>..] = [];
    var f_b: f64[>..] = [];
    for j in 0..n {
        for i in 0..d { pa[i] = rs_unit(r); }
        for i in 0..d { pb[i] = rs_unit(r); }
        let va, ca = rs_stat(pa[0..d]) { s => F(s) };
        if ca != OK { return 0, ca; }
        let vb, cb = rs_stat(pb[0..d]) { s => F(s) };
        if cb != OK { return 0, cb; }
        f_a.push(va);
        f_b.push(vb);
    }
    let m, v, cv = rs_sobol_center(f_a, f_b);
    if cv != OK { return 0, cv; }
    first.resize(d, 0.0);
    total.resize(d, 0.0);
    for j in 0..n {
        for i in 0..d { pa[i] = rs_unit(replay); }
        for i in 0..d { pb[i] = rs_unit(replay); }
        for k in 0..d {
            for i in 0..d { pt[i] = pa[i]; }
            pt[k] = pb[k];
            let ab, cab = rs_stat(pt[0..d]) { s => F(s) };
            if cab != OK { first.resize(0); total.resize(0); return 0, cab; }
            first[k] += (f_b[j] - m) * (ab - f_a[j]);
            total[k] += (f_a[j] - ab) * (f_a[j] - ab);
        }
    }
    for k in 0..d {
        first[k] = first[k] / n / v;
        total[k] = 0.5 * total[k] / n / v;
        if is_bad(first[k]) || is_bad(total[k]) {
            first.resize(0);
            total.resize(0);
            return 0, ERR_OVERFLOW;
        }
    }
    return d, OK;
}
