// The shape of a sample: moments, skewness, kurtosis and z-scores.
//
// Every sum runs over deviations that are scaled by a power of two near
// the largest one. The scaled sums are exact. The skewness and the
// kurtosis do not depend on the scale, so they never overflow.
namespace statistics;

import std;
import math;

// The moments about the mean of a checked sample, scaled by a power of
// two. The rounded mean m0 is not the true mean. The sums run about m0,
// where the deviations are exact for close values. The shift to the true
// mean is exact algebra: the shift is the mean of the deviations, which
// is tiny, so the correction terms stay small and the result keeps its
// digits even for a large common offset.
struct Raw { m: f64[9], t: f64, same: bool }

fn raw_about_mean(xs: const f64[:], kmax: i64) -> Raw, i64 {
    let n = xs.len;
    let m0, c = exact_mean(xs);
    if c != 0 { return Raw { .. }, c; }
    var top: f64 = 0.0;
    for x in xs {
        let d = abs(x - m0);
        if is_inf(d) { return Raw { .. }, 8; }
        if d > top { top = d; }
    }
    if top == 0.0 { return Raw { same: true, .. }, 0; }
    let t = pow2_floor(top);
    var r = Raw { t: t, .. };
    r.m[0] = 1.0;
    for j in 1..kmax + 1 {
        var a = acc_new();
        for x in xs {
            let d = (x - m0) / t;
            var p = d;
            for i in 1..j { p = p * d; }
            acc_add(a, p);
        }
        if a.fault { return Raw { .. }, 7; }
        r.m[j] = acc_total(a) / n;
    }
    return r, 0;
}

// The central moment of order k, scaled, from the moments about m0.
fn central_scaled(r: const Raw&, k: i64) -> f64 {
    let delta = r.m[1];
    var total: f64 = 0.0;
    var binom: f64 = 1.0;
    for j in 0..k + 1 {
        var p: f64 = 1.0;
        for i in 0..k - j { p = p * (-delta); }
        total += binom * p * r.m[j];
        binom = binom * (k - j) / (j + 1);
    }
    return total;
}

struct Shape { m2: f64, m3: f64, m4: f64 }

// The scaled central moments 2, 3 and 4. A constant sample gives the
// undefined fault.
fn shape_of(xs: const f64[:]) -> Shape, i64 {
    let r, c = raw_about_mean(xs, 4);
    if c != 0 { return Shape { .. }, c; }
    if r.same { return Shape { .. }, 7; }
    let s = Shape { m2: central_scaled(r, 2), m3: central_scaled(r, 3), m4: central_scaled(r, 4) };
    if s.m2 <= 0.0 { return Shape { .. }, 7; }
    return s, 0;
}

// The central moment of order k: the mean of (x - mean)^k.
fn central_moment(xs: const f64[:], order: i64) -> f64, i64 {
    let c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    if order < 1 || order > 8 { return nan(), 5; }
    if order == 1 { return 0.0, 0; }
    let r, c2 = raw_about_mean(xs, order);
    if c2 != 0 { return nan(), c2; }
    if r.same { return 0.0, 0; }
    var v = central_scaled(r, order);
    for k in 0..order {
        v = v * r.t;
        if is_inf(v) { return nan(), 8; }
    }
    return v, 0;
}

// The raw moment of order k: the mean of x^k.
fn moment(xs: const f64[:], order: i64) -> f64, i64 {
    var c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    if order < 1 || order > 8 { return nan(), 5; }
    var top: f64 = 0.0;
    for x in xs {
        let d = abs(x);
        if d > top { top = d; }
    }
    if top == 0.0 { return 0.0, 0; }
    let t = pow2_floor(top);
    var a = acc_new();
    for x in xs {
        let d = x / t;
        var p = d;
        for k in 1..order { p = p * d; }
        acc_add(a, p);
    }
    if a.fault { return nan(), 7; }
    var r = acc_total(a) / xs.len;
    for k in 0..order {
        r = r * t;
        if is_inf(r) { return nan(), 8; }
    }
    return r, 0;
}

// The skewness. The biased form is m3 / m2^1.5. The unbiased form needs
// three values.
fn skew(xs: const f64[:], bias: bool) -> f64, i64 {
    let need = if bias { 1 } else { 3 };
    let c = check_sample(xs, need);
    if c != 0 { return nan(), c; }
    let s, c2 = shape_of(xs);
    if c2 != 0 { return nan(), c2; }
    let g1 = s.m3 / (s.m2 * sqrt(s.m2));
    if bias { return finish(g1); }
    let n = xs.len;
    return finish(g1 * sqrt((n * (n - 1)) as f64) / (n - 2));
}

// The kurtosis. The Fisher form subtracts 3. The unbiased form needs four
// values.
fn kurtosis(xs: const f64[:], fisher: bool, bias: bool) -> f64, i64 {
    let need = if bias { 1 } else { 4 };
    let c = check_sample(xs, need);
    if c != 0 { return nan(), c; }
    let s, c2 = shape_of(xs);
    if c2 != 0 { return nan(), c2; }
    let pearson = s.m4 / (s.m2 * s.m2);
    let g2 = pearson - 3.0;
    if bias {
        if fisher { return finish(g2); }
        return finish(pearson);
    }
    let n = xs.len;
    let f = ((n - 1) as f64) / ((n - 2) * (n - 3)) * ((n + 1) * g2 + 6.0);
    if fisher { return finish(f); }
    return finish(f + 3.0);
}

// Writes the z-score of every value into out. The argument ddof is 0 or 1.
// Returns the count and the code. A fault gives 0.
fn zscore(out: f64[>..]&, xs: const f64[:], ddof: i64) -> i64, i64 {
    out.resize(0);
    var c = check_sample(xs, 1);
    if c != 0 { return 0, c; }
    if ddof < 0 || ddof > 1 { return 0, 5; }
    if xs.len < 2 { return 0, 2; }
    let sigma, c2 = stdev_core(xs, ddof);
    if c2 != 0 { return 0, c2; }
    if sigma == 0.0 { return 0, 7; }
    let m, c3 = exact_mean(xs);
    if c3 != 0 { return 0, c3; }
    for x in xs { out.push((x - m) / sigma); }
    return out.len, 0;
}
