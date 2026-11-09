// Order statistics: the extremes, the medians, the modes and the
// quantiles.
//
// A call that needs the order of the values sorts a copy of the sample.
// The copy belongs to the call, so the data of the caller never changes.
// The sort is the sort of the standard library. It is a pattern defeating
// quicksort, and a heapsort takes over when the partitions stay bad, so
// the worst case is n log n. A checked sample holds no NaN, so the order
// by < is a strict weak order.
namespace statistics;

import std;
import math;

// A sorted copy of a checked sample.
fn sort_copy(xs: const f64[:]) -> f64[>..] {
    var out: f64[>..] = [];
    for x in xs { out.push(x); }
    sort(out) { a, b => a < b };
    return out;
}

// The mean of two values a <= b that cannot overflow.
fn mid(a: f64, b: f64) -> f64 {
    if a < 0.0 && b > 0.0 { return (a + b) / 2.0; }
    return a + (b - a) / 2.0;
}

fn minimum(xs: const f64[:]) -> f64, i64 {
    let c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    var r = xs[0];
    for x in xs {
        if x < r { r = x; }
    }
    return r, 0;
}

fn maximum(xs: const f64[:]) -> f64, i64 {
    let c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    var r = xs[0];
    for x in xs {
        if x > r { r = x; }
    }
    return r, 0;
}

// The largest minus the smallest. A range that does not fit gives the
// overflow fault.
fn ptp(xs: const f64[:]) -> f64, i64 {
    let lo, c = minimum(xs);
    if c != 0 { return nan(), c; }
    let hi, c2 = maximum(xs);
    if c2 != 0 { return nan(), c2; }
    let r = hi - lo;
    if is_inf(r) { return nan(), 8; }
    return r, 0;
}

fn median(xs: const f64[:]) -> f64, i64 {
    let c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    let s = sort_copy(xs);
    let n = s.len;
    if n % 2 == 1 { return s[n / 2], 0; }
    return mid(s[n / 2 - 1], s[n / 2]), 0;
}

// The median of a sample that the caller has sorted.
fn median_sorted(s: const f64[:]) -> f64 {
    let n = s.len;
    if n == 0 { return nan(); }
    if n % 2 == 1 { return s[n / 2]; }
    return mid(s[n / 2 - 1], s[n / 2]);
}

fn median_low(xs: const f64[:]) -> f64, i64 {
    let c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    let s = sort_copy(xs);
    return s[(s.len - 1) / 2], 0;
}

fn median_high(xs: const f64[:]) -> f64, i64 {
    let c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    let s = sort_copy(xs);
    return s[s.len / 2], 0;
}

// The grouped median of Python. The width is the width of a class. The
// value at the middle index is the midpoint of the median class.
fn median_grouped(xs: const f64[:], width: f64) -> f64, i64 {
    let c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    if is_nan(width) { return nan(), 3; }
    if is_inf(width) || width <= 0.0 { return nan(), 5; }
    let s = sort_copy(xs);
    let n = s.len;
    let v = s[n / 2];
    var below: i64 = 0;
    var equal: i64 = 0;
    for x in s {
        if x < v { below += 1; }
        if x == v { equal += 1; }
    }
    let r = (v - width / 2.0) + width * (n / 2.0 - below) / equal;
    if is_inf(r) || is_nan(r) { return nan(), 8; }
    return r, 0;
}

struct Vi { v: f64, i: i64 }

// The runs of equal values of a sample. The run that wins has the most
// members. A tie goes to the run whose first member comes first in the
// sample. The first member of a run is the member with the lowest index.
// The output holds the first index of each run with the top count.
fn top_runs(firsts: i64[>..]&, xs: const f64[:]) {
    firsts.resize(0);
    var ps: Vi[>..] = [];
    for k in 0..xs.len { ps.push(Vi { v: xs[k], i: k }); }
    sort(ps) { a, b => a.v < b.v || (a.v == b.v && a.i < b.i) };
    var best: i64 = 0;
    var k: i64 = 0;
    while k < ps.len {
        var j = k;
        while j < ps.len && ps[j].v == ps[k].v { j += 1; }
        let count = j - k;
        if count > best {
            best = count;
            firsts.resize(0);
        }
        if count == best { firsts.push(ps[k].i); }
        k = j;
    }
}

fn mode(xs: const f64[:]) -> f64, i64 {
    let c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    var firsts: i64[>..] = [];
    top_runs(firsts, xs);
    var first = firsts[0];
    for f in firsts {
        if f < first { first = f; }
    }
    return xs[first], 0;
}

// Writes every most common value into out, in the order of the sample.
// Returns the count of values written and the code. A fault gives 0.
fn multimode(out: f64[>..]&, xs: const f64[:]) -> i64, i64 {
    out.resize(0);
    let c = check_sample(xs, 1);
    if c != 0 { return 0, c; }
    var firsts: i64[>..] = [];
    top_runs(firsts, xs);
    sort(firsts) { a, b => a < b };
    for f in firsts { out.push(xs[f]); }
    return out.len, 0;
}

// Linear interpolation that is exact at both ends.
fn lerp_exact(a: f64, b: f64, t: f64) -> f64 {
    let d = b - a;
    if is_inf(d) { return a * (1.0 - t) + b * t; }
    if t >= 0.5 { return b - d * (1.0 - t); }
    return a + d * t;
}

fn clamp_index(k: f64, n: i64) -> i64 {
    if k <= 0.0 { return 0; }
    if k >= n - 1 { return n - 1; }
    return k as i64;
}

// The quantile of a sorted sample at q, by one of the nine methods. The
// arguments are already checked.
fn quantile_sorted(s: const f64[:], q: f64, method: i64) -> f64 {
    let n = s.len;
    let nq = n * q;
    if method == 0 {
        let idx = nq - 1.0;
        let p = floor(idx);
        if idx - p == 0.0 { return s[clamp_index(p, n)]; }
        return s[clamp_index(p + 1.0, n)];
    }
    if method == 2 {
        let idx = nq - 1.0 - 0.5;
        let p = floor(idx);
        if idx - p == 0.0 && (p as i64) % 2 == 1 { return s[clamp_index(p, n)]; }
        return s[clamp_index(p + 1.0, n)];
    }
    var v: f64 = 0.0;
    var g: f64 = 0.0;
    var fixed_g = false;
    if method == 1 {
        v = nq - 1.0;
        if v - floor(v) == 0.0 { g = 0.5; } else { g = 1.0; }
        fixed_g = true;
    } else if method == 6 {
        v = (n - 1) * q;
    } else {
        var a: f64 = 0.0;
        var b: f64 = 0.0;
        if method == 3 { a = 0.0; b = 1.0; }
        if method == 4 { a = 0.5; b = 0.5; }
        if method == 5 { a = 0.0; b = 0.0; }
        if method == 7 { a = 1.0 / 3.0; b = 1.0 / 3.0; }
        if method == 8 { a = 0.375; b = 0.375; }
        v = nq + (a + q * (1.0 - a - b)) - 1.0;
    }
    if v >= n - 1 { return s[n - 1]; }
    if v < 0.0 { return s[0]; }
    let p = floor(v);
    if !fixed_g { g = v - p; }
    let j = p as i64;
    return lerp_exact(s[j], s[j + 1], g);
}

// The checks of the arguments of a quantile call.
fn check_q(q: f64, method: i64) -> i64 {
    if is_nan(q) { return 3; }
    if q < 0.0 || q > 1.0 { return 5; }
    if method < 0 || method > 8 { return 5; }
    return 0;
}

// The quantile of a sample at q, from 0 to 1, by one of the nine methods.
fn quantile(xs: const f64[:], q: f64, method: i64) -> f64, i64 {
    var c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    c = check_q(q, method);
    if c != 0 { return nan(), c; }
    let s = sort_copy(xs);
    return quantile_sorted(s[..], q, method), 0;
}

// The interquartile range: the quantile at 0.75 minus the quantile at 0.25.
fn iqr(xs: const f64[:], method: i64) -> f64, i64 {
    var c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    c = check_q(0.5, method);
    if c != 0 { return nan(), c; }
    let s = sort_copy(xs);
    let r = quantile_sorted(s[..], 0.75, method) - quantile_sorted(s[..], 0.25, method);
    if is_inf(r) { return nan(), 8; }
    return r, 0;
}

// Writes groups - 1 cut points into out. Returns the count and the code.
// The count of groups is from 2 to MAX_INPUT.
fn cut_points(out: f64[>..]&, xs: const f64[:], groups: i64, method: i64) -> i64, i64 {
    out.resize(0);
    var c = check_sample(xs, 1);
    if c != 0 { return 0, c; }
    c = check_q(0.5, method);
    if c != 0 { return 0, c; }
    if groups < 2 || groups > 4194304 { return 0, 5; }
    let s = sort_copy(xs);
    for i in 1..groups { out.push(quantile_sorted(s[..], (i as f64) / groups, method)); }
    return out.len, 0;
}
