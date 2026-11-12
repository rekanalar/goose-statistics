// The exact sum of Shewchuk, the method of math.fsum in Python.
//
// The accumulator keeps a short list of partial sums. The partials do not
// overlap, and their sum is exact. The list grows with the sample, and a
// double spans 2098 bit positions, so it stays short. A list that passes
// MAX_PARTIALS sets the fault flag and stops. That limit is a guard and no
// input reaches it.
namespace statistics;

import std;
import math;

// A sample with a value of 2^1000 or more is scaled by 2^-100 before it
// is summed, so that no partial overflows.
const SCALE_LIMIT: f64 = 1.0715086071862673e301;
const SCALE_DOWN: f64 = 7.888609052210118e-31;
const SCALE_UP: f64 = 1.2676506002282294e30;

// The partials do not overlap, so a double holds at most 2098 of them.
const MAX_PARTIALS: i64 = 2200;

struct Acc {
    n: i64,
    fault: bool,
    p: f64[>..],
}

fn acc_new() -> Acc {
    return Acc { .. };
}

// Add one value. The value must be finite.
fn acc_add(a: Acc&, v: f64) {
    var x = v;
    var i: i64 = 0;
    for j in 0..a.n {
        var y = a.p[j];
        if abs(x) < abs(y) {
            let t = x;
            x = y;
            y = t;
        }
        let hi = x + y;
        let lo = y - (hi - x);
        if lo != 0.0 {
            a.p[i] = lo;
            i += 1;
        }
        x = hi;
    }
    if x != 0.0 {
        if i >= MAX_PARTIALS {
            a.fault = true;
            return;
        }
        if i < a.p.len { a.p[i] = x; } else { a.p.push(x); }
        i += 1;
    }
    a.n = i;
}

// The rounded sum of the partials, with the rounding to even that makes
// the result the same for every order of the values.
fn acc_total(a: Acc&) -> f64 {
    var n = a.n;
    var hi: f64 = 0.0;
    if n > 0 {
        n -= 1;
        hi = a.p[n];
        var lo: f64 = 0.0;
        while n > 0 {
            let x = hi;
            n -= 1;
            let y = a.p[n];
            hi = x + y;
            lo = y - (hi - x);
            if lo != 0.0 { break; }
        }
        if n > 0 && ((lo < 0.0 && a.p[n - 1] < 0.0) || (lo > 0.0 && a.p[n - 1] > 0.0)) {
            let y = lo * 2.0;
            let x = hi + y;
            if y == x - hi { hi = x; }
        }
    }
    return hi;
}

// The scale for a sample: 1, or 2^-100 when a value is near the limit.
fn scale_of(xs: const f64[:]) -> f64 {
    var m: f64 = 0.0;
    for x in xs {
        let a = abs(x);
        if a > m { m = a; }
    }
    if m >= SCALE_LIMIT { return SCALE_DOWN; }
    return 1.0;
}

// The exact sum of a checked sample, divided by d, with one rounding of
// the sum. The division happens in the scaled range, so a mean that fits
// is never lost. The fault is 0, 8 for a result that does not fit, or 7
// when the partial array is full.
fn exact_sum_div(xs: const f64[:], d: f64) -> f64, i64 {
    let s = scale_of(xs);
    var a = acc_new();
    for x in xs { acc_add(a, x * s); }
    if a.fault { return nan(), 7; }
    var r = acc_total(a) / d;
    if s != 1.0 { r = r * SCALE_UP; }
    if is_inf(r) || is_nan(r) { return nan(), 8; }
    return r, 0;
}

fn nan() -> f64 { return 0.0 / 0.0; }

fn inf() -> f64 { return 1.0 / 0.0; }

// A result that is not finite means the true result does not fit.
fn finish(r: f64) -> f64, i64 {
    if is_inf(r) || is_nan(r) { return nan(), 8; }
    return r, 0;
}

// The mean of a checked sample. A sample of equal values has exactly that
// value as its mean, which the division of a rounded sum would not give.
fn exact_mean(xs: const f64[:]) -> f64, i64 {
    var same = xs[0] != 0.0;
    for x in xs {
        if x != xs[0] { same = false; break; }
    }
    if same { return xs[0], 0; }
    return exact_sum_div(xs, xs.len);
}

// ---- public calls -----------------------------------------------------

// The exact sum of a sample, rounded once.
fn sum(xs: const f64[:]) -> f64, i64 {
    let c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    return exact_sum_div(xs, 1.0);
}

// The mean of a sample. A mean that fits is never lost, even when the sum
// is near the limit of the type.
fn mean(xs: const f64[:]) -> f64, i64 {
    let c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    return exact_mean(xs);
}

// The same as mean. The name is for a reader who knows Python.
fn fmean(xs: const f64[:]) -> f64, i64 {
    return mean(xs);
}
