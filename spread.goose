// The spread of a sample: the variance family, the geometric and the
// harmonic mean.
//
// The variance uses the corrected two pass method. The first pass gets the
// mean with the exact sum. The second pass sums the scaled deviations and
// their squares with the exact sum. The deviations are scaled by a power
// of two near the largest one, so no square overflows or underflows. The
// correction term repairs the rounding of the mean.
namespace statistics;

import std;
import math;

// The power of two at or below m. The value m must be above zero.
fn pow2_floor(m: f64) -> f64 {
    return pow(2.0, floor(log2(m)));
}

// The scaled variance of a checked sample with the divisor n - ddof, and
// the scale t. The variance is the scaled value times t * t. The root of
// the scaled value times t is the deviation, which can fit when the
// variance itself underflows.
fn variance_scaled(xs: const f64[:], ddof: i64) -> f64, f64, i64 {
    let n = xs.len;
    let m, c = exact_mean(xs);
    if c != 0 { return nan(), 1.0, c; }
    var top: f64 = 0.0;
    for x in xs {
        let d = abs(x - m);
        if is_inf(d) { return nan(), 1.0, 8; }
        if d > top { top = d; }
    }
    if top == 0.0 { return 0.0, 1.0, 0; }
    let t = pow2_floor(top);
    var s1 = acc_new();
    var s2 = acc_new();
    for x in xs {
        let d = (x - m) / t;
        acc_add(s1, d);
        acc_add(s2, d * d);
    }
    if s1.fault || s2.fault { return nan(), 1.0, 7; }
    let a = acc_total(s1);
    var s = acc_total(s2) - a * a / n;
    if s < 0.0 { s = 0.0; }
    return s / (n - ddof), t, 0;
}

// The variance of a checked sample with the divisor n - ddof.
fn variance_core(xs: const f64[:], ddof: i64) -> f64, i64 {
    let v, t, c = variance_scaled(xs, ddof);
    if c != 0 { return nan(), c; }
    let r = v * t * t;
    if is_inf(r) || is_nan(r) { return nan(), 8; }
    return r, 0;
}

// The deviation of a checked sample with the divisor n - ddof.
fn stdev_core(xs: const f64[:], ddof: i64) -> f64, i64 {
    let v, t, c = variance_scaled(xs, ddof);
    if c != 0 { return nan(), c; }
    let r = sqrt(v) * t;
    if is_inf(r) || is_nan(r) { return nan(), 8; }
    return r, 0;
}

// The population variance. One value gives 0.
fn pvariance(xs: const f64[:]) -> f64, i64 {
    let c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    return variance_core(xs, 0);
}

// The sample variance. It needs two values.
fn variance(xs: const f64[:]) -> f64, i64 {
    let c = check_sample(xs, 2);
    if c != 0 { return nan(), c; }
    return variance_core(xs, 1);
}

fn pstdev(xs: const f64[:]) -> f64, i64 {
    let c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    return stdev_core(xs, 0);
}

fn stdev(xs: const f64[:]) -> f64, i64 {
    let c = check_sample(xs, 2);
    if c != 0 { return nan(), c; }
    return stdev_core(xs, 1);
}

// The standard error of the mean: the sample deviation over the root of n.
fn sem(xs: const f64[:]) -> f64, i64 {
    let s, c = stdev(xs);
    if c != 0 { return nan(), c; }
    return s / sqrt(xs.len as f64), 0;
}

// The coefficient of variation: the sample deviation over the mean. A mean
// of 0 gives the undefined fault.
fn variation(xs: const f64[:]) -> f64, i64 {
    let s, c = stdev(xs);
    if c != 0 { return nan(), c; }
    let m, c2 = mean(xs);
    if c2 != 0 { return nan(), c2; }
    if m == 0.0 { return nan(), 7; }
    return finish(s / m);
}

// The geometric mean. A zero gives 0. A negative value gives the domain
// fault.
fn geometric_mean(xs: const f64[:]) -> f64, i64 {
    let c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    var zero = false;
    for x in xs {
        if x < 0.0 { return nan(), 5; }
        if x == 0.0 { zero = true; }
    }
    if zero { return 0.0, 0; }
    // Each value is m * 2^e with m near 1. The logs of the m values stay
    // small, so their rounding error stays small. The exponents add as
    // whole numbers. A plain mean of logs would multiply the error of a log
    // by the size of the log, which is 700 for a value near 1e-300.
    var a = acc_new();
    var total: i64 = 0;
    for x in xs {
        let e = floor(log2(x));
        acc_add(a, log(x / pow(2.0, e)));
        total += e as i64;
    }
    if a.fault { return nan(), 7; }
    let n = xs.len;
    let q = div_euclid(total, n);
    let rest = total - q * n;
    let r = exp(acc_total(a) / n) * pow(2.0, rest as f64 / n) * pow(2.0, q as f64);
    if is_inf(r) || is_nan(r) { return nan(), 8; }
    return r, 0;
}

// The harmonic mean. A zero gives 0. A negative value gives the domain
// fault. A value so small that its reciprocal does not fit gives the
// overflow fault, unless the whole sample is small enough to scale.
fn harmonic_mean(xs: const f64[:]) -> f64, i64 {
    let c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    var lo: f64 = 1.7976931348623157e308;
    var hi: f64 = 0.0;
    for x in xs {
        if x < 0.0 { return nan(), 5; }
        if x == 0.0 { return 0.0, 0; }
        if x < lo { lo = x; }
        if x > hi { hi = x; }
    }
    var s: f64 = 1.0;
    if lo < 1e-300 && hi < 1e200 { s = 1.6069380442589903e60; }
    var a = acc_new();
    for x in xs {
        let r = 1.0 / (x * s);
        if is_inf(r) { return nan(), 8; }
        acc_add(a, r);
    }
    if a.fault { return nan(), 7; }
    let r = (xs.len / acc_total(a)) / s;
    if is_inf(r) || is_nan(r) { return nan(), 8; }
    return r, 0;
}
