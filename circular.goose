// The circular statistics: the mean, the deviation and the variance of
// a sample of angles, and the mean direction of a set of vectors.
namespace statistics;

import std;
import math;

// The ratio of a full turn to a radian.
const TWO_PI: f64 = 6.283185307179586;

// The check of the range of the angles. The high end must be above the
// low end, and both must be finite.
fn ci_range(low: f64, high: f64) -> i64 {
    let lc = check_scalar(low);
    if lc != OK { return lc; }
    let hc = check_scalar(high);
    if hc != OK { return hc; }
    if high <= low { return ERR_DOMAIN; }
    return OK;
}

// The sine and the cosine of the scaled angles of a sample. The scale
// maps the range of the angles to a full turn.
fn ci_sincos(xs: const f64[:], low: f64, high: f64) -> f64[>..], f64[>..] {
    let period = high - low;
    let scale = TWO_PI / period;
    var s: f64[>..] = [];
    var c: f64[>..] = [];
    for x in xs {
        let t = x * scale;
        s.push(sin(t));
        c.push(cos(t));
    }
    return s, c;
}

// Signature: circmean(xs, low, high) -> value, code
// The circular mean of a sample of angles. The range of the angles is
// from low to high, and the default is a full turn from 0. The result
// lies in the range. A sample whose resultant is 0 gives
// ERR_UNDEFINED, as in SciPy, which gives the arctangent of 0 over 0.
fn circmean(xs: const f64[:], low: f64, high: f64) -> f64, i64 {
    let c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    let rc = ci_range(low, high);
    if rc != 0 { return nan(), rc; }
    let s, co = ci_sincos(xs, low, high);
    let sin_sum = ds_sum(s);
    let cos_sum = ds_sum(co);
    // An angle that overflows when scaled gives a NaN sum.
    if is_bad(sin_sum) || is_bad(cos_sum) { return nan(), ERR_OVERFLOW; }
    if sin_sum == 0.0 && cos_sum == 0.0 { return nan(), ERR_UNDEFINED; }
    let period = high - low;
    let r = atan2(sin_sum, cos_sum);
    let v = r * (period / TWO_PI) - low;
    var out = v - period * floor(v / period);
    return out + low, OK;
}

// The sum of a list of values.
fn ds_sum(xs: const f64[:]) -> f64 {
    var acc = acc_new();
    for x in xs { acc_add(acc, x); }
    if acc.fault { return nan(); }
    return acc_total(acc);
}

// Signature: circvar(xs, low, high) -> value, code
// The circular variance of a sample of angles: 1 - R, where R is the
// mean resultant length. A sample whose resultant is 0 gives
// ERR_UNDEFINED.
fn circvar(xs: const f64[:], low: f64, high: f64) -> f64, i64 {
    let r, c = ci_resultant(xs, low, high);
    if c != 0 { return nan(), c; }
    return 1.0 - r, OK;
}

// Signature: circstd(xs, low, high, normalize) -> value, code
// The circular deviation of a sample of angles: the root of -2 log R,
// where R is the mean resultant length. The result is scaled by the
// range of the angles unless normalize is true.
fn circstd(xs: const f64[:], low: f64, high: f64, normalize: bool) -> f64, i64 {
    let r, c = ci_resultant(xs, low, high);
    if c != 0 { return nan(), c; }
    if r <= 0.0 { return nan(), ERR_UNDEFINED; }
    var rr = r;
    if rr > 1.0 { rr = 1.0; }
    var v = sqrt(-2.0 * log(rr));
    if !normalize { v = v * ((high - low) / TWO_PI); }
    return v, OK;
}

// The mean resultant length of a sample of angles. The value is at
// most 1.
fn ci_resultant(xs: const f64[:], low: f64, high: f64) -> f64, i64 {
    let c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    let rc = ci_range(low, high);
    if rc != 0 { return nan(), rc; }
    let s, co = ci_sincos(xs, low, high);
    let n = xs.len as f64;
    let sin_mean = ds_sum(s) / n;
    let cos_mean = ds_sum(co) / n;
    if is_bad(sin_mean) || is_bad(cos_mean) { return nan(), ERR_OVERFLOW; }
    let r = sqrt(sin_mean * sin_mean + cos_mean * cos_mean);
    if r > 1.0 { return 1.0, OK; }
    return r, OK;
}

// Signature: directional_stats(dir, xs, dims, normalize) -> resultant, code
// The mean direction of a set of vectors and the mean resultant length.
// The list xs holds the rows of the vectors one after the other, and
// dims is the number of the components of a vector. Each row is scaled
// to unit length when normalize is true. The list dir receives the mean
// direction, a unit vector. A row of length 0 gives ERR_UNDEFINED, as
// SciPy gives NaN. A mean of length 0 gives ERR_UNDEFINED.
fn directional_stats(dir: f64[>..]&, xs: const f64[:], dims: i64, normalize: bool) -> f64, i64 {
    dir.resize(0);
    let c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    if dims < 1 { return nan(), ERR_DOMAIN; }
    if xs.len % dims != 0 { return nan(), ERR_LENGTH; }
    let rows = xs.len / dims;
    var sums: f64[>..] = [];
    for i in 0..dims { sums.push(0.0); }
    var r: i64 = 0;
    while r < rows {
        var norm: f64 = 0.0;
        for i in 0..dims {
            let v = xs[r * dims + i];
            norm += v * v;
        }
        norm = sqrt(norm);
        if norm == 0.0 && normalize { return nan(), ERR_UNDEFINED; }
        for i in 0..dims {
            let v = xs[r * dims + i];
            if normalize { sums[i] += v / norm; } else { sums[i] += v; }
        }
        r += 1;
    }
    var res: f64 = 0.0;
    for i in 0..dims {
        let v = sums[i] / rows;
        sums[i] = v;
        res += v * v;
    }
    res = sqrt(res);
    if res == 0.0 { return nan(), ERR_UNDEFINED; }
    for i in 0..dims { dir.push(sums[i] / res); }
    return res, OK;
}
