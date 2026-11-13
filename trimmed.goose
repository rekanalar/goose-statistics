// The trimmed forms of a sample: the extreme of a range, the spread
// over a range and the trimmed sample itself.
namespace statistics;

import std;
import math;

// Signature: tmin(xs, low, inclusive) -> value, code
// The lowest value of a sample that lies at or above the limit low. The
// limit is inclusive when inclusive is true. A sample with no value
// above the limit gives ERR_UNDEFINED, as in SciPy, which raises.
fn tmin(xs: const f64[:], low: f64, inclusive: bool) -> f64, i64 {
    let c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    let sc = check_scalar(low);
    if sc != OK { return nan(), sc; }
    var best: f64 = nan();
    var found = false;
    for x in xs {
        if !((x > low) || (x == low && inclusive)) { continue; }
        if !found { best = x; found = true; }
        else if x < best { best = x; }
    }
    if !found { return nan(), ERR_UNDEFINED; }
    return best, OK;
}

// Signature: tmax(xs, high, inclusive) -> value, code
// The highest value of a sample that lies at or below the limit high.
// The limit is inclusive when inclusive is true. A sample with no value
// below the limit gives ERR_UNDEFINED.
fn tmax(xs: const f64[:], high: f64, inclusive: bool) -> f64, i64 {
    let c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    let sc = check_scalar(high);
    if sc != OK { return nan(), sc; }
    var best: f64 = nan();
    var found = false;
    for x in xs {
        if !((x < high) || (x == high && inclusive)) { continue; }
        if !found { best = x; found = true; }
        else if x > best { best = x; }
    }
    if !found { return nan(), ERR_UNDEFINED; }
    return best, OK;
}

// The values of a sample that lie inside the limits. The flags tell
// which side is on. The list is empty when no value lies inside.
fn tr_values(xs: const f64[:], low: f64, high: f64, use_low: bool, use_high: bool, lo_inc: bool, hi_inc: bool) -> f64[>..], i64 {
    var out: f64[>..] = [];
    for x in xs {
        if use_low {
            if x < low { continue; }
            if x == low && !lo_inc { continue; }
        }
        if use_high {
            if x > high { continue; }
            if x == high && !hi_inc { continue; }
        }
        out.push(x);
    }
    return out, OK;
}

// Signature: tvar(xs, low, high, use_low, use_high, lo_inc, hi_inc, ddof) -> value, code
// The variance of the values of a sample that lie inside the limits.
// The divisor is n - ddof, so the default ddof is 1. A sample with no
// value inside gives ERR_UNDEFINED. One value with ddof 1 gives
// ERR_TOO_SHORT.
fn tvar(xs: const f64[:], low: f64, high: f64, use_low: bool, use_high: bool, lo_inc: bool, hi_inc: bool, ddof: i64) -> f64, i64 {
    let c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    if ddof < 0 { return nan(), ERR_DOMAIN; }
    if use_low {
        let lc = check_scalar(low);
        if lc != OK { return nan(), lc; }
    }
    if use_high {
        let hc = check_scalar(high);
        if hc != OK { return nan(), hc; }
    }
    if use_low && use_high {
        let lc = check_limits(low, high);
        if lc != 0 { return nan(), lc; }
    }
    let keep, c2 = tr_values(xs, low, high, use_low, use_high, lo_inc, hi_inc);
    if c2 != 0 { return nan(), c2; }
    if keep.len == 0 { return nan(), ERR_UNDEFINED; }
    if keep.len <= ddof { return nan(), ERR_TOO_SHORT; }
    let v, c3 = variance_core(keep[..], ddof);
    if c3 != 0 { return nan(), c3; }
    return v, OK;
}

// Signature: tstd(xs, low, high, use_low, use_high, lo_inc, hi_inc, ddof) -> value, code
// The deviation of the values of a sample that lie inside the limits.
fn tstd(xs: const f64[:], low: f64, high: f64, use_low: bool, use_high: bool, lo_inc: bool, hi_inc: bool, ddof: i64) -> f64, i64 {
    let v, c = tvar(xs, low, high, use_low, use_high, lo_inc, hi_inc, ddof);
    if c != 0 { return nan(), c; }
    return sqrt(v), OK;
}

// Signature: tsem(xs, low, high, use_low, use_high, lo_inc, hi_inc, ddof) -> value, code
// The standard error of the mean of the values of a sample that lie
// inside the limits.
fn tsem(xs: const f64[:], low: f64, high: f64, use_low: bool, use_high: bool, lo_inc: bool, hi_inc: bool, ddof: i64) -> f64, i64 {
    let sd, c = tstd(xs, low, high, use_low, use_high, lo_inc, hi_inc, ddof);
    if c != 0 { return nan(), c; }
    let keep, c2 = tr_values(xs, low, high, use_low, use_high, lo_inc, hi_inc);
    if c2 != 0 { return nan(), c2; }
    return sd / sqrt(keep.len as f64), OK;
}

// Signature: trimboth(out, xs, proportion) -> length, code
// The sample without the given proportion of the values at each end.
// The cut count is int(proportion * n). The kept values are in
// ascending order, as the package does for every sorted result. A
// proportion of 0.5 or more gives ERR_DOMAIN, as in SciPy, which
// raises. A proportion below 0 gives ERR_DOMAIN.
fn trimboth(out: f64[>..]&, xs: const f64[:], proportion: f64) -> i64, i64 {
    out.resize(0);
    let c = check_sample(xs, 1);
    if c != 0 { return 0, c; }
    let pc = check_scalar(proportion);
    if pc != OK { return 0, pc; }
    if proportion < 0.0 || proportion >= 0.5 { return 0, ERR_DOMAIN; }
    let n = xs.len;
    let cut = (proportion * n) as i64;
    let upper = n - cut;
    if cut >= upper { return 0, ERR_DOMAIN; }
    let s = sort_copy(xs);
    for i in cut..upper { out.push(s[i]); }
    return out.len, OK;
}

// Signature: trim1(out, xs, proportion, tail) -> length, code
// The sample without the given proportion of the values at one end. The
// tail is TAIL_LEFT for the lower end and TAIL_RIGHT for the upper
// end. The kept values are in ascending order. A proportion of 1 or
// more gives ERR_DOMAIN. A proportion below 0 gives ERR_DOMAIN.
fn trim1(out: f64[>..]&, xs: const f64[:], proportion: f64, tail: i64) -> i64, i64 {
    out.resize(0);
    let c = check_sample(xs, 1);
    if c != 0 { return 0, c; }
    let pc = check_scalar(proportion);
    if pc != OK { return 0, pc; }
    if proportion < 0.0 || proportion >= 1.0 { return 0, ERR_DOMAIN; }
    if tail != TAIL_LEFT && tail != TAIL_RIGHT { return 0, ERR_DOMAIN; }
    let n = xs.len;
    let cut = (proportion * n) as i64;
    let s = sort_copy(xs);
    if tail == TAIL_RIGHT {
        let upper = n - cut;
        if upper <= 0 { return 0, ERR_DOMAIN; }
        for i in 0..upper { out.push(s[i]); }
    } else {
        let lower = cut;
        if lower >= n { return 0, ERR_DOMAIN; }
        for i in lower..n { out.push(s[i]); }
    }
    return out.len, OK;
}
