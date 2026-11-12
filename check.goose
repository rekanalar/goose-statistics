// The input checks that every call runs first, in a fixed order, so the
// same input always gives the same code.
namespace statistics;

import std;
import math;

// The checks on a sample. The order is the length bound, the empty
// sample, NaN, infinity and then the length the call needs.
fn check_sample(xs: const f64[:], need: i64) -> i64 {
    if xs.len > MAX_INPUT { return 9; }
    if xs.len == 0 { return 1; }
    for x in xs {
        if is_nan(x) { return 3; }
    }
    for x in xs {
        if is_inf(x) { return 4; }
    }
    if xs.len < need { return 2; }
    return 0;
}

// The checks of two samples that a call reads apart. Both lengths come
// first, so a long sample costs no scan. Each sample needs at least
// `need` values.
fn check_two(xs: const f64[:], ys: const f64[:], need: i64) -> i64 {
    if xs.len > MAX_INPUT || ys.len > MAX_INPUT { return 9; }
    var c = check_sample(xs, need);
    if c != 0 { return c; }
    return check_sample(ys, need);
}

// The values of two samples in one buffer. The first sample comes
// first.
fn concat2(out: f64[>..]&, xs: const f64[:], ys: const f64[:]) {
    out.resize(0);
    for x in xs { out.push(x); }
    for y in ys { out.push(y); }
}

// The check of a set of groups in one buffer. The sizes array gives the
// length of each group. A call needs at least `min_groups` groups and
// each group needs at least `min_size` values. The order is the length
// bound, the empty sample, the count of the groups, a negative size, a
// short group and then the length match.
fn check_groups_n(groups: const f64[:], sizes: const i64[:], min_groups: i64, min_size: i64) -> i64 {
    let c = check_sample(groups, 1);
    if c != 0 { return c; }
    if sizes.len < min_groups { return ERR_TOO_SHORT; }
    var total: i64 = 0;
    for s in sizes {
        if s < 0 { return ERR_DOMAIN; }
        if s < min_size { return ERR_TOO_SHORT; }
        total += s;
        if total > groups.len { return ERR_LENGTH; }
    }
    if total != groups.len { return ERR_LENGTH; }
    return OK;
}

// True when a value is a NaN or an infinity. The calls use it to reject
// a value that has no defined result.
fn is_bad(x: f64) -> bool {
    return is_nan(x) || is_inf(x);
}

// The check on a list of p-values. A value outside 0 to 1 gives the
// domain fault.
fn check_probabilities(ps: const f64[:]) -> i64 {
    for p in ps {
        if p < 0.0 || p > 1.0 { return 5; }
    }
    return 0;
}

// The checks on a pair of columns. The order is the length bound, the
// empty sample, the length match, NaN, infinity and then the length the
// call needs.
fn check_pair(xs: const f64[:], ys: const f64[:], need: i64) -> i64 {
    if xs.len > 4194304 || ys.len > 4194304 { return 9; }
    if xs.len == 0 || ys.len == 0 { return 1; }
    if xs.len != ys.len { return 6; }
    for x in xs {
        if is_nan(x) { return 3; }
    }
    for y in ys {
        if is_nan(y) { return 3; }
    }
    for x in xs {
        if is_inf(x) { return 4; }
    }
    for y in ys {
        if is_inf(y) { return 4; }
    }
    if xs.len < need { return 2; }
    return 0;
}

// The checks on a pair of columns and their weights. The weight checks
// come after the checks of the columns.
fn check_weighted_pair(xs: const f64[:], ys: const f64[:], ws: const f64[:], need: i64) -> i64 {
    let c = check_pair(xs, ys, 0);
    if c != 0 { return c; }
    return check_weights(ws, xs.len, need);
}

// The checks on a sample and its weights. The weight checks come after
// the sample checks and before the argument checks of the call.
fn check_weighted(xs: const f64[:], ws: const f64[:], need: i64) -> i64 {
    let sc = check_sample(xs, 0);
    if sc != 0 { return sc; }
    return check_weights(ws, xs.len, need);
}

// The checks of the weights of a weighted call. The order is the length
// match, NaN, infinity, a negative weight, an all zero column and the
// length the call needs.
fn check_weights(ws: const f64[:], n: i64, need: i64) -> i64 {
    if ws.len != n { return 6; }
    for w in ws {
        if is_nan(w) { return 3; }
    }
    for w in ws {
        if is_inf(w) { return 4; }
    }
    for w in ws {
        if w < 0.0 { return 5; }
    }
    var top: f64 = 0.0;
    for w in ws {
        if w > top { top = w; }
    }
    if top <= 0.0 { return 7; }
    if n < need { return 2; }
    return 0;
}

// The check of a confidence level: it lies strictly between 0 and 1. A
// NaN gives the NaN code and an infinity the infinity code, so every call
// that takes a level gives the same code for the same bad value.
fn check_level(level: f64) -> i64 {
    let c = check_scalar(level);
    if c != 0 { return c; }
    if level <= 0.0 || level >= 1.0 { return 5; }
    return 0;
}
