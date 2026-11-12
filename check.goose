// The input checks that every call runs first, in a fixed order, so the
// same input always gives the same code.
namespace statistics;

import std;
import math;

// The checks on a sample. The order is the length bound, the empty
// sample, NaN, infinity and then the length the call needs.
fn check_sample(xs: const f64[:], need: i64) -> i64 {
    if xs.len > 4194304 { return 9; }
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
