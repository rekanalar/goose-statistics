// The checks and the generators that the hostile suites share.
namespace suite;

import std;
import math;
import statistics from "../..";

// A column of weights of the given length. Some of them are hostile, so
// the weight checks of the package give a fault.
fn weights_of(st: u64&, n: i64) -> f64[>..] {
    var out: f64[>..] = [];
    var i: i64 = 0;
    while i < n {
        if rnd(st) % 4 == 0 { out.push(hostile_value(st)); }
        else { out.push((rnd(st) % 4) as f64); }
        i += 1;
    }
    return out;
}

// A list of edges that may or may not ascend.
fn edges_of(st: u64&) -> f64[>..] {
    var out: f64[>..] = [];
    var i: i64 = 0;
    while i < 4 {
        out.push((rnd(st) % 7) as f64 - 3.0);
        i += 1;
    }
    return out;
}

// True when a frequencies answer holds numbers or is a clean fault.
fn sound_freqs(f: statistics::Frequencies, c: i64) -> bool {
    if c == 0 { return !is_nan(f.lowerlimit) && !is_nan(f.binsize); }
    if c < 1 || c > 10 { return false; }
    return f.lowerlimit == 0.0 && f.binsize == 0.0;
}

// True when a summary answer holds numbers or is a clean fault.
fn sound_summary(s: statistics::Summary, c: i64) -> bool {
    if c == 0 { return !is_nan(s.mean) && !is_nan(s.variance); }
    if c < 1 || c > 10 { return false; }
    return s.n == 0 && s.mean == 0.0 && s.variance == 0.0;
}

// True when a value answer is a number or a clean fault.
fn sound_value(v: f64, c: i64) -> bool {
    if c == 0 { return !is_nan(v); }
    if c < 1 || c > 10 { return false; }
    return is_nan(v);
}

// True when a count answer matches the length of the output list.
fn sound_count(n: i64, c: i64, len: i64) -> bool {
    if c == 0 { return n == len; }
    if c < 1 || c > 10 { return false; }
    return n == 0 && len == 0;
}

// True when a table answer is a set of sizes or a clean fault.
fn sound_table(t: statistics::BinnedTable, c: i64) -> bool {
    if c == 0 { return t.dims >= 1 && t.cells >= 1; }
    if c < 1 || c > 10 { return false; }
    return t.cells == 0 && t.points == 0;
}
