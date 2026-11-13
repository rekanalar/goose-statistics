// Time tests of the contingency calls on the largest allowed table.
namespace suite;

import std;
import math;
import os;
import statistics from "../..";

fn run_time_contingency(n: i64) -> i64 {
    var bad: i64 = 0;
    var table: f64[>..] = [];
    let side = 1000;
    for i in 0..(side * side) { table.push(1.0 + (i % 5)); }
    var shape: i64[>..] = [side, side];
    var out: f64[>..] = [];
    var t0 = clock();
    let c1 = statistics::expected_freq(out, table[..], shape[..]);
    time_shape(bad, 1.0, c1, "expected_freq", t0);
    t0 = clock();
    let v2, c2 = statistics::association(table[..], side, side, statistics::ASSOC_CRAMER);
    time_shape(bad, v2 + 1.0, c2, "association", t0);
    t0 = clock();
    let c3 = statistics::margins(out, table[..], shape[..], 1);
    time_shape(bad, 1.0, c3, "margins", t0);
    return bad;
}
