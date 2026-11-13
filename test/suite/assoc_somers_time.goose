// Time tests of Somers D on large samples and on a table at its limit.
namespace suite;

import std;
import math;
import os;
import statistics from "../..";

fn run_time_assoc_somers(n: i64) -> i64 {
    var bad: i64 = 0;
    var st: u64 = 8123;
    let xs = time_sample(n, true, st);
    let ys = time_sample(n, true, st);
    let t0 = clock();
    let r, c = statistics::somersd(xs[..], ys[..]);
    time_shape(bad, r.statistic, c, "somersd", t0);
    var table: f64[>..] = [];
    for k in 0..statistics::MAX_TABLE_SIDE * statistics::MAX_TABLE_SIDE { table.push((k % 7) as f64); }
    let t1 = clock();
    let r2, c2 = statistics::somersd_table(table[..], statistics::MAX_TABLE_SIDE, statistics::MAX_TABLE_SIDE);
    time_shape(bad, r2.statistic, c2, "somersd_table", t1);
    return bad;
}
