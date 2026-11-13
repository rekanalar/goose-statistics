// The fault edges of the binned calls: the empty sample, the NaN, the
// infinity, the domain errors, the bounds and the input that stays
// unchanged.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_binned_tests() -> i64 {
    var bad: i64 = 0;
    var stats: f64[>..] = [];
    var edges: f64[>..] = [];
    var shape: f64[>..] = [];
    var binnumber: f64[>..] = [];
    var empty: f64[>..] = [];
    let empty_sample: f64[0] = [];
    let one: f64[1] = [2.0];
    let xs: f64[4] = [1.0, 2.0, 3.0, 4.0];
    let values: f64[4] = [10.0, 20.0, 30.0, 40.0];
    let nan_sample: f64[3] = [1.0, 0.0 / 0.0, 3.0];
    let inf_sample: f64[3] = [1.0, 1.0 / 0.0, 3.0];
    let given: f64[3] = [1.0, 2.0, 3.0];
    let bad_edges: f64[3] = [1.0, 1.0, 2.0];

    // binned_statistic
    let _t1, c1 = statistics::binned_statistic(stats, edges, shape, binnumber, empty, empty_sample[..], empty_sample[..], statistics::STAT_MEAN, 4);
    chk_code(bad, "binned_statistic empty", c1, statistics::ERR_EMPTY);
    let _t2, c2 = statistics::binned_statistic(stats, edges, shape, binnumber, empty, xs[..], values[..], 99, 4);
    chk_code(bad, "binned_statistic unknown statistic", c2, statistics::ERR_DOMAIN);
    let _t3, c3 = statistics::binned_statistic(stats, edges, shape, binnumber, empty, xs[..], values[..], statistics::STAT_MEAN, 0);
    chk_code(bad, "binned_statistic zero bins", c3, statistics::ERR_DOMAIN);
    let _t4, c4 = statistics::binned_statistic(stats, edges, shape, binnumber, empty, xs[..], values[..], statistics::STAT_MEAN, 1048577);
    chk_code(bad, "binned_statistic bins above the bound", c4, statistics::ERR_TOO_LONG);
    let _t5, c5 = statistics::binned_statistic(stats, edges, shape, binnumber, empty, nan_sample[..], nan_sample[..], statistics::STAT_MEAN, 4);
    chk_code(bad, "binned_statistic NaN", c5, statistics::ERR_NAN);
    let _t6, c6 = statistics::binned_statistic(stats, edges, shape, binnumber, empty, inf_sample[..], inf_sample[..], statistics::STAT_MEAN, 4);
    chk_code(bad, "binned_statistic infinity", c6, statistics::ERR_INFINITY);
    let _t7, c7 = statistics::binned_statistic(stats, edges, shape, binnumber, empty, xs[..], one[..], statistics::STAT_MEAN, 4);
    chk_code(bad, "binned_statistic length", c7, statistics::ERR_LENGTH);

    // binned_statistic_edges
    let _t8, c8 = statistics::binned_statistic_edges(stats, edges, shape, binnumber, empty, xs[..], values[..], statistics::STAT_MEAN, bad_edges[..]);
    chk_code(bad, "binned_statistic_edges edges that do not ascend", c8, statistics::ERR_DOMAIN);
    let _t9, c9 = statistics::binned_statistic_edges(stats, edges, shape, binnumber, empty, xs[..], values[..], statistics::STAT_MEAN, one[..]);
    chk_code(bad, "binned_statistic_edges one edge", c9, statistics::ERR_DOMAIN);

    // binned_statistic_2d and binned_statistic_dd
    let three_values: f64[3] = [10.0, 20.0, 30.0];
    let _t10, c10 = statistics::binned_statistic_2d(stats, edges, shape, binnumber, empty, xs[..], three_values[..], values[..], statistics::STAT_MEAN, 2, 2);
    chk_code(bad, "binned_statistic_2d length", c10, statistics::ERR_LENGTH);
    let _t11, c11 = statistics::binned_statistic_2d(stats, edges, shape, binnumber, empty, xs[..], xs[..], values[..], statistics::STAT_MEAN, 0, 2);
    chk_code(bad, "binned_statistic_2d zero bins", c11, statistics::ERR_DOMAIN);
    let _t12, c12 = statistics::binned_statistic_dd(stats, edges, shape, binnumber, empty, xs[..], 3, values[..], statistics::STAT_MEAN, 2);
    chk_code(bad, "binned_statistic_dd dimension mismatch", c12, statistics::ERR_LENGTH);
    let _t13, c13 = statistics::binned_statistic_dd(stats, edges, shape, binnumber, empty, xs[..], 0, values[..], statistics::STAT_MEAN, 2);
    chk_code(bad, "binned_statistic_dd zero dimensions", c13, statistics::ERR_DOMAIN);
    let _t14, c14 = statistics::binned_statistic_dd(stats, edges, shape, binnumber, empty, xs[..], 17, values[..], statistics::STAT_MEAN, 2);
    chk_code(bad, "binned_statistic_dd dimensions above the bound", c14, statistics::ERR_TOO_LONG);

    // More cells than the bound gives the too long fault.
    let _t15b, c15b = statistics::binned_statistic_2d(stats, edges, shape, binnumber, empty, xs[..], xs[..], values[..], statistics::STAT_MEAN, 2050, 2050);
    chk_code(bad, "binned_statistic_2d cells above the bound", c15b, statistics::ERR_TOO_LONG);

    // An empty bin gets the value 0 and a flag of 1.
    let sparse: f64[2] = [1.0, 4.0];
    let sparse_values: f64[2] = [10.0, 40.0];
    let _t15, c15 = statistics::binned_statistic(stats, edges, shape, binnumber, empty, sparse[..], sparse_values[..], statistics::STAT_MEAN, 4);
    var ok = c15 == 0 && stats.len == 4 && empty.len == 4;
    if ok {
        if stats[1] != 0.0 || empty[1] != 1.0 { ok = false; }
        if stats[0] != 10.0 || empty[0] != 0.0 { ok = false; }
    }
    check(bad, ok, "binned_statistic, an empty bin");

    // A cell count that would wrap gives the too long fault, never a table.
    let one_row_v: f64[1] = [5.0];
    let _t20, c20 = statistics::binned_statistic_dd(stats, edges, shape, binnumber, empty, xs[..], 4, one_row_v[..], statistics::STAT_COUNT, 1048576);
    chk_code(bad, "binned_statistic_dd cell count that wraps", c20, statistics::ERR_TOO_LONG);
    let _t21, c21 = statistics::binned_statistic_dd(stats, edges, shape, binnumber, empty, xs[..], 4, one_row_v[..], statistics::STAT_COUNT, 55109);
    chk_code(bad, "binned_statistic_dd cell count that wraps to a negative", c21, statistics::ERR_TOO_LONG);
    // A gap so small that the rounding scale overflows keeps the edge value in the last bin.
    let tiny_edges: f64[2] = [0.0, 1e-303];
    let tiny_x: f64[1] = [1e-303];
    let tiny_v: f64[1] = [1.0];
    let _t22, c22 = statistics::binned_statistic_edges(stats, edges, shape, binnumber, empty, tiny_x[..], tiny_v[..], statistics::STAT_COUNT, tiny_edges[..]);
    check(bad, c22 == 0 && stats.len == 1 && stats[0] == 1.0, "binned_statistic_edges, a tiny gap keeps the edge value");

    // The calls leave the list of the caller unchanged.
    let keep: f64[4] = [1.0, 2.0, 3.0, 4.0];
    let before: f64[4] = [1.0, 2.0, 3.0, 4.0];
    let _w1, _c16 = statistics::binned_statistic(stats, edges, shape, binnumber, empty, keep[..], values[..], statistics::STAT_MEAN, 2);
    let _w2, _c17 = statistics::binned_statistic_edges(stats, edges, shape, binnumber, empty, keep[..], values[..], statistics::STAT_MEAN, given[..]);
    var same = true;
    for i in 0..keep.len { if keep[i] != before[i] { same = false; } }
    check(bad, same, "the binned calls leave the sample unchanged");
    return bad;
}
