// The binned statistics: the statistic of the values of each bin of one
// dimension, of two dimensions and of many dimensions. The trimmed and the
// circular calls sit beside them.
//
// Run it with: goose run examples/binning.goose
import std;
import statistics from "..";

fn main() {
    let xs: f64[8] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
    let values: f64[8] = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0];
    var stats: f64[>..] = [];
    var edges: f64[>..] = [];
    var shape: f64[>..] = [];
    var binnumber: f64[>..] = [];
    var empty: f64[>..] = [];

    // The mean of the values of each bin of one dimension.
    let t1, c1 = statistics::binned_statistic(stats, edges, shape, binnumber, empty,
        xs[..], values[..], statistics::STAT_MEAN, 4);
    print("bins: ", stats, "\n");
    print("edges: ", edges, ", code ", c1, "\n");
    print("bin numbers: ", binnumber, "\n");
    print("empty cells: ", empty, ", cells ", t1.cells, ", points ", t1.points, "\n");

    // The statistic of the caller, with a trailing block.
    var stats2: f64[>..] = [];
    var edges2: f64[>..] = [];
    var shape2: f64[>..] = [];
    var binnumber2: f64[>..] = [];
    var empty2: f64[>..] = [];
    let t2, c2 = statistics::binned_statistic_with(stats2, edges2, shape2, binnumber2, empty2,
        xs[..], values[..], 2) { v => v[0] + v[v.len - 1] };
    print("statistic of the caller: ", stats2, ", code ", c2, "\n");

    // Two dimensions.
    var stats3: f64[>..] = [];
    var edges3: f64[>..] = [];
    var shape3: f64[>..] = [];
    var binnumber3: f64[>..] = [];
    var empty3: f64[>..] = [];
    let t3, c3 = statistics::binned_statistic_2d(stats3, edges3, shape3, binnumber3, empty3,
        xs[..], xs[..], values[..], statistics::STAT_SUM, 2, 2);
    print("two dimensions: ", stats3, ", cells ", t3.cells, ", code ", c3, "\n");

    // Many dimensions. The rows hold the coordinates of the points.
    let sample: f64[10] = [1.0, 1.0, 2.0, 2.0, 3.0, 1.0, 4.0, 2.0, 5.0, 1.0];
    let weights: f64[5] = [1.0, 2.0, 3.0, 4.0, 5.0];
    var stats4: f64[>..] = [];
    var edges4: f64[>..] = [];
    var shape4: f64[>..] = [];
    var binnumber4: f64[>..] = [];
    var empty4: f64[>..] = [];
    let t4, c4 = statistics::binned_statistic_dd(stats4, edges4, shape4, binnumber4, empty4,
        sample[..], 2, weights[..], statistics::STAT_COUNT, 2);
    print("many dimensions: ", stats4, ", dimensions ", t4.dims, ", code ", c4, "\n");

    // The trimmed forms of a sample.
    let sample2: f64[8] = [2.5, -1.0, 4.0, 3.25, 10.5, 0.5, 7.75, 3.0];
    var kept: f64[>..] = [];
    let _n, c5 = statistics::trimboth(kept, sample2[..], 0.1);
    print("trimboth: ", kept, ", code ", c5, "\n");
    let sd, c6 = statistics::tstd(sample2[..], 2.0, 8.0, true, true, true, true, 1);
    print("deviation of the values inside the limits: ", sd, ", code ", c6, "\n");

    // The circular mean of a sample of angles.
    let angles: f64[5] = [0.1, 0.5, 1.0, 2.0, 3.0];
    let cm, c7 = statistics::circmean(angles[..], 0.0, 6.283185307179586);
    let cs, c8 = statistics::circstd(angles[..], 0.0, 6.283185307179586, false);
    print("circular mean ", cm, ", circular deviation ", cs, ", codes ", c7, " ", c8, "\n");
}
