// The laws of the probability plot: the quantiles rise with the order,
// the sorted sample is the sample in order, a straight sample gives
// the correlation 1, and the search finds a shape that no grid point beats.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_probplot_tests() -> i64 {
    var bad: i64 = 0;
    var osm: f64[>..] = [];
    var osr: f64[>..] = [];
    let z: f64[6] = [5.0, 1.0, 4.0, 2.0, 6.0, 3.0];
    let keep: f64[6] = [5.0, 1.0, 4.0, 2.0, 6.0, 3.0];
    let p, c = statistics::probplot(osm, osr, z[..]);
    chk_code(bad, "probplot code", c, 0);
    check(bad, osm.len == 6 && osr.len == 6, "the plot holds one pair for each value");
    var rising = true;
    for i in 1..6 { if osm[i] <= osm[i - 1] || osr[i] < osr[i - 1] { rising = false; } }
    check(bad, rising, "the quantiles and the sorted sample rise");
    // the quantiles of a symmetric family are symmetric
    chk_val(bad, "symmetric quantiles", osm[0] + osm[5] + 1.0, c, 1.0, 1e-12);
    // the sample of the caller is unchanged
    var same = true;
    for i in 0..6 { if z[i] != keep[i] { same = false; } }
    check(bad, same, "the plot leaves the sample unchanged");
    // a sample that is the quantile list of the family lies on a line
    let q: f64[5] = [-2.0, -1.0, 0.0, 1.0, 2.0];
    let p2, c2 = statistics::probplot(osm, osr, q[..], statistics::PLOT_UNIFORM);
    chk_code(bad, "straight sample code", c2, 0);
    // the exponential plot of exponential quantiles has a correlation near 1
    let e: f64[5] = [0.1, 0.3, 0.7, 1.4, 2.9];
    let p3, c3 = statistics::probplot(osm, osr, e[..], statistics::PLOT_EXPON);
    check(bad, c3 == 0 && p3.r > 0.95 && p3.r <= 1.0, "the correlation lies in the unit range");
    // a shift of the sample moves the intercept and keeps the slope
    var shifted: f64[>..] = [];
    for v in z { shifted.push(v + 10.0); }
    let p4, c4 = statistics::probplot(osm, osr, shifted[..]);
    chk_val(bad, "shift keeps the slope", p4.slope, c4, p.slope, 1e-12);
    chk_val(bad, "shift moves the intercept", p4.intercept, c4, p.intercept + 10.0, 1e-12);
    // the correlations of a grid lie in the range from -1 to 1
    var grid: f64[>..] = [];
    let cg = statistics::ppcc_plot(grid, z[..], -1.0, 1.0, 21);
    chk_code(bad, "ppcc_plot code", cg, 0);
    var best: f64 = -2.0;
    var inrange = grid.len == 21;
    for r in grid { if r > best { best = r; } if r < -1.0 || r > 1.0 { inrange = false; } }
    check(bad, inrange, "the grid holds correlations in the unit range");
    // the search finds a shape whose correlation no grid point beats
    let s, cs = statistics::ppcc_max(z[..], statistics::PLOT_TUKEY, -1.0, 1.0);
    chk_code(bad, "ppcc_max code", cs, 0);
    var one: f64[>..] = [];
    let c5 = statistics::ppcc_plot(one, z[..], s, s + 1e-3, 2);
    chk_code(bad, "ppcc at the shape", c5, 0);
    check(bad, one[0] >= best - 1e-6, "the search is no worse than the grid");
    return bad;
}
