// The Siegel estimate of a line: a robust fit that uses the median of the
// slopes.
//
// The slope is the median of the medians of the slopes. The slope of one
// point is the median of the slopes from that point to the other points.
// The slope of the line is the median of those medians.
//
// The intercept of the hierarchical method is the median of the values of
// y minus the slope times x. The intercept of the separate method is the
// median of the medians of the intercepts of the slopes of each point.
//
// Every call returns the slope, the intercept and then a fault code. A
// fault gives NaN in both values.
namespace statistics;

import std;
import math;

// The slopes from one point to the other points, and the intercepts of
// those slopes. A pair with the same x gives no slope. The list of the
// slopes holds the slope from the point to each other point.
fn siegel_row(slopes: f64[>..]&, inter: f64[>..]&, ys: const f64[:], xs: const f64[:], i: i64) {
    slopes.resize(0);
    inter.resize(0);
    for j in 0..ys.len {
        if j == i { continue; }
        let dx = xs[i] - xs[j];
        if dx == 0.0 { continue; }
        let s = (ys[i] - ys[j]) / dx;
        slopes.push(s);
        inter.push(ys[i] - s * xs[i]);
    }
}

// The median of a list that is not sorted. The list is sorted first.
fn siegel_median(out: f64[>..]&, xs: const f64[:]) -> f64 {
    out.resize(0);
    for x in xs { out.push(x); }
    sort(out) { a, b => a < b };
    return median_sorted(out[..]);
}

// The Siegel estimate of the line of two samples. The method is the
// hierarchical method or the separate method.
//
// A sample with fewer than two values gives the too short fault. A first
// list whose values are all equal gives the undefined fault, because no
// pair gives a slope. A sample above MAX_SIEGEL_POINTS gives the too long
// fault, because the estimate reads every pair.
fn siegelslopes(ys: const f64[:], xs: const f64[:], method: i64 = SIEGEL_HIERARCHICAL) -> f64, f64, i64 {
    let c = check_pair(xs, ys, 2);
    if c != OK { return nan(), nan(), c; }
    if method != SIEGEL_HIERARCHICAL && method != SIEGEL_SEPARATE { return nan(), nan(), ERR_DOMAIN; }
    if xs.len > MAX_SIEGEL_POINTS { return nan(), nan(), ERR_TOO_LONG; }
    let n = xs.len;
    var slopes: f64[>..] = [];
    var inter: f64[>..] = [];
    var row_slope: f64[>..] = [];
    var row_inter: f64[>..] = [];
    var med_slope: f64[>..] = [];
    var med_inter: f64[>..] = [];
    for i in 0..n {
        siegel_row(slopes, inter, ys, xs, i);
        if slopes.len == 0 { return nan(), nan(), ERR_UNDEFINED; }
        let s = siegel_median(row_slope, slopes[..]);
        if is_nan(s) || is_inf(s) { return nan(), nan(), ERR_OVERFLOW; }
        med_slope.push(s);
        if method == SIEGEL_SEPARATE {
            let b = siegel_median(row_inter, inter[..]);
            if is_nan(b) || is_inf(b) { return nan(), nan(), ERR_OVERFLOW; }
            med_inter.push(b);
        }
    }
    let slope = siegel_median(row_slope, med_slope[..]);
    if is_nan(slope) || is_inf(slope) { return nan(), nan(), ERR_OVERFLOW; }
    var intercept: f64 = 0.0;
    if method == SIEGEL_SEPARATE {
        intercept = siegel_median(row_inter, med_inter[..]);
    } else {
        var res: f64[>..] = [];
        for i in 0..n { res.push(ys[i] - slope * xs[i]); }
        intercept = siegel_median(row_inter, res[..]);
    }
    if is_nan(intercept) || is_inf(intercept) { return nan(), nan(), ERR_OVERFLOW; }
    return slope, intercept, OK;
}
