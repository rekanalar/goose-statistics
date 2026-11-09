// The line fit of two samples and the Theil-Sen estimate.
//
// The line uses the exact sums of the cross product. The slope of the
// line through the origin uses the exact sums of the products. The
// standard errors follow the formulas of SciPy. The Theil-Sen estimate is
// the median of the slopes of all pairs.
namespace statistics;

import std;
import math;

// The least squares line of two samples.
struct Line { slope: f64, intercept: f64, r: f64, r2: f64 }

// The line of a fault. Every field holds NaN.
fn line_nan() -> Line {
    let v = nan();
    return Line { slope: v, intercept: v, r: v, r2: v };
}

// The line of two samples. With the argument proportional the intercept is
// zero and the slope runs through the origin. The line needs two points. A
// constant column of x gives the undefined fault.
fn linear_regression(xs: const f64[:], ys: const f64[:], proportional: bool) -> Line, i64 {
    let c = check_pair(xs, ys, 2);
    if c != 0 { return line_nan(), c; }
    let r, c2 = cross_core(xs, ys, true);
    if c2 != 0 { return line_nan(), c2; }
    let corr, c3 = correlation_value(r.n_xy, r.n_xx, r.n_yy);
    if c3 != 0 { return line_nan(), c3; }
    if proportional {
        if r.s_xx == 0.0 { return line_nan(), 7; }
        let slope = scale_pow2(r.s_xy / r.s_xx, r.exp_x - r.exp_y);
        if is_nan(slope) || is_inf(slope) { return line_nan(), 8; }
        let l = Line { slope: slope, intercept: 0.0, r: corr, r2: corr * corr };
        return l, 0;
    }
    if r.n_xx <= 0.0 { return line_nan(), 7; }
    let slope = scale_pow2(r.n_xy / r.n_xx, r.exp_x - r.exp_y);
    let icpt = scale_pow2(r.intercept, -r.exp_y);
    if is_nan(slope) || is_inf(slope) || is_nan(icpt) || is_inf(icpt) { return line_nan(), 8; }
    let l = Line {
        slope: slope,
        intercept: icpt,
        r: corr,
        r2: corr * corr,
    };
    return l, 0;
}

// The standard error of the slope and the standard error of the intercept.
// The call needs three points. A constant column gives the undefined fault.
fn regression_errors(xs: const f64[:], ys: const f64[:]) -> f64, f64, i64 {
    let c = check_pair(xs, ys, 3);
    if c != 0 { return nan(), nan(), c; }
    let r, c2 = cross_core(xs, ys, false);
    if c2 != 0 { return nan(), nan(), c2; }
    if r.n_xx <= 0.0 { return nan(), nan(), 7; }
    let corr, c3 = correlation_value(r.n_xy, r.n_xx, r.n_yy);
    if c3 != 0 { return nan(), nan(), c3; }
    let n = xs.len as f64;
    let ssx = r.n_xx / n;
    let ssy = r.n_yy / n;
    if ssx <= 0.0 { return nan(), nan(), 7; }
    let r2 = corr * corr;
    let se = sqrt((1.0 - r2) * ssy / ssx / (n - 2.0));
    let slope_err = scale_pow2(se, r.exp_x - r.exp_y);
    let intercept_err = scale_pow2(se * sqrt(r.s_xx / n), -r.exp_y);
    if is_inf(slope_err) || is_nan(slope_err) { return nan(), nan(), 8; }
    if is_inf(intercept_err) || is_nan(intercept_err) { return nan(), nan(), 8; }
    return slope_err, intercept_err, 0;
}

// The Theil-Sen estimate: the slope is the median of the slopes of all
// pairs with a different x. The intercept is the median of y minus the
// slope times the median of x. A sample above MAX_THEIL gives the too long
// fault. A constant x gives the undefined fault.
fn theil_sen(xs: const f64[:], ys: const f64[:]) -> f64, f64, i64 {
    let c = check_pair(xs, ys, 0);
    if c != 0 { return nan(), nan(), c; }
    if xs.len > MAX_THEIL { return nan(), nan(), 9; }
    if xs.len < 2 { return nan(), nan(), 2; }
    let n = xs.len;
    var slopes: f64[>..] = [];
    for i in 0..n {
        for j in i + 1..n {
            let dx = xs[j] - xs[i];
            if dx == 0.0 { continue; }
            slopes.push((ys[j] - ys[i]) / dx);
        }
    }
    if slopes.len == 0 { return nan(), nan(), 7; }
    sort(slopes) { a, b => a < b };
    let slope = median_sorted(slopes[..]);
    let mx, c2 = median(xs);
    if c2 != 0 { return nan(), nan(), c2; }
    let my, c3 = median(ys);
    if c3 != 0 { return nan(), nan(), c3; }
    let intercept = my - slope * mx;
    if is_inf(intercept) || is_nan(intercept) { return nan(), nan(), 8; }
    return slope, intercept, 0;
}
