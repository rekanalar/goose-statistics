// The trimmed and the weighted forms.
//
// A weight list is divided by its largest weight before any sum, so the
// sums of the weights and the products of weight and value cannot
// overflow.
namespace statistics;

import std;
import math;

// The mean of the values from lower to upper, both included. A limit may
// be an infinity. No value inside the limits gives the undefined fault.
fn tmean(xs: const f64[:], lower: f64, upper: f64) -> f64, i64 {
    let c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    if is_nan(lower) || is_nan(upper) { return nan(), 3; }
    var kept: f64[>..] = [];
    for x in xs {
        if x >= lower && x <= upper { kept.push(x); }
    }
    if kept.len == 0 { return nan(), 7; }
    return exact_sum_div(kept[..], kept.len);
}

// The mean after the cut of floor(n * proportion) values from each end.
// The proportion is from 0 up to, but not including, 0.5.
fn trim_mean(xs: const f64[:], proportion: f64) -> f64, i64 {
    let c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    if is_nan(proportion) { return nan(), 3; }
    if proportion < 0.0 || proportion >= 0.5 { return nan(), 5; }
    let s = sort_copy(xs);
    let cut = floor(xs.len * proportion) as i64;
    let rest = s[cut..s.len - cut];
    return exact_sum_div(rest, rest.len);
}

// The largest weight of a checked weight list.
fn top_weight(ws: const f64[:]) -> f64 {
    var top: f64 = 0.0;
    for w in ws {
        if w > top { top = w; }
    }
    return top;
}

// The divisor of the weights. Dividing rounds the weights, so the raw
// weights are kept whenever the products and the sums stay far from the
// limits of the type. Otherwise the weights are divided by the largest.
fn weight_top_for(ws: const f64[:], big: f64) -> f64 {
    let top = top_weight(ws);
    if top >= 1e-100 && top <= 1e100 && big <= 1e100 { return 1.0; }
    return top;
}

fn weight_divisor(xs: const f64[:], ws: const f64[:]) -> f64 {
    var big: f64 = 0.0;
    for x in xs {
        if abs(x) > big { big = abs(x); }
    }
    return weight_top_for(ws, big);
}

// The weighted mean: the sum of w * x over the sum of w.
fn weighted_mean(xs: const f64[:], ws: const f64[:]) -> f64, i64 {
    let c = check_weighted(xs, ws, 1);
    if c != 0 { return nan(), c; }
    return wmean_core(xs, ws, weight_divisor(xs, ws));
}

fn wmean_core(xs: const f64[:], ws: const f64[:], top: f64) -> f64, i64 {
    var wn = acc_new();
    var prod: f64[>..] = [];
    for k in 0..xs.len {
        let w = ws[k] / top;
        acc_add(wn, w);
        prod.push(w * xs[k]);
    }
    if wn.fault { return nan(), 7; }
    return exact_sum_div(prod[..], acc_total(wn));
}

// The weighted variance. The form is one of the WEIGHT_ constants.
fn weighted_variance(xs: const f64[:], ws: const f64[:], form: i64) -> f64, i64 {
    let c = check_weighted(xs, ws, 2);
    if c != 0 { return nan(), c; }
    if form < 0 || form > 2 { return nan(), 5; }
    let top = weight_divisor(xs, ws);
    let m, c2 = wmean_core(xs, ws, top);
    if c2 != 0 { return nan(), c2; }
    var dmax: f64 = 0.0;
    for x in xs {
        let d = abs(x - m);
        if is_inf(d) { return nan(), 8; }
        if d > dmax { dmax = d; }
    }
    var sw = acc_new();
    var sw2 = acc_new();
    for w in ws {
        let u = w / top;
        acc_add(sw, u);
        acc_add(sw2, u * u);
    }
    let wsum = acc_total(sw);
    let w2sum = acc_total(sw2);
    if dmax == 0.0 { return 0.0, 0; }
    let t = pow2_floor(dmax);
    var s1 = acc_new();
    var s2 = acc_new();
    for k in 0..xs.len {
        let u = ws[k] / top;
        let d = (xs[k] - m) / t;
        acc_add(s1, u * d);
        acc_add(s2, u * d * d);
    }
    if s1.fault || s2.fault { return nan(), 7; }
    let a = acc_total(s1);
    var s = acc_total(s2) - a * a / wsum;
    if s < 0.0 { s = 0.0; }
    var r: f64 = 0.0;
    if form == 0 {
        r = (s / wsum) * t * t;
    } else if form == 2 {
        let denom = wsum - w2sum / wsum;
        if denom <= 0.0 { return nan(), 7; }
        r = (s / denom) * t * t;
    } else {
        let total = top * wsum;
        if is_inf(total) {
            r = (s / wsum) * t * t;
        } else {
            let denom = total - 1.0;
            if denom <= 0.0 { return nan(), 7; }
            r = (s / denom) * top * t * t;
        }
    }
    if is_inf(r) || is_nan(r) { return nan(), 8; }
    return r, 0;
}

struct Xw { x: f64, w: f64 }

// The weighted median: the first value whose cumulative weight reaches
// half of the total. It does not average at a tie.
fn weighted_median(xs: const f64[:], ws: const f64[:]) -> f64, i64 {
    let c = check_weighted(xs, ws, 1);
    if c != 0 { return nan(), c; }
    var plain: f64 = 0.0;
    for w in ws { plain += w; }
    var div: f64 = 1.0;
    if is_inf(plain) { div = top_weight(ws); }
    var ps: Xw[>..] = [];
    var total: f64 = 0.0;
    for k in 0..xs.len {
        let w = ws[k] / div;
        total += w;
        ps.push(Xw { x: xs[k], w: w });
    }
    sort(ps) { a, b => a.x < b.x };
    let half = total / 2.0;
    var cum: f64 = 0.0;
    for p in ps {
        cum += p.w;
        if p.w > 0.0 && cum >= half { return p.x, 0; }
    }
    return ps[ps.len - 1].x, 0;
}
