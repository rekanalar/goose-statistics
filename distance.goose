// The distance calls of the package: the Wasserstein distance, the energy
// distance and the O'Brien transform.
//
// The Wasserstein distance and the energy distance compare two samples.
// The transform changes each group of a sample, so the groups have the
// same spread.
//
// Every call returns the result first and a fault code last. A fault gives
// NaN as the result.
namespace statistics;

import std;
import math;

// The value of a sample with its weight. The sort uses the value, and the
// weight follows the value.
fn dist_pairs(out: Xw[>..]&, xs: const f64[:], ws: const f64[:], weighted: bool) {
    out.resize(0);
    for k in 0..xs.len {
        var w: f64 = 1.0;
        if weighted { w = ws[k]; }
        out.push(Xw { x: xs[k], w: w });
    }
    sort(out) { a, b => a.x < b.x };
}

// The gap of the two cumulative functions over the support of the samples.
// The result is the sum of the p-th power of the gap times the length of
// each interval of the support.
//
// The samples are sorted. The sweep walks the values of both samples in
// order. The cumulative weight of a sample holds the weight of the values
// that are not larger than the current value. The gap at a value is the
// difference of the two shares.
fn dist_gap(u: const Xw[:], v: const Xw[:], p: i64) -> f64 {
    var tu: f64 = 0.0;
    for e in u { tu += e.w; }
    var tv: f64 = 0.0;
    for e in v { tv += e.w; }
    var i: i64 = 0;
    var j: i64 = 0;
    var cu: f64 = 0.0;
    var cv: f64 = 0.0;
    var sum: f64 = 0.0;
    while i < u.len || j < v.len {
        // The next value of the support of the two samples.
        var x: f64 = 0.0;
        if i < u.len { x = u[i].x; }
        if j < v.len && (i >= u.len || v[j].x < x) { x = v[j].x; }
        // The weights of the values that equal the current value.
        while i < u.len && u[i].x == x { cu += u[i].w; i += 1; }
        while j < v.len && v[j].x == x { cv += v[j].w; j += 1; }
        // The next value of the support. The gap of the cumulative
        // functions holds on the interval from the current value to the
        // next one.
        if i >= u.len && j >= v.len { break; }
        var x2: f64 = 0.0;
        if i < u.len { x2 = u[i].x; }
        if j < v.len && (i >= u.len || v[j].x < x2) { x2 = v[j].x; }
        let d = cu / tu - cv / tv;
        let gap = x2 - x;
        if p == 1 { sum += abs(d) * gap; }
        else { sum += d * d * gap; }
    }
    return sum;
}

// The shared body of the four distance calls. The power of the gap is 1
// for the Wasserstein distance and 2 for the energy distance.
fn dist_run(u: const f64[:], v: const f64[:], uw: const f64[:], vw: const f64[:], weighted: bool, p: i64) -> f64, i64 {
    let c = check_two(u, v, 1);
    if c != OK { return nan(), c; }
    if weighted {
        let cw = check_weights(uw, u.len, 1);
        if cw != OK { return nan(), cw; }
        let cv = check_weights(vw, v.len, 1);
        if cv != OK { return nan(), cv; }
    }
    var pu: Xw[>..] = [];
    var pv: Xw[>..] = [];
    dist_pairs(pu, u, uw, weighted);
    dist_pairs(pv, v, vw, weighted);
    if weighted {
        var su: f64 = 0.0;
        for e in pu { su += e.w; }
        var sv: f64 = 0.0;
        for e in pv { sv += e.w; }
        if is_inf(su) || is_inf(sv) { return nan(), ERR_OVERFLOW; }
    }
    let gap = dist_gap(pu[..], pv[..], p);
    var w = gap;
    if p == 2 { w = sqrt(2.0) * sqrt(gap); }
    if is_nan(w) { return nan(), ERR_UNDEFINED; }
    if is_inf(w) { return nan(), ERR_OVERFLOW; }
    return w, OK;
}

// The Wasserstein distance of two samples. The samples may hold a
// different count of values.
fn wasserstein_distance(u: const f64[:], v: const f64[:]) -> f64, i64 {
    return dist_run(u, v, [], [], false, 1);
}

// The Wasserstein distance of two samples with weights. The weight list of
// a sample must hold one weight for each value.
fn wasserstein_distance_weighted(u: const f64[:], v: const f64[:], uw: const f64[:], vw: const f64[:]) -> f64, i64 {
    return dist_run(u, v, uw, vw, true, 1);
}

// The energy distance of two samples. It is the square root of twice the
// sum of the squared gap of the cumulative functions.
fn energy_distance(u: const f64[:], v: const f64[:]) -> f64, i64 {
    return dist_run(u, v, [], [], false, 2);
}

// The energy distance of two samples with weights.
fn energy_distance_weighted(u: const f64[:], v: const f64[:], uw: const f64[:], vw: const f64[:]) -> f64, i64 {
    return dist_run(u, v, uw, vw, true, 2);
}

// The O'Brien transform of each group of a sample. The sizes array gives
// the length of each group. The transformed values go into out.
//
// The transform of a group of n values with the mean m and the sum of the
// squared deviations s2 is
//
//     r = ((n - 1.5) * n * (x - m)^2 - 0.5 * s2) / ((n - 1) * (n - 2))
//
// The mean of the transformed values is the variance of the group, so a
// test of the spread of the groups can use the transform.
fn obrientransform(out: f64[>..]&, groups: const f64[:], sizes: const i64[:]) -> i64 {
    out.resize(0);
    let c = check_groups_n(groups, sizes, 1, 3);
    if c != OK { return c; }
    var at: i64 = 0;
    for g in 0..sizes.len {
        let n = sizes[g];
        let part = groups[at..at + n];
        at += n;
        let m, cm = exact_mean(part);
        if cm != OK { return cm; }
        var sq: f64[>..] = [];
        for x in part { sq.push((x - m) * (x - m)); }
        let s2, cs = exact_sum_div(sq[..], 1.0);
        if cs != OK { return cs; }
        let nf: f64 = n;
        let den = (nf - 1.0) * (nf - 2.0);
        for x in part {
            let d = x - m;
            let r = ((nf - 1.5) * nf * d * d - 0.5 * s2) / den;
            if is_nan(r) || is_inf(r) { return ERR_OVERFLOW; }
            out.push(r);
        }
    }
    return OK;
}
