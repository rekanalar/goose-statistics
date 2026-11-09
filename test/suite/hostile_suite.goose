// Tests with hostile input. Every call must answer with a fault code
// that the package names, give NaN with every fault, leave the data of the
// caller alone, and finish in the time of an n log n sort.
namespace suite;

import std;
import math;
import os;
import statistics from "../..";

// One value from a pool of values that break naive code.
fn pick(st: u64&) -> f64 {
    let r = rnd(st) % 16;
    if r == 0 { return 0.0 / 0.0; }
    if r == 1 { return 1.0 / 0.0; }
    if r == 2 { return -1.0 / 0.0; }
    if r == 3 { return 1e308; }
    if r == 4 { return -1e308; }
    if r == 5 { return 5e-324; }
    if r == 6 { return 1e16; }
    if r == 7 { return -1e16; }
    if r == 8 { return 0.0; }
    if r == 9 { return -0.0; }
    if r == 10 { return 1.0; }
    if r == 11 { return 1e-310; }
    return ((rnd(st) % 2001) - 1000) / 8.0;
}

// The code that the sample alone fixes: 0, empty, NaN or infinity.
fn sample_code(xs: const f64[:]) -> i64 {
    if xs.len == 0 { return statistics::ERR_EMPTY; }
    for x in xs {
        if is_nan(x) { return statistics::ERR_NAN; }
    }
    for x in xs {
        if is_inf(x) { return statistics::ERR_INFINITY; }
    }
    return 0;
}

// A result and a code must agree: a code of 0 gives a finite number, and
// any other code gives a known fault and NaN. A sample fault fixes the code.
fn judge(bad: i64&, v: f64, c: i64, sc: i64, what: u8[:]) {
    var ok = true;
    if c < 0 || c > 9 { ok = false; }
    if c == 0 && (is_nan(v) || is_inf(v)) { ok = false; }
    if c != 0 && !is_nan(v) { ok = false; }
    if sc != 0 && c != sc { ok = false; }
    if !ok {
        bad += 1;
        if bad <= 12 { print("FAIL hostile ", what, " code ", c, " value ", v, " sample code ", sc, "\n"); }
    }
}

// The same for a call that gives a count: a fault gives a count of 0.
fn judge_count(bad: i64&, n: i64, c: i64, sc: i64, what: u8[:]) {
    var ok = c >= 0 && c <= 9 && n >= 0;
    if c != 0 && n != 0 { ok = false; }
    if sc != 0 && c != sc { ok = false; }
    if !ok {
        bad += 1;
        if bad <= 12 { print("FAIL hostile ", what, " code ", c, " count ", n, " sample code ", sc, "\n"); }
    }
}

fn same(a: f64, b: f64) -> bool {
    return a == b || (is_nan(a) && is_nan(b));
}

// 100,000 seeded inputs through every call.
fn hostile_random() -> i64 {
    var bad: i64 = 0;
    var st: u64 = 20261109;
    for iter in 0..100000 {
        var xs: f64[>..] = [];
        var ws: f64[>..] = [];
        let n = rnd(st) % 9;
        for k in 0..n { xs.push(pick(st)); }
        let wn = if rnd(st) % 8 == 0 { rnd(st) % 9 } else { n };
        for k in 0..wn { ws.push(pick(st)); }
        var before: f64[>..] = [];
        for x in xs { before.push(x); }
        let sc = sample_code(xs[..]);
        let q = pick(st) / 1000.0;
        let method = (rnd(st) % 12) - 2;
        let order = (rnd(st) % 12) - 2;
        let groups = (rnd(st) % 8) - 1;

        let a1, c1 = statistics::sum(xs[..]);
        judge(bad, a1, c1, sc, "sum");
        let a2, c2 = statistics::mean(xs[..]);
        judge(bad, a2, c2, sc, "mean");
        let a3, c3 = statistics::geometric_mean(xs[..]);
        judge(bad, a3, c3, sc, "geometric_mean");
        let a4, c4 = statistics::harmonic_mean(xs[..]);
        judge(bad, a4, c4, sc, "harmonic_mean");
        let a5, c5 = statistics::minimum(xs[..]);
        judge(bad, a5, c5, sc, "minimum");
        let a6, c6 = statistics::maximum(xs[..]);
        judge(bad, a6, c6, sc, "maximum");
        let a7, c7 = statistics::ptp(xs[..]);
        judge(bad, a7, c7, sc, "ptp");
        let a8, c8 = statistics::median(xs[..]);
        judge(bad, a8, c8, sc, "median");
        let a9, c9 = statistics::median_low(xs[..]);
        judge(bad, a9, c9, sc, "median_low");
        let a10, c10 = statistics::median_high(xs[..]);
        judge(bad, a10, c10, sc, "median_high");
        let a11, c11 = statistics::median_grouped(xs[..], q);
        judge(bad, a11, c11, sc, "median_grouped");
        let a12, c12 = statistics::mode(xs[..]);
        judge(bad, a12, c12, sc, "mode");
        var modes: f64[>..] = [];
        let k13, c13 = statistics::multimode(modes, xs[..]);
        judge_count(bad, k13, c13, sc, "multimode");
        let a14, c14 = statistics::pvariance(xs[..]);
        judge(bad, a14, c14, sc, "pvariance");
        let a15, c15 = statistics::variance(xs[..]);
        judge(bad, a15, c15, sc, "variance");
        let a16, c16 = statistics::pstdev(xs[..]);
        judge(bad, a16, c16, sc, "pstdev");
        let a17, c17 = statistics::stdev(xs[..]);
        judge(bad, a17, c17, sc, "stdev");
        let a18, c18 = statistics::sem(xs[..]);
        judge(bad, a18, c18, sc, "sem");
        let a19, c19 = statistics::variation(xs[..]);
        judge(bad, a19, c19, sc, "variation");
        let a20, c20 = statistics::iqr(xs[..], method);
        judge(bad, a20, c20, sc, "iqr");
        let a21, c21 = statistics::quantile(xs[..], q, method);
        judge(bad, a21, c21, sc, "quantile");
        var cuts: f64[>..] = [];
        let k22, c22 = statistics::cut_points(cuts, xs[..], groups, method);
        judge_count(bad, k22, c22, sc, "cut_points");
        let a23, c23 = statistics::moment(xs[..], order);
        judge(bad, a23, c23, sc, "moment");
        let a24, c24 = statistics::central_moment(xs[..], order);
        judge(bad, a24, c24, sc, "central_moment");
        let a25, c25 = statistics::skew(xs[..], rnd(st) % 2 == 0);
        judge(bad, a25, c25, sc, "skew");
        let a26, c26 = statistics::kurtosis(xs[..], rnd(st) % 2 == 0, rnd(st) % 2 == 0);
        judge(bad, a26, c26, sc, "kurtosis");
        var zs: f64[>..] = [];
        let k27, c27 = statistics::zscore(zs, xs[..], order);
        judge_count(bad, k27, c27, sc, "zscore");
        let a28, c28 = statistics::tmean(xs[..], q, q * 2.0);
        judge(bad, a28, c28, sc, "tmean");
        let a29, c29 = statistics::trim_mean(xs[..], q);
        judge(bad, a29, c29, sc, "trim_mean");
        let a30, c30 = statistics::weighted_mean(xs[..], ws[..]);
        judge(bad, a30, c30, sc, "weighted_mean");
        let a31, c31 = statistics::weighted_variance(xs[..], ws[..], method);
        judge(bad, a31, c31, sc, "weighted_variance");
        let a32, c32 = statistics::weighted_median(xs[..], ws[..]);
        judge(bad, a32, c32, sc, "weighted_median");

        // No call changes the data of the caller.
        for k in 0..n {
            if !same(xs[k], before[k]) {
                bad += 1;
                if bad <= 12 { print("FAIL hostile a call changed the sample\n"); }
            }
        }
    }
    return bad;
}

// The list of partials grows with a hostile sample, but it never reaches
// the guard. The sample holds 160 values spread over the exponent range,
// in many orders. The list must stay below the number of bit positions of
// a double.
fn hostile_partials() -> i64 {
    var bad: i64 = 0;
    var st: u64 = 99;
    for round in 0..200 {
        var vals: f64[>..] = [];
        for k in 0..160 { vals.push(pow(2.0, -1070.0 + 52.0 * (k % 40))); }
        for k in 0..160 {
            let j = rnd(st) % 160;
            let t = vals[k];
            vals[k] = vals[j];
            vals[j] = t;
        }
        var a = statistics::acc_new();
        for v in vals { statistics::acc_add(a, v); }
        if a.fault || a.n > 2098 {
            bad += 1;
            print("FAIL partials ", a.n, " fault ", a.fault, "\n");
        }
    }
    return bad;
}

// Fills a sample of the maximum length in one of the crafted shapes.
fn crafted(xs: f64[>..]&, shape: i64) {
    let n = statistics::MAX_INPUT;
    xs.resize(0);
    var st: u64 = 7;
    if shape == 0 { for i in 0..n { xs.push(i); } }
    if shape == 1 { for i in 0..n { xs.push(n - i); } }
    if shape == 2 { for i in 0..n { xs.push(1.0); } }
    if shape == 3 {
        for i in 0..n {
            if i < n / 2 { xs.push(i); } else { xs.push(n - i); }
        }
    }
    if shape == 4 {
        // The median of three killer of Musser: it forces the worst case on
        // a quicksort that takes the median of three as its pivot.
        for i in 0..n { xs.push(0.0); }
        let k = n / 2;
        var i = 1;
        while i <= k {
            if i % 2 == 1 {
                xs[i - 1] = i;
                xs[i] = k + i;
            }
            xs[k + i - 1] = 2 * i;
            i += 1;
        }
    }
    if shape == 5 { for i in 0..n { xs.push(rnd(st) % 1000); } }
    if shape == 6 { for i in 0..n { xs.push(rnd(st)); } }
}

// The sort of a crafted sample of the maximum length takes at most five
// times as long as the sort of a random sample of the same length.
fn hostile_sort_time() -> i64 {
    var bad: i64 = 0;
    var xs: f64[>..] = [];
    crafted(xs, 6);
    var t0 = clock();
    let _r, _c = statistics::median(xs[..]);
    let random_time = clock() - t0;
    let names: const u8[:][7] = ["sorted", "reversed", "all equal", "organ pipe", "median of three killer", "few distinct", "random"];
    for shape in 0..6 {
        crafted(xs, shape);
        t0 = clock();
        let _m, c = statistics::median(xs[..]);
        let spent = clock() - t0;
        if c != 0 || spent > 5.0 * random_time + 0.05 {
            bad += 1;
            print("FAIL sort time of ", names[shape], ": ", spent, " seconds, random ", random_time, " seconds\n");
        }
    }
    return bad;
}

fn run_hostile() -> i64 {
    var bad: i64 = 0;
    bad += hostile_random();
    bad += hostile_partials();
    bad += hostile_sort_time();
    return bad;
}
