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

// The code that a pair of columns alone fixes, in the order of the
// package: empty, length match, NaN, infinity, then the length the call
// needs.
fn pair_code(xs: const f64[:], ys: const f64[:], need: i64) -> i64 {
    if xs.len == 0 || ys.len == 0 { return statistics::ERR_EMPTY; }
    if xs.len != ys.len { return statistics::ERR_LENGTH; }
    for x in xs {
        if is_nan(x) { return statistics::ERR_NAN; }
    }
    for y in ys {
        if is_nan(y) { return statistics::ERR_NAN; }
    }
    for x in xs {
        if is_inf(x) { return statistics::ERR_INFINITY; }
    }
    for y in ys {
        if is_inf(y) { return statistics::ERR_INFINITY; }
    }
    if xs.len < need { return statistics::ERR_TOO_SHORT; }
    return 0;
}

// The code of a pair and its weights: the pair first, then the weight
// length, NaN, infinity, a negative weight, the largest weight and the
// length the call needs.
fn weighted_code(xs: const f64[:], ys: const f64[:], ws: const f64[:], need: i64) -> i64 {
    let c = pair_code(xs, ys, 0);
    if c != 0 { return c; }
    if ws.len != xs.len { return statistics::ERR_LENGTH; }
    for w in ws {
        if is_nan(w) { return statistics::ERR_NAN; }
    }
    for w in ws {
        if is_inf(w) { return statistics::ERR_INFINITY; }
    }
    for w in ws {
        if w < 0.0 { return statistics::ERR_DOMAIN; }
    }
    var top: f64 = 0.0;
    for w in ws {
        if w > top { top = w; }
    }
    if top <= 0.0 { return statistics::ERR_UNDEFINED; }
    if xs.len < need { return statistics::ERR_TOO_SHORT; }
    return 0;
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

        // The association calls take a pair of columns.
        var ys: f64[>..] = [];
        let yn = if rnd(st) % 8 == 0 { rnd(st) % 9 } else { n };
        for k in 0..yn { ys.push(pick(st)); }
        var before_y: f64[>..] = [];
        for y in ys { before_y.push(y); }
        let a33, c33 = statistics::pcovariance(xs[..], ys[..]);
        judge(bad, a33, c33, pair_code(xs[..], ys[..], 1), "pcovariance");
        let a34, c34 = statistics::covariance(xs[..], ys[..]);
        judge(bad, a34, c34, pair_code(xs[..], ys[..], 2), "covariance");
        let a35, c35 = statistics::correlation(xs[..], ys[..]);
        judge(bad, a35, c35, pair_code(xs[..], ys[..], 2), "correlation");
        let a36, c36 = statistics::spearman(xs[..], ys[..]);
        judge(bad, a36, c36, pair_code(xs[..], ys[..], 2), "spearman");
        let tau_b, c37 = statistics::kendall(xs[..], ys[..], statistics::KENDALL_B);
        judge(bad, tau_b, c37, pair_code(xs[..], ys[..], 2), "kendall b");
        let tau_c, c38 = statistics::kendall(xs[..], ys[..], (rnd(st) % 3) - 1);
        judge(bad, tau_c, c38, pair_code(xs[..], ys[..], 2), "kendall");
        let l39, c39 = statistics::linear_regression(xs[..], ys[..], rnd(st) % 2 == 0);
        judge(bad, l39.slope, c39, pair_code(xs[..], ys[..], 2), "linear_regression slope");
        judge(bad, l39.intercept, c39, pair_code(xs[..], ys[..], 2), "linear_regression intercept");
        judge(bad, l39.r, c39, pair_code(xs[..], ys[..], 2), "linear_regression r");
        judge(bad, l39.r2, c39, pair_code(xs[..], ys[..], 2), "linear_regression r2");
        let se, ie, c40 = statistics::regression_errors(xs[..], ys[..]);
        judge(bad, se, c40, pair_code(xs[..], ys[..], 3), "regression_errors slope");
        judge(bad, ie, c40, pair_code(xs[..], ys[..], 3), "regression_errors intercept");
        let ts41, ti41, c41 = statistics::theil_sen(xs[..], ys[..]);
        judge(bad, ts41, c41, pair_code(xs[..], ys[..], 2), "theil_sen slope");
        judge(bad, ti41, c41, pair_code(xs[..], ys[..], 2), "theil_sen intercept");
        let a42, c42 = statistics::weighted_covariance(xs[..], ys[..], ws[..], (rnd(st) % 5) - 1);
        judge(bad, a42, c42, weighted_code(xs[..], ys[..], ws[..], 1), "weighted_covariance");
        let a43, c43 = statistics::weighted_correlation(xs[..], ys[..], ws[..]);
        judge(bad, a43, c43, weighted_code(xs[..], ys[..], ws[..], 2), "weighted_correlation");
        var ranks: f64[>..] = [];
        var cells: f64[>..] = [];
        let k44, c44 = statistics::rankdata(ranks, xs[..], (rnd(st) % 7) - 1);
        judge_count(bad, k44, c44, sc, "rankdata");
        let ddx = (rnd(st) % 3) - 1;
        let k45, c45 = statistics::covariance_matrix(cells, xs[..], (rnd(st) % 5) - 1, ddx);
        judge_count(bad, k45, c45, sc, "covariance_matrix");
        let k46, c46 = statistics::correlation_matrix(cells, xs[..], (rnd(st) % 5) - 1);
        judge_count(bad, k46, c46, sc, "correlation_matrix");

        // No call changes the data of the caller.
        for k in 0..n {
            if !same(xs[k], before[k]) {
                bad += 1;
                if bad <= 12 { print("FAIL hostile a call changed the sample\n"); }
            }
        }
        for k in 0..yn {
            if !same(ys[k], before_y[k]) {
                bad += 1;
                if bad <= 12 { print("FAIL hostile a call changed the second sample\n"); }
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
    bad += hostile_kendall_time();
    bad += hostile_theil_time();
    return bad;
}

// Fills two columns of the given length in one of the crafted shapes.
fn crafted_pair(xs: f64[>..]&, ys: f64[>..]&, shape: i64, n: i64) {
    xs.resize(0);
    ys.resize(0);
    var st: u64 = 11;
    if shape == 0 {
        for i in 0..n { xs.push(i); ys.push(i); }
    }
    if shape == 1 {
        for i in 0..n { xs.push(i); ys.push(n - i); }
    }
    if shape == 2 {
        for i in 0..n { xs.push(1.0); ys.push(1.0); }
    }
    if shape == 3 {
        for i in 0..n {
            xs.push(i);
            if i < n / 2 { ys.push(i); } else { ys.push(n - i); }
        }
    }
    if shape == 4 {
        // The median of three killer: it forces the worst case on a
        // quicksort that takes the median of three as its pivot.
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
        for i in 0..n { ys.push(xs[i]); }
    }
    if shape == 5 {
        for i in 0..n { xs.push(rnd(st) % 1000); ys.push(rnd(st) % 1000); }
    }
}

// The Kendall count must run in the time of a sort. A crafted sample of
// the maximum length must not take much longer than a random sample of
// the same length.
fn hostile_kendall_time() -> i64 {
    var bad: i64 = 0;
    var xs: f64[>..] = [];
    var ys: f64[>..] = [];
    let n = statistics::MAX_INPUT / 16;
    crafted_pair(xs, ys, 5, n);
    var t0 = clock();
    let _t, _c = statistics::kendall(xs[..], ys[..], statistics::KENDALL_B);
    let random_time = clock() - t0;
    let names: const u8[:][6] = ["sorted", "reversed", "all equal", "organ pipe", "median of three killer", "random"];
    for shape in 0..5 {
        crafted_pair(xs, ys, shape, n);
        t0 = clock();
        let _t, c = statistics::kendall(xs[..], ys[..], statistics::KENDALL_B);
        let spent = clock() - t0;
        if c < 0 || c > 9 {
            bad += 1;
            print("FAIL kendall code ", c, "\n");
            continue;
        }
        if spent > 5.0 * random_time + 1.0 {
            bad += 1;
            print("FAIL kendall time of ", names[shape], ": ", spent, " seconds, random ", random_time, " seconds\n");
        }
    }
    return bad;
}

// The Theil-Sen call reaches its bound of 4,096 points. The call must
// finish, and it must stay below the guard of the pair count.
fn hostile_theil_time() -> i64 {
    var bad: i64 = 0;
    let n = statistics::MAX_THEIL;
    var xs: f64[>..] = [];
    var ys: f64[>..] = [];
    var st: u64 = 13;
    for i in 0..n {
        xs.push(i);
        ys.push((rnd(st) % 2001) - 1000);
    }
    let t0 = clock();
    let slope, intercept, c = statistics::theil_sen(xs[..], ys[..]);
    let spent = clock() - t0;
    if c != 0 || is_nan(slope) || is_nan(intercept) {
        bad += 1;
        print("FAIL theil sen at the bound, code ", c, "\n");
    }
    print("theil sen at ", n, " points took ", spent, " seconds\n");
    // One point above the bound gives the too long fault.
    xs.push(0.0);
    ys.push(0.0);
    let _s, _i, c2 = statistics::theil_sen(xs[..], ys[..]);
    if c2 != statistics::ERR_TOO_LONG {
        bad += 1;
        print("FAIL theil sen above the bound, code ", c2, "\n");
    }
    return bad;
}
