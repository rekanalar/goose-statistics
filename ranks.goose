// The ranks of a sample and the rank correlations.
//
// The ranks come from the sorted order of an index array. The sort is
// stable, so the equal values keep the order of the input. The Kendall
// count sorts the pairs by x and then by y. A bottom up merge sort counts
// the discordant pairs, so the count runs in time of n log n.
namespace statistics;

import std;
import math;

// ---- tie methods ------------------------------------------------------
// The tie methods are the TIE_ constants of statistics.goose.

// ---- Kendall variants --------------------------------------------------
// The Kendall variants are the KENDALL_ constants of statistics.goose.

// A pair of values.
struct Kp { x: f64, y: f64 }

// Writes the rank of every value into out. The ranks are in the order of
// the sample. The method is one of the TIE_ constants. Another value gives
// the domain fault. Returns the count and the code. A fault gives 0.
fn rankdata(out: f64[>..]&, xs: const f64[:], method: i64) -> i64, i64 {
    out.resize(0);
    let c = check_sample(xs, 1);
    if c != 0 { return 0, c; }
    if method < 0 || method > TIE_ORDINAL { return 0, 5; }
    let n = xs.len;
    var idx: Vi[>..] = [];
    for k in 0..n { idx.push(Vi { v: xs[k], i: k }); }
    stable_sort(idx) { a, b => a.v < b.v };
    for k in 0..n { out.push(0.0); }
    var dense: i64 = 0;
    var k: i64 = 0;
    while k < n {
        var j = k + 1;
        while j < n && idx[j].v == idx[k].v { j += 1; }
        let count = j - k;
        let first = (k + 1) as f64;
        for m in k..j {
            var rank = first;
            if method == TIE_AVERAGE { rank = first + (count - 1) / 2.0; }
            if method == TIE_MAX { rank = first + (count - 1); }
            if method == TIE_DENSE { rank = dense + 1; }
            if method == TIE_ORDINAL { rank = first + (m - k); }
            out[idx[m].i] = rank;
        }
        if method == TIE_DENSE { dense += 1; }
        k = j;
    }
    return out.len, 0;
}

// The Spearman coefficient: the Pearson correlation of the average ranks.
// A column of equal values gives the undefined fault.
fn spearman(xs: const f64[:], ys: const f64[:]) -> f64, i64 {
    let c = check_pair(xs, ys, 2);
    if c != 0 { return nan(), c; }
    var rx: f64[>..] = [0.0];
    var ry: f64[>..] = [0.0];
    let n1, c1 = rankdata(rx, xs, TIE_AVERAGE);
    if c1 != 0 { return nan(), c1; }
    let n2, c2 = rankdata(ry, ys, TIE_AVERAGE);
    if c2 != 0 { return nan(), c2; }
    return correlation(rx[..], ry[..]);
}

// The count of the discordant pairs. The merge sort sorts the y column
// and counts the swaps. The sort is bottom up, so no function recurses and
// no growable value crosses a recursion. The merge takes the left value
// when the two values are equal, so a tie in y is not a swap.
fn merge_count(a: f64[>..]&, tmp: f64[>..]&) -> i64 {
    let n = a.len;
    var width: i64 = 1;
    var swaps: i64 = 0;
    while width < n {
        var lo: i64 = 0;
        while lo < n {
            var mid = lo + width;
            if mid > n { mid = n; }
            var hi = lo + 2 * width;
            if hi > n { hi = n; }
            var i = lo;
            var j = mid;
            var k = lo;
            while i < mid && j < hi {
                if a[i] <= a[j] {
                    tmp[k] = a[i];
                    i += 1;
                } else {
                    tmp[k] = a[j];
                    j += 1;
                    swaps += mid - i;
                }
                k += 1;
            }
            while i < mid { tmp[k] = a[i]; i += 1; k += 1; }
            while j < hi { tmp[k] = a[j]; j += 1; k += 1; }
            lo = lo + 2 * width;
        }
        for k in 0..n { a[k] = tmp[k]; }
        width = width * 2;
    }
    return swaps;
}

// The tie counts of a pair. The count n1 holds the pairs tied in x, the
// count n2 the pairs tied in y and the count n3 the pairs tied in both.
// The counts ux and uy hold the distinct values of the two columns.
struct Ties { n1: i64, n2: i64, n3: i64, ux: i64, uy: i64 }

fn kendall_ties(ps: Kp[>..]&, ys: f64[>..]&) -> Ties {
    var t = Ties { .. };
    var k: i64 = 0;
    while k < ps.len {
        var j = k + 1;
        while j < ps.len && ps[j].x == ps[k].x { j += 1; }
        let c = j - k;
        t.n1 += c * (c - 1) / 2;
        t.ux += 1;
        k = j;
    }
    k = 0;
    while k < ps.len {
        var j = k + 1;
        while j < ps.len && ps[j].x == ps[k].x && ps[j].y == ps[k].y { j += 1; }
        let c = j - k;
        t.n3 += c * (c - 1) / 2;
        k = j;
    }
    k = 0;
    while k < ys.len {
        var j = k + 1;
        while j < ys.len && ys[j] == ys[k] { j += 1; }
        let c = j - k;
        t.n2 += c * (c - 1) / 2;
        t.uy += 1;
        k = j;
    }
    return t;
}

// The Kendall coefficient of two samples. The variant is one of the
// KENDALL_ constants. Another value gives the domain fault. A column of
// equal values gives the undefined fault.
fn kendall(xs: const f64[:], ys: const f64[:], variant: i64) -> f64, i64 {
    let c = check_pair(xs, ys, 2);
    if c != 0 { return nan(), c; }
    if variant < 0 || variant > 1 { return nan(), 5; }
    let n = xs.len;
    var ps: Kp[>..] = [];
    for k in 0..n { ps.push(Kp { x: xs[k], y: ys[k] }); }
    sort(ps) { a, b => a.x < b.x || (a.x == b.x && a.y < b.y) };
    var ys2: f64[>..] = [];
    var tmp: f64[>..] = [];
    for p in ps {
        ys2.push(p.y);
        tmp.push(p.y);
    }
    let swaps = merge_count(ys2, tmp);
    let t = kendall_ties(ps, ys2);
    let n0 = n * (n - 1) / 2;
    let cd = n0 - t.n1 - t.n2 + t.n3 - 2 * swaps;
    if variant == 0 {
        let a = n0 - t.n1;
        let b = n0 - t.n2;
        if a <= 0 || b <= 0 { return nan(), 7; }
        let ab: f64 = a * b;
        let den = sqrt(ab);
        if den == 0.0 { return nan(), 7; }
        return cd / den, 0;
    }
    let m = if t.ux < t.uy { t.ux } else { t.uy };
    if m < 2 { return nan(), 7; }
    let nf: f64 = n;
    let mf: f64 = m;
    let den = nf * nf * (mf - 1.0) / mf;
    if den == 0.0 { return nan(), 7; }
    return 2.0 * cd / den, 0;
}
