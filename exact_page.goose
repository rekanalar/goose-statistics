// The Page trend test of the package. Each call returns a TestResult and
// then a fault code. A fault gives NaN in every field of the result.
namespace statistics;

import std;
import math;

fn ex_page_row(out: f64[>..]&, k: i64) {
    out.resize(0);
    let a = k * (k + 1) * (k + 2) / 6;
    let b = k * (k + 1) * (2 * k + 1) / 6;
    let w = b + 1;
    let full = 1 << k;
    var pos: i64[>..] = [];
    var count: i64[>..] = [];
    for t in 0..k + 1 { count.push(0); }
    for mask in 0..full {
        var bits: i64 = 0;
        for i in 0..k { if (mask >> i) & 1 == 1 { bits += 1; } }
        pos.push(count[bits]);
        count[bits] += 1;
    }
    var cur: f64[>..] = [];
    for i in 0..w { cur.push(0.0); }
    cur[0] = 1.0;
    for t in 0..k {
        var nxt: f64[>..] = [];
        for i in 0..count[t + 1] * w { nxt.push(0.0); }
        for mask in 0..full {
            var bits: i64 = 0;
            for i in 0..k { if (mask >> i) & 1 == 1 { bits += 1; } }
            if bits != t { continue; }
            let row = pos[mask] * w;
            for v in 0..k {
                if (mask >> v) & 1 == 1 { continue; }
                let add = (t + 1) * (v + 1);
                let nrow = pos[mask | (1 << v)] * w;
                for s in 0..w - add {
                    let x = cur[row + s];
                    if x != 0.0 { nxt[nrow + s + add] += x; }
                }
            }
        }
        cur.resize(0);
        for i in 0..nxt.len { cur.push(nxt[i]); }
    }
    var fact: f64 = 1.0;
    for i in 2..k + 1 { fact = fact * i; }
    for l in a..b + 1 { out.push(cur[l] / fact); }
}

// Signature: page_trend_test(data, columns, ranked = false, predicted = [],
// method = PAGE_AUTO) -> TestResult, code. The data hold one subject for
// each row, row by row, with `columns` values in a row. An empty
// predicted list means the order 1 to columns, and otherwise it holds a
// permutation of 1 to columns. With ranked set, each row holds ranks from
// 1 to columns already. The statistic is L. The df field holds PAGE_EXACT
// or PAGE_ASYMPTOTIC, the method that gave the p-value. The exact form
// needs at most MAX_PAGE_COLUMNS columns. The p-value is the upper tail,
// and the exact form cuts L to a whole number, as SciPy does.
fn page_trend_test(data: const f64[:], columns: i64, ranked: bool = false, predicted: const f64[:] = [], method: i64 = PAGE_AUTO) -> TestResult, i64 {
    var c = check_sample(data, 1);
    if c != OK { return fault_result(c); }
    if method < PAGE_AUTO || method > PAGE_ASYMPTOTIC { return fault_result(ERR_DOMAIN); }
    if columns < 3 { return fault_result(ERR_TOO_SHORT); }
    if data.len % columns != 0 { return fault_result(ERR_LENGTH); }
    let m = data.len / columns;
    if m < 2 { return fault_result(ERR_TOO_SHORT); }
    let n = columns;
    var pred: f64[>..] = [];
    if predicted.len == 0 {
        for j in 0..n { pred.push(j + 1); }
    } else {
        if predicted.len != n { return fault_result(ERR_LENGTH); }
        var seen: i64[>..] = [];
        for j in 0..n { seen.push(0); }
        for v in predicted {
            if v < 1.0 || v > n || v != floor(v) { return fault_result(ERR_DOMAIN); }
            seen[v as i64 - 1] += 1;
        }
        for j in 0..n { if seen[j] != 1 { return fault_result(ERR_DOMAIN); } }
        for v in predicted { pred.push(v); }
    }
    var sums: f64[>..] = [];
    for j in 0..n { sums.push(0.0); }
    if ranked {
        for v in data { if v < 1.0 || v > n { return fault_result(ERR_DOMAIN); } }
        for i in 0..m { for j in 0..n { sums[j] += data[i * n + j]; } }
    } else {
        var row: f64[>..] = [];
        var rk: f64[>..] = [];
        for i in 0..m {
            row.resize(0);
            for j in 0..n { row.push(data[i * n + j]); }
            let cnt, cr = rankdata(rk, row[..], TIE_AVERAGE);
            if cr != OK { return fault_result(cr); }
            for j in 0..n { sums[j] += rk[j]; }
        }
    }
    var lacc = acc_new();
    for j in 0..n { acc_add(lacc, pred[j] * sums[j]); }
    let l = acc_total(lacc);
    var exact = method == PAGE_EXACT;
    if method == PAGE_AUTO { exact = !(n > 8 || (m > 12 && n > 3) || m > 20); }
    if !exact {
        let mf = m as f64;
        let nf = n as f64;
        let e0 = mf * nf * (nf + 1.0) * (nf + 1.0) / 4.0;
        let v0 = mf * nf * nf * (nf + 1.0) * (nf * nf - 1.0) / 144.0;
        let pv, cp = norm_sf((l - e0) / sqrt(v0));
        if cp != OK { return fault_result(cp); }
        return good_result(l, pv, PAGE_ASYMPTOTIC, nan());
    }
    if n > MAX_PAGE_COLUMNS { return fault_result(ERR_TOO_LONG); }
    let a = n * (n + 1) * (n + 2) / 6;
    let b = n * (n + 1) * (2 * n + 1) / 6;
    let w1 = (b - a + 1) as f64;
    if m * (b - a) + 1 > MAX_EXACT_CELLS || m * m * w1 * w1 / 2.0 > 500000000.0 { return fault_result(ERR_TOO_LONG); }
    var p1: f64[>..] = [];
    ex_page_row(p1, n);
    var cur: f64[>..] = [];
    for v in p1 { cur.push(v); }
    for r in 1..m {
        var nxt: f64[>..] = [];
        for i in 0..cur.len + p1.len - 1 { nxt.push(0.0); }
        for i in 0..cur.len {
            for j in 0..p1.len { nxt[i + j] += cur[i] * p1[j]; }
        }
        cur.resize(0);
        for v in nxt { cur.push(v); }
    }
    let first = trunc(l) as i64 - m * a;
    var tail = acc_new();
    var from = first;
    if from < 0 { from = 0; }
    for i in from..cur.len { acc_add(tail, cur[i]); }
    return good_result(l, acc_total(tail), PAGE_EXACT, nan());
}
