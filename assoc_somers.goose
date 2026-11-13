// Somers D of a table or of two samples, with the asymptotic p-value.
//
// D measures the order of Y given X. It is (P - Q) / (N^2 - sum of the
// squared row sums), where P and Q are twice the counts of the concordant
// and the discordant pairs of cells, and N is the total. The p-value uses
// the asymptotic null variance of the statistic: Z = (P - Q) / sqrt(4 S)
// with S = sum(A_ij (C_ij - D_ij)^2) - (P - Q)^2 / N.
//
// A table with a single row or a single column has no order and gives the
// undefined fault. A variance of zero with P different from Q gives the
// p-value 0, and a variance of zero with P equal to Q is undefined.
//
// The call on two samples builds the table of the distinct values, and a
// table that would hold more than MAX_TABLE_CELLS cells gives the too
// long fault.
namespace statistics;

import std;
import math;

// The statistic and the p-value of a checked table.
fn somers_core(table: const f64[:], rows: i64, cols: i64, alt: i64) -> CorrelationResult, i64 {
    var n: f64 = 0.0;
    var rs: f64[>..] = [];
    var cs: f64[>..] = [];
    for i in 0..rows { rs.push(0.0); }
    for j in 0..cols { cs.push(0.0); }
    var nonzero: i64 = 0;
    for i in 0..rows {
        for j in 0..cols {
            let v = table[i * cols + j];
            n += v;
            rs[i] += v;
            cs[j] += v;
            if v > 0.0 { nonzero += 1; }
        }
    }
    if nonzero < 2 { return corr_fault(ERR_DOMAIN), ERR_DOMAIN; }
    if rows < 2 || cols < 2 { return corr_fault(ERR_UNDEFINED), ERR_UNDEFINED; }
    // pre holds the sums of the cells up to a row and a column
    var pre: f64[>..] = [];
    for k in 0..(rows + 1) * (cols + 1) { pre.push(0.0); }
    for i in 0..rows {
        for j in 0..cols {
            pre[(i + 1) * (cols + 1) + j + 1] = table[i * cols + j] + pre[i * (cols + 1) + j + 1]
                + pre[(i + 1) * (cols + 1) + j] - pre[i * (cols + 1) + j];
        }
    }
    var rcum: f64[>..] = [0.0];
    var ccum: f64[>..] = [0.0];
    for i in 0..rows { rcum.push(rcum[i] + rs[i]); }
    for j in 0..cols { ccum.push(ccum[j] + cs[j]); }
    var p: f64 = 0.0;
    var q: f64 = 0.0;
    var s: f64 = 0.0;
    for i in 0..rows {
        for j in 0..cols {
            let a = table[i * cols + j];
            if a == 0.0 { continue; }
            let upleft = pre[i * (cols + 1) + j];
            let lowright = n - rcum[i + 1] - ccum[j + 1] + pre[(i + 1) * (cols + 1) + j + 1];
            let upright = rcum[i] - pre[i * (cols + 1) + j + 1];
            let lowleft = ccum[j] - pre[(i + 1) * (cols + 1) + j];
            let conc = upleft + lowright;
            let disc = upright + lowleft;
            p += a * conc;
            q += a * disc;
            s += a * (conc - disc) * (conc - disc);
        }
    }
    var sri2: f64 = 0.0;
    for r in rs { sri2 += r * r; }
    let den = n * n - sri2;
    if !(den > 0.0) { return corr_fault(ERR_UNDEFINED), ERR_UNDEFINED; }
    let d = (p - q) / den;
    let var_s = s - (p - q) * (p - q) / n;
    var pv: f64 = 0.0;
    if var_s <= 0.0 {
        if p == q { return corr_fault(ERR_UNDEFINED), ERR_UNDEFINED; }
        pv = 0.0;
        if alt == ALT_LESS && p > q { pv = 1.0; }
        if alt == ALT_GREATER && p < q { pv = 1.0; }
    } else {
        let z = (p - q) / sqrt(4.0 * var_s);
        let upper, cu = norm_sf(z);
        let lower, cl = norm_cdf(z);
        if cu != OK { return corr_fault(cu), cu; }
        if cl != OK { return corr_fault(cl), cl; }
        if alt == ALT_LESS { pv = lower; }
        else if alt == ALT_GREATER { pv = upper; }
        else { pv = 2.0 * min(lower, upper); }
    }
    if is_nan(d) || is_inf(d) { return corr_fault(ERR_OVERFLOW), ERR_OVERFLOW; }
    return CorrelationResult { statistic: d, pvalue: clip01(pv) }, OK;
}

// Somers D of a table of whole counts in row order. The table holds
// the counts of X in the rows and of Y in the columns.
fn somersd_table(table: const f64[:], rows: i64, cols: i64, alt: i64 = ALT_TWO_SIDED) -> CorrelationResult, i64 {
    let shape: i64[2] = [rows, cols];
    let c = check_nd_table(table, shape[..]);
    if c != OK { return corr_fault(c), c; }
    let ca = check_alt(alt);
    if ca != OK { return corr_fault(ca), ca; }
    for v in table { if v != floor(v) { return corr_fault(ERR_DOMAIN), ERR_DOMAIN; } }
    return somers_core(table, rows, cols, alt);
}

// The index of a value in a sorted list of distinct values.
fn somers_index(vals: const f64[:], v: f64) -> i64 {
    var lo: i64 = 0;
    var hi: i64 = vals.len - 1;
    while lo < hi {
        let mid = (lo + hi) / 2;
        if vals[mid] < v { lo = mid + 1; } else { hi = mid; }
    }
    return lo;
}

// The sorted distinct values of a sample.
fn somers_levels(out: f64[>..]&, xs: const f64[:]) {
    out.resize(0);
    let s = sort_copy(xs);
    for v in s { if out.len == 0 || out[out.len - 1] != v { out.push(v); } }
}

// Somers D of two samples: the order of ys given xs. A sample with a
// single value has no order and gives the undefined fault.
fn somersd(xs: const f64[:], ys: const f64[:], alt: i64 = ALT_TWO_SIDED) -> CorrelationResult, i64 {
    let c = check_pair(xs, ys, 2);
    if c != OK { return corr_fault(c), c; }
    let ca = check_alt(alt);
    if ca != OK { return corr_fault(ca), ca; }
    var lx: f64[>..] = [];
    var ly: f64[>..] = [];
    somers_levels(lx, xs);
    somers_levels(ly, ys);
    if lx.len * ly.len > MAX_TABLE_CELLS { return corr_fault(ERR_TOO_LONG), ERR_TOO_LONG; }
    var table: f64[>..] = [];
    for k in 0..lx.len * ly.len { table.push(0.0); }
    for k in 0..xs.len {
        let i = somers_index(lx[..], xs[k]);
        let j = somers_index(ly[..], ys[k]);
        table[i * ly.len + j] += 1.0;
    }
    return somers_core(table[..], lx.len, ly.len, alt);
}
