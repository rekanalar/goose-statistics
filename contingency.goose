// The contingency table calls of the package: the expected counts, the
// margins and the measures of association. A table is a flat list in row
// order, with the length of each axis in a second list.
//
// Each call returns a fault code last. A fault leaves no result in the
// output list.
namespace statistics;

import std;
import math;

// The check of a table with any number of axes. The shape lists the
// length of each axis, and its product equals the length of the table.
fn check_nd_table(table: const f64[:], shape: const i64[:]) -> i64 {
    if table.len > MAX_TABLE_CELLS { return ERR_TOO_LONG; }
    if table.len == 0 { return ERR_EMPTY; }
    if shape.len < 1 { return ERR_DOMAIN; }
    if shape.len > MAX_DIMS { return ERR_TOO_LONG; }
    var cells: i64 = 1;
    for s in shape {
        if s < 1 { return ERR_DOMAIN; }
        if s > MAX_TABLE_SIDE { return ERR_TOO_LONG; }
        if cells > MAX_TABLE_CELLS { return ERR_TOO_LONG; }
        cells *= s;
    }
    for v in table {
        if is_nan(v) { return ERR_NAN; }
    }
    for v in table {
        if is_inf(v) { return ERR_INFINITY; }
    }
    if cells != table.len { return ERR_LENGTH; }
    for v in table {
        if v < 0.0 || v > MAX_COUNT { return ERR_DOMAIN; }
    }
    return OK;
}

// The stride of an axis in a table of this shape: the number of cells
// that one step on the axis skips.
fn axis_stride(shape: const i64[:], axis: i64) -> i64 {
    var s: i64 = 1;
    var d = axis + 1;
    while d < shape.len { s *= shape[d]; d += 1; }
    return s;
}

// The sums of the table over all axes but one. The length of out is the
// length of the axis.
fn margins(out: f64[>..]&, table: const f64[:], shape: const i64[:], axis: i64) -> i64 {
    out.resize(0);
    let c = check_nd_table(table, shape);
    if c != OK { return c; }
    if axis < 0 || axis >= shape.len { return ERR_DOMAIN; }
    return margins_core(out, table, shape, axis);
}

// The sums of the table over all axes but one, without the check of the
// table. The length of out is the length of the axis.
fn margins_core(out: f64[>..]&, table: const f64[:], shape: const i64[:], axis: i64) -> i64 {
    out.resize(0);
    if axis < 0 || axis >= shape.len { return ERR_DOMAIN; }
    let stride = axis_stride(shape, axis);
    let side = shape[axis];
    // A compensated sum for each index: the total and the lost low part.
    var hi: f64[>..] = [];
    var lo: f64[>..] = [];
    for i in 0..side { hi.push(0.0); lo.push(0.0); }
    for k in 0..table.len {
        let j = (k / stride) % side;
        let s, e = two_sum(hi[j], table[k]);
        hi[j] = s;
        lo[j] += e;
    }
    for i in 0..side { out.push(hi[i] + lo[i]); }
    return OK;
}

// The expected counts of the table when the axes are independent. A cell
// is the total times the product of the shares of its margins.
fn expected_freq(out: f64[>..]&, table: const f64[:], shape: const i64[:]) -> i64 {
    out.resize(0);
    let c = check_nd_table(table, shape);
    if c != OK { return c; }
    var all = acc_new();
    for v in table { acc_add(all, v); }
    let total = acc_total(all);
    if total == 0.0 { return ERR_UNDEFINED; }
    var shares: f64[>..] = [];
    var offsets: i64[>..] = [];
    var strides: i64[>..] = [];
    for d in 0..shape.len {
        strides.push(axis_stride(shape, d));
        offsets.push(shares.len);
        var m: f64[>..] = [];
        let cm = margins_core(m, table, shape, d);
        if cm != OK { return cm; }
        for v in m { shares.push(v / total); }
    }
    for k in 0..table.len {
        var e = total;
        for d in 0..shape.len {
            let idx = (k / strides[d]) % shape[d];
            e *= shares[offsets[d] + idx];
        }
        out.push(e);
    }
    return OK;
}

// The statistic of a table against its expected counts, with the
// continuity correction of Yates for a table of one degree of freedom.
fn table_stat(table: const f64[:], rows: i64, cols: i64, correction: bool, lam: f64) -> f64, i64 {
    var rs: f64[>..] = [];
    var cs: f64[>..] = [];
    for i in 0..rows { rs.push(0.0); }
    for j in 0..cols { cs.push(0.0); }
    var all = acc_new();
    for i in 0..rows {
        var ra = acc_new();
        for j in 0..cols { acc_add(ra, table[i * cols + j]); }
        rs[i] = acc_total(ra);
        acc_add(all, rs[i]);
    }
    for j in 0..cols {
        var ca = acc_new();
        for i in 0..rows { acc_add(ca, table[i * cols + j]); }
        cs[j] = acc_total(ca);
    }
    let n = acc_total(all);
    let yates = correction && (rows - 1) * (cols - 1) == 1;
    var acc = acc_new();
    for i in 0..rows {
        for j in 0..cols {
            let e = rs[i] * cs[j] / n;
            var o = table[i * cols + j];
            if yates {
                let diff = e - o;
                let mag = min(0.5, abs(diff));
                if diff > 0.0 { o += mag; }
                if diff < 0.0 { o -= mag; }
            }
            let t = power_term(o, e, lam);
            if is_nan(t) || is_inf(t) { return nan(), ERR_UNDEFINED; }
            acc_add(acc, t);
        }
    }
    return acc_total(acc), OK;
}

// The association of the two variables of a table of whole counts: the
// method gives Cramer V, Tschuprow T or the contingency coefficient of
// Pearson. The statistic is the chi-square statistic, or the power
// divergence of the power lam. A table needs two rows and two columns.
fn association(table: const f64[:], rows: i64, cols: i64, method: i64 = ASSOC_CRAMER, correction: bool = false, lam: f64 = 1.0) -> f64, i64 {
    let c = check_table(table, rows, cols);
    if c != OK { return nan(), c; }
    if is_nan(lam) { return nan(), ERR_NAN; }
    if is_inf(lam) { return nan(), ERR_INFINITY; }
    if method < ASSOC_CRAMER || method > ASSOC_PEARSON { return nan(), ERR_DOMAIN; }
    for v in table {
        if v != floor(v) { return nan(), ERR_DOMAIN; }
    }
    if rows < 2 || cols < 2 { return nan(), ERR_TOO_SHORT; }
    let stat, cs = table_stat(table, rows, cols, correction, lam);
    if cs != OK { return nan(), cs; }
    var all = acc_new();
    for v in table { acc_add(all, v); }
    let phi2 = stat / acc_total(all);
    var value: f64 = 0.0;
    if method == ASSOC_CRAMER { value = phi2 / (min(rows, cols) - 1); }
    else if method == ASSOC_TSCHUPROW { value = phi2 / sqrt(((rows - 1) * (cols - 1)) as f64); }
    else { value = phi2 / (1.0 + phi2); }
    let r = sqrt(value);
    if is_nan(r) || is_inf(r) { return nan(), ERR_OVERFLOW; }
    return r, OK;
}
