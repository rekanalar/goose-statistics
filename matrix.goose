// The covariance matrix and the correlation matrix of several columns.
//
// The buffer of the data holds the columns one after the other. The result
// buffer holds the cells row by row. Each cell is the call for the pair of
// columns of that cell. The work of the call is the count of the cells
// times the count of the rows.
namespace statistics;

import std;
import math;

// The checks of a matrix call. The order is the length bound, the empty
// buffer, NaN, infinity, the argument, the count of columns, the length
// of the buffer, the length of a column, the work and the length the call
// needs.
fn check_matrix(data: const f64[:], ncols: i64, ddof: i64) -> i64 {
    if data.len > MAX_INPUT * MAX_COLUMNS { return 9; }
    if data.len == 0 { return 1; }
    for x in data {
        if is_nan(x) { return 3; }
    }
    for x in data {
        if is_inf(x) { return 4; }
    }
    if ddof < 0 || ddof > 1 { return 5; }
    if ncols > MAX_COLUMNS { return 9; }
    if ncols < 1 { return 5; }
    if data.len % ncols != 0 { return 6; }
    let rows = data.len / ncols;
    if rows > MAX_INPUT { return 9; }
    if ncols * ncols * rows > MAX_WORK { return 9; }
    if rows < 2 { return 2; }
    return 0;
}

// The covariance of a pair of columns with the divisor n * (n - ddof).
fn covariance_ddof(xs: const f64[:], ys: const f64[:], ddof: i64) -> f64, i64 {
    let r, c = cross_core(xs, ys, false);
    if c != 0 { return nan(), c; }
    let n = xs.len as f64;
    return covariance_value(r.n_xy, r.exp_x, r.exp_y, n * (n - ddof));
}

// Writes the covariance matrix into out, row by row. The argument ddof is
// 0 for the population form and 1 for the sample form. Returns the count
// of the cells and the code. A fault gives 0.
fn covariance_matrix(out: f64[>..]&, data: const f64[:], ncols: i64, ddof: i64) -> i64, i64 {
    out.resize(0);
    let c = check_matrix(data, ncols, ddof);
    if c != 0 { return 0, c; }
    let rows = data.len / ncols;
    var cells: f64[>..] = [];
    for k in 0..ncols * ncols { cells.push(0.0); }
    for i in 0..ncols {
        for j in i..ncols {
            let a = data[i * rows..i * rows + rows];
            let b = data[j * rows..j * rows + rows];
            let v, c2 = covariance_ddof(a, b, ddof);
            if c2 != 0 { return 0, c2; }
            cells[i * ncols + j] = v;
            cells[j * ncols + i] = v;
        }
    }
    for v in cells { out.push(v); }
    return out.len, 0;
}

// Writes the correlation matrix into out, row by row. A constant column
// gives the undefined fault. The call writes no cell. Returns the count of
// the cells and the code. A fault gives 0.
fn correlation_matrix(out: f64[>..]&, data: const f64[:], ncols: i64) -> i64, i64 {
    out.resize(0);
    let c = check_matrix(data, ncols, 0);
    if c != 0 { return 0, c; }
    let rows = data.len / ncols;
    var cells: f64[>..] = [];
    for k in 0..ncols * ncols { cells.push(0.0); }
    for i in 0..ncols {
        for j in i..ncols {
            let a = data[i * rows..i * rows + rows];
            let b = data[j * rows..j * rows + rows];
            let r, c2 = cross_core(a, b, false);
            if c2 != 0 { return 0, c2; }
            let v, c3 = correlation_value(r.n_xy, r.n_xx, r.n_yy);
            if c3 != 0 { return 0, c3; }
            cells[i * ncols + j] = v;
            cells[j * ncols + i] = v;
        }
    }
    for v in cells { out.push(v); }
    return out.len, 0;
}
