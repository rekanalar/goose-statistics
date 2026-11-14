// The representations of a covariance matrix: the factories of the
// Covariance class of SciPy 1.18.1 (from_diagonal, from_precision,
// from_cholesky and from_eigendecomposition) and its calls (covariance,
// log_pdet, rank, whiten and colorize).
//
// A Goose call cannot hold a stored object, so a representation is a kind
// and one flat buffer, and every call takes both. The layout of the buffer
// depends on the kind:
//   COV_DIAGONAL     the d diagonal values;
//   COV_PRECISION    the d*d precision matrix, in row order;
//   COV_CHOLESKY     the d*d lower factor, in row order;
//   COV_EIGEN        the d eigenvalues, then the d*d eigenvectors in row
//                    order, one eigenvector in each column.
// A diagonal and an eigen form may hold a zero or a negative value. Such
// a value drops out of log_pdet and rank, and whiten treats its inverse
// root as 0, as SciPy does. The precision and the Cholesky forms must be
// of full rank.
//
// Every call returns a fault code, or the result and a fault code.
namespace statistics;

import std;
import math;

const COV_DIAGONAL: i64 = 1;
const COV_PRECISION: i64 = 2;
const COV_CHOLESKY: i64 = 3;
const COV_EIGEN: i64 = 4;
const COV_MAX_DIM: i64 = 256;

// The check of a buffer for a kind and a dimension.
fn covariance_check(kind: i64, data: const f64[:], d: i64) -> i64 {
    if d < 1 || d > COV_MAX_DIM { return ERR_DOMAIN; }
    var need: i64 = 0;
    if kind == COV_DIAGONAL { need = d; }
    else if kind == COV_PRECISION || kind == COV_CHOLESKY { need = d * d; }
    else if kind == COV_EIGEN { need = d + d * d; }
    else { return ERR_DOMAIN; }
    if data.len != need { return ERR_LENGTH; }
    for v in data {
        if is_nan(v) { return ERR_NAN; }
        if is_inf(v) { return ERR_INFINITY; }
    }
    if kind == COV_CHOLESKY {
        for i in 0..d {
            if data[i * d + i] <= 0.0 { return ERR_DOMAIN; }
            for j in i + 1..d {
                if data[i * d + j] != 0.0 { return ERR_DOMAIN; }
            }
        }
    }
    if kind == COV_PRECISION {
        for i in 0..d {
            for j in 0..i {
                if abs(data[i * d + j] - data[j * d + i]) > 1e-9 * (abs(data[i * d + j]) + 1.0) {
                    return ERR_DOMAIN;
                }
            }
        }
    }
    return OK;
}

fn cov_check_x(x: const f64[:], d: i64) -> i64 {
    return mv_check_flat(x, d);
}

// The log of the pseudo-determinant. A value that is 0 or less drops out.
fn covariance_log_pdet(kind: i64, data: const f64[:], d: i64) -> f64, i64 {
    let c = covariance_check(kind, data, d);
    if c != OK { return nan(), c; }
    var s: f64 = 0.0;
    if kind == COV_DIAGONAL || kind == COV_EIGEN {
        for i in 0..d { if data[i] > 0.0 { s = s + log(data[i]); } }
        return s, OK;
    }
    if kind == COV_CHOLESKY {
        for i in 0..d { s = s + 2.0 * log(data[i * d + i]); }
        return s, OK;
    }
    var f: f64[>..] = [];
    let cf = mv_cholesky(f, data, d);
    if cf != OK { return nan(), cf; }
    for i in 0..d { s = s - 2.0 * log(f[i * d + i]); }
    return s, OK;
}

// The rank of the matrix.
fn covariance_rank(kind: i64, data: const f64[:], d: i64) -> i64, i64 {
    let c = covariance_check(kind, data, d);
    if c != OK { return 0, c; }
    if kind == COV_DIAGONAL || kind == COV_EIGEN {
        var r: i64 = 0;
        for i in 0..d { if data[i] > 0.0 { r += 1; } }
        return r, OK;
    }
    if kind == COV_PRECISION {
        var f: f64[>..] = [];
        let cf = mv_cholesky(f, data, d);
        if cf != OK { return 0, cf; }
    }
    return d, OK;
}

// The covariance matrix, in row order.
fn covariance_matrix(out: f64[>..]&, kind: i64, data: const f64[:], d: i64) -> i64 {
    out.resize(0);
    let c = covariance_check(kind, data, d);
    if c != OK { return c; }
    if kind == COV_DIAGONAL {
        for i in 0..d {
            for j in 0..d {
                if i == j { out.push(data[i]); } else { out.push(0.0); }
            }
        }
        return OK;
    }
    if kind == COV_CHOLESKY {
        for i in 0..d {
            for j in 0..d {
                var s: f64 = 0.0;
                for k in 0..d { s = s + data[i * d + k] * data[j * d + k]; }
                out.push(s);
            }
        }
        return OK;
    }
    if kind == COV_EIGEN {
        for i in 0..d {
            for j in 0..d {
                var s: f64 = 0.0;
                for k in 0..d { s = s + data[d + i * d + k] * data[k] * data[d + j * d + k]; }
                out.push(s);
            }
        }
        return OK;
    }
    // The inverse of the precision matrix, one column at a time.
    var f: f64[>..] = [];
    let cf = mv_cholesky(f, data, d);
    if cf != OK { return cf; }
    var cols: f64[>..] = [];
    for j in 0..d {
        var e: f64[>..] = [];
        for i in 0..d { if i == j { e.push(1.0); } else { e.push(0.0); } }
        var y: f64[>..] = [];
        mv_solve_lower(f[..], e[..], d, y);
        var z: f64[>..] = [];
        mv_solve_lower_transpose(f[..], y[..], d, z);
        for k in 0..d { cols.push(z[k]); }
    }
    for i in 0..d {
        for j in 0..d { out.push(cols[j * d + i]); }
    }
    return OK;
}

// The whitened point: the point times the inverse of a factor of the
// covariance. The factor depends on the kind.
fn covariance_whiten(out: f64[>..]&, kind: i64, data: const f64[:], d: i64, x: const f64[:]) -> i64 {
    out.resize(0);
    let c = first_code(covariance_check(kind, data, d), cov_check_x(x, d));
    if c != OK { return c; }
    if kind == COV_DIAGONAL {
        for i in 0..d {
            if data[i] > 0.0 { out.push(x[i] / sqrt(data[i])); } else { out.push(0.0); }
        }
        return OK;
    }
    if kind == COV_CHOLESKY {
        mv_solve_lower(data, x, d, out);
        return OK;
    }
    if kind == COV_EIGEN {
        for j in 0..d {
            var s: f64 = 0.0;
            if data[j] > 0.0 {
                for i in 0..d { s = s + x[i] * data[d + i * d + j]; }
                s = s / sqrt(data[j]);
            }
            out.push(s);
        }
        return OK;
    }
    var f: f64[>..] = [];
    let cf = mv_cholesky(f, data, d);
    if cf != OK { return cf; }
    for j in 0..d {
        var s: f64 = 0.0;
        for i in 0..d { s = s + x[i] * f[i * d + j]; }
        out.push(s);
    }
    return OK;
}

// The coloured point: the inverse of the whitening.
fn covariance_colorize(out: f64[>..]&, kind: i64, data: const f64[:], d: i64, x: const f64[:]) -> i64 {
    out.resize(0);
    let c = first_code(covariance_check(kind, data, d), cov_check_x(x, d));
    if c != OK { return c; }
    if kind == COV_DIAGONAL {
        for i in 0..d {
            if data[i] < 0.0 { return ERR_DOMAIN; }
            out.push(x[i] * sqrt(data[i]));
        }
        return OK;
    }
    if kind == COV_CHOLESKY {
        for j in 0..d {
            var s: f64 = 0.0;
            for i in 0..d { s = s + x[i] * data[j * d + i]; }
            out.push(s);
        }
        return OK;
    }
    if kind == COV_EIGEN {
        for j in 0..d {
            var s: f64 = 0.0;
            for i in 0..d {
                if data[i] < 0.0 { return ERR_DOMAIN; }
                s = s + x[i] * data[d + j * d + i] * sqrt(data[i]);
            }
            out.push(s);
        }
        return OK;
    }
    // The solve of the transpose of the factor of the precision matrix.
    var f: f64[>..] = [];
    let cf = mv_cholesky(f, data, d);
    if cf != OK { return cf; }
    var z: f64[>..] = [];
    mv_solve_lower_transpose(f[..], x, d, z);
    for v in z { out.push(v); }
    return OK;
}
