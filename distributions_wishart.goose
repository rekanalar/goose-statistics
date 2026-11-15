// The multivariate families of SciPy 1.18.1 that need a matrix of a
// scale: wishart, invwishart and normal_inverse_gamma.
//
// Every family writes its result into a flat buffer and gives the count
// of the values and the code. The matrices hold their values in row
// order. The draws take a generator, so one seed gives one answer.
namespace statistics;

import std;
import math;

// ---- wishart --------------------------------------------------------------

// The draws of the law. The buffer holds size*d*d values, matrix by
// matrix, row by row. The matrix is symmetric, so only its lower part
// holds values. The call writes the whole matrix.
fn wishart_rvs(out: f64[>..]&, shape: const f64[:], df: f64, d: i64, size: i64, r: rng&) -> i64 {
    out.resize(0);
    let c1 = mv_check_matrix(shape, d);
    if c1 != OK { return c1; }
    let c2 = check_scalar(df);
    if c2 != OK { return c2; }
    if df < d { return ERR_DOMAIN; }
    if size < 0 || size > MV_MAX_DRAWS { return ERR_DOMAIN; }
    var l: f64[>..] = [];
    let cl = mv_cholesky(l, shape, d);
    if cl != OK { return cl; }
    for s in 0..size {
        // The lower triangular factor of the draw comes from the Bartlett
        // form. The part below the diagonal of each row is a normal draw.
        // The diagonal is the root of a chi law over the degrees of
        // freedom that are left.
        var f: f64[>..] = [];
        for i in 0..d * d { f.push(0.0); }
        for i in 0..d {
            for j in 0..i + 1 {
                if j < i {
                    f[i * d + j] = rs_normal(r);
                } else {
                    // The degrees of freedom that are left, at least one.
                    let g, gc = mv_chisq_draw(df - (i), r);
                    if gc != OK { return gc; }
                    f[i * d + i] = sqrt(g);
                }
            }
        }
        // The draw is L A A' L', where L is the factor of the shape.
        var a: f64[>..] = [];
        for i in 0..d * d { a.push(0.0); }
        for i in 0..d {
            for j in 0..d {
                var v: f64 = 0.0;
                for k in 0..i + 1 { v = v + l[i * d + k] * f[k * d + j]; }
                a[i * d + j] = v;
            }
        }
        var at: f64[>..] = [];
        mv_transpose(at, a, d, d);
        var p: f64[>..] = [];
        mv_matmul(p, a, at, d, d, d);
        for v in p { out.push(v); }
    }
    return OK;
}

// The mean of the law, row by row. It holds d*d values.
fn wishart_mean(out: f64[>..]&, shape: const f64[:], df: f64, d: i64) -> i64 {
    out.resize(0);
    let c1 = mv_check_matrix(shape, d);
    if c1 != OK { return c1; }
    let c2 = check_scalar(df);
    if c2 != OK { return c2; }
    if df < d { return ERR_DOMAIN; }
    for i in 0..d {
        for j in 0..d {
            if i == j { out.push(df * shape[i * d + i]); } else { out.push(0.0); }
        }
    }
    return OK;
}

// ---- invwishart -----------------------------------------------------------

// The draws of the law. The buffer holds size*d*d values, matrix by
// matrix, row by row. The law is the inverse of a Wishart draw.
fn invwishart_rvs(out: f64[>..]&, shape: const f64[:], df: f64, d: i64, size: i64, r: rng&) -> i64 {
    out.resize(0);
    let c1 = mv_check_matrix(shape, d);
    if c1 != OK { return c1; }
    let c2 = check_scalar(df);
    if c2 != OK { return c2; }
    if df <= d { return ERR_DOMAIN; }
    if size < 0 || size > MV_MAX_DRAWS { return ERR_DOMAIN; }
    var w: f64[>..] = [];
    let cw = wishart_rvs(w, shape, df, d, size, r);
    if cw != OK { return cw; }
    for s in 0..size {
        // The inverse of the draw, by the solve of each column.
        var m: f64[>..] = [];
        for i in 0..d * d { m.push(w[s * d * d + i]); }
        var l: f64[>..] = [];
        let cl = mv_cholesky(l, m, d);
        if cl != OK { return cl; }
        for j in 0..d {
            var e: f64[>..] = [];
            for i in 0..d { if i == j { e.push(1.0); } else { e.push(0.0); } }
            // The column of the inverse solves W y = e_j. With W = L L' the
            // solve runs forward on L and back on L'.
            var x: f64[>..] = [];
            mv_solve_lower(l, e, d, x);
            var y: f64[>..] = [];
            mv_solve_lower_transpose(l, x, d, y);
            for i in 0..d { out.push(y[i]); }
        }
    }
    return OK;
}

// ---- normal_inverse_gamma -------------------------------------------------

// The draws of the law. The buffer holds size*(d + 1) values: the
// dimensions of the draw and then the shape of the law.
fn normal_inverse_gamma_rvs(out: f64[>..]&, mu: const f64[:], lam: f64, a: f64, b: f64, d: i64, size: i64, r: rng&) -> i64 {
    out.resize(0);
    let c1 = mv_check_vector(mu, d);
    if c1 != OK { return c1; }
    let c2 = check_shape(lam);
    if c2 != OK { return c2; }
    let c3 = check_shape(a);
    if c3 != OK { return c3; }
    let c4 = check_shape(b);
    if c4 != OK { return c4; }
    if size < 0 || size > MV_MAX_DRAWS { return ERR_DOMAIN; }
    for s in 0..size {
        // The shape of the law is an inverse gamma draw.
        let g, gc2 = mv_igamma_draw(a, b, r);
        if gc2 != OK { return gc2; }
        for i in 0..d {
            let z = rs_normal(r);
            out.push(mu[i] + z * sqrt(lam / g));
        }
        out.push(g);
    }
    return OK;
}

// A draw of the inverse gamma law. The value is the scale over a draw of
// the gamma law of the shape a.
fn mv_igamma_draw(a: f64, b: f64, r: rng&) -> f64, i64 {
    let g, c = mv_gamma_draw(a, r);
    if c != OK { return 0.0, c; }
    if g <= 0.0 { return 0.0, ERR_OVERFLOW; }
    return b / g, OK;
}

