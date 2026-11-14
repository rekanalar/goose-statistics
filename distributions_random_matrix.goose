// Seven more multivariate families of SciPy 1.18.1: vonmises_fisher,
// uniform_direction, ortho_group, special_ortho_group, unitary_group,
// random_correlation and random_table.
//
// The draws take a generator, so one seed gives one answer. Every call
// writes its result into a flat buffer and gives the code. The matrices
// hold their values in row order. A unitary matrix holds its real parts
// first and then its imaginary parts, each in row order.
namespace statistics;

import std;
import math;


// ---- uniform_direction ----------------------------------------------------

// The draws of the law. A draw is a normal vector with its length scaled
// to 1. The buffer holds size*dim values, row by row.
fn uniform_direction_rvs(out: f64[>..]&, dim: i64, size: i64, r: rng&) -> i64 {
    out.resize(0);
    if dim < 1 || dim > MV_MAX_DIM { return ERR_DOMAIN; }
    if size < 0 || size > MV_MAX_DRAWS { return ERR_DOMAIN; }
    for s in 0..size {
        var v: f64[>..] = [];
        var n2: f64 = 0.0;
        for i in 0..dim {
            let z = rs_normal(r);
            v.push(z);
            n2 = n2 + z * z;
        }
        let n = sqrt(n2);
        for i in 0..dim { out.push(v[i] / n); }
    }
    return OK;
}

// The log of the density of the law. The density is the reciprocal of the
// area of the unit sphere, so it does not depend on the point.
fn uniform_direction_logpdf(dim: i64) -> f64, i64 {
    if dim < 1 || dim > MV_MAX_DIM { return nan(), ERR_DOMAIN; }
    let d = dim as f64;
    // The area of the unit sphere is 2 pi^(d/2) / Gamma(d/2).
    let lg, sg, cg = lgamma(0.5 * d);
    if cg != OK { return nan(), cg; }
    return lg - log(2.0) - 0.5 * d * log(PI), OK;
}

// ---- ortho_group, special_ortho_group -------------------------------------

// The factor of the Gram-Schmidt step: the draw of an orthogonal matrix
// from the Haar law. The matrix comes from the QR factor of a normal
// matrix, with the sign of the diagonal of R fixed to plus.
fn rm_haar(out: f64[>..]&, d: i64, r: rng&) {
    out.resize(0);
    for i in 0..d * d { out.push(0.0); }
    // The columns of q are built one at a time from the columns of a
    // normal matrix.
    for j in 0..d {
        var col: f64[>..] = [];
        for i in 0..d { col.push(rs_normal(r)); }
        for k in 0..j {
            var dot: f64 = 0.0;
            for i in 0..d { dot = dot + col[i] * out[i * d + k]; }
            for i in 0..d { col[i] = col[i] - dot * out[i * d + k]; }
        }
        var n2: f64 = 0.0;
        for i in 0..d { n2 = n2 + col[i] * col[i]; }
        let n = sqrt(n2);
        for i in 0..d { out[i * d + j] = col[i] / n; }
    }
}

// The determinant of a matrix, by the elimination with row swaps.
fn rm_det(a: const f64[:], d: i64) -> f64 {
    var m: f64[>..] = [];
    for v in a { m.push(v); }
    var det: f64 = 1.0;
    for c in 0..d {
        var piv = c;
        for i in c + 1..d {
            if abs(m[i * d + c]) > abs(m[piv * d + c]) { piv = i; }
        }
        if m[piv * d + c] == 0.0 { return 0.0; }
        if piv != c {
            for j in 0..d {
                let t = m[c * d + j];
                m[c * d + j] = m[piv * d + j];
                m[piv * d + j] = t;
            }
            det = -det;
        }
        det = det * m[c * d + c];
        for i in c + 1..d {
            let f = m[i * d + c] / m[c * d + c];
            for j in c..d { m[i * d + j] = m[i * d + j] - f * m[c * d + j]; }
        }
    }
    return det;
}

// The draws of the law. The buffer holds size*dim*dim values, matrix by
// matrix, row by row.
fn ortho_group_rvs(out: f64[>..]&, dim: i64, size: i64, r: rng&) -> i64 {
    out.resize(0);
    if dim < 1 || dim > MV_MAX_DIM { return ERR_DOMAIN; }
    if size < 0 || size > MV_MAX_DRAWS { return ERR_DOMAIN; }
    for s in 0..size {
        var q: f64[>..] = [];
        rm_haar(q, dim, r);
        for v in q { out.push(v); }
    }
    return OK;
}

// The draws of the law. The determinant of every matrix is 1.
fn special_ortho_group_rvs(out: f64[>..]&, dim: i64, size: i64, r: rng&) -> i64 {
    out.resize(0);
    if dim < 2 || dim > MV_MAX_DIM { return ERR_DOMAIN; }
    if size < 0 || size > MV_MAX_DRAWS { return ERR_DOMAIN; }
    for s in 0..size {
        var q: f64[>..] = [];
        rm_haar(q, dim, r);
        // A matrix with a determinant of -1 changes the sign of its first
        // column.
        if rm_det(q[..], dim) < 0.0 {
            for i in 0..dim { q[i * dim] = -q[i * dim]; }
        }
        for v in q { out.push(v); }
    }
    return OK;
}

// ---- unitary_group --------------------------------------------------------

// The draws of the law. A matrix takes 2*dim*dim values: the real parts
// then the imaginary parts, each in row order. The matrix comes from the
// Gram-Schmidt step on the columns of a complex normal matrix.
fn unitary_group_rvs(out: f64[>..]&, dim: i64, size: i64, r: rng&) -> i64 {
    out.resize(0);
    if dim < 1 || dim > MV_MAX_DIM { return ERR_DOMAIN; }
    if size < 0 || size > MV_MAX_DRAWS { return ERR_DOMAIN; }
    for s in 0..size {
        var qr: f64[>..] = [];
        var qi: f64[>..] = [];
        for i in 0..dim * dim {
            qr.push(0.0);
            qi.push(0.0);
        }
        for j in 0..dim {
            var cr: f64[>..] = [];
            var ci: f64[>..] = [];
            for i in 0..dim {
                cr.push(rs_normal(r));
                ci.push(rs_normal(r));
            }
            for k in 0..j {
                // The inner product of the column with an earlier column.
                var dr: f64 = 0.0;
                var di: f64 = 0.0;
                for i in 0..dim {
                    dr = dr + qr[i * dim + k] * cr[i] + qi[i * dim + k] * ci[i];
                    di = di + qr[i * dim + k] * ci[i] - qi[i * dim + k] * cr[i];
                }
                for i in 0..dim {
                    cr[i] = cr[i] - (dr * qr[i * dim + k] - di * qi[i * dim + k]);
                    ci[i] = ci[i] - (dr * qi[i * dim + k] + di * qr[i * dim + k]);
                }
            }
            var n2: f64 = 0.0;
            for i in 0..dim { n2 = n2 + cr[i] * cr[i] + ci[i] * ci[i]; }
            let n = sqrt(n2);
            for i in 0..dim {
                qr[i * dim + j] = cr[i] / n;
                qi[i * dim + j] = ci[i] / n;
            }
        }
        for v in qr { out.push(v); }
        for v in qi { out.push(v); }
    }
    return OK;
}

// ---- random_correlation ---------------------------------------------------

// The draw of a correlation matrix with the given eigenvalues. The
// eigenvalues must be at least 0 and must sum to their count. The matrix
// comes from the product Q L Q', where Q is a Haar matrix. The call
// then rotates the matrix until its diagonal is 1, as SciPy does.
fn random_correlation_rvs(out: f64[>..]&, eigs: const f64[:], r: rng&) -> i64 {
    out.resize(0);
    let d = eigs.len;
    if d < 1 || d > MV_MAX_DIM { return ERR_DOMAIN; }
    var total: f64 = 0.0;
    for v in eigs {
        if is_nan(v) { return ERR_NAN; }
        if is_inf(v) { return ERR_INFINITY; }
        if v < 0.0 { return ERR_DOMAIN; }
        total = total + v;
    }
    if abs(total - d as f64) > 1e-9 * d as f64 { return ERR_DOMAIN; }
    var q: f64[>..] = [];
    rm_haar(q, d, r);
    var m: f64[>..] = [];
    for i in 0..d {
        for j in 0..d {
            var v: f64 = 0.0;
            for k in 0..d { v = v + q[i * d + k] * eigs[k] * q[j * d + k]; }
            m.push(v);
        }
    }
    // The Givens rotations drive the diagonal to 1, one row at a time, in
    // the order of SciPy.
    for i in 0..d - 1 {
        let aii = m[i * d + i];
        if aii == 1.0 { continue; }
        var j = i + 1;
        if aii > 1.0 {
            while j < d - 1 && !(m[j * d + j] < 1.0) { j += 1; }
        } else {
            while j < d - 1 && !(m[j * d + j] > 1.0) { j += 1; }
        }
        let ajj = m[j * d + j];
        let aij = m[i * d + j];
        var co: f64 = 0.0;
        var si: f64 = 1.0;
        if ajj != 1.0 {
            let dd = sqrt(max(aij * aij - (aii - 1.0) * (ajj - 1.0), 0.0));
            // The sign of the root follows the sign of the off-diagonal.
            var sd = dd;
            if aij < 0.0 { sd = -dd; }
            let t = (aij + sd) / (ajj - 1.0);
            co = 1.0 / sqrt(1.0 + t * t);
            si = co * t;
        }
        for k in 0..d {
            let x = m[i * d + k];
            let y = m[j * d + k];
            m[i * d + k] = co * x - si * y;
            m[j * d + k] = si * x + co * y;
        }
        for k in 0..d {
            let x = m[k * d + i];
            let y = m[k * d + j];
            m[k * d + i] = co * x - si * y;
            m[k * d + j] = si * x + co * y;
        }
    }
    for i in 0..d {
        if abs(m[i * d + i] - 1.0) > 1e-7 { return ERR_FIT; }
    }
    for v in m { out.push(v); }
    return OK;
}

// ---- vonmises_fisher ------------------------------------------------------

// The log of the normalising factor of the law of the dimension dim and
// the concentration kappa: (dim/2 - 1) log(kappa) - (dim/2) log(2 pi) -
// log(I_{dim/2 - 1}(kappa)). The Bessel function comes from the scaled
// form of the package.
fn vmf_lognorm(dim: i64, kappa: f64) -> f64, i64 {
    let nu = 0.5 * dim as f64 - 1.0;
    let lb = vmf_log_besseli(nu, kappa);
    if is_nan(lb) { return nan(), ERR_FIT; }
    return nu * log(kappa) - 0.5 * dim as f64 * log(2.0 * PI) - lb, OK;
}

// The log of the modified Bessel function of the first kind and the
// order nu, by its power series for a small argument and by the
// asymptotic form for a large one.
fn vmf_log_besseli(nu: f64, x: f64) -> f64 {
    if x < 30.0 {
        let q = 0.25 * x * x;
        var term = exp(nu * log(0.5 * x) - fam_lngamma(nu + 1.0));
        var total = term;
        var k: f64 = 1.0;
        while k < 500.0 {
            term = term * q / (k * (k + nu));
            total = total + term;
            if term < 1e-17 * total { break; }
            k = k + 1.0;
        }
        return log(total);
    }
    // The asymptotic series e^x / sqrt(2 pi x) times a sum.
    let mu = 4.0 * nu * nu;
    var term: f64 = 1.0;
    var total: f64 = 1.0;
    var k: f64 = 1.0;
    while k < 60.0 {
        let f = (mu - (2.0 * k - 1.0) * (2.0 * k - 1.0)) / (k * 8.0 * x);
        let next = -term * f;
        if abs(next) >= abs(term) { break; }
        term = next;
        total = total + term;
        if abs(term) < 1e-17 * abs(total) { break; }
        k = k + 1.0;
    }
    return x - 0.5 * log(2.0 * PI * x) + log(total);
}

// The log of the density of the law at a unit vector.
fn vonmises_fisher_logpdf(x: const f64[:], mu: const f64[:], kappa: f64) -> f64, i64 {
    let d = mu.len;
    if d < 2 || d > MV_MAX_DIM { return nan(), ERR_DOMAIN; }
    if x.len != d { return nan(), ERR_LENGTH; }
    let ck = check_shape_nonneg(kappa);
    if ck != OK { return nan(), ck; }
    var nx: f64 = 0.0;
    var nm: f64 = 0.0;
    var dot: f64 = 0.0;
    for i in 0..d {
        if is_nan(x[i]) || is_nan(mu[i]) { return nan(), ERR_NAN; }
        if is_inf(x[i]) || is_inf(mu[i]) { return nan(), ERR_INFINITY; }
        nx = nx + x[i] * x[i];
        nm = nm + mu[i] * mu[i];
        dot = dot + x[i] * mu[i];
    }
    if abs(nx - 1.0) > 1e-8 || abs(nm - 1.0) > 1e-8 { return nan(), ERR_DOMAIN; }
    if kappa == 0.0 {
        // The law is the uniform law of the sphere.
        return uniform_direction_logpdf(d);
    }
    let ln, cn = vmf_lognorm(d, kappa);
    if cn != OK { return nan(), cn; }
    return ln + kappa * dot, OK;
}

fn vonmises_fisher_pdf(x: const f64[:], mu: const f64[:], kappa: f64) -> f64, i64 {
    let l, c = vonmises_fisher_logpdf(x, mu, kappa);
    if c != OK { return nan(), c; }
    return exp(l), OK;
}

// The draws of the law, by the rejection form of Wood. The buffer holds
// size*dim values, row by row.
fn vonmises_fisher_rvs(out: f64[>..]&, mu: const f64[:], kappa: f64, size: i64, r: rng&) -> i64 {
    out.resize(0);
    let d = mu.len;
    if d < 2 || d > MV_MAX_DIM { return ERR_DOMAIN; }
    let ck = check_shape_nonneg(kappa);
    if ck != OK { return ck; }
    if size < 0 || size > MV_MAX_DRAWS { return ERR_DOMAIN; }
    var nm: f64 = 0.0;
    for v in mu { nm = nm + v * v; }
    if abs(nm - 1.0) > 1e-8 { return ERR_DOMAIN; }
    let dm1 = (d - 1) as f64;
    // The constants of the sampler.
    let b = dm1 / (2.0 * kappa + sqrt(4.0 * kappa * kappa + dm1 * dm1));
    let x0 = (1.0 - b) / (1.0 + b);
    let c = kappa * x0 + dm1 * log(1.0 - x0 * x0);
    for s in 0..size {
        // The component along the mean comes from the rejection step.
        var w: f64 = 0.0;
        if kappa == 0.0 {
            w = 2.0 * rs_unit(r) - 1.0;
        } else {
            var tries: i64 = 0;
            var done = false;
            while !done && tries < 10000 {
                let z = rs_unit(r);
                let u = rs_unit(r);
                w = (1.0 - (1.0 + b) * z) / (1.0 - (1.0 - b) * z);
                if kappa * w + dm1 * log(1.0 - x0 * w) - c >= log(u) { done = true; }
                tries += 1;
            }
            if !done { return ERR_FIT; }
        }
        // The part across the mean is a unit vector of a normal draw.
        var v: f64[>..] = [];
        var n2: f64 = 0.0;
        for i in 0..d {
            let z = rs_normal(r);
            v.push(z);
            n2 = n2 + z * z;
        }
        // The projection onto the plane across the mean.
        var dot: f64 = 0.0;
        for i in 0..d { dot = dot + v[i] * mu[i]; }
        var perp: f64[>..] = [];
        var p2: f64 = 0.0;
        for i in 0..d {
            let t = v[i] - dot * mu[i];
            perp.push(t);
            p2 = p2 + t * t;
        }
        let pn = sqrt(p2);
        let sc = sqrt(max(0.0, 1.0 - w * w));
        for i in 0..d { out.push(w * mu[i] + sc * perp[i] / pn); }
    }
    return OK;
}

// ---- random_table ---------------------------------------------------------

// The draws of the law: a table with the given row sums and column sums,
// from the law of Patefield. The buffer holds size*r*c values, table by
// table, row by row. The call draws each cell from the hypergeometric
// law of what is left.
fn random_table_rvs(out: f64[>..]&, row: const f64[:], col: const f64[:], size: i64, g: rng&) -> i64 {
    out.resize(0);
    let nr = row.len;
    let nc = col.len;
    if nr < 1 || nc < 1 || nr > MV_MAX_DIM || nc > MV_MAX_DIM { return ERR_DOMAIN; }
    if size < 0 || size > MV_MAX_DRAWS { return ERR_DOMAIN; }
    var sr: f64 = 0.0;
    var sc: f64 = 0.0;
    for v in row {
        if is_nan(v) { return ERR_NAN; }
        if v < 0.0 || v != floor(v) { return ERR_DOMAIN; }
        sr = sr + v;
    }
    for v in col {
        if is_nan(v) { return ERR_NAN; }
        if v < 0.0 || v != floor(v) { return ERR_DOMAIN; }
        sc = sc + v;
    }
    if sr != sc { return ERR_DOMAIN; }
    if sr > 1000000.0 { return ERR_DOMAIN; }
    for s in 0..size {
        var cl: f64[>..] = [];
        for v in col { cl.push(v); }
        for i in 0..nr {
            var rleft = row[i];
            for j in 0..nc {
                if j == nc - 1 {
                    out.push(rleft);
                    cl[j] = cl[j] - rleft;
                } else {
                    // The cell is a hypergeometric draw: rleft draws from
                    // the columns that are left, with cl[j] marked.
                    var pool: f64 = 0.0;
                    for k in j..nc { pool = pool + cl[k]; }
                    let cell = rm_hyper_draw(pool, cl[j], rleft, g);
                    out.push(cell);
                    cl[j] = cl[j] - cell;
                    rleft = rleft - cell;
                }
            }
        }
    }
    return OK;
}

// A draw of the hypergeometric law: draws of the count n from a pool of
// the size m with k marked items. The draw takes the items one by one.
fn rm_hyper_draw(m: f64, k: f64, n: f64, g: rng&) -> f64 {
    var left = m;
    var marked = k;
    var c: f64 = 0.0;
    var i: f64 = 0.0;
    while i < n && left > 0.0 {
        if rs_unit(g) * left < marked {
            c = c + 1.0;
            marked = marked - 1.0;
        }
        left = left - 1.0;
        i = i + 1.0;
    }
    return c;
}
