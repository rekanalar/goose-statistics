// The multivariate families of SciPy 1.18.1 whose counts follow the
// law of a draw: multinomial, multivariate_hypergeom and dirichlet.
//
// Every family writes its result into a flat buffer and gives the count
// of the values and the code. The draws take a generator, so one seed
// gives one answer.
namespace statistics;

import std;
import math;

// ---- multinomial ----------------------------------------------------------

// The mass of a vector of counts.
fn multinomial_pmf(x: const f64[:], n: i64, p: const f64[:]) -> f64, i64 {
    let c1 = mv_check_vector(x, p.len);
    if c1 != OK { return nan(), c1; }
    if n < 0 { return nan(), ERR_DOMAIN; }
    var total: f64 = 0.0;
    for v in x {
        if v < 0.0 || v != floor(v) { return nan(), ERR_DOMAIN; }
        total = total + v;
    }
    if total != n { return nan(), ERR_DOMAIN; }
    var sp: f64 = 0.0;
    for v in p {
        if is_nan(v) { return nan(), ERR_NAN; }
        if v < 0.0 || v > 1.0 { return nan(), ERR_DOMAIN; }
        sp = sp + v;
    }
    if abs(sp - 1.0) > 1e-9 { return nan(), ERR_DOMAIN; }
    // The mass is n! times the product of p_i to the count, over the
    // factorials of the counts.
    let ln, sgn, cl = lgamma(n + 1.0);
    if cl != OK { return nan(), cl; }
    var s = ln;
    for i in 0..x.len {
        if x[i] > 0.0 {
            let t, st, ct = lgamma(x[i] + 1.0);
            if ct != OK { return nan(), ct; }
            s = s - t + x[i] * log(p[i]);
        }
    }
    return exp(s), OK;
}

// The mean of each category.
fn multinomial_mean(out: f64[>..]&, n: i64, p: const f64[:]) -> i64 {
    out.resize(0);
    for v in p {
        if v < 0.0 || v > 1.0 { return ERR_DOMAIN; }
        out.push(n * v);
    }
    return OK;
}

// The covariance matrix of the law, row by row.
fn multinomial_cov(out: f64[>..]&, n: i64, p: const f64[:]) -> i64 {
    out.resize(0);
    let k = p.len;
    for i in 0..k {
        for j in 0..k {
            if i == j { out.push(n * p[i] * (1.0 - p[i])); }
            else { out.push(-n as f64 * p[i] * p[j]); }
        }
    }
    return OK;
}

// The draws of the law. The buffer holds size*k values, row by row.
fn multinomial_rvs(out: f64[>..]&, n: i64, p: const f64[:], size: i64, r: rng&) -> i64 {
    out.resize(0);
    let k = p.len;
    if k == 0 || k > MV_MAX_DIM { return ERR_DOMAIN; }
    if n < 0 || n > MV_MAX_DRAWS { return ERR_DOMAIN; }
    if size < 0 || size > MV_MAX_DRAWS { return ERR_DOMAIN; }
    var total: f64 = 0.0;
    for v in p {
        if v < 0.0 || v > 1.0 { return ERR_DOMAIN; }
        total = total + v;
    }
    if abs(total - 1.0) > 1e-9 { return ERR_DOMAIN; }
    for s in 0..size {
        // The counts follow the law of a sum of Bernoulli draws. The
        // count of each category is a draw of a binomial law.
        var left = n;
        for i in 0..k - 1 {
            let v, c = binom_ppf(rs_unit(r), left, p[i]);
            if c != OK { return c; }
            let cnt = v as i64;
            out.push(cnt);
            left = left - cnt;
        }
        out.push(left);
    }
    return OK;
}

// ---- multivariate_hypergeom -----------------------------------------------

// The mass of a vector of counts. The first count holds the total.
fn multivariate_hypergeom_pmf(x: const f64[:], n: const f64[:], N: f64) -> f64, i64 {
    let c1 = mv_check_vector(x, n.len);
    if c1 != OK { return nan(), c1; }
    let c2 = check_scalar(N);
    if c2 != OK { return nan(), c2; }
    if N < 0.0 || N != floor(N) { return nan(), ERR_DOMAIN; }
    var total: f64 = 0.0;
    for v in x {
        if v < 0.0 || v != floor(v) { return nan(), ERR_DOMAIN; }
        total = total + v;
    }
    if total != N { return nan(), ERR_DOMAIN; }
    var nt: f64 = 0.0;
    for v in n {
        if v < 0.0 || v != floor(v) { return nan(), ERR_DOMAIN; }
        nt = nt + v;
    }
    if nt < N { return nan(), ERR_DOMAIN; }
    // The mass is the product of C(n_i, x_i) over C(sum n, N). A count
    // above its pool has mass 0.
    var s = 0.0;
    for i in 0..x.len {
        if x[i] > n[i] { return 0.0, OK; }
        let v, cc = lchoose(n[i], x[i]);
        if cc != OK { return nan(), cc; }
        s = s + v;
    }
    let t, ct = lchoose(nt, N);
    if ct != OK { return nan(), ct; }
    return exp(s - t), OK;
}

// The mean of each category.
fn multivariate_hypergeom_mean(out: f64[>..]&, n: const f64[:], N: f64) -> i64 {
    out.resize(0);
    var nt: f64 = 0.0;
    for v in n {
        if v < 0.0 { return ERR_DOMAIN; }
        nt = nt + v;
    }
    if nt <= 0.0 { return ERR_DOMAIN; }
    for v in n { out.push(N * v / nt); }
    return OK;
}

// ---- dirichlet ------------------------------------------------------------

// The log of the density of the law at a point of the simplex.
fn dirichlet_logpdf(x: const f64[:], alpha: const f64[:]) -> f64, i64 {
    let c = mv_check_vector(x, alpha.len);
    if c != OK { return nan(), c; }
    var a0: f64 = 0.0;
    for v in alpha {
        if is_nan(v) { return nan(), ERR_NAN; }
        if v <= 0.0 { return nan(), ERR_DOMAIN; }
        a0 = a0 + v;
    }
    var s: f64 = 0.0;
    for v in x {
        if v < 0.0 || v > 1.0 { return nan(), ERR_DOMAIN; }
        s = s + v;
    }
    if abs(s - 1.0) > 1e-9 { return nan(), ERR_DOMAIN; }
    let lg, sg, cg = lgamma(a0);
    if cg != OK { return nan(), cg; }
    var l = lg;
    for i in 0..alpha.len {
        let la, sa, ca = lgamma(alpha[i]);
        if ca != OK { return nan(), ca; }
        l = l - la + (alpha[i] - 1.0) * log(x[i]);
    }
    return l, OK;
}

fn dirichlet_pdf(x: const f64[:], alpha: const f64[:]) -> f64, i64 {
    let l, c = dirichlet_logpdf(x, alpha);
    if c != OK { return nan(), c; }
    return exp(l), OK;
}

// The mean of each category.
fn dirichlet_mean(out: f64[>..]&, alpha: const f64[:]) -> i64 {
    out.resize(0);
    var a0: f64 = 0.0;
    for v in alpha {
        if v <= 0.0 { return ERR_DOMAIN; }
        a0 = a0 + v;
    }
    for v in alpha { out.push(v / a0); }
    return OK;
}

// The variance of each category.
fn dirichlet_var(out: f64[>..]&, alpha: const f64[:]) -> i64 {
    out.resize(0);
    var a0: f64 = 0.0;
    for v in alpha {
        if v <= 0.0 { return ERR_DOMAIN; }
        a0 = a0 + v;
    }
    for v in alpha { out.push(v * (a0 - v) / (a0 * a0 * (a0 + 1.0))); }
    return OK;
}

// The covariance matrix of the law, row by row.
fn dirichlet_cov(out: f64[>..]&, alpha: const f64[:]) -> i64 {
    out.resize(0);
    var a0: f64 = 0.0;
    for v in alpha {
        if v <= 0.0 { return ERR_DOMAIN; }
        a0 = a0 + v;
    }
    let k = alpha.len;
    for i in 0..k {
        for j in 0..k {
            if i == j {
                out.push(alpha[i] * (a0 - alpha[i]) / (a0 * a0 * (a0 + 1.0)));
            } else {
                out.push(-alpha[i] * alpha[j] / (a0 * a0 * (a0 + 1.0)));
            }
        }
    }
    return OK;
}

// The draws of the law. The buffer holds size*k values, row by row. The
// draw uses the law of the ratios of the gamma draws.
fn dirichlet_rvs(out: f64[>..]&, alpha: const f64[:], size: i64, r: rng&) -> i64 {
    out.resize(0);
    let k = alpha.len;
    if k == 0 || k > MV_MAX_DIM { return ERR_DOMAIN; }
    if size < 0 || size > MV_MAX_DRAWS { return ERR_DOMAIN; }
    for v in alpha {
        if v <= 0.0 { return ERR_DOMAIN; }
    }
    for s in 0..size {
        var total: f64 = 0.0;
        var g: f64[>..] = [];
        for i in 0..k {
            let v = mv_gamma_draw(alpha[i], r);
            g.push(v);
            total = total + v;
        }
        for i in 0..k { out.push(g[i] / total); }
    }
    return OK;
}

// A draw of the gamma law of the shape a, by the inverse of its
// cumulative function. The scale is 1.
fn mv_gamma_draw(a: f64, r: rng&) -> f64, i64 {
    let v, c = gammadist_ppf(rs_unit(r), a, 1.0);
    if c != OK { return 0.0, c; }
    return v, OK;
}
