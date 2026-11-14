// Five continuous families of the search kind: loguniform, irwinhall,
// gausshyper, geninvgauss and genhyperbolic. The names and the
// parameters follow SciPy 1.18.1.
//
// Every family has three calls: pdf, cdf and ppf. The quantile comes
// from a bracketed root of the cumulative function. Every call returns
// the result first and a fault code last. A value outside the support
// gives a density of 0.
//
// The cumulative function of gausshyper, geninvgauss and genhyperbolic
// comes from the quadrature of the density. The tail side of the
// quadrature keeps the digits of a small tail.
namespace statistics;

import std;
import math;

const IRWINHALL_MAX_N: f64 = 1000.0;   // the largest number of terms

// ---- loguniform: the law of a uniform law on a logarithm ------------------

fn loguniform_check(a: f64, b: f64) -> i64 {
    let c1 = check_shape(a);
    if c1 != OK { return c1; }
    let c2 = check_shape(b);
    if c2 != OK { return c2; }
    if b <= a { return ERR_DOMAIN; }
    return OK;
}

fn loguniform_pdf(x: f64, a: f64, b: f64) -> f64, i64 {
    let k = loguniform_check(a, b);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x < a || x > b { return 0.0, OK; }
    return 1.0 / (x * log(b / a)), OK;
}

fn loguniform_cdf(x: f64, a: f64, b: f64) -> f64, i64 {
    let k = loguniform_check(a, b);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= a { return 0.0, OK; }
    if x >= b { return 1.0, OK; }
    return log(x / a) / log(b / a), OK;
}

fn loguniform_ppf(p: f64, a: f64, b: f64) -> f64, i64 {
    let k = loguniform_check(a, b);
    if k != OK { return nan(), k; }
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    if p == 0.0 { return a, OK; }
    if p == 1.0 { return b, OK; }
    return exp(log(a) + p * log(b / a)), OK;
}

// ---- irwinhall: the sum of n uniform laws ---------------------------------

fn irwinhall_check(n: f64) -> i64 {
    let c = check_whole_shape(n);
    if c != OK { return c; }
    if n > IRWINHALL_MAX_N { return ERR_DOMAIN; }
    return OK;
}

// The density and the cumulative function at x, 0 <= x <= n. The cardinal
// B-spline of the order m has the values B_m(f + j) for j from 0 to m - 1,
// where f is the part of x below 1. The recursion
// B_m(t) = (t B_(m-1)(t) + (m - t) B_(m-1)(t - 1)) / (m - 1) adds only
// positive terms, so no digits cancel. The density is the spline of the
// order n. The cumulative function is the sum of the spline of the order
// n + 1 over the whole steps below x.
fn irwinhall_both(x: f64, n: i64) -> f64, f64 {
    let fl = floor(x) as i64;
    let f = x - floor(x);
    var v: f64[>..] = [];
    v.push(1.0);
    var dens: f64 = 0.0;
    if n == 1 { dens = 1.0; }
    for m in 2..n + 2 {
        var next: f64[>..] = [];
        for j in 0..m {
            var a: f64 = 0.0;
            var b: f64 = 0.0;
            if j < m - 1 { a = (f + j) * v[j]; }
            if j > 0 { b = (m - f - j) * v[j - 1]; }
            next.push((a + b) / (m - 1) as f64);
        }
        v = copy(next);
        if m == n && fl < v.len { dens = v[fl]; }
    }
    var cdf: f64 = 0.0;
    for j in 0..min(fl + 1, v.len) { cdf += v[j]; }
    return dens, cdf;
}

fn irwinhall_pdf(x: f64, n: f64) -> f64, i64 {
    let k = irwinhall_check(n);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 || x >= n { return 0.0, OK; }
    let d, c = irwinhall_both(x, n as i64);
    return d, OK;
}

fn irwinhall_cdf(x: f64, n: f64) -> f64, i64 {
    let k = irwinhall_check(n);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    if x >= n { return 1.0, OK; }
    // The upper half uses the symmetry of the law.
    if x > 0.5 * n {
        let d, c = irwinhall_both(n - x, n as i64);
        return 1.0 - c, OK;
    }
    let d, c = irwinhall_both(x, n as i64);
    return c, OK;
}

fn irwinhall_cdf_value(x: f64, n: f64) -> f64 {
    let v, c = irwinhall_cdf(x, n);
    if c != OK { return nan(); }
    return v;
}

fn irwinhall_ppf(p: f64, n: f64) -> f64, i64 {
    let k = irwinhall_check(n);
    if k != OK { return nan(), k; }
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    if p == 0.0 { return 0.0, OK; }
    if p == 1.0 { return n, OK; }
    let x, cx = brent_root(0.0, n, 1e-14, 4e-16) { v => irwinhall_cdf_value(v, n) - p };
    return x, cx;
}

// ---- gausshyper: the Gauss hypergeometric law -----------------------------

fn gausshyper_check(a: f64, b: f64, c: f64, z: f64) -> i64 {
    let c1 = check_shape(a);
    if c1 != OK { return c1; }
    let c2 = check_shape(b);
    if c2 != OK { return c2; }
    let c3 = check_scalar(c);
    if c3 != OK { return c3; }
    let c4 = check_scalar(z);
    if c4 != OK { return c4; }
    if z <= -1.0 { return ERR_DOMAIN; }
    return OK;
}

// The density without its constant.
fn gausshyper_raw(x: f64, a: f64, b: f64, c: f64, z: f64) -> f64 {
    if x <= 0.0 || x >= 1.0 { return 0.0; }
    return exp((a - 1.0) * log(x) + (b - 1.0) * log1p_core(-x) - c * log1p_core(z * x));
}

// The constant of the density is the integral of the raw density over 0
// to 1, which equals B(a, b) 2F1(c, a; a + b; -z).
fn gausshyper_norm(a: f64, b: f64, c: f64, z: f64) -> f64 {
    return fam_quad_graded(0.0, 1.0) { x => gausshyper_raw(x, a, b, c, z) };
}

fn gausshyper_pdf(x: f64, a: f64, b: f64, c: f64, z: f64) -> f64, i64 {
    let k = gausshyper_check(a, b, c, z);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 || x >= 1.0 { return 0.0, OK; }
    let nrm = gausshyper_norm(a, b, c, z);
    if !(nrm > 0.0) || is_inf(nrm) { return nan(), ERR_OVERFLOW; }
    return gausshyper_raw(x, a, b, c, z) / nrm, OK;
}

fn gausshyper_cdf(x: f64, a: f64, b: f64, c: f64, z: f64) -> f64, i64 {
    let k = gausshyper_check(a, b, c, z);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    if x >= 1.0 { return 1.0, OK; }
    let nrm = gausshyper_norm(a, b, c, z);
    if !(nrm > 0.0) || is_inf(nrm) { return nan(), ERR_OVERFLOW; }
    var v: f64 = 0.0;
    if x <= 0.5 {
        v = fam_quad_graded(0.0, x) { t => gausshyper_raw(t, a, b, c, z) } / nrm;
    } else {
        v = 1.0 - fam_quad_graded(x, 1.0) { t => gausshyper_raw(t, a, b, c, z) } / nrm;
    }
    if v < 0.0 { return 0.0, OK; }
    if v > 1.0 { return 1.0, OK; }
    return v, OK;
}

fn gausshyper_ppf(p: f64, a: f64, b: f64, c: f64, z: f64) -> f64, i64 {
    let k = gausshyper_check(a, b, c, z);
    if k != OK { return nan(), k; }
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    if p == 0.0 { return 0.0, OK; }
    if p == 1.0 { return 1.0, OK; }
    let nrm = gausshyper_norm(a, b, c, z);
    if !(nrm > 0.0) || is_inf(nrm) { return nan(), ERR_OVERFLOW; }
    let x, cx = brent_root(0.0, 1.0, 1e-14, 4e-16) { v => gausshyper_cdf_value(v, a, b, c, z, nrm) - p };
    return x, cx;
}

// The cumulative function with the constant given, without the code.
fn gausshyper_cdf_value(x: f64, a: f64, b: f64, c: f64, z: f64, nrm: f64) -> f64 {
    if x <= 0.0 { return 0.0; }
    if x >= 1.0 { return 1.0; }
    if x <= 0.5 {
        return fam_quad_graded(0.0, x) { t => gausshyper_raw(t, a, b, c, z) } / nrm;
    }
    return 1.0 - fam_quad_graded(x, 1.0) { t => gausshyper_raw(t, a, b, c, z) } / nrm;
}

// ---- geninvgauss: the generalized inverse Gaussian law --------------------

fn geninvgauss_check(p: f64, b: f64) -> i64 {
    let c1 = check_scalar(p);
    if c1 != OK { return c1; }
    return check_shape(b);
}

// The log of the density, without the checks. The value of kve is the
// scaled Bessel function K_p(b) exp(b), which the caller gives.
fn geninvgauss_logdens(x: f64, p: f64, b: f64, lkve: f64) -> f64 {
    if x <= 0.0 { return -inf(); }
    return (p - 1.0) * log(x) - 0.5 * b * (x + 1.0 / x) + b - log(2.0) - lkve;
}

fn geninvgauss_lkve(p: f64, b: f64) -> f64 {
    var bx: f64[>..] = [];
    var bw: f64[>..] = [];
    fam_bessel_nodes(bx, bw);
    return log(fam_kve(p, b, bx[..], bw[..]));
}

fn geninvgauss_pdf(x: f64, p: f64, b: f64) -> f64, i64 {
    let k = geninvgauss_check(p, b);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    let l = geninvgauss_logdens(x, p, b, geninvgauss_lkve(p, b));
    if is_nan(l) { return nan(), ERR_FIT; }
    return exp(l), OK;
}

// The mode of the law. The quasi-density x^(p-1) exp(-b (x + 1/x) / 2)
// has its peak where (p - 1) / x - b (1 - 1/x^2) / 2 is 0.
fn geninvgauss_mode(p: f64, b: f64) -> f64 {
    return ((p - 1.0) + sqrt((p - 1.0) * (p - 1.0) + b * b)) / b;
}

fn geninvgauss_cdf_value(x: f64, p: f64, b: f64, lkve: f64) -> f64 {
    if x <= 0.0 { return 0.0; }
    let m = geninvgauss_mode(p, b);
    var v: f64 = 0.0;
    if x <= m {
        v = fam_quad(0.0, x) { t => exp(geninvgauss_logdens(t, p, b, lkve)) };
    } else {
        v = 1.0 - fam_quad_up(x) { t => exp(geninvgauss_logdens(t, p, b, lkve)) };
    }
    if v < 0.0 { return 0.0; }
    if v > 1.0 { return 1.0; }
    return v;
}

fn geninvgauss_cdf(x: f64, p: f64, b: f64) -> f64, i64 {
    let k = geninvgauss_check(p, b);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    return geninvgauss_cdf_value(x, p, b, geninvgauss_lkve(p, b)), OK;
}

fn geninvgauss_ppf(p: f64, pp: f64, b: f64) -> f64, i64 {
    let k = geninvgauss_check(pp, b);
    if k != OK { return nan(), k; }
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    if p == 0.0 { return 0.0, OK; }
    if p == 1.0 { return inf(), OK; }
    let lk = geninvgauss_lkve(pp, b);
    return fam_search_ppf(p, 0.0, geninvgauss_mode(pp, b)) { v => geninvgauss_cdf_value(v, pp, b, lk) };
}

// ---- genhyperbolic: the generalized hyperbolic law ------------------------

fn genhyperbolic_check(p: f64, a: f64, b: f64) -> i64 {
    let c1 = check_scalar(p);
    if c1 != OK { return c1; }
    let c2 = check_shape(a);
    if c2 != OK { return c2; }
    let c3 = check_scalar(b);
    if c3 != OK { return c3; }
    if p >= 0.0 && abs(b) >= a { return ERR_DOMAIN; }
    if p < 0.0 && abs(b) > a { return ERR_DOMAIN; }
    return OK;
}

// The log of the constant of the density. At |b| = a the Bessel function
// of the constant has no value, and the limit of the ratio gives the
// constant. The limit exists for p below 0.
fn genhyperbolic_logc(p: f64, a: f64, b: f64, bx: const f64[:], bw: const f64[:]) -> f64 {
    let d2 = (a + b) * (a - b);
    let base = -0.5 * log(2.0 * PI) - (p - 0.5) * log(a);
    if d2 <= 0.0 {
        let lg, sg, cg = lgamma(-p);
        return base + (1.0 - (-p)) * log(2.0) - lg + 0.0 * p;
    }
    let d = sqrt(d2);
    return base + 0.5 * p * log(d2) - (log(fam_kve(p, d, bx, bw)) - d);
}

fn genhyperbolic_logdens(x: f64, p: f64, a: f64, b: f64, lc: f64, bx: const f64[:], bw: const f64[:]) -> f64 {
    let sq = sqrt(1.0 + x * x);
    let z = a * sq;
    return lc + b * x + log(fam_kve(p - 0.5, z, bx, bw)) - z - (0.5 - p) * log(sq);
}

fn genhyperbolic_pdf(x: f64, p: f64, a: f64, b: f64) -> f64, i64 {
    let k = genhyperbolic_check(p, a, b);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    var bx: f64[>..] = [];
    var bw: f64[>..] = [];
    fam_bessel_nodes(bx, bw);
    let lc = genhyperbolic_logc(p, a, b, bx[..], bw[..]);
    let l = genhyperbolic_logdens(x, p, a, b, lc, bx[..], bw[..]);
    if is_nan(l) { return nan(), ERR_FIT; }
    return exp(l), OK;
}

// The mean of the law, which splits the quadrature. At |b| = a the split
// is 0.
fn genhyperbolic_split(p: f64, a: f64, b: f64, bx: const f64[:], bw: const f64[:]) -> f64 {
    let d2 = (a + b) * (a - b);
    if d2 <= 0.0 { return 0.0; }
    let d = sqrt(d2);
    return b / d * fam_kve(p + 1.0, d, bx, bw) / fam_kve(p, d, bx, bw);
}

fn genhyperbolic_cdf_value(x: f64, p: f64, a: f64, b: f64, lc: f64, split: f64, bx: const f64[:], bw: const f64[:]) -> f64 {
    var v: f64 = 0.0;
    if x <= split {
        v = fam_quad_down(x) { t => exp(genhyperbolic_logdens(t, p, a, b, lc, bx, bw)) };
    } else {
        v = 1.0 - fam_quad_up(x) { t => exp(genhyperbolic_logdens(t, p, a, b, lc, bx, bw)) };
    }
    if v < 0.0 { return 0.0; }
    if v > 1.0 { return 1.0; }
    return v;
}

fn genhyperbolic_cdf(x: f64, p: f64, a: f64, b: f64) -> f64, i64 {
    let k = genhyperbolic_check(p, a, b);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    var bx: f64[>..] = [];
    var bw: f64[>..] = [];
    fam_bessel_nodes(bx, bw);
    let lc = genhyperbolic_logc(p, a, b, bx[..], bw[..]);
    let split = genhyperbolic_split(p, a, b, bx[..], bw[..]);
    return genhyperbolic_cdf_value(x, p, a, b, lc, split, bx[..], bw[..]), OK;
}

fn genhyperbolic_ppf(q: f64, p: f64, a: f64, b: f64) -> f64, i64 {
    let k = genhyperbolic_check(p, a, b);
    if k != OK { return nan(), k; }
    let pc = check_prob(q);
    if pc != OK { return nan(), pc; }
    if q == 0.0 { return -inf(), OK; }
    if q == 1.0 { return inf(), OK; }
    var bx: f64[>..] = [];
    var bw: f64[>..] = [];
    fam_bessel_nodes(bx, bw);
    let lc = genhyperbolic_logc(p, a, b, bx[..], bw[..]);
    let split = genhyperbolic_split(p, a, b, bx[..], bw[..]);
    // The bracket grows on both sides of the split.
    var lo = split - 1.0;
    var hi = split + 1.0;
    var steps: i64 = 0;
    while steps < MAX_FAMILY_SEARCH && genhyperbolic_cdf_value(lo, p, a, b, lc, split, bx[..], bw[..]) > q {
        lo = split - 2.0 * (split - lo);
        steps += 1;
    }
    while steps < MAX_FAMILY_SEARCH && genhyperbolic_cdf_value(hi, p, a, b, lc, split, bx[..], bw[..]) < q {
        hi = split + 2.0 * (hi - split);
        steps += 1;
    }
    if steps >= MAX_FAMILY_SEARCH { return nan(), ERR_FIT; }
    let x, cx = brent_root(lo, hi, 1e-14, 4e-16) { v => genhyperbolic_cdf_value(v, p, a, b, lc, split, bx[..], bw[..]) - q };
    return x, cx;
}
