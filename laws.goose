// Laws built from a family: the truncated law, the law of an order
// statistic and the Dirichlet multinomial law.
//
// A family is a number, LAW_NORMAL and the others, with a location and a
// scale. The uniform family reads its two ends instead. Every call
// returns the result first and a fault code last.
//
// The truncated law keeps the mass of the family between two bounds. The
// bounds may be infinite. A range with no mass that a float can hold
// gives the undefined fault. The mean and the variance come from a
// Gauss-Legendre rule on 64 panels. An infinite bound is cut where the
// tail mass of the family falls below 1e-16. A bound of a Cauchy law that
// is infinite gives the undefined fault, because the mean does not exist.
//
// The law of the r-th order statistic of n draws has the cdf
// I_F(x)(r, n - r + 1). The rank is from 1 to n, and n is at most
// MAX_ORDER_N.
namespace statistics;

import std;
import math;

const LAW_NODES: i64 = 16;                 // nodes of the rule on one panel
const LAW_PANELS: i64 = 32;                // panels on each side of the location
const LAW_TAIL: f64 = 1e-18;               // lower tail mass left out of a moment
const LAW_TOP: f64 = 0.9999999999999999;   // the upper probability of a moment

// ---- the logistic family, with the others behind one switch --------

fn logi_pdf(x: f64, a: f64, b: f64) -> f64, i64 {
    let c = loc_code(x, a, b);
    if c != OK { return nan(), c; }
    let e = exp(-abs((x - a) / b));
    return loc_result(e / (b * (1.0 + e) * (1.0 + e)));
}

fn logi_lower(z: f64) -> f64 {
    if z >= 0.0 { return 1.0 / (1.0 + exp(-z)); }
    let e = exp(z);
    return e / (1.0 + e);
}

fn logi_cdf(x: f64, a: f64, b: f64) -> f64, i64 {
    let c = loc_code(x, a, b);
    if c != OK { return nan(), c; }
    return logi_lower((x - a) / b), OK;
}

fn logi_sf(x: f64, a: f64, b: f64) -> f64, i64 {
    let c = loc_code(x, a, b);
    if c != OK { return nan(), c; }
    return logi_lower(-(x - a) / b), OK;
}

fn logi_ppf(p: f64, a: f64, b: f64) -> f64, i64 {
    let c = loc_code_p(p, a, b);
    if c != OK { return nan(), c; }
    if p == 0.0 || p == 1.0 { return nan(), ERR_DOMAIN; }
    return loc_result(a + b * (log(p) - log(1.0 - p)));
}

fn law_pdf(fam: i64, a: f64, b: f64, x: f64) -> f64, i64 {
    if fam == LAW_NORMAL { return norm_pdf(x, a, b); }
    if fam == LAW_UNIFORM { return unif_pdf(x, a, b); }
    if fam == LAW_LOGISTIC { return logi_pdf(x, a, b); }
    if fam == LAW_CAUCHY { return cauchy_pdf(x, a, b); }
    if fam == LAW_LAPLACE { return laplace_pdf(x, a, b); }
    if fam == LAW_RAYLEIGH { return rayleigh_pdf(x, a, b); }
    if fam == LAW_GUMBEL_R { return gumbel_r_pdf(x, a, b); }
    if fam == LAW_GUMBEL_L { return gumbel_l_pdf(x, a, b); }
    return nan(), ERR_DOMAIN;
}

fn law_cdf(fam: i64, a: f64, b: f64, x: f64) -> f64, i64 {
    if fam == LAW_NORMAL { return norm_cdf(x, a, b); }
    if fam == LAW_UNIFORM { return unif_cdf(x, a, b); }
    if fam == LAW_LOGISTIC { return logi_cdf(x, a, b); }
    if fam == LAW_CAUCHY { return cauchy_cdf(x, a, b); }
    if fam == LAW_LAPLACE { return laplace_cdf(x, a, b); }
    if fam == LAW_RAYLEIGH { return rayleigh_cdf(x, a, b); }
    if fam == LAW_GUMBEL_R { return gumbel_r_cdf(x, a, b); }
    if fam == LAW_GUMBEL_L { return gumbel_l_cdf(x, a, b); }
    return nan(), ERR_DOMAIN;
}

fn law_sf(fam: i64, a: f64, b: f64, x: f64) -> f64, i64 {
    if fam == LAW_NORMAL { return norm_sf(x, a, b); }
    if fam == LAW_UNIFORM { return unif_sf(x, a, b); }
    if fam == LAW_LOGISTIC { return logi_sf(x, a, b); }
    if fam == LAW_CAUCHY { return cauchy_sf(x, a, b); }
    if fam == LAW_LAPLACE { return laplace_sf(x, a, b); }
    if fam == LAW_RAYLEIGH { return rayleigh_sf(x, a, b); }
    if fam == LAW_GUMBEL_R { return gumbel_r_sf(x, a, b); }
    if fam == LAW_GUMBEL_L { return gumbel_l_sf(x, a, b); }
    return nan(), ERR_DOMAIN;
}

fn law_ppf(fam: i64, a: f64, b: f64, p: f64) -> f64, i64 {
    if fam == LAW_NORMAL { return norm_ppf(p, a, b); }
    if fam == LAW_UNIFORM { return unif_ppf(p, a, b); }
    if fam == LAW_LOGISTIC { return logi_ppf(p, a, b); }
    if fam == LAW_CAUCHY { return cauchy_ppf(p, a, b); }
    if fam == LAW_LAPLACE { return laplace_ppf(p, a, b); }
    if fam == LAW_RAYLEIGH { return rayleigh_ppf(p, a, b); }
    if fam == LAW_GUMBEL_R { return gumbel_r_ppf(p, a, b); }
    if fam == LAW_GUMBEL_L { return gumbel_l_ppf(p, a, b); }
    return nan(), ERR_DOMAIN;
}

// The end of the support of a family that is not infinite.
fn law_floor(fam: i64, a: f64) -> f64 {
    if fam == LAW_UNIFORM || fam == LAW_RAYLEIGH { return a; }
    return -inf();
}

fn law_ceil(fam: i64, b: f64) -> f64 {
    if fam == LAW_UNIFORM { return b; }
    return inf();
}

// The code of the family and of its two parameters.
fn law_check(fam: i64, a: f64, b: f64) -> i64 {
    let _v, c = law_cdf(fam, a, b, a);
    return c;
}

// ---- the integral of a moment --------------------------------------

// The density of the law that the moment reads. Mode 0 is the truncated
// law with the mass p1. Mode 1 is the order statistic of rank p1 in p2
// draws, and p3 is the log of the beta function of its law.
fn law_dens(mode: i64, fam: i64, a: f64, b: f64, x: f64, p1: f64, p2: f64, p3: f64) -> f64 {
    let f, cf = law_pdf(fam, a, b, x);
    if cf != OK || f == 0.0 { return 0.0; }
    if mode == 0 { return f / p1; }
    let lo, c1 = law_cdf(fam, a, b, x);
    let hi, c2 = law_sf(fam, a, b, x);
    var t: f64 = 0.0;
    if p1 > 1.0 {
        if lo <= 0.0 { return 0.0; }
        t += (p1 - 1.0) * log(lo);
    }
    if p2 > p1 {
        if hi <= 0.0 { return 0.0; }
        t += (p2 - p1) * log(hi);
    }
    return f * exp(t - p3);
}

// The integral of (x - center)^order times the density, from lo to hi.
// The location splits the range, so that a kink of the density sits on
// the edge of a panel.
fn law_integral(gx: const f64[:], gw: const f64[:], mode: i64, fam: i64, a: f64, b: f64, lo: f64, hi: f64, p1: f64, p2: f64, p3: f64, center: f64, order: i64) -> f64 {
    var cuts: f64[>..] = [lo];
    if a > lo && a < hi { cuts.push(a); }
    cuts.push(hi);
    var total: f64 = 0.0;
    for s in 0..cuts.len - 1 {
        let left = cuts[s];
        let width = (cuts[s + 1] - left) / LAW_PANELS;
        for k in 0..LAW_PANELS {
            let lo_k = left + width * k;
            for j in 0..LAW_NODES {
                let x = lo_k + width * (gx[j] + 1.0) / 2.0;
                let d = x - center;
                var term = law_dens(mode, fam, a, b, x, p1, p2, p3);
                if order == 1 { term *= x; }
                else if order == 2 { term *= d * d; }
                total += term * gw[j] * width / 2.0;
            }
        }
    }
    return total;
}

// The finite range of an integral: the bounds, cut at the tail of the
// family where a bound is infinite. A range that is empty gives false.
fn law_range(fam: i64, a: f64, b: f64, lb: f64, ub: f64, lo: f64&, hi: f64&) -> bool {
    lo = lb;
    hi = ub;
    if is_inf(lb) {
        let cut, c = law_ppf(fam, a, b, LAW_TAIL);
        if c != OK { return false; }
        lo = cut;
    }
    if is_inf(ub) {
        let cut, c = law_ppf(fam, a, b, LAW_TOP);
        if c != OK { return false; }
        hi = cut;
    }
    let fl = law_floor(fam, a);
    let ce = law_ceil(fam, b);
    if lo < fl { lo = fl; }
    if hi > ce { hi = ce; }
    return lo < hi;
}

// The mean and the variance of a law from its integral.
fn law_moments(mode: i64, fam: i64, a: f64, b: f64, lo: f64, hi: f64, p1: f64, p2: f64, p3: f64, want_var: bool) -> f64, i64 {
    var gx: f64[>..] = [];
    var gw: f64[>..] = [];
    gauss_legendre(gx, gw, LAW_NODES);
    let m = law_integral(gx[..], gw[..], mode, fam, a, b, lo, hi, p1, p2, p3, 0.0, 1);
    var r = m;
    if want_var { r = law_integral(gx[..], gw[..], mode, fam, a, b, lo, hi, p1, p2, p3, m, 2); }
    if is_nan(r) || is_inf(r) { return nan(), ERR_OVERFLOW; }
    return r, OK;
}

// ---- the truncated law ---------------------------------------------

// The cdf and the survival function of the family at the two bounds, and
// the mass between them. An infinite bound has the cdf 0 or 1.
fn trunc_prep(fam: i64, a: f64, b: f64, lb: f64, ub: f64, fa: f64&, sfa: f64&, mass: f64&) -> i64 {
    if is_nan(lb) || is_nan(ub) { return ERR_NAN; }
    let c = law_check(fam, a, b);
    if c != OK { return c; }
    if lb >= ub { return ERR_DOMAIN; }
    fa = 0.0;
    sfa = 1.0;
    if !is_inf(lb) {
        let f, c1 = law_cdf(fam, a, b, lb);
        let s, c2 = law_sf(fam, a, b, lb);
        if c1 != OK { return c1; }
        if c2 != OK { return c2; }
        fa = f;
        sfa = s;
    }
    var fb: f64 = 1.0;
    var sfb: f64 = 0.0;
    if !is_inf(ub) {
        let f, c1 = law_cdf(fam, a, b, ub);
        let s, c2 = law_sf(fam, a, b, ub);
        if c1 != OK { return c1; }
        if c2 != OK { return c2; }
        fb = f;
        sfb = s;
    }
    if fa > 0.5 { mass = sfa - sfb; } else { mass = fb - fa; }
    if !(mass > 0.0) { return ERR_UNDEFINED; }
    return OK;
}

fn truncate_pdf(x: f64, fam: i64, a: f64, b: f64, lb: f64 = -inf(), ub: f64 = inf()) -> f64, i64 {
    var fa: f64 = 0.0;
    var sfa: f64 = 0.0;
    var mass: f64 = 0.0;
    let c = trunc_prep(fam, a, b, lb, ub, fa, sfa, mass);
    if c != OK { return nan(), c; }
    if is_nan(x) { return nan(), ERR_NAN; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    if x <= lb || x >= ub { return 0.0, OK; }
    let f, cf = law_pdf(fam, a, b, x);
    if cf != OK { return nan(), cf; }
    return loc_result(f / mass);
}

fn truncate_cdf(x: f64, fam: i64, a: f64, b: f64, lb: f64 = -inf(), ub: f64 = inf()) -> f64, i64 {
    var fa: f64 = 0.0;
    var sfa: f64 = 0.0;
    var mass: f64 = 0.0;
    let c = trunc_prep(fam, a, b, lb, ub, fa, sfa, mass);
    if c != OK { return nan(), c; }
    if is_nan(x) { return nan(), ERR_NAN; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    if x <= lb { return 0.0, OK; }
    if x >= ub { return 1.0, OK; }
    var v: f64 = 0.0;
    if fa > 0.5 {
        let s, cs = law_sf(fam, a, b, x);
        if cs != OK { return nan(), cs; }
        v = (sfa - s) / mass;
    } else {
        let f, cf = law_cdf(fam, a, b, x);
        if cf != OK { return nan(), cf; }
        v = (f - fa) / mass;
    }
    return min(1.0, max(0.0, v)), OK;
}

fn truncate_ppf(p: f64, fam: i64, a: f64, b: f64, lb: f64 = -inf(), ub: f64 = inf()) -> f64, i64 {
    var fa: f64 = 0.0;
    var sfa: f64 = 0.0;
    var mass: f64 = 0.0;
    let c = trunc_prep(fam, a, b, lb, ub, fa, sfa, mass);
    if c != OK { return nan(), c; }
    if is_nan(p) { return nan(), ERR_NAN; }
    if is_inf(p) { return nan(), ERR_INFINITY; }
    if p < 0.0 || p > 1.0 { return nan(), ERR_DOMAIN; }
    if p == 0.0 { if is_inf(lb) { return nan(), ERR_DOMAIN; } return lb, OK; }
    if p == 1.0 { if is_inf(ub) { return nan(), ERR_DOMAIN; } return ub, OK; }
    var target = fa + p * mass;
    if fa > 0.5 { target = 1.0 - (sfa - p * mass); }
    let v, cv = law_ppf(fam, a, b, min(target, LAW_TOP));
    if cv != OK { return nan(), cv; }
    return min(ub, max(lb, v)), OK;
}

// The mean and the variance of the truncated law.
fn trunc_moment(fam: i64, a: f64, b: f64, lb: f64, ub: f64, want_var: bool) -> f64, i64 {
    var fa: f64 = 0.0;
    var sfa: f64 = 0.0;
    var mass: f64 = 0.0;
    let c = trunc_prep(fam, a, b, lb, ub, fa, sfa, mass);
    if c != OK { return nan(), c; }
    if fam == LAW_CAUCHY && (is_inf(lb) || is_inf(ub)) { return nan(), ERR_UNDEFINED; }
    var lo: f64 = 0.0;
    var hi: f64 = 0.0;
    if !law_range(fam, a, b, lb, ub, lo, hi) { return nan(), ERR_UNDEFINED; }
    return law_moments(0, fam, a, b, lo, hi, mass, 0.0, 0.0, want_var);
}

fn truncate_mean(fam: i64, a: f64, b: f64, lb: f64 = -inf(), ub: f64 = inf()) -> f64, i64 {
    return trunc_moment(fam, a, b, lb, ub, false);
}

fn truncate_var(fam: i64, a: f64, b: f64, lb: f64 = -inf(), ub: f64 = inf()) -> f64, i64 {
    return trunc_moment(fam, a, b, lb, ub, true);
}

// ---- the order statistic -------------------------------------------

// The code of a family, a rank and a size.
fn order_check(fam: i64, a: f64, b: f64, r: i64, n: i64) -> i64 {
    let c = law_check(fam, a, b);
    if c != OK { return c; }
    if n < 1 || r < 1 || r > n { return ERR_DOMAIN; }
    if n > MAX_ORDER_N { return ERR_TOO_LONG; }
    return OK;
}

fn order_statistic_pdf(x: f64, fam: i64, a: f64, b: f64, r: i64, n: i64) -> f64, i64 {
    let c = order_check(fam, a, b, r, n);
    if c != OK { return nan(), c; }
    if is_nan(x) { return nan(), ERR_NAN; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    let lb, cb = lbeta(r, n - r + 1);
    if cb != OK { return nan(), cb; }
    return loc_result(law_dens(1, fam, a, b, x, r, n, lb));
}

fn order_statistic_cdf(x: f64, fam: i64, a: f64, b: f64, r: i64, n: i64) -> f64, i64 {
    let c = order_check(fam, a, b, r, n);
    if c != OK { return nan(), c; }
    let f, cf = law_cdf(fam, a, b, x);
    if cf != OK { return nan(), cf; }
    return betainc(r, n - r + 1, f);
}

fn order_statistic_ppf(p: f64, fam: i64, a: f64, b: f64, r: i64, n: i64) -> f64, i64 {
    let c = order_check(fam, a, b, r, n);
    if c != OK { return nan(), c; }
    let q, cq = betaincinv(r, n - r + 1, p);
    if cq != OK { return nan(), cq; }
    return law_ppf(fam, a, b, q);
}

fn order_moment(fam: i64, a: f64, b: f64, r: i64, n: i64, want_var: bool) -> f64, i64 {
    let c = order_check(fam, a, b, r, n);
    if c != OK { return nan(), c; }
    if fam == LAW_CAUCHY { return nan(), ERR_UNDEFINED; }
    let lb, cb = lbeta(r, n - r + 1);
    if cb != OK { return nan(), cb; }
    var lo: f64 = 0.0;
    var hi: f64 = 0.0;
    if !law_range(fam, a, b, -inf(), inf(), lo, hi) { return nan(), ERR_UNDEFINED; }
    return law_moments(1, fam, a, b, lo, hi, r, n, lb, want_var);
}

fn order_statistic_mean(fam: i64, a: f64, b: f64, r: i64, n: i64) -> f64, i64 {
    return order_moment(fam, a, b, r, n, false);
}

fn order_statistic_var(fam: i64, a: f64, b: f64, r: i64, n: i64) -> f64, i64 {
    return order_moment(fam, a, b, r, n, true);
}

// ---- the Dirichlet multinomial law ---------------------------------

// The code of the counts, the concentrations and the number of draws.
// The check order is: bound, empty, length, NaN, infinity, domain.
fn dirmult_check(x: const f64[:], alpha: const f64[:], n: i64, need_counts: bool) -> i64 {
    if alpha.len > MAX_DIRMULT_K { return ERR_TOO_LONG; }
    if alpha.len == 0 { return ERR_EMPTY; }
    if need_counts && x.len != alpha.len { return ERR_LENGTH; }
    var total: f64 = 0.0;
    for v in alpha { if is_nan(v) { return ERR_NAN; } }
    for v in alpha { if is_inf(v) { return ERR_INFINITY; } }
    if need_counts {
        for v in x { if is_nan(v) { return ERR_NAN; } }
        for v in x { if is_inf(v) { return ERR_INFINITY; } }
    }
    for v in alpha { if v <= 0.0 { return ERR_DOMAIN; } }
    if n < 0 { return ERR_DOMAIN; }
    if need_counts {
        for v in x {
            if v < 0.0 || v != floor(v) { return ERR_DOMAIN; }
            total += v;
        }
        if total != n { return ERR_DOMAIN; }
    }
    return OK;
}

fn dirichlet_multinomial_logpmf(x: const f64[:], alpha: const f64[:], n: i64) -> f64, i64 {
    let c = dirmult_check(x, alpha, n, true);
    if c != OK { return nan(), c; }
    if n == 0 { return 0.0, OK; }
    var a0: f64 = 0.0;
    for v in alpha { a0 += v; }
    let head, ch = lbeta(a0, n);
    if ch != OK { return nan(), ch; }
    var s = log(n as f64) + head;
    for i in 0..x.len {
        if x[i] > 0.0 {
            let t, ct = lbeta(alpha[i], x[i]);
            if ct != OK { return nan(), ct; }
            s -= log(x[i]) + t;
        }
    }
    return loc_result(s);
}

fn dirichlet_multinomial_pmf(x: const f64[:], alpha: const f64[:], n: i64) -> f64, i64 {
    let l, c = dirichlet_multinomial_logpmf(x, alpha, n);
    if c != OK { return nan(), c; }
    return exp(l), OK;
}

// The mean of each category.
fn dirichlet_multinomial_mean(out: f64[>..]&, alpha: const f64[:], n: i64) -> i64 {
    out.resize(0);
    let none: f64[0] = [];
    let c = dirmult_check(none[..], alpha, n, false);
    if c != OK { return c; }
    var a0: f64 = 0.0;
    for v in alpha { a0 += v; }
    if is_inf(a0) { return ERR_OVERFLOW; }
    for v in alpha { out.push(n * v / a0); }
    return OK;
}

// The variance of each category.
fn dirichlet_multinomial_var(out: f64[>..]&, alpha: const f64[:], n: i64) -> i64 {
    out.resize(0);
    let none: f64[0] = [];
    let c = dirmult_check(none[..], alpha, n, false);
    if c != OK { return c; }
    var a0: f64 = 0.0;
    for v in alpha { a0 += v; }
    let k = (n + a0) / (a0 * a0 * (a0 + 1.0));
    if is_inf(a0) || is_nan(k) || k == 0.0 { return ERR_OVERFLOW; }
    for v in alpha { out.push(n * v * (a0 - v) * k); }
    return OK;
}

// The covariance matrix, row by row. The count of categories is at most
// MAX_DIRMULT_K, so the matrix holds at most MAX_DIRMULT_K squared values.
fn dirichlet_multinomial_cov(out: f64[>..]&, alpha: const f64[:], n: i64) -> i64 {
    out.resize(0);
    let none: f64[0] = [];
    let c = dirmult_check(none[..], alpha, n, false);
    if c != OK { return c; }
    var a0: f64 = 0.0;
    for v in alpha { a0 += v; }
    let k = (n + a0) / (a0 * a0 * (a0 + 1.0));
    if is_inf(a0) || is_nan(k) || k == 0.0 { return ERR_OVERFLOW; }
    for i in 0..alpha.len {
        for j in 0..alpha.len {
            if i == j { out.push(n * alpha[i] * (a0 - alpha[i]) * k); }
            else { out.push(-n * alpha[i] * alpha[j] * k); }
        }
    }
    return OK;
}
