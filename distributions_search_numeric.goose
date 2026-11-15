// The shared parts of the families whose quantile needs a search.
//
// The search of a quantile is a bracketed root. The bracket grows until
// the cumulative function of the law passes the probability, and the
// root comes from the Brent method of the package. The search stops
// after MAX_FAMILY_SEARCH steps and gives the fit fault.
namespace statistics;

import std;
import math;

const MAX_FAMILY_SEARCH: i64 = 200;   // steps of the search of a quantile

// The quantile of a law whose cumulative function is F. The lower bound
// is the lower bound of the support and the upper bound grows until the
// cumulative function passes the probability. The block gives the
// cumulative function of a point.
fn fam_search_ppf<F>(p: f64, lo: f64, hi0: f64) -> f64, i64 {
    var hi = hi0;
    var steps: i64 = 0;
    var found = false;
    while steps < MAX_FAMILY_SEARCH {
        let c = F(hi);
        if is_nan(c) { return nan(), ERR_NAN; }
        if c >= p { found = true; break; }
        hi = hi * 2.0 + 1.0;
        steps += 1;
    }
    if !found { return nan(), ERR_FIT; }
    // The root sits above the lower bound and under the upper bound. The
    // lower bound may hold a cumulative function above the probability
    // when the law has no bound there. The search grows the bracket
    // downward until the root lies inside it.
    let flo0 = F(lo);
    if is_nan(flo0) { return nan(), ERR_NAN; }
    if flo0 >= p {
        var l = lo;
        var fl = flo0;
        var steps: i64 = 0;
        while steps < MAX_FAMILY_SEARCH && fl >= p {
            l = l * 2.0 - 1.0;
            let f2 = F(l);
            if is_nan(f2) { return nan(), ERR_NAN; }
            fl = f2;
            steps += 1;
        }
        if fl >= p { return nan(), ERR_FIT; }
        let x, cx = brent_root(l, hi, 1e-14, 1e-15) { v => F(v) - p };
        if cx != OK { return nan(), cx; }
        return x, OK;
    }
    let flo = flo0;
    let x, c = brent_root(lo, hi, 1e-14, 1e-15) { v => F(v) - p };
    if c != OK { return nan(), c; }
    return x, OK;
}

// The quantile of a law whose survival function is F. The upper bound
// grows until the survival function falls under the probability.
fn fam_search_ppf_sf<F>(p: f64, lo: f64, hi0: f64) -> f64, i64 {
    var hi = hi0;
    var steps: i64 = 0;
    var found = false;
    while steps < MAX_FAMILY_SEARCH {
        let c = F(hi);
        if is_nan(c) { return nan(), ERR_NAN; }
        if c <= p { found = true; break; }
        hi = hi * 2.0 + 1.0;
        steps += 1;
    }
    if !found { return nan(), ERR_FIT; }
    let x, c = brent_root(lo, hi, 1e-14, 1e-15) { v => p - F(v) };
    if c != OK { return nan(), c; }
    return x, OK;
}

// The check of a shape that must be a whole number. The Erlang law and
// the Irwin-Hall law want a whole number of terms.
fn check_whole_shape(a: f64) -> i64 {
    let c = check_shape(a);
    if c != OK { return c; }
    if a != floor(a) { return ERR_DOMAIN; }
    return OK;
}

// The check of a shape that must be 0 or more.
fn check_shape_nonneg(a: f64) -> i64 {
    if is_nan(a) { return ERR_NAN; }
    if is_inf(a) { return ERR_INFINITY; }
    if a < 0.0 { return ERR_DOMAIN; }
    return OK;
}

// ---- the quadrature of the families ---------------------------------------

const FAM_QUAD_NODES: i64 = 16;       // nodes of one panel
const FAM_QUAD_START: i64 = 8;        // panels of the first estimate
const FAM_QUAD_DOUBLINGS: i64 = 9;    // the panels double at most this often
const FAM_QUAD_TOL: f64 = 1e-14;      // relative change that ends the doubling
const FAM_GRADE_LEVELS: i64 = 40;     // levels of the graded panels toward an end

// The sum of one estimate: n equal panels of the rule on lo to hi.
fn fam_quad_panels<F>(lo: f64, hi: f64, n: i64, gx: const f64[:], gw: const f64[:]) -> f64 {
    let width = (hi - lo) / n as f64;
    var total: f64 = 0.0;
    for k in 0..n {
        let base = lo + width * k as f64;
        for j in 0..gx.len {
            let v = F(base + width * (gx[j] + 1.0) / 2.0);
            if !is_nan(v) && !is_inf(v) { total += v * gw[j] * width / 2.0; }
        }
    }
    return total;
}

// The integral of F on lo to hi. The panels double until two estimates
// agree to FAM_QUAD_TOL, or FAM_QUAD_DOUBLINGS doublings pass.
fn fam_quad<F>(lo: f64, hi: f64) -> f64 {
    var gx: f64[>..] = [];
    var gw: f64[>..] = [];
    gauss_legendre(gx, gw, FAM_QUAD_NODES);
    var n = FAM_QUAD_START;
    var prev = fam_quad_panels(lo, hi, n, gx[..], gw[..]) { x => F(x) };
    for d in 0..FAM_QUAD_DOUBLINGS {
        n = n * 2;
        let cur = fam_quad_panels(lo, hi, n, gx[..], gw[..]) { x => F(x) };
        // A zero estimate may hide a narrow peak, so it never ends the loop.
        if cur != 0.0 && abs(cur - prev) <= FAM_QUAD_TOL * abs(cur) { return cur; }
        prev = cur;
    }
    return prev;
}

// The integral of F from lo to infinity, by the map x = lo + t / (1 - t).
fn fam_quad_up<F>(lo: f64) -> f64 {
    return fam_quad(0.0, 1.0) { t => F(lo + t / (1.0 - t)) / ((1.0 - t) * (1.0 - t)) };
}

// The integral of F from minus infinity to hi.
fn fam_quad_down<F>(hi: f64) -> f64 {
    return fam_quad(0.0, 1.0) { t => F(hi - t / (1.0 - t)) / ((1.0 - t) * (1.0 - t)) };
}

// The integral of F on lo to hi with the panels graded toward both ends.
// A density with a power of the distance at an end gives a fast sum. The
// left panels cover the distance w/2^(k+1) to w/2^k from lo, and the
// right panels do the same from hi. The panels meet at the middle.
fn fam_quad_graded<F>(lo: f64, hi: f64) -> f64 {
    var gx: f64[>..] = [];
    var gw: f64[>..] = [];
    gauss_legendre(gx, gw, FAM_QUAD_NODES);
    let w = hi - lo;
    var total: f64 = 0.0;
    var d: f64 = 0.25;
    for lev in 0..FAM_GRADE_LEVELS {
        total += fam_quad_panels(lo + w * d, lo + w * 2.0 * d, 1, gx[..], gw[..]) { x => F(x) };
        total += fam_quad_panels(hi - w * 2.0 * d, hi - w * d, 1, gx[..], gw[..]) { x => F(x) };
        d = d * 0.5;
    }
    return total;
}

// The scaled modified Bessel function of the second kind, exp(x) K_v(x),
// for x above 0. The nodes come from the caller, so a loop builds them
// once.
fn fam_kve(v: f64, x: f64, bx: const f64[:], bw: const f64[:]) -> f64 {
    return ex_bessel_ke(v, x, bx, bw);
}

// The nodes of the rule of the Bessel function.
fn fam_bessel_nodes(bx: f64[>..]&, bw: f64[>..]&) {
    gauss_legendre(bx, bw, EX_BESSEL_NODES);
}
