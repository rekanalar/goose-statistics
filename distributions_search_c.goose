// Five continuous families of the search kind: norminvgauss,
// recipinvgauss, rel_breitwigner, studentized_range and vonmises. The
// names and the parameters follow SciPy 1.18.1.
//
// Every family has three calls: pdf, cdf and ppf. The quantile comes
// from a bracketed root of the cumulative function. Every call returns
// the result first and a fault code last.
//
// The law of vonmises is a law of an angle, so its cumulative function
// grows by one at every turn of pi. SciPy 1.18.1 gives the same values.
namespace statistics;

import std;
import math;

// ---- norminvgauss: the normal inverse Gaussian law ------------------------

fn norminvgauss_check(a: f64, b: f64) -> i64 {
    let c1 = check_shape(a);
    if c1 != OK { return c1; }
    let c2 = check_scalar(b);
    if c2 != OK { return c2; }
    if abs(b) >= a { return ERR_DOMAIN; }
    return OK;
}

// The log of the density of the normal inverse Gaussian law. The value of
// g is the square root of (a + b) (a - b).
fn norminvgauss_logdens(x: f64, a: f64, b: f64, g: f64) -> f64 {
    let s = sqrt(1.0 + x * x);
    var bx: f64[>..] = [];
    var bw: f64[>..] = [];
    fam_bessel_nodes(bx, bw);
    // The scaled Bessel function of the first kind, K1(a s) exp(a s).
    let lk = log(fam_kve(1.0, a * s, bx[..], bw[..]));
    if is_nan(lk) { return nan(); }
    return log(a / PI) + lk - a * s + b * x + g - log(s);
}

// The density of the law, without the checks.
fn norminvgauss_dens(x: f64, a: f64, b: f64, g: f64) -> f64 {
    let l = norminvgauss_logdens(x, a, b, g);
    if is_nan(l) { return nan(); }
    if l < -700.0 { return 0.0; }
    return exp(l);
}

fn norminvgauss_pdf(x: f64, a: f64, b: f64) -> f64, i64 {
    let k = norminvgauss_check(a, b);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    let g = sqrt((a + b) * (a - b));
    let v = norminvgauss_dens(x, a, b, g);
    if is_nan(v) { return nan(), ERR_FIT; }
    return v, OK;
}

// The cumulative function of the law, without the code. The integral runs
// from minus infinity to the point.
fn norminvgauss_cdf_value(x: f64, a: f64, b: f64, g: f64) -> f64 {
    let mean = b / g;
    var v: f64 = 0.0;
    if x <= mean {
        v = fam_quad_down(x) { t => norminvgauss_dens(t, a, b, g) };
    } else {
        v = 1.0 - fam_quad_up(x) { t => norminvgauss_dens(t, a, b, g) };
    }
    if v < 0.0 { return 0.0; }
    if v > 1.0 { return 1.0; }
    return v;
}

fn norminvgauss_cdf(x: f64, a: f64, b: f64) -> f64, i64 {
    let k = norminvgauss_check(a, b);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    let g = sqrt((a + b) * (a - b));
    let v = norminvgauss_cdf_value(x, a, b, g);
    if is_nan(v) { return nan(), ERR_FIT; }
    return v, OK;
}

fn norminvgauss_ppf(q: f64, a: f64, b: f64) -> f64, i64 {
    let k = norminvgauss_check(a, b);
    if k != OK { return nan(), k; }
    let pc = check_prob(q);
    if pc != OK { return nan(), pc; }
    if q == 0.0 { return -inf(), OK; }
    if q == 1.0 { return inf(), OK; }
    let g = sqrt((a + b) * (a - b));
    let mean = b / g;
    // The bracket grows on both sides of the mean.
    var lo = mean - 1.0;
    var hi = mean + 1.0;
    var steps: i64 = 0;
    while steps < MAX_FAMILY_SEARCH && norminvgauss_cdf_value(lo, a, b, g) > q {
        lo = mean - 2.0 * (mean - lo);
        steps += 1;
    }
    while steps < MAX_FAMILY_SEARCH && norminvgauss_cdf_value(hi, a, b, g) < q {
        hi = mean + 2.0 * (hi - mean);
        steps += 1;
    }
    if steps >= MAX_FAMILY_SEARCH { return nan(), ERR_FIT; }
    let x, cx = brent_root(lo, hi, 1e-13, 4e-16) { v => norminvgauss_cdf_value(v, a, b, g) - q };
    return x, cx;
}

// ---- recipinvgauss: the reciprocal of the inverse Gaussian law -------------

fn recipinvgauss_check(mu: f64) -> i64 {
    return check_shape(mu);
}

fn recipinvgauss_pdf(x: f64, mu: f64) -> f64, i64 {
    let k = recipinvgauss_check(mu);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    let d = 1.0 - mu * x;
    let l = -d * d / (2.0 * x * mu * mu) - 0.5 * log(2.0 * PI * x);
    return exp(l), OK;
}

fn recipinvgauss_cdf(x: f64, mu: f64) -> f64, i64 {
    let k = recipinvgauss_check(mu);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    let t1 = 1.0 / mu - x;
    let t2 = 1.0 / mu + x;
    let isq = 1.0 / sqrt(x);
    let n1, c1 = norm_cdf(-isq * t1);
    if c1 != OK { return nan(), c1; }
    let n2, c2 = norm_cdf(-isq * t2);
    if c2 != OK { return nan(), c2; }
    let e = exp(2.0 / mu);
    if is_inf(e) { return nan(), ERR_OVERFLOW; }
    let v = n1 - e * n2;
    if v < 0.0 { return 0.0, OK; }
    if v > 1.0 { return 1.0, OK; }
    return v, OK;
}

fn recipinvgauss_cdf_value(x: f64, mu: f64) -> f64 {
    let v, c = recipinvgauss_cdf(x, mu);
    if c != OK { return nan(); }
    return v;
}

fn recipinvgauss_ppf(q: f64, mu: f64) -> f64, i64 {
    let k = recipinvgauss_check(mu);
    if k != OK { return nan(), k; }
    let pc = check_prob(q);
    if pc != OK { return nan(), pc; }
    if q == 0.0 { return 0.0, OK; }
    return fam_search_ppf(q, 0.0, 1.0) { v => recipinvgauss_cdf_value(v, mu) };
}

// ---- rel_breitwigner: the relative Breit-Wigner law ------------------------

fn rel_breitwigner_check(rho: f64) -> i64 {
    return check_shape(rho);
}

// The constant of the density.
fn rel_breitwigner_c(rho: f64) -> f64 {
    let r2 = 1.0 + 1.0 / (rho * rho);
    return sqrt(2.0 * r2 / (1.0 + sqrt(r2))) * 2.0 / PI;
}

fn rel_breitwigner_pdf(x: f64, rho: f64) -> f64, i64 {
    let k = rel_breitwigner_check(rho);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    let q = ((x - rho) * (x + rho) / rho);
    return rel_breitwigner_c(rho) / (q * q + 1.0), OK;
}

// The cumulative function of the law. The value comes from the tangent
// of a complex angle. The law writes that angle with a square root and
// an arctangent of a complex number, so the call splits the complex
// numbers into their real and their imaginary parts.
fn rel_breitwigner_cdf(x: f64, rho: f64) -> f64, i64 {
    let k = rel_breitwigner_check(rho);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    // The constant of the cumulative function.
    let c = sqrt(2.0 / (1.0 + sqrt(1.0 + 1.0 / (rho * rho)))) / PI;
    // The square root of -rho * (rho + i). Its real part and its
    // imaginary part follow from the size of the number.
    let r = rho * sqrt(rho * rho + 1.0);
    let wr = sqrt((r - rho * rho) / 2.0);
    let wi = -sqrt((r + rho * rho) / 2.0);
    // The point divided by that root.
    let pr = x * wr / r;
    let pi = -x * wi / r;
    // The arctangent of a complex number. Its real part and its
    // imaginary part follow from the real and the imaginary parts.
    let ar = 0.5 * atan2(2.0 * pr, 1.0 - pr * pr - pi * pi);
    let ai = 0.25 * log((pr * pr + (pi + 1.0) * (pi + 1.0)) / (pr * pr + (pi - 1.0) * (pi - 1.0)));
    // The square root of -1 + i / rho.
    let q = sqrt(1.0 + 1.0 / (rho * rho));
    let sr = sqrt((q - 1.0) / 2.0);
    let si = sqrt((q + 1.0) / 2.0);
    // The imaginary part of the product of the two numbers.
    let im = sr * ai + si * ar;
    let v = c * 2.0 * im;
    if is_nan(v) { return nan(), ERR_FIT; }
    if v < 0.0 { return 0.0, OK; }
    if v > 1.0 { return 1.0, OK; }
    return v, OK;
}

fn rel_breitwigner_cdf_value(x: f64, rho: f64) -> f64 {
    let v, c = rel_breitwigner_cdf(x, rho);
    if c != OK { return nan(); }
    return v;
}

fn rel_breitwigner_ppf(q: f64, rho: f64) -> f64, i64 {
    let k = rel_breitwigner_check(rho);
    if k != OK { return nan(), k; }
    let pc = check_prob(q);
    if pc != OK { return nan(), pc; }
    if q == 0.0 { return 0.0, OK; }
    return fam_search_ppf(q, 0.0, rho) { v => rel_breitwigner_cdf_value(v, rho) };
}

// ---- studentized_range: the studentized range law -------------------------

const SR_Z_NODES: i64 = 10;          // nodes of one panel of the inner rule
const SR_S_NODES: i64 = 10;          // nodes of one panel of the outer rule
const SR_Z_REACH: f64 = 12.0;        // the reach of the inner rule
const SR_S_REACH: f64 = 12.0;        // the reach of the outer rule

fn studentized_range_check(k: f64, df: f64) -> i64 {
    let c1 = check_shape(k);
    if c1 != OK { return c1; }
    let c2 = check_shape(df);
    if c2 != OK { return c2; }
    // The law needs two samples at least.
    if k < 2.0 { return ERR_DOMAIN; }
    return OK;
}

// The nodes and the weights of the inner rule. The rule runs over the
// whole line, and its panels grow away from the middle. The panels hold
// the weight of the normal law, so the rule integrates the product of a
// smooth function and the normal density.
fn sr_z_rule(zs: f64[>..]&, zw: f64[>..]&) {
    zs.resize(0);
    zw.resize(0);
    var gx: f64[>..] = [];
    var gw: f64[>..] = [];
    gauss_legendre(gx, gw, SR_Z_NODES);
    // The edges of the panels on the positive side. The width doubles at
    // every step.
    var edges: f64[>..] = [];
    edges.push(0.0);
    var d = 0.125;
    while edges[edges.len - 1] < SR_Z_REACH {
        let last = edges[edges.len - 1] + d;
        edges.push(min(last, SR_Z_REACH));
        d = d * 2.0;
    }
    // The panels on the negative side, in the reverse order.
    var n = edges.len;
    for i in 0..n - 1 {
        let hi = -edges[i];
        let lo = -edges[i + 1];
        sr_z_panel(zs, zw, gx, gw, lo, hi);
    }
    for i in 0..n - 1 {
        let lo = edges[i];
        let hi = edges[i + 1];
        sr_z_panel(zs, zw, gx, gw, lo, hi);
    }
}

// The nodes and the weights of one panel of the inner rule.
fn sr_z_panel(zs: f64[>..]&, zw: f64[>..]&, gx: const f64[:], gw: const f64[:], lo: f64, hi: f64) {
    let w = hi - lo;
    for j in 0..gx.len {
        zs.push(lo + w * (gx[j] + 1.0) / 2.0);
        zw.push(gw[j] * w / 2.0);
    }
}

// The nodes and the weights of the outer rule. The rule runs from 0 to
// the reach, and its panels grow away from 0. The panels hold the weight
// of the divisor, so the rule integrates the product of a smooth
// function and a power of the point.
fn sr_s_rule(sp: f64[>..]&, sw: f64[>..]&) {
    sp.resize(0);
    sw.resize(0);
    var gx: f64[>..] = [];
    var gw: f64[>..] = [];
    gauss_legendre(gx, gw, SR_S_NODES);
    var edges: f64[>..] = [];
    edges.push(0.0);
    var d = 0.05;
    while edges[edges.len - 1] < SR_S_REACH {
        let last = edges[edges.len - 1] + d;
        edges.push(min(last, SR_S_REACH));
        d = min(d * 1.5, 1.0);
    }
    for i in 0..edges.len - 1 {
        let lo = edges[i];
        let hi = edges[i + 1];
        let w = hi - lo;
        for j in 0..gx.len {
            sp.push(lo + w * (gx[j] + 1.0) / 2.0);
            sw.push(gw[j] * w / 2.0);
        }
    }
}

// The logarithm of the constant of the cumulative function.
fn sr_log_cdf_const(k: f64, df: f64) -> f64 {
    let lg, sg, cg = lgamma(0.5 * df);
    if cg != OK { return nan(); }
    return log(k) + 0.5 * df * log(df) - lg - (0.5 * df - 1.0) * log(2.0);
}

// The cumulative function of the law, without the code. The law is the
// range of k normal values divided by the square root of a chi-square
// law. The integral runs over the divisor and over the point of the
// normal law.
fn studentized_range_cdf_value(x: f64, k: f64, df: f64) -> f64 {
    let lc = sr_log_cdf_const(k, df);
    if is_nan(lc) { return nan(); }
    var zs: f64[>..] = [];
    var zw: f64[>..] = [];
    sr_z_rule(zs, zw);
    var sp: f64[>..] = [];
    var sw: f64[>..] = [];
    sr_s_rule(sp, sw);
    var total: f64 = 0.0;
    for i in 0..sp.len {
        let s = sp[i];
        if s <= 0.0 { continue; }
        // The inner integral over the point of the normal law.
        var inner: f64 = 0.0;
        for j in 0..zs.len {
            let z = zs[j];
            let d = norm_cdf_diff(z + x * s, z);
            if d <= 0.0 { continue; }
            let ph, cph = norm_pdf(z);
            if cph != OK { continue; }
            inner += zw[j] * ph * pow(d, k - 1.0);
        }
        total += sw[i] * exp(lc + (df - 1.0) * log(s) - 0.5 * df * s * s) * inner;
    }
    if total < 0.0 { return 0.0; }
    if total > 1.0 { return 1.0; }
    return total;
}

// The difference of two cumulative functions of the normal law. The form
// keeps the digits of a small difference.
fn norm_cdf_diff(a: f64, b: f64) -> f64 {
    let na, ca = norm_cdf(a);
    if ca != OK { return nan(); }
    let nb, cb = norm_cdf(b);
    if cb != OK { return nan(); }
    return na - nb;
}

fn studentized_range_pdf(x: f64, k: f64, df: f64) -> f64, i64 {
    let c = studentized_range_check(k, df);
    if c != OK { return nan(), c; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    // The density is the derivative of the cumulative function. The
    // integral runs over the divisor and over the point of the normal
    // law.
    let lg, sg, cg = lgamma(0.5 * df);
    if cg != OK { return nan(), cg; }
    let lc = log(k) + log(k - 1.0) + 0.5 * df * log(df) - lg - (0.5 * df - 1.0) * log(2.0);
    var zs: f64[>..] = [];
    var zw: f64[>..] = [];
    sr_z_rule(zs, zw);
    var sp: f64[>..] = [];
    var sw: f64[>..] = [];
    sr_s_rule(sp, sw);
    var total: f64 = 0.0;
    for i in 0..sp.len {
        let s = sp[i];
        if s <= 0.0 { continue; }
        var inner: f64 = 0.0;
        for j in 0..zs.len {
            let z = zs[j];
            let d = norm_cdf_diff(z + x * s, z);
            if d <= 0.0 { continue; }
            let p1, c1 = norm_pdf(z);
            if c1 != OK { continue; }
            let p2, c2 = norm_pdf(z + x * s);
            if c2 != OK { continue; }
            inner += zw[j] * p1 * p2 * pow(d, k - 2.0);
        }
        total += sw[i] * exp(lc + df * log(s) - 0.5 * df * s * s) * inner;
    }
    if is_nan(total) { return nan(), ERR_FIT; }
    return total, OK;
}

fn studentized_range_cdf(x: f64, k: f64, df: f64) -> f64, i64 {
    let c = studentized_range_check(k, df);
    if c != OK { return nan(), c; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    let v = studentized_range_cdf_value(x, k, df);
    if is_nan(v) { return nan(), ERR_FIT; }
    return v, OK;
}


fn studentized_range_ppf(q: f64, k: f64, df: f64) -> f64, i64 {
    let c = studentized_range_check(k, df);
    if c != OK { return nan(), c; }
    let pc = check_prob(q);
    if pc != OK { return nan(), pc; }
    if q == 0.0 { return 0.0, OK; }
    if q == 1.0 { return inf(), OK; }
    return fam_search_ppf(q, 0.0, k + 4.0) { v => studentized_range_cdf_value(v, k, df) };
}

// ---- vonmises: the von Mises law -------------------------------------------

fn vonmises_check(kappa: f64) -> i64 {
    return check_shape_nonneg(kappa);
}

// The scaled modified Bessel function of the first kind, exp(-x) I0(x).
fn vonmises_i0e(x: f64) -> f64 {
    if x < 30.0 {
        let q = 0.25 * x * x;
        var term: f64 = 1.0;
        var total: f64 = 1.0;
        var k: f64 = 1.0;
        while k < 400.0 {
            term = term * q / (k * k);
            total += term;
            if term < 1e-17 * total { break; }
            k += 1.0;
        }
        return exp(-x) * total;
    }
    var term: f64 = 1.0;
    var total: f64 = 1.0;
    var k: f64 = 1.0;
    while k < 60.0 {
        let f = (2.0 * k - 1.0) / (8.0 * x * k);
        let next = term * f * (2.0 * k - 1.0);
        if abs(next) >= abs(term) { break; }
        term = next;
        total += term;
        if abs(term) < 1e-17 * abs(total) { break; }
        k += 1.0;
    }
    return total / sqrt(2.0 * PI * x);
}

fn vonmises_pdf(x: f64, kappa: f64) -> f64, i64 {
    let k = vonmises_check(kappa);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    let c = cos(x) - 1.0;
    return exp(kappa * c) / (2.0 * PI * vonmises_i0e(kappa)), OK;
}

// The cumulative function of the law. The law is a law of an angle, so
// the value grows by one at every turn of pi. The block gives the
// integral of the density of the law from zero to a point.
fn vonmises_cdf_value(x: f64, kappa: f64) -> f64 {
    // The number of full turns before the point.
    let turns = floor((x + PI) / (2.0 * PI));
    // The point inside the turn.
    let y = x - 2.0 * PI * turns;
    // The law is symmetric about zero, so half of the mass lies below
    // the zero angle.
    let v = fam_quad(0.0, y) { t => exp(kappa * (cos(t) - 1.0)) / (2.0 * PI * vonmises_i0e(kappa)) };
    return turns + 0.5 + v;
}

fn vonmises_cdf(x: f64, kappa: f64) -> f64, i64 {
    let k = vonmises_check(kappa);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    let v = vonmises_cdf_value(x, kappa);
    if is_nan(v) { return nan(), ERR_FIT; }
    return v, OK;
}

fn vonmises_ppf(q: f64, kappa: f64) -> f64, i64 {
    let k = vonmises_check(kappa);
    if k != OK { return nan(), k; }
    let pc = check_prob(q);
    if pc != OK { return nan(), pc; }
    if q == 0.0 { return -inf(), OK; }
    if q == 1.0 { return inf(), OK; }
    // The law repeats every turn, so the quantile of a probability sits
    // in the first turn.
    let f = q - floor(q);
    if f <= 0.0 { return inf(), OK; }
    // The law is symmetric about the zero angle, so half of the mass lies
    // below it. The answer of that point is exact.
    if f == 0.5 { return 0.0, OK; }
    let x, cx = brent_root(-PI, PI, 1e-14, 4e-16) { v => vonmises_cdf_value(v, kappa) - f };
    return x, cx;
}
