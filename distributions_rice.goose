// The Rice law. The name and the shape follow SciPy: rice(b), with the
// scale 1 and the support 0 to infinity.
//
// The density is x exp(-(x^2 + b^2) / 2) I0(x b), where I0 is the modified
// Bessel function of the first kind. The cumulative function is
// 1 - Q1(b, x), where Q1 is the Marcum function. The sum below writes the
// Marcum function as a Poisson mixture of regularized incomplete gamma
// functions, which is the form of the noncentral chi-square law with 2
// degrees of freedom.
//
// The shape must be 0 or more. A shape above RICE_MAX_B gives the domain
// fault, because the sum of the mixture needs a term count that grows
// with the shape.
namespace statistics;

import std;
import math;

const RICE_MAX_B: f64 = 3000.0;
const RICE_SERIES_LIMIT: f64 = 30.0;

// The check of the shape of the Rice law.
fn rice_check(b: f64) -> i64 {
    if is_nan(b) { return ERR_NAN; }
    if is_inf(b) { return ERR_INFINITY; }
    if b < 0.0 || b > RICE_MAX_B { return ERR_DOMAIN; }
    return OK;
}

// exp(-z) I0(z) for z of 0 or more. The power series serves a small z.
// All its terms are positive, so no digits cancel. The asymptotic series
// serves a large z.
fn rice_i0e(z: f64) -> f64 {
    if z < RICE_SERIES_LIMIT {
        let q = 0.25 * z * z;
        var term: f64 = 1.0;
        var total: f64 = 1.0;
        var k: f64 = 1.0;
        while k < 400.0 {
            term = term * q / (k * k);
            total += term;
            if term < 1e-17 * total { break; }
            k += 1.0;
        }
        return exp(-z) * total;
    }
    var term: f64 = 1.0;
    var total: f64 = 1.0;
    var k: f64 = 1.0;
    while k < 60.0 {
        let f = (2.0 * k - 1.0) / (8.0 * z * k);
        let next = term * f * (2.0 * k - 1.0);
        if abs(next) >= abs(term) { break; }
        term = next;
        total += term;
        if abs(term) < 1e-17 * abs(total) { break; }
        k += 1.0;
    }
    return total / sqrt(2.0 * PI * z);
}

// The density of the Rice law, without the checks.
fn rice_dens(x: f64, b: f64) -> f64 {
    if x <= 0.0 { return 0.0; }
    let d = x - b;
    return x * exp(-0.5 * d * d) * rice_i0e(x * b);
}

// The Marcum sum. The result is the cumulative function when upper is
// false and the survival function when upper is true. The sum runs from
// the mode of the Poisson weights outward, and it stops when the weight
// falls under 1e-18.
//
// The incomplete gamma function of the count above the mode follows from
// the function of the count below by one step, and the same holds below
// the mode. One call of the incomplete gamma function then serves the
// whole sum.
fn rice_sum(x: f64, b: f64, upper: bool) -> f64, i64 {
    let half = 0.5 * x * x;
    let lam = 0.5 * b * b;
    let j0 = floor(lam);
    let a0 = j0 + 1.0;
    var g: f64 = 0.0;
    var cc: i64 = OK;
    if upper {
        g, cc = gammaincc(a0, half);
    } else {
        g, cc = gammainc(a0, half);
    }
    if cc != OK { return nan(), cc; }
    let lg, sg, cg = lgamma(a0 + 1.0);
    if cg != OK { return nan(), cg; }
    // e holds half^a exp(-half) / gamma(a + 1), the step of the sum.
    let e0 = exp(a0 * log(half) - half - lg);
    let g0 = g;
    var e = e0;
    var a = a0;
    let llam = log(lam);
    let lcut = log(1e-18);
    var total: f64 = 0.0;
    // The count at the mode.
    let lg0, sg0, cg0 = lgamma(j0 + 1.0);
    if cg0 != OK { return nan(), cg0; }
    total += exp(-lam + j0 * llam - lg0) * g;
    // The counts above the mode.
    var j = j0 + 1.0;
    while true {
        let lgj, sgj, cgj = lgamma(j + 1.0);
        if cgj != OK { return nan(), cgj; }
        let lw = -lam + j * llam - lgj;
        if lw < lcut { break; }
        if upper { g = g + e; } else { g = g - e; }
        e = e * half / (a + 1.0);
        a = a + 1.0;
        total += exp(lw) * g;
        j += 1.0;
    }
    // The counts below the mode.
    e = e0;
    g = g0;
    a = a0;
    j = j0 - 1.0;
    while j >= 0.0 {
        let lgj, sgj, cgj = lgamma(j + 1.0);
        if cgj != OK { return nan(), cgj; }
        let lw = -lam + j * llam - lgj;
        if lw < lcut { break; }
        if half > 0.0 { e = e * a / half; }
        a = a - 1.0;
        if upper { g = g - e; } else { g = g + e; }
        total += exp(lw) * g;
        j -= 1.0;
    }
    if total > 1.0 { total = 1.0; }
    return total, OK;
}

fn rice_pdf(x: f64, b: f64) -> f64, i64 {
    let k = rice_check(b);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    return rice_dens(x, b), OK;
}

fn rice_cdf(x: f64, b: f64) -> f64, i64 {
    let k = rice_check(b);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    // A shape so small that its square drops below the resolution of a
    // double gives the law of Rayleigh.
    if 0.5 * b * b == 0.0 {
        let m, cm = expm1(-0.5 * x * x);
        if cm != OK { return nan(), cm; }
        return -m, OK;
    }
    return rice_sum(x, b, false);
}

fn rice_sf(x: f64, b: f64) -> f64, i64 {
    let k = rice_check(b);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 1.0, OK; }
    if 0.5 * b * b == 0.0 { return exp(-0.5 * x * x), OK; }
    return rice_sum(x, b, true);
}

fn rice_ppf(p: f64, b: f64) -> f64, i64 {
    let k = rice_check(b);
    if k != OK { return nan(), k; }
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(p, 0.0, inf());
    if e { return ev, ec; }
    // The root of the cumulative function for the lower half and of the
    // survival function for the upper half keeps the digits of the tail.
    if p <= 0.5 {
        return fam_search_ppf(p, 0.0, b + 8.0) { v => rice_cdf_value(v, b) };
    }
    return fam_search_ppf_sf(1.0 - p, 0.0, b + 8.0) { v => rice_sf_value(v, b) };
}

// The cumulative function of the law, without the code. The search of
// the quantile takes this value.
fn rice_cdf_value(x: f64, b: f64) -> f64 {
    let v, c = rice_cdf(x, b);
    if c != OK { return nan(); }
    return v;
}

// The survival function of the law, without the code. The search of the
// quantile takes this value.
fn rice_sf_value(x: f64, b: f64) -> f64 {
    let v, c = rice_sf(x, b);
    if c != OK { return nan(); }
    return v;
}

fn rice_mean(b: f64) -> f64, i64 {
    let k = rice_check(b);
    if k != OK { return nan(), k; }
    let m1, m2, m3, m4, ck = fam_raws(1, 1.0) { x => rice_dens(x, b) };
    if ck != OK { return nan(), ck; }
    return m1, OK;
}

fn rice_var(b: f64) -> f64, i64 {
    let k = rice_check(b);
    if k != OK { return nan(), k; }
    let m1, m2, m3, m4, ck = fam_raws(2, 1.0) { x => rice_dens(x, b) };
    if ck != OK { return nan(), ck; }
    return fam_var(m1, m2), OK;
}

fn rice_skew(b: f64) -> f64, i64 {
    let k = rice_check(b);
    if k != OK { return nan(), k; }
    let m1, m2, m3, m4, ck = fam_raws(3, 1.0) { x => rice_dens(x, b) };
    if ck != OK { return nan(), ck; }
    return fam_skew(m1, m2, m3), OK;
}

fn rice_kurt(b: f64) -> f64, i64 {
    let k = rice_check(b);
    if k != OK { return nan(), k; }
    let m1, m2, m3, m4, ck = fam_raws(4, 1.0) { x => rice_dens(x, b) };
    if ck != OK { return nan(), ck; }
    return fam_kurt(m1, m2, m3, m4), OK;
}
