// The Landau law as SciPy 1.18.1 gives it, with a location of 0 and a
// scale of 1. The standard density is the integral over t from 0 to
// infinity of exp(-t log t - u t) sin(pi t) / pi. SciPy takes the scale
// of the library Boost: its density at x is pi / 2 times the standard
// density at u = pi x / 2 + log(pi / 2). The survival function has the
// same integral with a divisor t.
//
// The law has no finite mean. The left tail falls faster than any
// exponential, so a point below LANDAU_LEFT gives a density of 0.
//
// The call has three forms: pdf, cdf and ppf. Every call returns the
// result first and a fault code last.
namespace statistics;

import std;
import math;

const LANDAU_LEFT: f64 = -3.5;          // below this standard point the call gives 0
const LANDAU_X_LEFT: f64 = -2.5;        // the same point in the scale of the call
const LANDAU_PANELS: i64 = 60000;       // the bound on the panels of one integral
const LANDAU_CUT: f64 = -50.0;          // the exponent where the integral stops

// The exponent of the integrand.
fn landau_expo(t: f64, x: f64) -> f64 {
    if t <= 0.0 { return 0.0; }
    return -t * log(t) - x * t;
}

// The integral of the integrand over t. The flag tells whether the
// integrand holds the divisor t. The panels have a width of 1 and the
// integral stops once the exponent falls under the cut.
fn landau_integral(x: f64, with_div: bool) -> f64, i64 {
    var total: f64 = 0.0;
    var lo: f64 = 0.0;
    var n: i64 = 0;
    while n < LANDAU_PANELS {
        let hi = lo + 1.0;
        total += fam_quad(lo, hi) { t => landau_term(t, x, with_div) };
        // The exponent falls for t above its peak. It falls below the cut
        // soon after.
        if landau_expo(hi, x) < LANDAU_CUT && hi > exp(-x - 1.0) { return total, OK; }
        lo = hi;
        n += 1;
    }
    return nan(), ERR_FIT;
}

// One term of the integrand.
fn landau_term(t: f64, x: f64, with_div: bool) -> f64 {
    if t <= 0.0 {
        // The limit of sin(pi t) / t is pi.
        if with_div { return PI; }
        return 0.0;
    }
    let e = exp(landau_expo(t, x)) * sin(PI * t);
    if with_div { return e / t; }
    return e;
}

// The standard point of a point of the call.
fn landau_std(x: f64) -> f64 { return 0.5 * PI * x + log(0.5 * PI); }

fn landau_pdf(x: f64) -> f64, i64 {
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    let u = landau_std(x);
    if u < LANDAU_LEFT { return 0.0, OK; }
    let v, c = landau_integral(u, false);
    if c != OK { return nan(), c; }
    return max(0.5 * v, 0.0), OK;
}

fn landau_sf_value(x: f64) -> f64 {
    let u = landau_std(x);
    if u < LANDAU_LEFT { return 1.0; }
    let v, c = landau_integral(u, true);
    if c != OK { return nan(); }
    return min(max(v / PI, 0.0), 1.0);
}

fn landau_cdf(x: f64) -> f64, i64 {
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if landau_std(x) < LANDAU_LEFT { return 0.0, OK; }
    let sf = landau_sf_value(x);
    if is_nan(sf) { return nan(), ERR_FIT; }
    return 1.0 - sf, OK;
}

fn landau_cdf_value(x: f64) -> f64 {
    let v, c = landau_cdf(x);
    if c != OK { return nan(); }
    return v;
}

fn landau_ppf(q: f64) -> f64, i64 {
    let k = check_prob(q);
    if k != OK { return nan(), k; }
    if q == 0.0 { return -inf(), OK; }
    if q == 1.0 { return inf(), OK; }
    return fam_search_ppf(q, LANDAU_X_LEFT, 1.0) { v => landau_cdf_value(v) };
}
