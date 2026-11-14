// The stable law of the Lévy kind.
//
// The name and the parameters follow SciPy 1.18.1. SciPy uses the
// parameterization S1, where the variable has a location.
//
// The density and the cumulative function come from the integral of
// Nolan, which runs over an angle. The law needs the shape alpha and
// the skewness beta. The shape must be above 0 and at most 2. The
// skewness must be between -1 and 1.
//
// Every call returns the result first and a fault code last.
namespace statistics;

import std;
import math;

const LEVY_NODES: i64 = 40;           // nodes of one panel of the integral
const LEVY_ALPHA_TOL: f64 = 0.005;    // a shape near one becomes one
const LEVY_X_TOL: f64 = 0.005;        // a point near zeta becomes zeta
const LEVY_ROOT_STEPS: i64 = 60;      // steps of the search of a crossing

// The constants of the integral of Nolan.
struct LevyN {
    alpha: f64,
    zeta: f64,
    xi: f64,
    c1: f64,
    c2: f64,
    c3: f64,
    one: bool,
    zeta_prefactor: f64,
    alpha_exp: f64,
    alpha_xi: f64,
    zeta_offset: f64,
    two_beta_div_pi: f64,
    pi_div_two_beta: f64,
    x0_div_term: f64,
}

fn levy_check(alpha: f64, beta: f64) -> i64 {
    let c1 = check_scalar(alpha);
    if c1 != OK { return c1; }
    let c2 = check_scalar(beta);
    if c2 != OK { return c2; }
    if alpha <= 0.0 || alpha > 2.0 { return ERR_DOMAIN; }
    if abs(beta) > 1.0 { return ERR_DOMAIN; }
    return OK;
}

// The constants of the integral of Nolan. The point is the point of the
// form S0, where the variable has no location.
fn levy_nolan(alpha: f64, beta: f64, x0: f64) -> LevyN {
    let n = LevyN {
        alpha: alpha,
        zeta: 0.0,
        xi: 0.0,
        c1: 0.0,
        c2: 0.0,
        c3: 0.0,
        one: false,
        zeta_prefactor: 0.0,
        alpha_exp: 0.0,
        alpha_xi: 0.0,
        zeta_offset: 0.0,
        two_beta_div_pi: 0.0,
        pi_div_two_beta: 0.0,
        x0_div_term: 0.0
    };
    n.zeta = -beta * tan(0.5 * PI * alpha);
    if alpha == 1.0 {
        n.one = true;
        n.xi = 0.5 * PI;
        n.two_beta_div_pi = beta * 2.0 / PI;
        n.pi_div_two_beta = 0.5 * PI / beta;
        n.x0_div_term = x0 / n.two_beta_div_pi;
        n.c1 = 0.0;
        n.c2 = 0.5 / abs(beta);
        n.c3 = 1.0 / PI;
    } else {
        n.xi = atan(-n.zeta) / alpha;
        n.zeta_prefactor = pow(pow(n.zeta, 2.0) + 1.0, -1.0 / (2.0 * (alpha - 1.0)));
        n.alpha_exp = alpha / (alpha - 1.0);
        n.alpha_xi = atan(-n.zeta);
        n.zeta_offset = x0 - n.zeta;
        if alpha < 1.0 {
            n.c1 = 0.5 - n.xi / PI;
            n.c3 = 1.0 / PI;
        } else {
            n.c1 = 1.0;
            n.c3 = -1.0 / PI;
        }
        n.c2 = alpha / (PI * abs(alpha - 1.0) * (x0 - n.zeta));
    }
    return n;
}

// The helper of the integral of Nolan. The helper gives the argument of
// the exponential. The shape one has its own form.
fn levy_g(n: LevyN, theta: f64) -> f64 {
    if n.one {
        if theta == -n.xi { return 0.0; }
        if theta == 0.5 * PI { return inf(); }
        let c = cos(theta);
        if c == 0.0 { return inf(); }
        let v = (1.0 + theta * n.two_beta_div_pi)
            * exp((n.pi_div_two_beta + theta) * tan(theta) - n.x0_div_term) / c;
        if is_nan(v) { return 0.0; }
        return v;
    }
    if theta == -n.xi {
        if n.alpha < 1.0 { return 0.0; }
        return inf();
    }
    if theta == 0.5 * PI {
        if n.alpha < 1.0 { return inf(); }
        return 0.0;
    }
    let c = cos(theta);
    if c == 0.0 { return inf(); }
    let s = sin(n.alpha_xi + n.alpha * theta);
    if s == 0.0 {
        // The lower end of the range. The helper goes to 0 for a shape
        // below one and to infinity for a shape above one.
        if n.alpha < 1.0 { return 0.0; }
        return inf();
    }
    let base = c / s * n.zeta_offset;
    if base <= 0.0 { return inf(); }
    let p = n.zeta_prefactor * pow(base, n.alpha_exp);
    let r = cos(n.alpha_xi + (n.alpha - 1.0) * theta) / c;
    let v = p * r;
    if is_nan(v) { return 0.0; }
    if v < 0.0 { return 0.0; }
    return v;
}

// The value of the integrand of the law at an angle. The density of the
// law has the helper and the exponential of the helper in it.
fn levy_integrand(n: LevyN, want_pdf: bool, theta: f64) -> f64 {
    let g = levy_g(n, theta);
    if want_pdf {
        if g > 0.0 && !is_inf(g) { return g * exp(-g); }
        return 0.0;
    }
    if is_inf(g) { return 0.0; }
    return exp(-g);
}

// The integral of the law over one panel.
fn levy_panel(n: LevyN, want_pdf: bool, a: f64, b: f64, gx: const f64[:], gw: const f64[:]) -> f64 {
    let w = b - a;
    var total: f64 = 0.0;
    for j in 0..gx.len {
        let v = levy_integrand(n, want_pdf, a + w * (gx[j] + 1.0) / 2.0);
        if !is_nan(v) && !is_inf(v) { total += v * gw[j] * w / 2.0; }
    }
    return total;
}

// The heights of the helper where the integral is cut. The helper rises
// from 0 to infinity, or falls from infinity to 0, over the range of
// the angle. A cut at a height keeps the integrand of the panel near
// one. The heights cover the small values, where the integrand changes
// fast, and the large values, where the integrand is tiny.
fn levy_levels() -> f64[>..] {
    var out: f64[>..] = [];
    let nice: f64[12] = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 25.0, 100.0, 400.0];
    for i in 0..12 { out.push(nice[i]); }
    // The geometric heights from one hundred quadrillion to a hundred
    // million, with a ratio of two.
    var v = 1e-16;
    while v <= 1e8 {
        out.push(v);
        v = v * 2.0;
    }
    return out;
}

// The point of the law, snapped to zeta when it lies near it. The bound
// uses the shape of the law.
fn levy_snap(x: f64, alpha: f64, zeta: f64) -> f64 {
    let bound = LEVY_X_TOL * pow(alpha, 1.0 / alpha);
    if abs(x - zeta) < bound { return zeta; }
    return x;
}

// The angle where the helper reaches the height. The search fails when
// the helper does not reach the height, and the caller then leaves the
// angle out. The helper has one turning point at most in the range, so
// the search finds the crossing.
fn levy_cross(n: LevyN, h: f64) -> f64, i64 {
    var a = -n.xi;
    var b = 0.5 * PI;
    var fa = levy_g(n, a) - h;
    var fb = levy_g(n, b) - h;
    if is_nan(fa) || is_nan(fb) { return nan(), ERR_NAN; }
    if fa == 0.0 { return a, OK; }
    if fb == 0.0 { return b, OK; }
    if (fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0) { return nan(), ERR_DOMAIN; }
    var i: i64 = 0;
    while i < LEVY_ROOT_STEPS {
        let m = 0.5 * (a + b);
        let fm = levy_g(n, m) - h;
        if fm == 0.0 { return m, OK; }
        if is_nan(fm) { return nan(), ERR_NAN; }
        if (fa > 0.0 && fm > 0.0) || (fa < 0.0 && fm < 0.0) {
            a = m;
            fa = fm;
        } else {
            b = m;
        }
        i += 1;
    }
    return 0.5 * (a + b), OK;
}

// The angles that bracket the peak of the density. The peak sits where
// the helper equals one. The bracket runs from the angle where the
// helper is a half to the angle where it is two. The two angles come in
// either order, because the helper rises for a shape below one and
// falls for a shape above one.
fn levy_peak_bracket(n: LevyN) -> f64, f64, i64 {
    let lo, c1 = levy_cross(n, 0.5);
    if c1 != OK { return nan(), nan(), c1; }
    let hi, c2 = levy_cross(n, 2.0);
    if c2 != OK { return nan(), nan(), c2; }
    if lo < hi { return lo, hi, OK; }
    return hi, lo, OK;
}

// The integral of the law over the angle.
//
// The panels end at the angles where the helper reaches the heights of
// the list. The panel of the density near the height of one carries the
// peak. The peak itself gets its own panels, so that a narrow peak
// keeps its digits.
fn levy_integrate(n: LevyN, want_pdf: bool) -> f64 {
    var gx: f64[>..] = [];
    var gw: f64[>..] = [];
    gauss_legendre(gx, gw, LEVY_NODES);
    let levels = levy_levels();
    var pts: f64[>..] = [];
    pts.push(-n.xi);
    pts.push(0.5 * PI);
    for i in 0..levels.len {
        let t, c = levy_cross(n, levels[i]);
        if c == OK { pts.push(t); }
    }
    if want_pdf {
        // The peak of the density gets its own cut points. The bracket
        // of the peak is split into sixteen parts. The integrand of the
        // density is narrow, and a rule that never lands on it loses
        // digits.
        let blo, bhi, cp = levy_peak_bracket(n);
        if cp == OK {
            for i in 0..17 {
                pts.push(blo + (bhi - blo) * (i as f64) / 16.0);
            }
        }
    }
    // The points must rise. The count of the points is small.
    for i in 1..pts.len {
        let v = pts[i];
        var j = i;
        while j > 0 && pts[j - 1] > v {
            pts[j] = pts[j - 1];
            j -= 1;
        }
        pts[j] = v;
    }
    var total: f64 = 0.0;
    for k in 0..pts.len - 1 {
        let a = pts[k];
        let b = pts[k + 1];
        if b - a <= 0.0 { continue; }
        total += levy_panel(n, want_pdf, a, b, gx[..], gw[..]);
    }
    return total;
}

// The point and the shape of the law in the form of Nolan. The point
// comes from the form S1, where the variable has a location. A shape
// near one becomes one and a point near zeta becomes zeta, as
// SciPy 1.18.1 does.
fn levy_prepare(x: f64, alpha: f64, beta: f64) -> f64, f64 {
    let zeta = -beta * tan(0.5 * PI * alpha);
    var x0 = x;
    if alpha != 1.0 { x0 = x + zeta; }
    var a = alpha;
    if abs(alpha - 1.0) < LEVY_ALPHA_TOL { a = 1.0; }
    return levy_snap(x0, a, zeta), a;
}

// The density of the law in the form of Nolan, without the known cases.
recursive fn levy_pdf_core(x0: f64, alpha: f64, beta: f64) -> f64, i64 {
    // The shape two gives the normal law.
    if alpha == 2.0 {
        let v, c = norm_pdf(x0 / sqrt(2.0));
        if c != OK { return nan(), c; }
        return v / sqrt(2.0), OK;
    }
    let n = levy_nolan(alpha, beta, x0);
    let x = levy_snap(x0, alpha, n.zeta);
    if x == n.zeta {
        let l = fam_lngamma(1.0 + 1.0 / alpha);
        if is_nan(l) { return nan(), ERR_NAN; }
        let d = pow(n.zeta * n.zeta + 1.0, 1.0 / (2.0 * alpha));
        return exp(l) * cos(n.xi) / PI / d, OK;
    }
    // A point under zeta gives the law of the point with the other sign.
    if x < n.zeta {
        return levy_pdf_core(-x, alpha, -beta);
    }
    // The range of the angle is empty.
    if abs(-n.xi - 0.5 * PI) < 1e-14 { return 0.0, OK; }
    let intg = levy_integrate(n, true);
    return n.c2 * intg, OK;
}

// The cumulative function of the law in the form of Nolan, without the
// known cases.
recursive fn levy_cdf_core(x0: f64, alpha: f64, beta: f64) -> f64, i64 {
    if alpha == 2.0 {
        let v, c = norm_cdf(x0 / sqrt(2.0));
        if c != OK { return nan(), c; }
        return v, OK;
    }
    let n = levy_nolan(alpha, beta, x0);
    let x = levy_snap(x0, alpha, n.zeta);
    // A shape of one with a negative skewness, or a point under zeta,
    // gives the law of the point with the other sign.
    if (alpha == 1.0 && beta < 0.0) || x < n.zeta {
        let v, c = levy_cdf_core(-x, alpha, -beta);
        if c != OK { return nan(), c; }
        return 1.0 - v, OK;
    }
    if x == n.zeta { return 0.5 - n.xi / PI, OK; }
    if abs(-n.xi - 0.5 * PI) < 1e-14 { return n.c1, OK; }
    let intg = levy_integrate(n, false);
    return n.c1 + n.c3 * intg, OK;
}

// The density of the law in the form of Nolan, with the known cases.
fn levy_pdf_z0(x0: f64, alpha: f64, beta: f64) -> f64, i64 {
    if alpha == 0.5 && beta == 1.0 {
        // The law of Lévy.
        let xx = x0 + 1.0;
        if xx <= 0.0 { return 0.0, OK; }
        return exp(-0.5 / xx) / (sqrt(2.0 * PI) * xx * sqrt(xx)), OK;
    }
    if alpha == 1.0 && beta == 0.0 {
        // The law of Cauchy.
        return 1.0 / (1.0 + x0 * x0) / PI, OK;
    }
    return levy_pdf_core(x0, alpha, beta);
}

// The cumulative function of the law in the form of Nolan, with the
// known cases.
fn levy_cdf_z0(x0: f64, alpha: f64, beta: f64) -> f64, i64 {
    if alpha == 0.5 && beta == 1.0 {
        // The law of Lévy.
        let xx = x0 + 1.0;
        if xx <= 0.0 { return 0.0, OK; }
        let v, c = erfc(sqrt(0.5 / xx));
        if c != OK { return nan(), c; }
        return v, OK;
    }
    if alpha == 1.0 && beta == 0.0 {
        // The law of Cauchy.
        return 0.5 + atan(x0) / PI, OK;
    }
    return levy_cdf_core(x0, alpha, beta);
}

fn levy_stable_pdf(x: f64, alpha: f64, beta: f64) -> f64, i64 {
    let k = levy_check(alpha, beta);
    if k != OK { return nan(), k; }
    if is_nan(x) { return nan(), ERR_NAN; }
    if is_inf(x) { return 0.0, OK; }
    let x0, a = levy_prepare(x, alpha, beta);
    return levy_pdf_z0(x0, a, beta);
}

fn levy_stable_cdf(x: f64, alpha: f64, beta: f64) -> f64, i64 {
    let k = levy_check(alpha, beta);
    if k != OK { return nan(), k; }
    if is_nan(x) { return nan(), ERR_NAN; }
    if is_inf(x) {
        if x < 0.0 { return 0.0, OK; }
        return 1.0, OK;
    }
    let x0, a = levy_prepare(x, alpha, beta);
    return levy_cdf_z0(x0, a, beta);
}

fn levy_stable_cdf_value(x: f64, alpha: f64, beta: f64) -> f64 {
    let v, c = levy_stable_cdf(x, alpha, beta);
    if c != OK { return nan(); }
    return v;
}

fn levy_stable_ppf(q: f64, alpha: f64, beta: f64) -> f64, i64 {
    let k = levy_check(alpha, beta);
    if k != OK { return nan(), k; }
    let pc = check_prob(q);
    if pc != OK { return nan(), pc; }
    if q == 0.0 { return -inf(), OK; }
    if q == 1.0 { return inf(), OK; }
    return fam_search_ppf(q, -40.0, 0.0) { v => levy_stable_cdf_value(v, alpha, beta) };
}
