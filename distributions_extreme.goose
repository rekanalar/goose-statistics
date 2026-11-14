// The generalised extreme value law.
//
// The name and the parameters follow SciPy.
//
// The law has eight calls: pdf, cdf, sf, ppf, mean, var, skew and kurt.
// The kurtosis is the excess kurtosis. Every call returns the result
// first and a fault code last. A value outside the support gives a density
// of 0. The quantile of 0 or 1 gives the bound of the support.
namespace statistics;

import std;
import math;

// ---- genextreme: the generalised extreme value law ------------------------

// The machine epsilon of a double number. The cutoffs of the skewness and
// the kurtosis are powers of it, so the law of the Gumbel takes over when
// the shape is small.
const FAM_EPS_POW_029: f64 = 8.709635899560812e-05;
const FAM_EPS_POW_023: f64 = 0.0006025595860743576;

fn genextreme_pdf(x: f64, c: f64) -> f64, i64 {
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    let k = check_scalar(c);
    if k != OK { return nan(), k; }
    if c == 0.0 {
        // The law is then the Gumbel law.
        return exp(-x - exp(-x)), OK;
    }
    let cx = c * x;
    if cx >= 1.0 { return 0.0, OK; }
    let l, cl = log1p(-cx);
    if cl != OK { return nan(), cl; }
    let t = l / c;
    let lp = -exp(t) + (1.0 / c - 1.0) * l;
    if is_inf(lp) { return 0.0, OK; }
    return exp(lp), OK;
}

fn genextreme_cdf(x: f64, c: f64) -> f64, i64 {
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    let k = check_scalar(c);
    if k != OK { return nan(), k; }
    if c == 0.0 {
        let e = exp(-x);
        return exp(-e), OK;
    }
    let cx = c * x;
    if cx >= 1.0 {
        if c > 0.0 { return 1.0, OK; }
        return 0.0, OK;
    }
    let l, cl = log1p(-cx);
    if cl != OK { return nan(), cl; }
    return exp(-exp(l / c)), OK;
}

fn genextreme_sf(x: f64, c: f64) -> f64, i64 {
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    let k = check_scalar(c);
    if k != OK { return nan(), k; }
    if c == 0.0 {
        let e = exp(-x);
        let d, cd = expm1(-e);
        if cd != OK { return nan(), cd; }
        return -d, OK;
    }
    let cx = c * x;
    if cx >= 1.0 {
        if c > 0.0 { return 0.0, OK; }
        return 1.0, OK;
    }
    let l, cl = log1p(-cx);
    if cl != OK { return nan(), cl; }
    let d, cd = expm1(-exp(l / c));
    if cd != OK { return nan(), cd; }
    return -d, OK;
}

fn genextreme_ppf(q: f64, c: f64) -> f64, i64 {
    let pc = check_prob(q);
    if pc != OK { return nan(), pc; }
    let k = check_scalar(c);
    if k != OK { return nan(), k; }
    let e, ev, ec = fam_end(q, -inf(), inf());
    if e { return ev, ec; }
    let l = log(q);
    let x = -log(-l);
    if c == 0.0 { return x, OK; }
    let d, cd = expm1(-c * x);
    if cd != OK { return nan(), cd; }
    return -d / c, OK;
}

// The gamma function of the moments of the generalised extreme value law.
fn genextreme_gamma(n: f64, c: f64) -> f64 {
    let a = n * c + 1.0;
    if a <= 0.0 { return nan(); }
    return exp(fam_lngamma(a));
}

fn genextreme_mean(c: f64) -> f64, i64 {
    let k = check_scalar(c);
    if k != OK { return nan(), k; }
    if c < -1.0 { return nan(), ERR_UNDEFINED; }
    if abs(c) < 1e-14 { return EULER_GAMMA, OK; }
    let l = fam_lngamma(c + 1.0);
    if is_nan(l) { return nan(), ERR_UNDEFINED; }
    let d, cd = expm1(l);
    if cd != OK { return nan(), cd; }
    return -d / c, OK;
}

fn genextreme_var(c: f64) -> f64, i64 {
    let k = check_scalar(c);
    if k != OK { return nan(), k; }
    if c < -0.5 { return nan(), ERR_UNDEFINED; }
    if abs(c) < 1e-7 {
        let g = genextreme_gamma(1.0, c);
        let v = g * g * PI * PI / 6.0;
        if is_inf(v) { return nan(), ERR_OVERFLOW; }
        return v, OK;
    }
    let g1 = genextreme_gamma(1.0, c);
    if is_inf(g1) { return nan(), ERR_OVERFLOW; }
    let l = fam_lngamma(2.0 * c + 1.0) - 2.0 * fam_lngamma(c + 1.0);
    if is_nan(l) { return nan(), ERR_UNDEFINED; }
    let d, cd = expm1(l);
    if cd != OK { return nan(), cd; }
    let v = g1 * g1 * d / (c * c);
    if is_inf(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

fn genextreme_skew(c: f64) -> f64, i64 {
    let k = check_scalar(c);
    if k != OK { return nan(), k; }
    if abs(c) <= FAM_EPS_POW_029 { return 12.0 * sqrt(6.0) * ZETA3 / (PI * PI * PI), OK; }
    if c < -1.0 / 3.0 { return nan(), ERR_UNDEFINED; }
    let g1 = genextreme_gamma(1.0, c);
    let g2 = genextreme_gamma(2.0, c);
    let g3 = genextreme_gamma(3.0, c);
    if is_nan(g1) || is_nan(g2) || is_nan(g3) { return nan(), ERR_UNDEFINED; }
    if is_inf(g1) || is_inf(g2) || is_inf(g3) { return nan(), ERR_OVERFLOW; }
    let gm = genextreme_var_minus(c);
    if is_nan(gm) { return nan(), ERR_UNDEFINED; }
    if is_inf(gm) { return nan(), ERR_OVERFLOW; }
    let s = if c > 0.0 { 1.0 } else { -1.0 };
    let r = s * (-g3 + (g2 + 2.0 * gm) * g1) / (gm * sqrt(gm));
    if is_nan(r) { return nan(), ERR_OVERFLOW; }
    if is_inf(r) { return nan(), ERR_OVERFLOW; }
    return r, OK;
}

fn genextreme_kurt(c: f64) -> f64, i64 {
    let k = check_scalar(c);
    if k != OK { return nan(), k; }
    if abs(c) <= FAM_EPS_POW_023 { return 2.4, OK; }
    if c < -0.25 { return nan(), ERR_UNDEFINED; }
    let g1 = genextreme_gamma(1.0, c);
    let g2 = genextreme_gamma(2.0, c);
    let g3 = genextreme_gamma(3.0, c);
    let g4 = genextreme_gamma(4.0, c);
    if is_nan(g1) || is_nan(g2) || is_nan(g3) || is_nan(g4) { return nan(), ERR_UNDEFINED; }
    if is_inf(g1) || is_inf(g2) || is_inf(g3) || is_inf(g4) { return nan(), ERR_OVERFLOW; }
    let gm = genextreme_var_minus(c);
    if is_nan(gm) { return nan(), ERR_UNDEFINED; }
    if is_inf(gm) { return nan(), ERR_OVERFLOW; }
    let t = (g4 + (-4.0 * g3 + 3.0 * (g2 + gm) * g1) * g1) / (gm * gm) - 3.0;
    if is_nan(t) { return nan(), ERR_OVERFLOW; }
    if is_inf(t) { return nan(), ERR_OVERFLOW; }
    return t, OK;
}

// The variance without the mean square: the difference of the second raw
// moment and the square of the first one.
fn genextreme_var_minus(c: f64) -> f64 {
    if abs(c) < 1e-7 { return c * c * PI * PI / 6.0; }
    let g1 = genextreme_gamma(1.0, c);
    let g2 = genextreme_gamma(2.0, c);
    if is_nan(g1) || is_nan(g2) { return nan(); }
    return g2 - g1 * g1;
}
