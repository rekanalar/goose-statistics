// Fourteen continuous families with a closed form: alpha, fatiguelife,
// levy, levy_l, logistic, loglaplace, lomax, pareto, powerlaw, reciprocal,
// rdist, skewcauchy, wrapcauchy and kappa4. The names and the parameters
// follow SciPy 1.18.1, with a location of 0 and a scale of 1.
//
// Every family has three calls: pdf, cdf and ppf. Every call returns the
// result first and a fault code last. A value outside the support gives
// a density of 0 and a cumulative function of 0 or 1.
namespace statistics;

import std;
import math;

// The normal cumulative function and its inverse, without the code. The
// callers have checked their input already.
fn tl_ncdf(x: f64) -> f64 {
    let v, c = norm_cdf(x);
    if c != OK { return nan(); }
    return v;
}

fn tl_nppf(p: f64) -> f64 {
    let v, c = norm_ppf(p);
    if c != OK { return nan(); }
    return v;
}

// The result of a call: a NaN gives the fit fault.
fn tl_out(v: f64) -> f64, i64 {
    if is_nan(v) { return nan(), ERR_FIT; }
    return v, OK;
}

// The checks that every call of a family makes on the point and on the
// probability.
fn tl_px(x: f64, k: i64) -> i64 {
    if k != OK { return k; }
    return fam_x(x);
}

fn tl_pq(q: f64, k: i64) -> i64 {
    if k != OK { return k; }
    return check_prob(q);
}

// ---- alpha ----------------------------------------------------------------

fn alpha_pdf(x: f64, a: f64) -> f64, i64 {
    let k = tl_px(x, check_shape(a));
    if k != OK { return nan(), k; }
    if x <= 0.0 { return 0.0, OK; }
    let p, c = norm_pdf(a - 1.0 / x);
    if c != OK { return nan(), c; }
    return tl_out(p / (x * x * tl_ncdf(a)));
}

fn alpha_cdf(x: f64, a: f64) -> f64, i64 {
    let k = tl_px(x, check_shape(a));
    if k != OK { return nan(), k; }
    if x <= 0.0 { return 0.0, OK; }
    return tl_out(tl_ncdf(a - 1.0 / x) / tl_ncdf(a));
}

fn alpha_ppf(q: f64, a: f64) -> f64, i64 {
    let k = tl_pq(q, check_shape(a));
    if k != OK { return nan(), k; }
    if q == 0.0 { return 0.0, OK; }
    if q == 1.0 { return inf(), OK; }
    return tl_out(1.0 / (a - tl_nppf(q * tl_ncdf(a))));
}

// ---- fatiguelife ----------------------------------------------------------

fn fatiguelife_pdf(x: f64, c: f64) -> f64, i64 {
    let k = tl_px(x, check_shape(c));
    if k != OK { return nan(), k; }
    if x <= 0.0 { return 0.0, OK; }
    let l = log(x + 1.0) - log(2.0 * c) - 0.5 * log(2.0 * PI * x * x * x)
        - (x - 1.0) * (x - 1.0) / (2.0 * x * c * c);
    return tl_out(exp(l));
}

fn fatiguelife_cdf(x: f64, c: f64) -> f64, i64 {
    let k = tl_px(x, check_shape(c));
    if k != OK { return nan(), k; }
    if x <= 0.0 { return 0.0, OK; }
    let s = sqrt(x);
    return tl_out(tl_ncdf((s - 1.0 / s) / c));
}

fn fatiguelife_ppf(q: f64, c: f64) -> f64, i64 {
    let k = tl_pq(q, check_shape(c));
    if k != OK { return nan(), k; }
    if q == 0.0 { return 0.0, OK; }
    if q == 1.0 { return inf(), OK; }
    let t = c * tl_nppf(q);
    let r = t + sqrt(t * t + 4.0);
    return tl_out(0.25 * r * r);
}

// ---- levy: the stable law of a shape one half -----------------------------

fn levy_pdf(x: f64) -> f64, i64 {
    let k = fam_x(x);
    if k != OK { return nan(), k; }
    if x <= 0.0 { return 0.0, OK; }
    return tl_out(exp(-1.0 / (2.0 * x) - 0.5 * log(2.0 * PI * x * x * x)));
}

fn levy_cdf(x: f64) -> f64, i64 {
    let k = fam_x(x);
    if k != OK { return nan(), k; }
    if x <= 0.0 { return 0.0, OK; }
    // The upper tail of the normal law keeps the digits of a small value.
    return tl_out(2.0 * tl_ncdf(-1.0 / sqrt(x)));
}

fn levy_ppf(q: f64) -> f64, i64 {
    let k = check_prob(q);
    if k != OK { return nan(), k; }
    if q == 0.0 { return 0.0, OK; }
    if q == 1.0 { return inf(), OK; }
    let v = tl_nppf(1.0 - 0.5 * q);
    return tl_out(1.0 / (v * v));
}

// ---- levy_l: the mirror of levy -------------------------------------------

fn levy_l_pdf(x: f64) -> f64, i64 {
    let k = fam_x(x);
    if k != OK { return nan(), k; }
    if x >= 0.0 { return 0.0, OK; }
    return levy_pdf(-x);
}

fn levy_l_cdf(x: f64) -> f64, i64 {
    let k = fam_x(x);
    if k != OK { return nan(), k; }
    if x >= 0.0 { return 1.0, OK; }
    // The survival function of the mirror, from the upper tail.
    return tl_out(1.0 - 2.0 * tl_ncdf(-1.0 / sqrt(-x)));
}

fn levy_l_ppf(q: f64) -> f64, i64 {
    let k = check_prob(q);
    if k != OK { return nan(), k; }
    if q == 0.0 { return -inf(), OK; }
    if q == 1.0 { return 0.0, OK; }
    let v = tl_nppf(0.5 * (q + 1.0));
    return tl_out(-1.0 / (v * v));
}

// ---- logistic -------------------------------------------------------------

fn logistic_pdf(x: f64) -> f64, i64 {
    let k = fam_x(x);
    if k != OK { return nan(), k; }
    // The form with the absolute value never overflows.
    let e = exp(-abs(x));
    return tl_out(e / ((1.0 + e) * (1.0 + e)));
}

fn logistic_cdf(x: f64) -> f64, i64 {
    let k = fam_x(x);
    if k != OK { return nan(), k; }
    if x >= 0.0 { return tl_out(1.0 / (1.0 + exp(-x))); }
    let e = exp(x);
    return tl_out(e / (1.0 + e));
}

fn logistic_ppf(q: f64) -> f64, i64 {
    let k = check_prob(q);
    if k != OK { return nan(), k; }
    if q == 0.0 { return -inf(), OK; }
    if q == 1.0 { return inf(), OK; }
    return tl_out(log(q) - log1p_raw(-q));
}

// The log of 1 + x, without the code. The callers pass x above -1.
fn log1p_raw(x: f64) -> f64 {
    let v, c = log1p(x);
    if c != OK { return nan(); }
    return v;
}

fn expm1_raw(x: f64) -> f64 {
    let v, c = expm1(x);
    if c != OK { return nan(); }
    return v;
}

// ---- loglaplace -----------------------------------------------------------

fn loglaplace_pdf(x: f64, c: f64) -> f64, i64 {
    let k = tl_px(x, check_shape(c));
    if k != OK { return nan(), k; }
    if x <= 0.0 { return 0.0, OK; }
    if x < 1.0 { return tl_out(0.5 * c * pow(x, c - 1.0)); }
    return tl_out(0.5 * c * pow(x, -c - 1.0));
}

fn loglaplace_cdf(x: f64, c: f64) -> f64, i64 {
    let k = tl_px(x, check_shape(c));
    if k != OK { return nan(), k; }
    if x <= 0.0 { return 0.0, OK; }
    if x < 1.0 { return tl_out(0.5 * pow(x, c)); }
    return tl_out(1.0 - 0.5 * pow(x, -c));
}

fn loglaplace_ppf(q: f64, c: f64) -> f64, i64 {
    let k = tl_pq(q, check_shape(c));
    if k != OK { return nan(), k; }
    if q == 0.0 { return 0.0, OK; }
    if q == 1.0 { return inf(), OK; }
    if q < 0.5 { return tl_out(pow(2.0 * q, 1.0 / c)); }
    return tl_out(pow(2.0 * (1.0 - q), -1.0 / c));
}

// ---- lomax ----------------------------------------------------------------

fn lomax_pdf(x: f64, c: f64) -> f64, i64 {
    let k = tl_px(x, check_shape(c));
    if k != OK { return nan(), k; }
    if x < 0.0 { return 0.0, OK; }
    return tl_out(c * exp(-(c + 1.0) * log1p_raw(x)));
}

fn lomax_cdf(x: f64, c: f64) -> f64, i64 {
    let k = tl_px(x, check_shape(c));
    if k != OK { return nan(), k; }
    if x <= 0.0 { return 0.0, OK; }
    return tl_out(-expm1_raw(-c * log1p_raw(x)));
}

fn lomax_ppf(q: f64, c: f64) -> f64, i64 {
    let k = tl_pq(q, check_shape(c));
    if k != OK { return nan(), k; }
    if q == 1.0 { return inf(), OK; }
    return tl_out(expm1_raw(-log1p_raw(-q) / c));
}

// ---- pareto ---------------------------------------------------------------

fn pareto_pdf(x: f64, b: f64) -> f64, i64 {
    let k = tl_px(x, check_shape(b));
    if k != OK { return nan(), k; }
    if x < 1.0 { return 0.0, OK; }
    return tl_out(b * pow(x, -b - 1.0));
}

fn pareto_cdf(x: f64, b: f64) -> f64, i64 {
    let k = tl_px(x, check_shape(b));
    if k != OK { return nan(), k; }
    if x <= 1.0 { return 0.0, OK; }
    return tl_out(-expm1_raw(-b * log(x)));
}

fn pareto_ppf(q: f64, b: f64) -> f64, i64 {
    let k = tl_pq(q, check_shape(b));
    if k != OK { return nan(), k; }
    if q == 1.0 { return inf(), OK; }
    return tl_out(exp(-log1p_raw(-q) / b));
}

// ---- powerlaw -------------------------------------------------------------

fn powerlaw_pdf(x: f64, a: f64) -> f64, i64 {
    let k = tl_px(x, check_shape(a));
    if k != OK { return nan(), k; }
    if x < 0.0 || x > 1.0 { return 0.0, OK; }
    return tl_out(a * pow(x, a - 1.0));
}

fn powerlaw_cdf(x: f64, a: f64) -> f64, i64 {
    let k = tl_px(x, check_shape(a));
    if k != OK { return nan(), k; }
    if x <= 0.0 { return 0.0, OK; }
    if x >= 1.0 { return 1.0, OK; }
    return tl_out(pow(x, a));
}

fn powerlaw_ppf(q: f64, a: f64) -> f64, i64 {
    let k = tl_pq(q, check_shape(a));
    if k != OK { return nan(), k; }
    return tl_out(pow(q, 1.0 / a));
}

// ---- reciprocal: the same law as loguniform -------------------------------

fn reciprocal_pdf(x: f64, a: f64, b: f64) -> f64, i64 { return loguniform_pdf(x, a, b); }
fn reciprocal_cdf(x: f64, a: f64, b: f64) -> f64, i64 { return loguniform_cdf(x, a, b); }
fn reciprocal_ppf(q: f64, a: f64, b: f64) -> f64, i64 { return loguniform_ppf(q, a, b); }

// ---- rdist ----------------------------------------------------------------

// The log of the beta function of c / 2 and c / 2.
fn rdist_lbeta(c: f64) -> f64 {
    return 2.0 * fam_lngamma(0.5 * c) - fam_lngamma(c);
}

fn rdist_pdf(x: f64, c: f64) -> f64, i64 {
    let k = tl_px(x, check_shape(c));
    if k != OK { return nan(), k; }
    if x <= -1.0 || x >= 1.0 {
        // The density at the bounds is finite only for a shape of 2 or more.
        if abs(x) == 1.0 && c == 2.0 { return tl_out(exp(-rdist_lbeta(c) - (c - 1.0) * log(2.0))); }
        return 0.0, OK;
    }
    let l = (0.5 * c - 1.0) * log1p_raw(-x * x) - rdist_lbeta(c) - (c - 1.0) * log(2.0);
    return tl_out(exp(l));
}

fn rdist_cdf(x: f64, c: f64) -> f64, i64 {
    let k = tl_px(x, check_shape(c));
    if k != OK { return nan(), k; }
    if x <= -1.0 { return 0.0, OK; }
    if x >= 1.0 { return 1.0, OK; }
    return betainc(0.5 * c, 0.5 * c, 0.5 * (x + 1.0));
}

fn rdist_ppf(q: f64, c: f64) -> f64, i64 {
    let k = tl_pq(q, check_shape(c));
    if k != OK { return nan(), k; }
    if q == 0.0 { return -1.0, OK; }
    if q == 1.0 { return 1.0, OK; }
    let v, cv = betaincinv(0.5 * c, 0.5 * c, q);
    if cv != OK { return nan(), cv; }
    return 2.0 * v - 1.0, OK;
}

// ---- skewcauchy -----------------------------------------------------------

fn skewcauchy_check(a: f64) -> i64 {
    let k = check_scalar(a);
    if k != OK { return k; }
    if abs(a) >= 1.0 { return ERR_DOMAIN; }
    return OK;
}

fn skewcauchy_pdf(x: f64, a: f64) -> f64, i64 {
    let k = tl_px(x, skewcauchy_check(a));
    if k != OK { return nan(), k; }
    var s = 1.0 + a;
    if x < 0.0 { s = 1.0 - a; }
    return tl_out(1.0 / (PI * (x * x / (s * s) + 1.0)));
}

fn skewcauchy_cdf(x: f64, a: f64) -> f64, i64 {
    let k = tl_px(x, skewcauchy_check(a));
    if k != OK { return nan(), k; }
    if x <= 0.0 { return tl_out(0.5 * (1.0 - a) + (1.0 - a) / PI * atan(x / (1.0 - a))); }
    return tl_out(0.5 * (1.0 - a) + (1.0 + a) / PI * atan(x / (1.0 + a)));
}

fn skewcauchy_ppf(q: f64, a: f64) -> f64, i64 {
    let k = tl_pq(q, skewcauchy_check(a));
    if k != OK { return nan(), k; }
    if q == 0.0 { return -inf(), OK; }
    if q == 1.0 { return inf(), OK; }
    let mid = 0.5 * (1.0 - a);
    if q < mid { return tl_out(tan(PI / (1.0 - a) * (q - mid)) * (1.0 - a)); }
    return tl_out(tan(PI / (1.0 + a) * (q - mid)) * (1.0 + a));
}

// ---- wrapcauchy -----------------------------------------------------------

fn wrapcauchy_check(c: f64) -> i64 {
    let k = check_scalar(c);
    if k != OK { return k; }
    if c <= 0.0 || c >= 1.0 { return ERR_DOMAIN; }
    return OK;
}

fn wrapcauchy_pdf(x: f64, c: f64) -> f64, i64 {
    let k = tl_px(x, wrapcauchy_check(c));
    if k != OK { return nan(), k; }
    if x < 0.0 || x > 2.0 * PI { return 0.0, OK; }
    return tl_out((1.0 - c * c) / (2.0 * PI * (1.0 + c * c - 2.0 * c * cos(x))));
}

fn wrapcauchy_cdf(x: f64, c: f64) -> f64, i64 {
    let k = tl_px(x, wrapcauchy_check(c));
    if k != OK { return nan(), k; }
    if x <= 0.0 { return 0.0, OK; }
    if x >= 2.0 * PI { return 1.0, OK; }
    let r = (1.0 + c) / (1.0 - c);
    if x < PI { return tl_out(atan(r * tan(0.5 * x)) / PI); }
    if x == PI { return 0.5, OK; }
    return tl_out(1.0 - atan(r * tan(0.5 * (2.0 * PI - x))) / PI);
}

fn wrapcauchy_ppf(q: f64, c: f64) -> f64, i64 {
    let k = tl_pq(q, wrapcauchy_check(c));
    if k != OK { return nan(), k; }
    if q == 0.0 { return 0.0, OK; }
    if q == 1.0 { return 2.0 * PI, OK; }
    if q == 0.5 { return PI, OK; }
    let v = (1.0 - c) / (1.0 + c);
    if q < 0.5 { return tl_out(2.0 * atan(v * tan(PI * q))); }
    return tl_out(2.0 * PI - 2.0 * atan(v * tan(PI * (1.0 - q))));
}

// ---- kappa4 ---------------------------------------------------------------

fn kappa4_check(h: f64, k: f64) -> i64 {
    let c1 = check_scalar(h);
    if c1 != OK { return c1; }
    return check_scalar(k);
}

// The value 1 - k x, and the value 1 - h (1 - k x)^(1 / k). The flag
// tells whether the point lies inside the support.
fn kappa4_parts(x: f64, h: f64, k: f64) -> f64, f64, bool {
    if k == 0.0 {
        let t = 1.0;
        var u = 1.0;
        if h != 0.0 { u = 1.0 - h * exp(-x); } else { u = exp(-x); }
        return t, u, true;
    }
    let y = 1.0 - k * x;
    if y <= 0.0 { return y, 0.0, false; }
    let p = pow(y, 1.0 / k);
    if h == 0.0 { return y, p, true; }
    let u = 1.0 - h * p;
    return y, u, u > 0.0;
}

fn kappa4_pdf(x: f64, h: f64, k: f64) -> f64, i64 {
    let c = tl_px(x, kappa4_check(h, k));
    if c != OK { return nan(), c; }
    let y, u, ok = kappa4_parts(x, h, k);
    if !ok { return 0.0, OK; }
    var l: f64 = 0.0;
    if k == 0.0 {
        if h == 0.0 { l = -x - exp(-x); } else { l = -x + (1.0 / h - 1.0) * log(u); }
    } else if h == 0.0 {
        l = (1.0 / k - 1.0) * log(y) - u;
    } else {
        l = (1.0 / k - 1.0) * log(y) + (1.0 / h - 1.0) * log(u);
    }
    return tl_out(exp(l));
}

fn kappa4_cdf(x: f64, h: f64, k: f64) -> f64, i64 {
    let c = tl_px(x, kappa4_check(h, k));
    if c != OK { return nan(), c; }
    let y, u, ok = kappa4_parts(x, h, k);
    if !ok {
        // Above the upper bound of the support the value is 1.
        if k > 0.0 && y <= 0.0 { return 1.0, OK; }
        return 0.0, OK;
    }
    if h == 0.0 { return tl_out(exp(-u)); }
    return tl_out(exp(log(u) / h));
}

fn kappa4_ppf(q: f64, h: f64, k: f64) -> f64, i64 {
    let c = tl_pq(q, kappa4_check(h, k));
    if c != OK { return nan(), c; }
    if q == 0.0 || q == 1.0 {
        // The bounds of the support follow the limits of the formula.
        if k == 0.0 && h == 0.0 { if q == 0.0 { return -inf(), OK; } return inf(), OK; }
    }
    if h != 0.0 && k != 0.0 { return tl_out(1.0 / k * (1.0 - pow((1.0 - pow(q, h)) / h, k))); }
    if h == 0.0 && k != 0.0 { return tl_out(1.0 / k * (1.0 - pow(-log(q), k))); }
    if h != 0.0 && k == 0.0 { return tl_out(-log((1.0 - pow(q, h)) / h)); }
    return tl_out(-log(-log(q)));
}
