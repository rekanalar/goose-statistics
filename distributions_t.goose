// The Student t and the F families of the package. Each family has a
// density, a cumulative function, a survival function, a quantile
// function and the four moments.
//
// Every call returns the result first and a fault code last. The checks
// run in one order: NaN, infinity, domain, pole.
namespace statistics;

import std;
import math;

fn tdist_pdf(x: f64, df: f64 = 1.0) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    let c = chk_df(df);
    if c == ERR_NAN { return nan(), c; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    if c != OK { return nan(), c; }
    // 1 / (sqrt(df) B(df/2, 1/2)) (1 + x^2/df)^(-(df+1)/2)
    let q = dd_div(dd_mul(dd_from(x), dd_from(x)), dd_from(df));
    // x^2 overflows above 1e154, and the density falls below the smallest subnormal
    if is_inf(q.hi) { return 0.0, OK; }
    let l1p = dd_log(dd_add_d(q, 1.0));
    var l = dd_neg(dd_mul_d(l1p, 0.5 * (df + 1.0)));
    l = dd_sub(l, dd_lnbeta(0.5 * df, 0.5));
    l = dd_sub(l, dd_mul_d(dd_log_d(df), 0.5));
    return exp_or_overflow(l);
}

// The pair (lower, upper) of the Student t tails at x. The tail beyond
// |x| is half the incomplete beta at df/(df+x^2) for |x| above the root of
// df, and the complement of the central part below it. Both forms keep
// the digits of the small tail, because each passes the exact argument.
fn tdist_tails(x: f64, df: f64) -> f64, f64, i64 {
    let ax = abs(x);
    let a2 = dd_mul(dd_from(ax), dd_from(ax));
    let den = dd_add_d(a2, df);
    var tail: f64 = 0.0;
    if ax * ax < df {
        let y = dd_div(a2, den);
        let t, c = betainc_upper(0.5, 0.5 * df, y.hi);
        if c != OK { return nan(), nan(), c; }
        tail = 0.5 * beta_step(t, 0.5, 0.5 * df, y, -1.0);
    } else {
        let w = dd_div(dd_from(df), den);
        let t, c = betainc(0.5 * df, 0.5, w.hi);
        if c != OK { return nan(), nan(), c; }
        tail = 0.5 * beta_step(t, 0.5 * df, 0.5, w, 1.0);
    }
    if x >= 0.0 { return 1.0 - tail, tail, OK; }
    return tail, 1.0 - tail, OK;
}

fn tdist_cdf(x: f64, df: f64 = 1.0) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    let c = chk_df(df);
    if c == ERR_NAN { return nan(), c; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    if c != OK { return nan(), c; }
    let lo, up, c2 = tdist_tails(x, df);
    return lo, c2;
}

fn tdist_sf(x: f64, df: f64 = 1.0) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    let c = chk_df(df);
    if c == ERR_NAN { return nan(), c; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    if c != OK { return nan(), c; }
    let lo, up, c2 = tdist_tails(x, df);
    return up, c2;
}

fn tdist_ppf(p: f64, df: f64 = 1.0) -> f64, i64 {
    if is_nan(p) { return nan(), ERR_NAN; }
    let c = chk_df(df);
    if c != OK { return nan(), c; }
    let c2 = chk_open(p);
    if c2 != OK { return nan(), c2; }
    if p == 0.5 { return 0.0, OK; }
    let tail = min(p, 1.0 - p);
    let d = abs(2.0 * p - 1.0);
    var t: f64 = 0.0;
    if d < 0.5 {
        let y, c3 = betaincinv(0.5, 0.5 * df, d);
        if c3 != OK { return nan(), c3; }
        t = sqrt(df * y / (1.0 - y));
    } else {
        let w, c3 = betaincinv(0.5 * df, 0.5, 2.0 * tail);
        if c3 != OK { return nan(), c3; }
        t = sqrt(df * (1.0 - w) / w);
    }
    if is_inf(t) || is_nan(t) { return nan(), ERR_OVERFLOW; }
    if p < 0.5 { t = -t; }
    return tdist_polish(t, p, df), OK;
}

// One Newton step on the cumulative function removes the error of the
// inverse of the incomplete beta. The step needs a density above 0 and it
// never moves the value by more than a part in a million.
fn tdist_polish(t: f64, p: f64, df: f64) -> f64 {
    let pdf, c1 = tdist_pdf(t, df);
    if c1 != OK || pdf <= 0.0 { return t; }
    var step: f64 = 0.0;
    if p < 0.5 {
        let f, c2 = tdist_cdf(t, df);
        if c2 != OK { return t; }
        step = (f - p) / pdf;
    } else {
        let f, c2 = tdist_sf(t, df);
        if c2 != OK { return t; }
        step = -(f - (1.0 - p)) / pdf;
    }
    if abs(step) > 1e-6 * abs(t) { return t; }
    return t - step;
}

fn tdist_mean(df: f64 = 1.0) -> f64, i64 {
    let c = chk_df(df);
    if c != OK { return nan(), c; }
    if df <= 1.0 { return nan(), ERR_UNDEFINED; }
    return 0.0, OK;
}

fn tdist_var(df: f64 = 1.0) -> f64, i64 {
    let c = chk_df(df);
    if c != OK { return nan(), c; }
    if df <= 2.0 { return nan(), ERR_UNDEFINED; }
    return df / (df - 2.0), OK;
}

fn tdist_skew(df: f64 = 1.0) -> f64, i64 {
    let c = chk_df(df);
    if c != OK { return nan(), c; }
    if df <= 3.0 { return nan(), ERR_UNDEFINED; }
    return 0.0, OK;
}

fn tdist_kurt(df: f64 = 1.0) -> f64, i64 {
    let c = chk_df(df);
    if c != OK { return nan(), c; }
    if df <= 4.0 { return nan(), ERR_UNDEFINED; }
    return 6.0 / (df - 4.0), OK;
}

// ---- the F family ----------------------------------------------------------------------

fn fdist_check(d1: f64, d2: f64) -> i64 {
    if is_nan(d1) || is_nan(d2) { return ERR_NAN; }
    return first_code(chk_df(d1), chk_df(d2));
}

fn fdist_pdf(x: f64, d1: f64 = 1.0, d2: f64 = 1.0) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    let c = fdist_check(d1, d2);
    if c == ERR_NAN { return nan(), c; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    if c != OK { return nan(), c; }
    if x < 0.0 { return 0.0, OK; }
    if x == 0.0 {
        if d1 < 2.0 { return nan(), ERR_OVERFLOW; }
        if d1 == 2.0 { return 1.0, OK; }
        return 0.0, OK;
    }
    // x^(d1/2-1) (d1/d2)^(d1/2) (1 + d1 x/d2)^(-(d1+d2)/2) / B(d1/2, d2/2)
    let r = dd_div(dd_mul(dd_from(d1), dd_from(x)), dd_from(d2));
    var l = dd_mul(dd_log_d(x), dd_from(0.5 * d1 - 1.0));
    l = dd_add(l, dd_mul_d(dd_sub(dd_log_d(d1), dd_log_d(d2)), 0.5 * d1));
    l = dd_sub(l, dd_mul_d(dd_log(dd_add_d(r, 1.0)), 0.5 * (d1 + d2)));
    l = dd_sub(l, dd_lnbeta(0.5 * d1, 0.5 * d2));
    return exp_or_overflow(l);
}

// The pair u = d1 x / (d1 x + d2) and the pair w = d2 / (d1 x + d2).
fn fdist_ratios(x: f64, d1: f64, d2: f64) -> DD, DD {
    let h, l = two_prod(d1, x);
    let t = DD { hi: h, lo: l };
    let s = dd_add_d(t, d2);
    return dd_div(t, s), dd_div(dd_from(d2), s);
}

fn fdist_cdf(x: f64, d1: f64 = 1.0, d2: f64 = 1.0) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    let c = fdist_check(d1, d2);
    if c == ERR_NAN { return nan(), c; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    if c != OK { return nan(), c; }
    if x <= 0.0 { return 0.0, OK; }
    let t = d1 * x;
    if is_inf(t) { return 1.0, OK; }
    let u, w = fdist_ratios(x, d1, d2);
    if t <= d2 {
        let v, c2 = betainc(0.5 * d1, 0.5 * d2, u.hi);
        if c2 != OK { return nan(), c2; }
        return beta_step(v, 0.5 * d1, 0.5 * d2, u, 1.0), OK;
    }
    let v, c2 = betainc_upper(0.5 * d2, 0.5 * d1, w.hi);
    if c2 != OK { return nan(), c2; }
    return beta_step(v, 0.5 * d2, 0.5 * d1, w, -1.0), OK;
}

fn fdist_sf(x: f64, d1: f64 = 1.0, d2: f64 = 1.0) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    let c = fdist_check(d1, d2);
    if c == ERR_NAN { return nan(), c; }
    if is_inf(x) { return nan(), ERR_INFINITY; }
    if c != OK { return nan(), c; }
    if x <= 0.0 { return 1.0, OK; }
    let t = d1 * x;
    if is_inf(t) { return 0.0, OK; }
    let u, w = fdist_ratios(x, d1, d2);
    if t >= d2 {
        let v, c2 = betainc(0.5 * d2, 0.5 * d1, w.hi);
        if c2 != OK { return nan(), c2; }
        return beta_step(v, 0.5 * d2, 0.5 * d1, w, 1.0), OK;
    }
    let v, c2 = betainc_upper(0.5 * d1, 0.5 * d2, u.hi);
    if c2 != OK { return nan(), c2; }
    return beta_step(v, 0.5 * d1, 0.5 * d2, u, -1.0), OK;
}

fn fdist_ppf(p: f64, d1: f64 = 1.0, d2: f64 = 1.0) -> f64, i64 {
    if is_nan(p) { return nan(), ERR_NAN; }
    let c = fdist_check(d1, d2);
    if c != OK { return nan(), c; }
    let c2 = chk_open(p);
    if c2 != OK { return nan(), c2; }
    var r: f64 = 0.0;
    if p <= 0.5 {
        let u, c3 = betaincinv(0.5 * d1, 0.5 * d2, p);
        if c3 != OK { return nan(), c3; }
        r = d2 * u / (d1 * (1.0 - u));
    } else {
        let v, c3 = betaincinv(0.5 * d2, 0.5 * d1, 1.0 - p);
        if c3 != OK { return nan(), c3; }
        r = d2 * (1.0 - v) / (d1 * v);
    }
    if is_inf(r) { return nan(), ERR_OVERFLOW; }
    return fdist_polish(r, p, d1, d2), OK;
}

// One Newton step on the cumulative function, as for the Student t.
fn fdist_polish(x: f64, p: f64, d1: f64, d2: f64) -> f64 {
    let pdf, c1 = fdist_pdf(x, d1, d2);
    if c1 != OK || pdf <= 0.0 || is_inf(pdf) { return x; }
    var step: f64 = 0.0;
    if p <= 0.5 {
        let f, c2 = fdist_cdf(x, d1, d2);
        if c2 != OK { return x; }
        step = (f - p) / pdf;
    } else {
        let f, c2 = fdist_sf(x, d1, d2);
        if c2 != OK { return x; }
        step = -(f - (1.0 - p)) / pdf;
    }
    if abs(step) > 1e-6 * abs(x) { return x; }
    return x - step;
}

fn fdist_mean(d1: f64 = 1.0, d2: f64 = 1.0) -> f64, i64 {
    let c = fdist_check(d1, d2);
    if c != OK { return nan(), c; }
    if d2 <= 2.0 { return nan(), ERR_UNDEFINED; }
    return d2 / (d2 - 2.0), OK;
}

fn fdist_var(d1: f64 = 1.0, d2: f64 = 1.0) -> f64, i64 {
    let c = fdist_check(d1, d2);
    if c != OK { return nan(), c; }
    if d2 <= 4.0 { return nan(), ERR_UNDEFINED; }
    let r = d2 / (d2 - 2.0);
    let v = 2.0 * r * r * ((d1 + d2 - 2.0) / d1) / (d2 - 4.0);
    if is_inf(v) || is_nan(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

fn fdist_skew(d1: f64 = 1.0, d2: f64 = 1.0) -> f64, i64 {
    let c = fdist_check(d1, d2);
    if c != OK { return nan(), c; }
    if d2 <= 6.0 { return nan(), ERR_UNDEFINED; }
    let v = (2.0 * d1 + d2 - 2.0) / (d2 - 6.0) * sqrt(8.0 * (d2 - 4.0) / d1 / (d1 + d2 - 2.0));
    if is_inf(v) || is_nan(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

fn fdist_kurt(d1: f64 = 1.0, d2: f64 = 1.0) -> f64, i64 {
    let c = fdist_check(d1, d2);
    if c != OK { return nan(), c; }
    if d2 <= 8.0 { return nan(), ERR_UNDEFINED; }
    let a = (5.0 * d2 - 22.0) / (d2 - 8.0);
    let b = ((d2 - 4.0) / d1) * ((d2 - 2.0) / (d2 - 8.0)) * ((d2 - 2.0) / (d1 + d2 - 2.0));
    let v = 12.0 / (d2 - 6.0) * (a + b);
    if is_inf(v) || is_nan(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}
