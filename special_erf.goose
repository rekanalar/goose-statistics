// The error function family of the package, and the accurate log1p and
// expm1.
//
// Every call returns the result first and a fault code last. A fault
// gives NaN as the result.
namespace statistics;

import std;
import math;

fn erf_series(x: f64) -> f64 {
    return x * horner(ERF_SERIES, x * x);
}

// q(x) = exp(x^2) erfc(x), for 0.5 <= x < 28. One fit per range, and the
// asymptotic series above 7.
fn scaled_q(x: f64) -> f64 {
    if x < 1.5 {
        return clenshaw(H1_CHEB, 2.0 * x - 2.0);
    }
    if x < 3.0 {
        return clenshaw(H2_CHEB, (4.0 * x - 9.0) / 3.0);
    }
    if x < 7.0 {
        return clenshaw(G_CHEB, (2.0 * x - 10.0) / 4.0) / x * INV_SQRT_PI;
    }
    return horner(ASYM_COEF, -0.5 / (x * x)) / x * INV_SQRT_PI;
}

// The erfc core for 0.5 <= x < 28. The split of x*x keeps the digits of
// the exponential, so the far tail holds its relative accuracy.
fn erfc_pos(x: f64) -> f64 {
    let hi, lo = two_prod(x, x);
    var q = scaled_q(x);
    q = q - q * lo;
    return q * exp(-hi);
}

// The error function of one value.
fn erf(x: f64) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    if is_inf(x) {
        if x > 0.0 { return 1.0, OK; }
        return -1.0, OK;
    }
    let ax = abs(x);
    var r: f64 = 0.0;
    if ax <= 0.5 {
        r = erf_series(ax);
    } else if ax < 28.0 {
        r = 1.0 - erfc_pos(ax);
    } else {
        r = 1.0;
    }
    if x < 0.0 { r = -r; }
    return r, OK;
}

// The complementary error function. The call never forms one minus the
// erf above 0.5, where that form loses every digit.
fn erfc(x: f64) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    if is_inf(x) {
        if x > 0.0 { return 0.0, OK; }
        return 2.0, OK;
    }
    if abs(x) <= 0.5 {
        return 1.0 - erf_series(x), OK;
    }
    if x > 27.25 { return 0.0, OK; }
    if x < -27.25 { return 2.0, OK; }
    if x > 0.0 { return erfc_pos(x), OK; }
    return 2.0 - erfc_pos(-x), OK;
}

// The start value of the inverse, from the form of Winitzki. The input
// is log(1 - y*y), never above 0.
fn winitzki(lg: f64) -> f64 {
    let b = 2.0 / (PI * 0.147) + 0.5 * lg;
    let eps = -lg / 0.147;
    var d: f64 = 0.0;
    if b > 0.0 {
        d = eps / (sqrt(b * b + eps) + b);
    } else {
        d = sqrt(b * b + eps) - b;
    }
    return sqrt(d);
}

// The inverse core for a probability p in (0, 0.5]. Above 1e-3 the steps
// refine the erfc. Below 1e-3 the steps refine the log of the erfc, so
// the far tail never forms a subnormal value.
fn erfcinv_core(p: f64) -> f64 {
    var x = winitzki(log(p * (2.0 - p)));
    if p >= 1e-3 {
        var i: i64 = 0;
        while i < 3 {
            let g = erfc(x) - p;
            let d = -TWO_OVER_SQRT_PI * exp(-x * x);
            let u = g / d;
            x = x - u / (1.0 + x * u);
            i += 1;
        }
        return x;
    }
    let lp = log(p);
    var i: i64 = 0;
    while i < 3 {
        let hi, lo = two_prod(x, x);
        let q = scaled_q(x);
        let h = (-hi - lp) + (log(q) - lo);
        let hp = -TWO_OVER_SQRT_PI / q;
        let u = h / hp;
        x = x - u / (1.0 + u * (x + 0.5 * hp));
        i += 1;
    }
    return x;
}

// The inverse erf for a value in (0, 0.5].
fn erfinv_core(y: f64) -> f64 {
    if y < 1e-8 {
        return HALF_SQRT_PI * y * (1.0 + PI_OVER_12 * y * y);
    }
    var x = winitzki(log1p_core(-y * y));
    var i: i64 = 0;
    while i < 3 {
        let g = erf_core(x) - y;
        let d = TWO_OVER_SQRT_PI * exp(-x * x);
        let u = g / d;
        x = x - u / (1.0 + x * u);
        i += 1;
    }
    return x;
}

// The core of the erf without the checks, for the steps above.
fn erf_core(x: f64) -> f64 {
    let ax = abs(x);
    var r: f64 = 0.0;
    if ax <= 0.5 {
        r = erf_series(ax);
    } else if ax < 28.0 {
        r = 1.0 - erfc_pos(ax);
    } else {
        r = 1.0;
    }
    if x < 0.0 { r = -r; }
    return r;
}

// The inverse of the error function. The range is the open interval from
// -1 to 1. The end points give the overflow code.
fn erfinv(y: f64) -> f64, i64 {
    let c = check_scalar(y);
    if c != 0 { return nan(), c; }
    if y > 1.0 || y < -1.0 { return nan(), ERR_DOMAIN; }
    let a = abs(y);
    if a == 1.0 { return nan(), ERR_OVERFLOW; }
    if a == 0.0 { return y, OK; }
    var r: f64 = 0.0;
    if a <= 0.5 {
        r = erfinv_core(a);
    } else {
        r = erfcinv_core(1.0 - a);
    }
    if y < 0.0 { r = -r; }
    return r, OK;
}

// The inverse of the erfc. The range runs from 0 to 2. The end points
// give the overflow code.
fn erfcinv(p: f64) -> f64, i64 {
    let c = check_scalar(p);
    if c != 0 { return nan(), c; }
    if p < 0.0 || p > 2.0 { return nan(), ERR_DOMAIN; }
    if p == 0.0 || p == 2.0 { return nan(), ERR_OVERFLOW; }
    if p <= 0.5 { return erfcinv_core(p), OK; }
    if p < 1.5 {
        let t = 1.0 - p;
        if t > 0.0 { return erfinv_core(t), OK; }
        if t < 0.0 { return -erfinv_core(-t), OK; }
        return 0.0, OK;
    }
    return -erfcinv_core(2.0 - p), OK;
}

// ---- the accurate log1p and expm1 ------------------------------------

// The log of one plus x. The two sum holds 1 + x exactly, so the small x
// keeps its digits.
fn log1p_core(x: f64) -> f64 {
    let u, e = two_sum(1.0, x);
    return log(u) + e / u;
}

fn log1p(x: f64) -> f64, i64 {
    let c = check_scalar(x);
    if c != 0 { return nan(), c; }
    if x < -1.0 { return nan(), ERR_DOMAIN; }
    if x == -1.0 { return nan(), ERR_OVERFLOW; }
    let r = log1p_core(x);
    if is_inf(r) { return nan(), ERR_OVERFLOW; }
    return r, OK;
}

fn expm1_core(x: f64) -> f64 {
    let k = floor(x * INV_LN2 + 0.5);
    let r = (x - k * LN2_HI) - k * LN2_LO;
    let p = r + r * r * horner(EXPM1_SERIES, r);
    if k == 0.0 { return p; }
    let s, e = two_sum(1.0, p);
    let ik = k as i64;
    if k > 0.0 {
        let t = (s * pow2_int(ik - 1)) * 2.0 - 1.0;
        return t + (e * pow2_int(ik - 1)) * 2.0;
    }
    return s * pow2_int(ik) - 1.0 + e * pow2_int(ik);
}

fn expm1(x: f64) -> f64, i64 {
    let c = check_scalar(x);
    if c != 0 { return nan(), c; }
    if x > 709.782712893384 { return nan(), ERR_OVERFLOW; }
    if x < -40.0 { return -1.0, OK; }
    let r = expm1_core(x);
    if is_inf(r) { return nan(), ERR_OVERFLOW; }
    return r, OK;
}
