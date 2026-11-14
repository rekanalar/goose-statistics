// Four continuous families of a power or an exponent: powernorm,
// powerlognorm, exponpow and exponweib. The names and the parameters
// follow SciPy.
//
// Every family has eight calls: pdf, cdf, sf, ppf, mean, var, skew and
// kurt. The kurtosis is the excess kurtosis. Every call returns the result
// first and a fault code last. A value outside the support gives a density
// of 0. The quantile of 0 or 1 gives the bound of the support.
//
// The moments of these families have no closed form. They come from the
// rule of the file distributions_numeric.goose. That rule integrates the
// density of the law. SciPy integrates the quantile function instead, and
// that route loses the accuracy of the third and the fourth raw moment.
// The README names each family and the size of the difference.
namespace statistics;

import std;
import math;

// The density of the power of the normal law, without the checks.
fn powernorm_dens(x: f64, c: f64) -> f64 {
    let pv, pc = norm_cdf(-x);
    let p = pv;
    if p == 0.0 { return 0.0; }
    return c * exp(-0.5 * x * x) / SQRT_2PI * pow(p, c - 1.0);
}

// ---- powernorm: the power of the normal law -----------------------------

fn powernorm_pdf(x: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    return powernorm_dens(x, c), OK;
}

fn powernorm_cdf(x: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    let p, cp = norm_sf(x);
    if cp != OK { return nan(), cp; }
    if p == 0.0 { return 1.0, OK; }
    let m, cm = expm1(c * log(p));
    if cm != OK { return nan(), cm; }
    return -m, OK;
}

fn powernorm_sf(x: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    let p, cp = norm_sf(x);
    if cp != OK { return nan(), cp; }
    return pow(p, c), OK;
}

fn powernorm_ppf(p: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(p, -inf(), inf());
    if e { return ev, ec; }
    let r = -norm_std_ppf(pow(1.0 - p, 1.0 / c));
    return fam_finite(r);
}

fn powernorm_mean(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let m1, m2, m3, m4, ck = fam_raws_line(1) { x => powernorm_dens(x, c) };
    if ck != OK { return nan(), ck; }
    if is_nan(m1) || is_inf(m1) { return nan(), ERR_OVERFLOW; }
    return m1, OK;
}

fn powernorm_var(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let m1, m2, m3, m4, ck = fam_raws_line(2) { x => powernorm_dens(x, c) };
    if ck != OK { return nan(), ck; }
    if is_nan(m2) || is_inf(m2) { return nan(), ERR_OVERFLOW; }
    return fam_var(m1, m2), OK;
}

fn powernorm_skew(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let m1, m2, m3, m4, ck = fam_raws_line(3) { x => powernorm_dens(x, c) };
    if ck != OK { return nan(), ck; }
    if is_nan(m3) || is_inf(m3) { return nan(), ERR_OVERFLOW; }
    return fam_skew(m1, m2, m3), OK;
}

fn powernorm_kurt(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let m1, m2, m3, m4, ck = fam_raws_line(4) { x => powernorm_dens(x, c) };
    if ck != OK { return nan(), ck; }
    if is_nan(m4) || is_inf(m4) { return nan(), ERR_OVERFLOW; }
    return fam_kurt(m1, m2, m3, m4), OK;
}

// The density of the power of the log-normal law, without the checks.
fn powerlognorm_dens(x: f64, c: f64, s: f64) -> f64 {
    if x <= 0.0 { return 0.0; }
    let z = log(x) / s;
    let pv, pc = norm_cdf(-z);
    let p = pv;
    if p == 0.0 { return 0.0; }
    return exp(-0.5 * z * z) / (SQRT_2PI * s * x) * c * pow(p, c - 1.0);
}

// ---- powerlognorm: the power of the log-normal law ----------------------

fn powerlognorm_pdf(x: f64, c: f64, s: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let ks = check_shape(s);
    if ks != OK { return nan(), ks; }
    let v = fam_x(x);
    if v != OK { return nan(), v; }
    if x <= 0.0 { return 0.0, OK; }
    return powerlognorm_dens(x, c, s), OK;
}

fn powerlognorm_cdf(x: f64, c: f64, s: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let ks = check_shape(s);
    if ks != OK { return nan(), ks; }
    let v = fam_x(x);
    if v != OK { return nan(), v; }
    if x <= 0.0 { return 0.0, OK; }
    let p, cp = norm_sf(log(x) / s);
    if cp != OK { return nan(), cp; }
    if p == 0.0 { return 1.0, OK; }
    let m, cm = expm1(c * log(p));
    if cm != OK { return nan(), cm; }
    return -m, OK;
}

fn powerlognorm_sf(x: f64, c: f64, s: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let ks = check_shape(s);
    if ks != OK { return nan(), ks; }
    let v = fam_x(x);
    if v != OK { return nan(), v; }
    if x <= 0.0 { return 1.0, OK; }
    let p, cp = norm_sf(log(x) / s);
    if cp != OK { return nan(), cp; }
    return pow(p, c), OK;
}

fn powerlognorm_ppf(p: f64, c: f64, s: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let ks = check_shape(s);
    if ks != OK { return nan(), ks; }
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(p, 0.0, inf());
    if e { return ev, ec; }
    return exp(-norm_std_ppf(pow(1.0 - p, 1.0 / c)) * s), OK;
}

fn powerlognorm_mean(c: f64, s: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let ks = check_shape(s);
    if ks != OK { return nan(), ks; }
    let m1, m2, m3, m4, ck = fam_raws(1, s) { x => powerlognorm_dens(x, c, s) };
    if ck != OK { return nan(), ck; }
    if is_nan(m1) || is_inf(m1) { return nan(), ERR_OVERFLOW; }
    return m1, OK;
}

fn powerlognorm_var(c: f64, s: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let ks = check_shape(s);
    if ks != OK { return nan(), ks; }
    let m1, m2, m3, m4, ck = fam_raws(2, s) { x => powerlognorm_dens(x, c, s) };
    if ck != OK { return nan(), ck; }
    if is_nan(m2) || is_inf(m2) { return nan(), ERR_OVERFLOW; }
    return fam_var(m1, m2), OK;
}

fn powerlognorm_skew(c: f64, s: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let ks = check_shape(s);
    if ks != OK { return nan(), ks; }
    let m1, m2, m3, m4, ck = fam_raws(3, s) { x => powerlognorm_dens(x, c, s) };
    if ck != OK { return nan(), ck; }
    if is_nan(m3) || is_inf(m3) { return nan(), ERR_OVERFLOW; }
    return fam_skew(m1, m2, m3), OK;
}

fn powerlognorm_kurt(c: f64, s: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let ks = check_shape(s);
    if ks != OK { return nan(), ks; }
    let m1, m2, m3, m4, ck = fam_raws(4, s) { x => powerlognorm_dens(x, c, s) };
    if ck != OK { return nan(), ck; }
    if is_nan(m4) || is_inf(m4) { return nan(), ERR_OVERFLOW; }
    return fam_kurt(m1, m2, m3, m4), OK;
}

// The density of the generalized exponential law, without the checks.
fn exponpow_dens(x: f64, b: f64) -> f64 {
    if x <= 0.0 { return 0.0; }
    let xb = pow(x, b);
    return b * pow(x, b - 1.0) * exp(1.0 + xb - exp(xb));
}

// ---- exponpow: the generalized exponential law --------------------------

fn exponpow_pdf(x: f64, b: f64) -> f64, i64 {
    let k = check_shape(b);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    return exponpow_dens(x, b), OK;
}

fn exponpow_cdf(x: f64, b: f64) -> f64, i64 {
    let k = check_shape(b);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    let xb = pow(x, b);
    // The exponential of x^b grows fast. A value of x^b above 40 gives a
    // cumulative function of 1 to the precision of a double.
    if xb > 40.0 { return 1.0, OK; }
    let m, cm = expm1(xb);
    if cm != OK { return nan(), cm; }
    let n, cn = expm1(-m);
    if cn != OK { return nan(), cn; }
    return -n, OK;
}

fn exponpow_sf(x: f64, b: f64) -> f64, i64 {
    let k = check_shape(b);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 1.0, OK; }
    let xb = pow(x, b);
    // The exponential of x^b grows fast. A value of x^b above 40 gives a
    // survival function of 0 to the precision of a double.
    if xb > 40.0 { return 0.0, OK; }
    let m, cm = expm1(xb);
    if cm != OK { return nan(), cm; }
    return exp(-m), OK;
}

fn exponpow_ppf(p: f64, b: f64) -> f64, i64 {
    let k = check_shape(b);
    if k != OK { return nan(), k; }
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(p, 0.0, inf());
    if e { return ev, ec; }
    let m, cm = log1p(-p);
    if cm != OK { return nan(), cm; }
    let n, cn = log1p(-m);
    if cn != OK { return nan(), cn; }
    let r = pow(n, 1.0 / b);
    return fam_finite(r);
}

fn exponpow_mean(b: f64) -> f64, i64 {
    let k = check_shape(b);
    if k != OK { return nan(), k; }
    let m1, m2, m3, m4, ck = fam_raws(1, 1.0) { x => exponpow_dens(x, b) };
    if ck != OK { return nan(), ck; }
    if is_nan(m1) || is_inf(m1) { return nan(), ERR_OVERFLOW; }
    return m1, OK;
}

fn exponpow_var(b: f64) -> f64, i64 {
    let k = check_shape(b);
    if k != OK { return nan(), k; }
    let m1, m2, m3, m4, ck = fam_raws(2, 1.0) { x => exponpow_dens(x, b) };
    if ck != OK { return nan(), ck; }
    if is_nan(m2) || is_inf(m2) { return nan(), ERR_OVERFLOW; }
    return fam_var(m1, m2), OK;
}

fn exponpow_skew(b: f64) -> f64, i64 {
    let k = check_shape(b);
    if k != OK { return nan(), k; }
    let m1, m2, m3, m4, ck = fam_raws(3, 1.0) { x => exponpow_dens(x, b) };
    if ck != OK { return nan(), ck; }
    if is_nan(m3) || is_inf(m3) { return nan(), ERR_OVERFLOW; }
    return fam_skew(m1, m2, m3), OK;
}

fn exponpow_kurt(b: f64) -> f64, i64 {
    let k = check_shape(b);
    if k != OK { return nan(), k; }
    let m1, m2, m3, m4, ck = fam_raws(4, 1.0) { x => exponpow_dens(x, b) };
    if ck != OK { return nan(), ck; }
    if is_nan(m4) || is_inf(m4) { return nan(), ERR_OVERFLOW; }
    return fam_kurt(m1, m2, m3, m4), OK;
}

// The density of the Weibull law of the exponential kind, without the
// checks.
fn exponweib_dens(x: f64, a: f64, c: f64) -> f64 {
    if x <= 0.0 { return 0.0; }
    let xc = pow(x, c);
    let mv, mc = expm1(-xc);
    let m = -mv;
    return a * c * pow(m, a - 1.0) * exp(-xc) * pow(x, c - 1.0);
}

// ---- exponweib: the Weibull law of the exponential kind -----------------

fn exponweib_pdf(x: f64, a: f64, c: f64) -> f64, i64 {
    let k = check_shape(a);
    if k != OK { return nan(), k; }
    let kc = check_shape(c);
    if kc != OK { return nan(), kc; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    return exponweib_dens(x, a, c), OK;
}

fn exponweib_cdf(x: f64, a: f64, c: f64) -> f64, i64 {
    let k = check_shape(a);
    if k != OK { return nan(), k; }
    let kc = check_shape(c);
    if kc != OK { return nan(), kc; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    let m, cm = expm1(-pow(x, c));
    if cm != OK { return nan(), cm; }
    return pow(-m, a), OK;
}

fn exponweib_sf(x: f64, a: f64, c: f64) -> f64, i64 {
    let k = check_shape(a);
    if k != OK { return nan(), k; }
    let kc = check_shape(c);
    if kc != OK { return nan(), kc; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 1.0, OK; }
    let m, cm = expm1(-pow(x, c));
    if cm != OK { return nan(), cm; }
    let n, cn = expm1(a * log(-m));
    if cn != OK { return nan(), cn; }
    return -n, OK;
}

fn exponweib_ppf(p: f64, a: f64, c: f64) -> f64, i64 {
    let k = check_shape(a);
    if k != OK { return nan(), k; }
    let kc = check_shape(c);
    if kc != OK { return nan(), kc; }
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(p, 0.0, inf());
    if e { return ev, ec; }
    let m, cm = log1p(-pow(p, 1.0 / a));
    if cm != OK { return nan(), cm; }
    let r = pow(-m, 1.0 / c);
    return fam_finite(r);
}

fn exponweib_mean(a: f64, c: f64) -> f64, i64 {
    let k = check_shape(a);
    if k != OK { return nan(), k; }
    let kc = check_shape(c);
    if kc != OK { return nan(), kc; }
    let m1, m2, m3, m4, ck = fam_raws(1, 1.0) { x => exponweib_dens(x, a, c) };
    if ck != OK { return nan(), ck; }
    if is_nan(m1) || is_inf(m1) { return nan(), ERR_OVERFLOW; }
    return m1, OK;
}

fn exponweib_var(a: f64, c: f64) -> f64, i64 {
    let k = check_shape(a);
    if k != OK { return nan(), k; }
    let kc = check_shape(c);
    if kc != OK { return nan(), kc; }
    let m1, m2, m3, m4, ck = fam_raws(2, 1.0) { x => exponweib_dens(x, a, c) };
    if ck != OK { return nan(), ck; }
    if is_nan(m2) || is_inf(m2) { return nan(), ERR_OVERFLOW; }
    return fam_var(m1, m2), OK;
}

fn exponweib_skew(a: f64, c: f64) -> f64, i64 {
    let k = check_shape(a);
    if k != OK { return nan(), k; }
    let kc = check_shape(c);
    if kc != OK { return nan(), kc; }
    let m1, m2, m3, m4, ck = fam_raws(3, 1.0) { x => exponweib_dens(x, a, c) };
    if ck != OK { return nan(), ck; }
    if is_nan(m3) || is_inf(m3) { return nan(), ERR_OVERFLOW; }
    return fam_skew(m1, m2, m3), OK;
}

fn exponweib_kurt(a: f64, c: f64) -> f64, i64 {
    let k = check_shape(a);
    if k != OK { return nan(), k; }
    let kc = check_shape(c);
    if kc != OK { return nan(), kc; }
    let m1, m2, m3, m4, ck = fam_raws(4, 1.0) { x => exponweib_dens(x, a, c) };
    if ck != OK { return nan(), ck; }
    if is_nan(m4) || is_inf(m4) { return nan(), ERR_OVERFLOW; }
    return fam_kurt(m1, m2, m3, m4), OK;
}
