// Three continuous families of the Johnson kind and of the skew t law:
// jf_skew_t, johnsonsb and johnsonsu. The names and the parameters follow
// SciPy.
//
// Every family has eight calls: pdf, cdf, sf, ppf, mean, var, skew and
// kurt. The kurtosis is the excess kurtosis. Every call returns the result
// first and a fault code last. A value outside the support gives a density
// of 0. The quantile of 0 or 1 gives the bound of the support.
//
// The law of the skew t law of the Johnson kind has a moment of order n
// when the shapes a and b are both above n over 2.
namespace statistics;

import std;
import math;

// ---- jf_skew_t: the skew t law of the Johnson kind ------------------------

fn jf_skew_t_pdf(x: f64, a: f64, b: f64) -> f64, i64 {
    let k = jf_skew_t_check(a, b);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    let r = sqrt(a + b + x * x);
    let d1 = pow(1.0 + x / r, a + 0.5);
    let d2 = pow(1.0 - x / r, b + 0.5);
    let bt, cb = beta(a, b);
    if cb != OK { return nan(), cb; }
    let c = pow(2.0, a + b - 1.0) * bt * sqrt(a + b);
    return d1 * d2 / c, OK;
}

fn jf_skew_t_cdf(x: f64, a: f64, b: f64) -> f64, i64 {
    let k = jf_skew_t_check(a, b);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    let r = sqrt(a + b + x * x);
    let y = (1.0 + x / r) / 2.0;
    let v, cv = betainc(a, b, y);
    if cv != OK { return nan(), cv; }
    return v, OK;
}

fn jf_skew_t_sf(x: f64, a: f64, b: f64) -> f64, i64 {
    let k = jf_skew_t_check(a, b);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    let r = sqrt(a + b + x * x);
    let y = (1.0 + x / r) / 2.0;
    let v, cv = betainc_upper(a, b, y);
    if cv != OK { return nan(), cv; }
    return v, OK;
}

fn jf_skew_t_ppf(q: f64, a: f64, b: f64) -> f64, i64 {
    let k = jf_skew_t_check(a, b);
    if k != OK { return nan(), k; }
    let pc = check_prob(q);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(q, -inf(), inf());
    if e { return ev, ec; }
    let d1, cd1 = betaincinv(a, b, q);
    if cd1 != OK { return nan(), cd1; }
    let d2 = (2.0 * d1 - 1.0) * sqrt(a + b);
    let d3 = 2.0 * sqrt(d1 * (1.0 - d1));
    return d2 / d3, OK;
}

// The check of the shapes of the skew t law of the Johnson kind. Both
// shapes must be above 0.
fn jf_skew_t_check(a: f64, b: f64) -> i64 {
    let ka = check_shape(a);
    if ka != OK { return ka; }
    let kb = check_shape(b);
    if kb != OK { return kb; }
    return OK;
}

// The raw moment of order n of the skew t law of the Johnson kind. It comes
// from a sum of beta functions. The moment exists when both shapes are
// above n over 2.
fn jf_skew_t_raw(n: i64, a: f64, b: f64) -> f64, i64 {
    let half = n as f64 / 2.0;
    if a <= half || b <= half { return nan(), ERR_UNDEFINED; }
    let bt, cb = beta(a, b);
    if cb != OK { return nan(), cb; }
    var sum: f64 = 0.0;
    for i in 0..n + 1 {
        let d, cd = beta(a + half - i, b - half + i);
        if cd != OK { return nan(), cd; }
        // The binomial coefficient of n over i.
        var comb: f64 = 1.0;
        for k in 0..i { comb = comb * (n - k) / (k + 1); }
        let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
        sum += comb * sign * d;
    }
    let num = pow(a + b, half);
    return num / (pow(2.0, n as f64) * bt) * sum, OK;
}

fn jf_skew_t_mean(a: f64, b: f64) -> f64, i64 {
    let k = jf_skew_t_check(a, b);
    if k != OK { return nan(), k; }
    return jf_skew_t_raw(1, a, b);
}

fn jf_skew_t_var(a: f64, b: f64) -> f64, i64 {
    let k = jf_skew_t_check(a, b);
    if k != OK { return nan(), k; }
    let m1, c1 = jf_skew_t_raw(1, a, b);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = jf_skew_t_raw(2, a, b);
    if c2 != OK { return nan(), c2; }
    return m2 - m1 * m1, OK;
}

fn jf_skew_t_skew(a: f64, b: f64) -> f64, i64 {
    let k = jf_skew_t_check(a, b);
    if k != OK { return nan(), k; }
    let m1, c1 = jf_skew_t_raw(1, a, b);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = jf_skew_t_raw(2, a, b);
    if c2 != OK { return nan(), c2; }
    let m3, c3 = jf_skew_t_raw(3, a, b);
    if c3 != OK { return nan(), c3; }
    return fam_skew_raw(m1, m2, m3), OK;
}

fn jf_skew_t_kurt(a: f64, b: f64) -> f64, i64 {
    let k = jf_skew_t_check(a, b);
    if k != OK { return nan(), k; }
    let m1, c1 = jf_skew_t_raw(1, a, b);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = jf_skew_t_raw(2, a, b);
    if c2 != OK { return nan(), c2; }
    let m3, c3 = jf_skew_t_raw(3, a, b);
    if c3 != OK { return nan(), c3; }
    let m4, c4 = jf_skew_t_raw(4, a, b);
    if c4 != OK { return nan(), c4; }
    return fam_kurt_raw(m1, m2, m3, m4), OK;
}

// ---- johnsonsb: the Johnson law of the B kind -----------------------------

fn johnsonsb_pdf(x: f64, a: f64, b: f64) -> f64, i64 {
    let k = johnson_check(a, b);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 || x >= 1.0 { return 0.0, OK; }
    let l = log(x / (1.0 - x));
    let p, cp = norm_pdf(a + b * l);
    if cp != OK { return nan(), cp; }
    return b * p / (x * (1.0 - x)), OK;
}

fn johnsonsb_cdf(x: f64, a: f64, b: f64) -> f64, i64 {
    let k = johnson_check(a, b);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    if x >= 1.0 { return 1.0, OK; }
    let l = log(x / (1.0 - x));
    let v, cv = norm_cdf(a + b * l);
    if cv != OK { return nan(), cv; }
    return v, OK;
}

fn johnsonsb_sf(x: f64, a: f64, b: f64) -> f64, i64 {
    let k = johnson_check(a, b);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 1.0, OK; }
    if x >= 1.0 { return 0.0, OK; }
    let l = log(x / (1.0 - x));
    let v, cv = norm_sf(a + b * l);
    if cv != OK { return nan(), cv; }
    return v, OK;
}

fn johnsonsb_ppf(q: f64, a: f64, b: f64) -> f64, i64 {
    let k = johnson_check(a, b);
    if k != OK { return nan(), k; }
    let pc = check_prob(q);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(q, 0.0, 1.0);
    if e { return ev, ec; }
    let z = norm_std_ppf(q);
    let t = (z - a) / b;
    // The logistic function keeps its digits on both sides of 0.
    if t >= 0.0 {
        let e = exp(-t);
        return 1.0 / (1.0 + e), OK;
    }
    let e = exp(t);
    return e / (1.0 + e), OK;
}

// The check of the shapes of the Johnson laws. The shape b must be above 0.
// The shape a may be any finite value.
fn johnson_check(a: f64, b: f64) -> i64 {
    let ka = check_scalar(a);
    if ka != OK { return ka; }
    let kb = check_shape(b);
    if kb != OK { return kb; }
    return OK;
}

// The raw moments up to the order count of the Johnson law of the B
// kind. They come from an integration of the density of the law. One
// pass over the grid serves every order: the density at a node comes
// once, and the powers of the node come from multiplications.
fn johnsonsb_raws(a: f64, b: f64, count: i64) -> f64, f64, f64, f64, i64 {
    var gx: f64[>..] = [];
    var gw: f64[>..] = [];
    gauss_legendre(gx, gw, FAM_NODES);
    var m1: f64 = 0.0;
    var m2: f64 = 0.0;
    var m3: f64 = 0.0;
    var m4: f64 = 0.0;
    let panels: i64 = 40;
    let width = 1.0 / panels as f64;
    for k in 0..panels {
        let left = width * k;
        for j in 0..FAM_NODES {
            let x = left + width * (gx[j] + 1.0) / 2.0;
            if x <= 0.0 || x >= 1.0 { continue; }
            let l = log(x / (1.0 - x));
            let p, cp = norm_pdf(a + b * l);
            if cp != OK { return nan(), nan(), nan(), nan(), cp; }
            let d = b * p / (x * (1.0 - x));
            let wj = d * gw[j] * width / 2.0;
            let x2 = x * x;
            if count >= 1 { m1 += wj * x; }
            if count >= 2 { m2 += wj * x2; }
            if count >= 3 { m3 += wj * x2 * x; }
            if count >= 4 { m4 += wj * x2 * x2; }
        }
    }
    if count >= 1 && (is_nan(m1) || is_inf(m1)) { return nan(), nan(), nan(), nan(), ERR_OVERFLOW; }
    if count >= 2 && (is_nan(m2) || is_inf(m2)) { return nan(), nan(), nan(), nan(), ERR_OVERFLOW; }
    if count >= 3 && (is_nan(m3) || is_inf(m3)) { return nan(), nan(), nan(), nan(), ERR_OVERFLOW; }
    if count >= 4 && (is_nan(m4) || is_inf(m4)) { return nan(), nan(), nan(), nan(), ERR_OVERFLOW; }
    return m1, m2, m3, m4, OK;
}

// The raw moment of order n of the Johnson law of the B kind. It comes
// from an integration of the density of the law.
fn johnsonsb_raw(n: i64, a: f64, b: f64) -> f64, i64 {
    let m1, m2, m3, m4, ck = johnsonsb_raws(a, b, n);
    if ck != OK { return nan(), ck; }
    if n == 1 { return m1, OK; }
    if n == 2 { return m2, OK; }
    if n == 3 { return m3, OK; }
    return m4, OK;
}

fn johnsonsb_mean(a: f64, b: f64) -> f64, i64 {
    let k = johnson_check(a, b);
    if k != OK { return nan(), k; }
    let m1, m2, m3, m4, ck = johnsonsb_raws(a, b, 1);
    if ck != OK { return nan(), ck; }
    return m1, OK;
}

fn johnsonsb_var(a: f64, b: f64) -> f64, i64 {
    let k = johnson_check(a, b);
    if k != OK { return nan(), k; }
    let m1, m2, m3, m4, ck = johnsonsb_raws(a, b, 2);
    if ck != OK { return nan(), ck; }
    return m2 - m1 * m1, OK;
}

fn johnsonsb_skew(a: f64, b: f64) -> f64, i64 {
    let k = johnson_check(a, b);
    if k != OK { return nan(), k; }
    let m1, m2, m3, m4, ck = johnsonsb_raws(a, b, 3);
    if ck != OK { return nan(), ck; }
    return fam_skew_raw(m1, m2, m3), OK;
}

fn johnsonsb_kurt(a: f64, b: f64) -> f64, i64 {
    let k = johnson_check(a, b);
    if k != OK { return nan(), k; }
    let m1, m2, m3, m4, ck = johnsonsb_raws(a, b, 4);
    if ck != OK { return nan(), ck; }
    return fam_kurt_raw(m1, m2, m3, m4), OK;
}

// The hyperbolic sine. The standard library of the older compiler holds
// no call for it. The form below keeps its digits on both sides of the
// line and does not overflow for a large argument.
fn johnson_sinh(x: f64) -> f64 {
    if x > 40.0 { return 0.5 * exp(x); }
    if x < -40.0 { return -0.5 * exp(-x); }
    return 0.5 * (exp(x) - exp(-x));
}

// The inverse of the hyperbolic sine. The standard library holds the
// hyperbolic sine but not its inverse. The form below keeps its digits on
// both sides of the line.
fn johnson_asinh(x: f64) -> f64 {
    if x >= 0.0 {
        return log(x + sqrt(x * x + 1.0));
    }
    return -log(-x + sqrt(x * x + 1.0));
}

// ---- johnsonsu: the Johnson law of the U kind -----------------------------

fn johnsonsu_pdf(x: f64, a: f64, b: f64) -> f64, i64 {
    let k = johnson_check(a, b);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    let h = johnson_asinh(x);
    let p, cp = norm_pdf(a + b * h);
    if cp != OK { return nan(), cp; }
    return b * p / sqrt(x * x + 1.0), OK;
}

fn johnsonsu_cdf(x: f64, a: f64, b: f64) -> f64, i64 {
    let k = johnson_check(a, b);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    let h = johnson_asinh(x);
    let v, cv = norm_cdf(a + b * h);
    if cv != OK { return nan(), cv; }
    return v, OK;
}

fn johnsonsu_sf(x: f64, a: f64, b: f64) -> f64, i64 {
    let k = johnson_check(a, b);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    let h = johnson_asinh(x);
    let v, cv = norm_sf(a + b * h);
    if cv != OK { return nan(), cv; }
    return v, OK;
}

fn johnsonsu_ppf(q: f64, a: f64, b: f64) -> f64, i64 {
    let k = johnson_check(a, b);
    if k != OK { return nan(), k; }
    let pc = check_prob(q);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(q, -inf(), inf());
    if e { return ev, ec; }
    let z = norm_std_ppf(q);
    return johnson_sinh((z - a) / b), OK;
}

// The raw moment of order n of the Johnson law of the U kind. The law is
// the sine of the hyperbolic of a normal variable. Its moments come from the
// binomial sum of the powers of the hyperbolic sine.
fn johnsonsu_raw(n: i64, a: f64, b: f64) -> f64, i64 {
    // The mean and the variance of the normal variable of the law.
    let m = -a / b;
    let s2 = 1.0 / (b * b);
    var total: f64 = 0.0;
    for k in 0..n + 1 {
        // The binomial coefficient of n over k.
        var comb: f64 = 1.0;
        for j in 0..k { comb = comb * (n - j) / (j + 1); }
        let sign = if k % 2 == 0 { 1.0 } else { -1.0 };
        let d = (n - 2 * k) as f64;
        let term = comb * sign * exp(d * m + d * d * s2 / 2.0);
        if is_inf(term) { return nan(), ERR_OVERFLOW; }
        total += term;
    }
    return total * pow(2.0, -(n as f64)), OK;
}

fn johnsonsu_mean(a: f64, b: f64) -> f64, i64 {
    let k = johnson_check(a, b);
    if k != OK { return nan(), k; }
    let m, cm = johnsonsu_raw(1, a, b);
    if cm != OK { return nan(), cm; }
    return m, OK;
}

fn johnsonsu_var(a: f64, b: f64) -> f64, i64 {
    let k = johnson_check(a, b);
    if k != OK { return nan(), k; }
    let m1, c1 = johnsonsu_raw(1, a, b);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = johnsonsu_raw(2, a, b);
    if c2 != OK { return nan(), c2; }
    return m2 - m1 * m1, OK;
}

fn johnsonsu_skew(a: f64, b: f64) -> f64, i64 {
    let k = johnson_check(a, b);
    if k != OK { return nan(), k; }
    let m1, c1 = johnsonsu_raw(1, a, b);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = johnsonsu_raw(2, a, b);
    if c2 != OK { return nan(), c2; }
    let m3, c3 = johnsonsu_raw(3, a, b);
    if c3 != OK { return nan(), c3; }
    return fam_skew_raw(m1, m2, m3), OK;
}

fn johnsonsu_kurt(a: f64, b: f64) -> f64, i64 {
    let k = johnson_check(a, b);
    if k != OK { return nan(), k; }
    let m1, c1 = johnsonsu_raw(1, a, b);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = johnsonsu_raw(2, a, b);
    if c2 != OK { return nan(), c2; }
    let m3, c3 = johnsonsu_raw(3, a, b);
    if c3 != OK { return nan(), c3; }
    let m4, c4 = johnsonsu_raw(4, a, b);
    if c4 != OK { return nan(), c4; }
    return fam_kurt_raw(m1, m2, m3, m4), OK;
}
