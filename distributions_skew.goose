// Two continuous families of a skew shape: the asymmetric Laplace law and
// the skew normal law.
//
// The names and the parameters follow SciPy.
//
// Every family has eight calls: pdf, cdf, sf, ppf, mean, var, skew and
// kurt. The kurtosis is the excess kurtosis. Every call returns the result
// first and a fault code last. A value outside the support gives a density
// of 0. The quantile of 0 or 1 gives the bound of the support.
namespace statistics;

import std;
import math;

// The cumulative function of the normal law.
fn fam_ndtr(z: f64) -> f64 {
    let e, c = erfc(-z * INV_SQRT2);
    if c != OK { return nan(); }
    return 0.5 * e;
}

// The T function of Owen. It gives the probability that a normal variable
// is above h while another one lies between 0 and a. The identity of Owen
// turns a large a into a small one, so the integral stays short.
fn fam_owens_t(h: f64, a: f64) -> f64 {
    if a == 0.0 { return 0.0; }
    var sign = 1.0;
    var aa = a;
    if aa < 0.0 { aa = -aa; sign = -1.0; }
    if abs(h) < 1e-12 { return sign * atan(aa) / (2.0 * PI); }
    if aa <= 1.0 { return sign * fam_owens_integral(h, aa); }
    // The identity of Owen turns a large a into a small one. Its right
    // side is the sum of the two T values, so the form with the upper
    // tails keeps the digits. The other form would take the difference of
    // two values that both stand at 1 in the far tail.
    var extra: f64;
    if h > 0.0 {
        let ph = fam_ndtr(-h);
        let pa = fam_ndtr(-aa * h);
        extra = 0.5 * (ph + pa) - ph * pa;
    } else {
        let ph = fam_ndtr(h);
        let pa = fam_ndtr(aa * h);
        extra = 0.5 * (ph + pa) - ph * pa;
    }
    return sign * (extra - fam_owens_integral(aa * h, 1.0 / aa));
}

// The integral of the T function of Owen over the range from 0 to a. The
// variable s = abs(h) * t keeps the peak of the integrand at the same
// width, whatever the height h is.
fn fam_owens_integral(h: f64, a: f64) -> f64 {
    var gx: f64[>..] = [];
    var gw: f64[>..] = [];
    gauss_legendre(gx, gw, FAM_NODES);
    let panels: i64 = 32;
    let ah = abs(h);
    var total: f64 = 0.0;
    if ah >= 1.0 {
        let top = min(ah * a, 12.0);
        let width = top / panels as f64;
        for k in 0..panels {
            let left = width * k as f64;
            for j in 0..FAM_NODES {
                let s = left + width * (gx[j] + 1.0) / 2.0;
                let t = s / ah;
                let d = 1.0 + t * t;
                total += exp(-0.5 * (h * h + s * s)) / d * gw[j] * width / 2.0;
            }
        }
        return total / (ah * 2.0 * PI);
    }
    let width = a / panels as f64;
    for k in 0..panels {
        let left = width * k as f64;
        for j in 0..FAM_NODES {
            let t = left + width * (gx[j] + 1.0) / 2.0;
            let d = 1.0 + t * t;
            total += exp(-0.5 * h * h * d) / d * gw[j] * width / 2.0;
        }
    }
    return total / (2.0 * PI);
}

// ---- laplace_asymmetric: the asymmetric Laplace law ------------------------

fn laplace_asymmetric_pdf(x: f64, kappa: f64) -> f64, i64 {
    let c = check_shape(kappa);
    if c != OK { return nan(), c; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    let ki = 1.0 / kappa;
    var l = x * ki;
    if x >= 0.0 { l = -x * kappa; }
    return exp(l - log(kappa + ki)), OK;
}

fn laplace_asymmetric_cdf(x: f64, kappa: f64) -> f64, i64 {
    let c = check_shape(kappa);
    if c != OK { return nan(), c; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    let ki = 1.0 / kappa;
    let tot = kappa + ki;
    if x >= 0.0 { return 1.0 - exp(-x * kappa) * ki / tot, OK; }
    return exp(x * ki) * kappa / tot, OK;
}

fn laplace_asymmetric_sf(x: f64, kappa: f64) -> f64, i64 {
    let c = check_shape(kappa);
    if c != OK { return nan(), c; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    let ki = 1.0 / kappa;
    let tot = kappa + ki;
    if x >= 0.0 { return exp(-x * kappa) * ki / tot, OK; }
    return 1.0 - exp(x * ki) * kappa / tot, OK;
}

fn laplace_asymmetric_ppf(q: f64, kappa: f64) -> f64, i64 {
    let c = check_shape(kappa);
    if c != OK { return nan(), c; }
    let pc = check_prob(q);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(q, -inf(), inf());
    if e { return ev, ec; }
    let ki = 1.0 / kappa;
    let tot = kappa + ki;
    if q >= kappa / tot {
        return -log((1.0 - q) * tot * kappa) * ki, OK;
    }
    return log(q * tot / kappa) * kappa, OK;
}

fn laplace_asymmetric_mean(kappa: f64) -> f64, i64 {
    let c = check_shape(kappa);
    if c != OK { return nan(), c; }
    return 1.0 / kappa - kappa, OK;
}

fn laplace_asymmetric_var(kappa: f64) -> f64, i64 {
    let c = check_shape(kappa);
    if c != OK { return nan(), c; }
    let ki = 1.0 / kappa;
    return ki * ki + kappa * kappa, OK;
}

fn laplace_asymmetric_skew(kappa: f64) -> f64, i64 {
    let c = check_shape(kappa);
    if c != OK { return nan(), c; }
    let k2 = kappa * kappa;
    let k4 = k2 * k2;
    let k6 = k4 * k2;
    let k8 = k4 * k4;
    return 2.0 * (1.0 - k6) / pow(1.0 + k4, 1.5), OK;
}

fn laplace_asymmetric_kurt(kappa: f64) -> f64, i64 {
    let c = check_shape(kappa);
    if c != OK { return nan(), c; }
    let k2 = kappa * kappa;
    let k4 = k2 * k2;
    let k8 = k4 * k4;
    return 6.0 * (1.0 + k8) / ((1.0 + k4) * (1.0 + k4)), OK;
}

// ---- skewnorm: the skew normal law ----------------------------------------

// Above this product of the shape and the distance the tail of the skew
// normal law is thin. There the difference of the parts of the law loses
// every digit, so the tail gets its own integral.
const FAM_SKEW_TAIL: f64 = 4.0;

// The thin tail of the skew normal law: the mass of the law below a
// negative point, when the shape is positive. The integral runs over the
// complementary error function, which stays positive, so the tail keeps
// every digit.
fn skewnorm_left_tail(w: f64, a: f64) -> f64 {
    let ww = -w;
    let b = 1.0 + a * a;
    var gx: f64[>..] = [];
    var gw: f64[>..] = [];
    gauss_legendre(gx, gw, FAM_NODES);
    let panels: i64 = 32;
    let width: f64 = 10.0 / panels as f64;
    var total: f64 = 0.0;
    for k in 0..panels {
        let left = width * k as f64;
        for j in 0..FAM_NODES {
            let u = left + width * (gx[j] + 1.0) / 2.0;
            let z = u * u / b;
            let s = sqrt(ww * ww + 2.0 * z);
            let e, ce = erfc(a * s * INV_SQRT2);
            if ce != OK { return nan(); }
            total += 0.5 * e * exp(-z) / s * (2.0 * u / b) * gw[j] * width / 2.0;
        }
    }
    return 2.0 / sqrt(2.0 * PI) * exp(-ww * ww / 2.0) * total;
}

fn skewnorm_pdf(x: f64, a: f64) -> f64, i64 {
    let c = check_scalar(x);
    if c != OK { return nan(), c; }
    let k = check_scalar(a);
    if k != OK { return nan(), k; }
    let p, cp = norm_pdf(x);
    if cp != OK { return nan(), cp; }
    return 2.0 * p * fam_ndtr(a * x), OK;
}

fn skewnorm_cdf(x: f64, a: f64) -> f64, i64 {
    let c = check_scalar(x);
    if c != OK { return nan(), c; }
    let k = check_scalar(a);
    if k != OK { return nan(), k; }
    if a == 0.0 { return norm_cdf(x); }
    if a > 0.0 && x < 0.0 && a * (-x) > FAM_SKEW_TAIL {
        return skewnorm_left_tail(x, a), OK;
    }
    let v = fam_ndtr(x) - 2.0 * fam_owens_t(x, a);
    if is_nan(v) { return nan(), ERR_NAN; }
    return min(1.0, max(0.0, v)), OK;
}

// The cumulative function of the skew normal law as one number. The search
// for the quantile needs it.
fn skewnorm_cdf_value(x: f64, a: f64) -> f64 {
    let v, c = skewnorm_cdf(x, a);
    if c != OK { return nan(); }
    return v;
}

// The survival function of the skew normal law as one number. The search
// for the quantile needs it.
fn skewnorm_sf_value(x: f64, a: f64) -> f64 {
    let v, c = skewnorm_sf(x, a);
    if c != OK { return nan(); }
    return v;
}

fn skewnorm_sf(x: f64, a: f64) -> f64, i64 {
    let c = check_scalar(x);
    if c != OK { return nan(), c; }
    let k = check_scalar(a);
    if k != OK { return nan(), k; }
    if a == 0.0 { return norm_sf(x); }
    if a < 0.0 && x > 0.0 && (-a) * x > FAM_SKEW_TAIL {
        return skewnorm_left_tail(-x, -a), OK;
    }
    let v = fam_ndtr(-x) + 2.0 * fam_owens_t(x, a);
    if is_nan(v) { return nan(), ERR_NAN; }
    return min(1.0, max(0.0, v)), OK;
}

fn skewnorm_ppf(q: f64, a: f64) -> f64, i64 {
    let pc = check_prob(q);
    if pc != OK { return nan(), pc; }
    let k = check_scalar(a);
    if k != OK { return nan(), k; }
    let e, ev, ec = fam_end(q, -inf(), inf());
    if e { return ev, ec; }
    // A target above one half is searched on the survival function, so the
    // small value of that function keeps its digits.
    if q > 0.5 { return skewnorm_root_sf(a, 1.0 - q); }
    return skewnorm_root_cdf(a, q);
}

// A root of the cumulative function of the skew normal law at the target.
// The function rises from 0 to 1, so the search walks out from the origin
// in both directions until it passes the target. The walk gives a point
// below the target and a point above it.
fn skewnorm_root_cdf(a: f64, q: f64) -> f64, i64 {
    var hi: f64 = 0.0;
    var d: f64 = 1.0;
    var iter: i64 = 0;
    while skewnorm_cdf_value(hi, a) < q && iter < 200 {
        hi = hi + d;
        d = d * 2.0;
        iter += 1;
    }
    var lo: f64 = 0.0;
    d = 1.0;
    iter = 0;
    while skewnorm_cdf_value(lo, a) > q && iter < 200 {
        lo = lo - d;
        d = d * 2.0;
        iter += 1;
    }
    let x, cx = brent_root(lo, hi, 1e-14, 4e-16) { v => skewnorm_cdf_value(v, a) - q };
    if cx != OK { return nan(), cx; }
    return x, OK;
}

// A root of the survival function of the skew normal law at the target.
// The reflection of the law turns the survival function into a cumulative
// function, so the search runs on the shape with the other sign.
fn skewnorm_root_sf(a: f64, q: f64) -> f64, i64 {
    let x, cx = skewnorm_root_cdf(-a, q);
    if cx != OK { return nan(), cx; }
    return -x, OK;
}

// The mean of the skew normal law. It is the only moment that the shape
// changes.
fn skewnorm_mu(a: f64) -> f64 {
    return sqrt(2.0 / PI) * a / sqrt(1.0 + a * a);
}

fn skewnorm_mean(a: f64) -> f64, i64 {
    let c = check_scalar(a);
    if c != OK { return nan(), c; }
    return skewnorm_mu(a), OK;
}

fn skewnorm_var(a: f64) -> f64, i64 {
    let c = check_scalar(a);
    if c != OK { return nan(), c; }
    let m = skewnorm_mu(a);
    return 1.0 - m * m, OK;
}

fn skewnorm_skew(a: f64) -> f64, i64 {
    let c = check_scalar(a);
    if c != OK { return nan(), c; }
    let m = skewnorm_mu(a);
    let v = 1.0 - m * m;
    if v <= 0.0 { return nan(), ERR_UNDEFINED; }
    return (4.0 - PI) / 2.0 * pow(m / sqrt(v), 3.0), OK;
}

fn skewnorm_kurt(a: f64) -> f64, i64 {
    let c = check_scalar(a);
    if c != OK { return nan(), c; }
    let m = skewnorm_mu(a);
    let v = 1.0 - m * m;
    if v <= 0.0 { return nan(), ERR_UNDEFINED; }
    let m2 = m * m;
    return 2.0 * (PI - 3.0) * m2 * m2 / (v * v), OK;
}
