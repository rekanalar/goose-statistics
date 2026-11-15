// Four continuous families that come from the gamma law: chi, dgamma,
// invgamma and loggamma. The names and the parameters follow SciPy.
//
// Every family has eight calls: pdf, cdf, sf, ppf, mean, var, skew and
// kurt. The kurtosis is the excess kurtosis. Every call returns the result
// first and a fault code last. A value outside the support gives a density
// of 0. The quantile of 0 or 1 gives the bound of the support.
//
// The moments of the law of the logarithm of the gamma law use the
// Hurwitz zeta function of the file distributions_growth.goose. The inverse gamma
// law has a moment of order n only when the shape is above n.
namespace statistics;

import std;
import math;

// ---- chi: the chi law ---------------------------------------------------

fn chi_pdf(x: f64, df: f64) -> f64, i64 {
    let k = check_shape(df);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x < 0.0 { return 0.0, OK; }
    if x == 0.0 {
        if df > 1.0 { return 0.0, OK; }
        if df < 1.0 { return nan(), ERR_OVERFLOW; }
        // The exponent of x is zero here, so the density has a limit.
        let lg0, sg0, cl0 = lgamma(df / 2.0);
        if cl0 != OK { return nan(), cl0; }
        return exp(-(df / 2.0 - 1.0) * log(2.0) - lg0), OK;
    }
    let lg, sg, cl = lgamma(df / 2.0);
    if cl != OK { return nan(), cl; }
    let l = (df - 1.0) * log(x) - x * x / 2.0 - (df / 2.0 - 1.0) * log(2.0) - lg;
    return exp(l), OK;
}

fn chi_cdf(x: f64, df: f64) -> f64, i64 {
    let k = check_shape(df);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    let v, cv = gammainc(df / 2.0, x * x / 2.0);
    if cv != OK { return nan(), cv; }
    return v, OK;
}

fn chi_sf(x: f64, df: f64) -> f64, i64 {
    let k = check_shape(df);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 1.0, OK; }
    let v, cv = gammaincc(df / 2.0, x * x / 2.0);
    if cv != OK { return nan(), cv; }
    return v, OK;
}

fn chi_ppf(q: f64, df: f64) -> f64, i64 {
    let k = check_shape(df);
    if k != OK { return nan(), k; }
    let pc = check_prob(q);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(q, 0.0, inf());
    if e { return ev, ec; }
    let v, cv = gammaincinv(df / 2.0, q);
    if cv != OK { return nan(), cv; }
    return sqrt(2.0 * v), OK;
}

// The mean of the chi law: the square root of 2 times the ratio of the
// gamma function of (df + 1) / 2 and of df / 2.
fn chi_mu(df: f64) -> f64, i64 {
    let a, sa, ca = lgamma(df / 2.0 + 0.5);
    if ca != OK { return nan(), ca; }
    let b, sb, cb = lgamma(df / 2.0);
    if cb != OK { return nan(), cb; }
    return SQRT2 * exp(a - b), OK;
}

fn chi_mean(df: f64) -> f64, i64 {
    let k = check_shape(df);
    if k != OK { return nan(), k; }
    return chi_mu(df);
}

fn chi_var(df: f64) -> f64, i64 {
    let k = check_shape(df);
    if k != OK { return nan(), k; }
    let mu, c = chi_mu(df);
    if c != OK { return nan(), c; }
    return df - mu * mu, OK;
}

fn chi_skew(df: f64) -> f64, i64 {
    let k = check_shape(df);
    if k != OK { return nan(), k; }
    let mu, c = chi_mu(df);
    if c != OK { return nan(), c; }
    let v = df - mu * mu;
    return (2.0 * mu * mu * mu + mu * (1.0 - 2.0 * df)) / (v * sqrt(v)), OK;
}

fn chi_kurt(df: f64) -> f64, i64 {
    let k = check_shape(df);
    if k != OK { return nan(), k; }
    let mu, c = chi_mu(df);
    if c != OK { return nan(), c; }
    let v = df - mu * mu;
    let m2 = mu * mu;
    let g = 2.0 * df * (1.0 - df) - 6.0 * m2 * m2 + 4.0 * m2 * (2.0 * df - 1.0);
    return g / (v * v), OK;
}

// ---- dgamma: the double gamma law ---------------------------------------

fn dgamma_pdf(x: f64, a: f64) -> f64, i64 {
    let k = check_shape(a);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    let ax = abs(x);
    if ax == 0.0 {
        if a > 1.0 { return 0.0, OK; }
        if a < 1.0 { return nan(), ERR_OVERFLOW; }
        // The exponent of x is zero here, so the density has a limit.
        let lg0, sg0, cl0 = lgamma(a);
        if cl0 != OK { return nan(), cl0; }
        return exp(-lg0) / 2.0, OK;
    }
    let lg, sg, cl = lgamma(a);
    if cl != OK { return nan(), cl; }
    return exp((a - 1.0) * log(ax) - ax - lg) / 2.0, OK;
}

fn dgamma_cdf(x: f64, a: f64) -> f64, i64 {
    let k = check_shape(a);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x > 0.0 {
        let v, cv = gammainc(a, x);
        if cv != OK { return nan(), cv; }
        return 0.5 + 0.5 * v, OK;
    }
    let v, cv = gammaincc(a, -x);
    if cv != OK { return nan(), cv; }
    return 0.5 * v, OK;
}

fn dgamma_sf(x: f64, a: f64) -> f64, i64 {
    let v, k = dgamma_cdf(-x, a);
    if k != OK { return nan(), k; }
    return v, OK;
}

fn dgamma_ppf(q: f64, a: f64) -> f64, i64 {
    let k = check_shape(a);
    if k != OK { return nan(), k; }
    let pc = check_prob(q);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(q, -inf(), inf());
    if e { return ev, ec; }
    if q > 0.5 {
        let v, cv = gammaincinv(a, 2.0 * q - 1.0);
        if cv != OK { return nan(), cv; }
        return v, OK;
    }
    let v, cv = gammainccinv(a, 2.0 * q);
    if cv != OK { return nan(), cv; }
    return -v, OK;
}

fn dgamma_mean(a: f64) -> f64, i64 {
    let k = check_shape(a);
    if k != OK { return nan(), k; }
    // The law is symmetric around 0, so its mean is 0.
    return 0.0, OK;
}

fn dgamma_var(a: f64) -> f64, i64 {
    let k = check_shape(a);
    if k != OK { return nan(), k; }
    return a * (a + 1.0), OK;
}

fn dgamma_skew(a: f64) -> f64, i64 {
    let k = check_shape(a);
    if k != OK { return nan(), k; }
    // The law is symmetric around 0, so its skewness is 0.
    return 0.0, OK;
}

fn dgamma_kurt(a: f64) -> f64, i64 {
    let k = check_shape(a);
    if k != OK { return nan(), k; }
    return (a + 2.0) * (a + 3.0) / (a * (a + 1.0)) - 3.0, OK;
}

// ---- invgamma: the inverse gamma law ------------------------------------

fn invgamma_pdf(x: f64, a: f64) -> f64, i64 {
    let k = check_shape(a);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    let lg, sg, cl = lgamma(a);
    if cl != OK { return nan(), cl; }
    return exp(-(a + 1.0) * log(x) - lg - 1.0 / x), OK;
}

fn invgamma_cdf(x: f64, a: f64) -> f64, i64 {
    let k = check_shape(a);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    let v, cv = gammaincc(a, 1.0 / x);
    if cv != OK { return nan(), cv; }
    return v, OK;
}

fn invgamma_sf(x: f64, a: f64) -> f64, i64 {
    let k = check_shape(a);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 1.0, OK; }
    let v, cv = gammainc(a, 1.0 / x);
    if cv != OK { return nan(), cv; }
    return v, OK;
}

fn invgamma_ppf(q: f64, a: f64) -> f64, i64 {
    let k = check_shape(a);
    if k != OK { return nan(), k; }
    let pc = check_prob(q);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(q, 0.0, inf());
    if e { return ev, ec; }
    let v, cv = gammainccinv(a, q);
    if cv != OK { return nan(), cv; }
    return 1.0 / v, OK;
}

fn invgamma_mean(a: f64) -> f64, i64 {
    let k = check_shape(a);
    if k != OK { return nan(), k; }
    if a <= 1.0 { return nan(), ERR_UNDEFINED; }
    return 1.0 / (a - 1.0), OK;
}

fn invgamma_var(a: f64) -> f64, i64 {
    let k = check_shape(a);
    if k != OK { return nan(), k; }
    if a <= 2.0 { return nan(), ERR_UNDEFINED; }
    return 1.0 / ((a - 1.0) * (a - 1.0) * (a - 2.0)), OK;
}

fn invgamma_skew(a: f64) -> f64, i64 {
    let k = check_shape(a);
    if k != OK { return nan(), k; }
    if a <= 3.0 { return nan(), ERR_UNDEFINED; }
    return 4.0 * sqrt(a - 2.0) / (a - 3.0), OK;
}

fn invgamma_kurt(a: f64) -> f64, i64 {
    let k = check_shape(a);
    if k != OK { return nan(), k; }
    if a <= 4.0 { return nan(), ERR_UNDEFINED; }
    return 6.0 * (5.0 * a - 11.0) / ((a - 3.0) * (a - 4.0)), OK;
}

// ---- loggamma: the law of the logarithm of the gamma law ------------------

fn loggamma_pdf(x: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    let lg, sg, cl = lgamma(c);
    if cl != OK { return nan(), cl; }
    return exp(c * x - exp(x) - lg), OK;
}

// The law holds a tail for a very low point: the cumulative function is
// the first term of the series of the incomplete gamma function there.
fn loggamma_cdf(x: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x < -700.0 {
        let lg, sg, cl = lgamma(c + 1.0);
        if cl != OK { return nan(), cl; }
        return exp(c * x - lg), OK;
    }
    let v, cv = gammainc(c, exp(x));
    if cv != OK { return nan(), cv; }
    return v, OK;
}

fn loggamma_sf(x: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x < -700.0 {
        let lg, sg, cl = lgamma(c + 1.0);
        if cl != OK { return nan(), cl; }
        return -expm1(c * x - lg), OK;
    }
    let v, cv = gammaincc(c, exp(x));
    if cv != OK { return nan(), cv; }
    return v, OK;
}

fn loggamma_ppf(q: f64, c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let pc = check_prob(q);
    if pc != OK { return nan(), pc; }
    let e, ev, ec = fam_end(q, -inf(), inf());
    if e { return ev, ec; }
    let g, cg = gammaincinv(c, q);
    if cg != OK { return nan(), cg; }
    if g < 1.0e-300 {
        // The inverse of the first term of the series of the incomplete
        // gamma function.
        let lg, sg, cl = lgamma(c + 1.0);
        if cl != OK { return nan(), cl; }
        return (log(q) + lg) / c, OK;
    }
    return log(g), OK;
}

fn loggamma_mean(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let p, cp = digamma(c);
    if cp != OK { return nan(), cp; }
    return p, OK;
}

// The variance is the trigamma function. The Hurwitz zeta function of the
// order 2 gives it.
fn loggamma_var(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    return fam_finite(genlogistic_zeta(2, c));
}

fn loggamma_skew(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let v = genlogistic_zeta(2, c);
    return fam_finite(-2.0 * genlogistic_zeta(3, c) / (v * sqrt(v)));
}

fn loggamma_kurt(c: f64) -> f64, i64 {
    let k = check_shape(c);
    if k != OK { return nan(), k; }
    let v = genlogistic_zeta(2, c);
    return fam_finite(6.0 * genlogistic_zeta(4, c) / (v * v));
}
