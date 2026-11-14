// The distribution objects of SciPy 1.18.1: Normal, Uniform, Logistic,
// Binomial and Mixture.
//
// A Goose block cannot be stored, so an object is a plain struct that
// holds a kind and two numbers. One set of calls serves every kind. The
// calls follow the method names of SciPy: pdf, cdf, ccdf, icdf, iccdf,
// logpdf, logcdf, logccdf, ilogcdf, ilogccdf, mean, median, mode,
// variance, standard_deviation, skewness, kurtosis, entropy, logentropy,
// moment and support. The Binomial kind has a mass in place of a density.
// The kurtosis is the plain kurtosis, as in SciPy: the normal law has 3.
//
// Every call returns the result first and a fault code last.
namespace statistics;

import std;
import math;

const DS_NORMAL: i64 = 1;
const DS_UNIFORM: i64 = 2;
const DS_LOGISTIC: i64 = 3;
const DS_BINOMIAL: i64 = 4;
const DS_MAX_ORDER: i64 = 64;   // the highest order of a raw moment

struct Dist {
    kind: i64,
    a: f64,     // mu, the low end, 0 or n
    b: f64,     // sigma, the high end, 1 or p
}

// ---- the constructors -----------------------------------------------------

fn Normal(mu: f64 = 0.0, sigma: f64 = 1.0) -> Dist, i64 {
    let c = first_code(check_scalar(mu), check_scalar(sigma));
    if c != OK { return Dist { kind: 0, a: 0.0, b: 0.0 }, c; }
    if sigma <= 0.0 { return Dist { kind: 0, a: 0.0, b: 0.0 }, ERR_DOMAIN; }
    return Dist { kind: DS_NORMAL, a: mu, b: sigma }, OK;
}

fn Uniform(a: f64 = 0.0, b: f64 = 1.0) -> Dist, i64 {
    let w, c = unif_check(a, b);
    if c != OK { return Dist { kind: 0, a: 0.0, b: 0.0 }, c; }
    return Dist { kind: DS_UNIFORM, a: a, b: b }, OK;
}

fn Logistic() -> Dist, i64 {
    return Dist { kind: DS_LOGISTIC, a: 0.0, b: 1.0 }, OK;
}

fn Binomial(n: f64, p: f64) -> Dist, i64 {
    let c = first_code(check_scalar(n), check_scalar(p));
    if c != OK { return Dist { kind: 0, a: 0.0, b: 0.0 }, c; }
    if n < 0.0 || n != floor(n) || p < 0.0 || p > 1.0 {
        return Dist { kind: 0, a: 0.0, b: 0.0 }, ERR_DOMAIN;
    }
    return Dist { kind: DS_BINOMIAL, a: n, b: p }, OK;
}

// ---- the support ----------------------------------------------------------

fn dist_support(d: Dist) -> f64, f64, i64 {
    if d.kind == DS_NORMAL || d.kind == DS_LOGISTIC { return -inf(), inf(), OK; }
    if d.kind == DS_UNIFORM { return d.a, d.b, OK; }
    if d.kind == DS_BINOMIAL { return 0.0, d.a, OK; }
    return nan(), nan(), ERR_DOMAIN;
}

// ---- the density, the mass and the cumulative functions -------------------

fn dist_pdf(d: Dist, x: f64) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    if d.kind == DS_NORMAL { return norm_pdf(x, d.a, d.b); }
    if d.kind == DS_UNIFORM { return unif_pdf(x, d.a, d.b); }
    if d.kind == DS_LOGISTIC { return logistic_pdf(x); }
    if d.kind == DS_BINOMIAL {
        // A point under the support of the law has no mass.
        if x < 0.0 || x != floor(x) { return 0.0, OK; }
        return binom_pmf(x, d.a, d.b);
    }
    return nan(), ERR_DOMAIN;
}

fn dist_cdf(d: Dist, x: f64) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    if d.kind == DS_NORMAL { return norm_cdf(x, d.a, d.b); }
    if d.kind == DS_UNIFORM { return unif_cdf(x, d.a, d.b); }
    if d.kind == DS_LOGISTIC { return logistic_cdf(x); }
    if d.kind == DS_BINOMIAL {
        // A point under the support of the law has no mass.
        if x < 0.0 { return 0.0, OK; }
        // The cumulative function of SciPy 1.18.1 is the continuous
        // extension of the law: the regularized incomplete beta function
        // I(1 - p) of the pair (n - x, x + 1).
        if x >= d.a { return 1.0, OK; }
        let v, c = betainc(d.a - x, x + 1.0, 1.0 - d.b);
        if c != OK { return nan(), c; }
        return v, OK;
    }
    return nan(), ERR_DOMAIN;
}

fn dist_ccdf(d: Dist, x: f64) -> f64, i64 {
    if is_nan(x) { return nan(), ERR_NAN; }
    if d.kind == DS_NORMAL { return norm_sf(x, d.a, d.b); }
    if d.kind == DS_UNIFORM { return unif_sf(x, d.a, d.b); }
    if d.kind == DS_LOGISTIC { return logistic_cdf(-x); }
    if d.kind == DS_BINOMIAL {
        // A point under the support of the law has all the mass.
        if x < 0.0 { return 1.0, OK; }
        let v, c = dist_cdf(d, x);
        if c != OK { return nan(), c; }
        return 1.0 - v, OK;
    }
    return nan(), ERR_DOMAIN;
}

fn dist_icdf(d: Dist, p: f64) -> f64, i64 {
    let c = check_prob(p);
    if c != OK { return nan(), c; }
    let lo, hi, cs = dist_support(d);
    if cs != OK { return nan(), cs; }
    if p == 0.0 { return lo, OK; }
    if p == 1.0 { return hi, OK; }
    if d.kind == DS_NORMAL { return norm_ppf(p, d.a, d.b); }
    if d.kind == DS_UNIFORM { return unif_ppf(p, d.a, d.b); }
    if d.kind == DS_LOGISTIC { return logistic_ppf(p); }
    return binom_ppf(p, d.a, d.b);
}

fn dist_iccdf(d: Dist, p: f64) -> f64, i64 {
    let c = check_prob(p);
    if c != OK { return nan(), c; }
    return dist_icdf(d, 1.0 - p);
}

// ---- the logarithms -------------------------------------------------------

fn dist_logpdf(d: Dist, x: f64) -> f64, i64 {
    if d.kind == DS_NORMAL {
        if is_nan(x) { return nan(), ERR_NAN; }
        let z = (x - d.a) / d.b;
        return -0.5 * z * z - log(d.b) - 0.5 * log(2.0 * PI), OK;
    }
    let v, c = dist_pdf(d, x);
    if c != OK { return nan(), c; }
    return log(v), OK;
}

fn dist_logcdf(d: Dist, x: f64) -> f64, i64 {
    if d.kind == DS_NORMAL {
        if is_nan(x) { return nan(), ERR_NAN; }
        return fam_log_ndtr((x - d.a) / d.b), OK;
    }
    let v, c = dist_cdf(d, x);
    if c != OK { return nan(), c; }
    return log(v), OK;
}

fn dist_logccdf(d: Dist, x: f64) -> f64, i64 {
    if d.kind == DS_NORMAL {
        if is_nan(x) { return nan(), ERR_NAN; }
        return fam_log_ndtr(-(x - d.a) / d.b), OK;
    }
    let v, c = dist_ccdf(d, x);
    if c != OK { return nan(), c; }
    return log(v), OK;
}

fn dist_check_log(l: f64) -> i64 {
    if is_nan(l) { return ERR_NAN; }
    if l > 0.0 { return ERR_DOMAIN; }
    return OK;
}

fn dist_ilogcdf(d: Dist, l: f64) -> f64, i64 {
    let c = dist_check_log(l);
    if c != OK { return nan(), c; }
    return dist_icdf(d, exp(l));
}

fn dist_ilogccdf(d: Dist, l: f64) -> f64, i64 {
    let c = dist_check_log(l);
    if c != OK { return nan(), c; }
    return dist_iccdf(d, exp(l));
}

// ---- the moments ----------------------------------------------------------

fn dist_mean(d: Dist) -> f64, i64 {
    if d.kind == DS_NORMAL { return d.a, OK; }
    if d.kind == DS_UNIFORM { return 0.5 * (d.a + d.b), OK; }
    if d.kind == DS_LOGISTIC { return 0.0, OK; }
    if d.kind == DS_BINOMIAL { return d.a * d.b, OK; }
    return nan(), ERR_DOMAIN;
}

fn dist_variance(d: Dist) -> f64, i64 {
    if d.kind == DS_NORMAL { return d.b * d.b, OK; }
    if d.kind == DS_UNIFORM { return (d.b - d.a) * (d.b - d.a) / 12.0, OK; }
    if d.kind == DS_LOGISTIC { return PI * PI / 3.0, OK; }
    if d.kind == DS_BINOMIAL { return d.a * d.b * (1.0 - d.b), OK; }
    return nan(), ERR_DOMAIN;
}

fn dist_standard_deviation(d: Dist) -> f64, i64 {
    let v, c = dist_variance(d);
    if c != OK { return nan(), c; }
    return sqrt(v), OK;
}

fn dist_median(d: Dist) -> f64, i64 {
    return dist_icdf(d, 0.5);
}

fn dist_mode(d: Dist) -> f64, i64 {
    if d.kind == DS_NORMAL { return d.a, OK; }
    if d.kind == DS_UNIFORM { return 0.5 * (d.a + d.b), OK; }
    if d.kind == DS_LOGISTIC { return 0.0, OK; }
    if d.kind == DS_BINOMIAL { return min(floor((d.a + 1.0) * d.b), d.a), OK; }
    return nan(), ERR_DOMAIN;
}

fn dist_skewness(d: Dist) -> f64, i64 {
    if d.kind == DS_NORMAL || d.kind == DS_UNIFORM || d.kind == DS_LOGISTIC { return 0.0, OK; }
    if d.kind == DS_BINOMIAL {
        let v = d.a * d.b * (1.0 - d.b);
        if v == 0.0 { return nan(), ERR_UNDEFINED; }
        return (1.0 - 2.0 * d.b) / sqrt(v), OK;
    }
    return nan(), ERR_DOMAIN;
}

fn dist_kurtosis(d: Dist) -> f64, i64 {
    if d.kind == DS_NORMAL { return 3.0, OK; }
    if d.kind == DS_UNIFORM { return 1.8, OK; }
    if d.kind == DS_LOGISTIC { return 4.2, OK; }
    if d.kind == DS_BINOMIAL {
        let v = d.a * d.b * (1.0 - d.b);
        if v == 0.0 { return nan(), ERR_UNDEFINED; }
        return 3.0 + (1.0 - 6.0 * d.b * (1.0 - d.b)) / v, OK;
    }
    return nan(), ERR_DOMAIN;
}

fn dist_entropy(d: Dist) -> f64, i64 {
    if d.kind == DS_NORMAL { return 0.5 * log(2.0 * PI * exp(1.0) * d.b * d.b), OK; }
    if d.kind == DS_UNIFORM { return log(d.b - d.a), OK; }
    if d.kind == DS_LOGISTIC { return 2.0, OK; }
    if d.kind == DS_BINOMIAL {
        var s: f64 = 0.0;
        var k: f64 = 0.0;
        while k <= d.a {
            let m, c = binom_pmf(k, d.a, d.b);
            if c != OK { return nan(), c; }
            if m > 0.0 { s = s - m * log(m); }
            k = k + 1.0;
        }
        return s, OK;
    }
    return nan(), ERR_DOMAIN;
}

// The log of the entropy. A negative entropy gives the undefined fault.
fn dist_logentropy(d: Dist) -> f64, i64 {
    let h, c = dist_entropy(d);
    if c != OK { return nan(), c; }
    if h <= 0.0 { return nan(), ERR_UNDEFINED; }
    return log(h), OK;
}

// The raw moment of a whole order. The order is at most DS_MAX_ORDER.
// The raw moments of the logistic law. The moments of an odd order are
// zero, because the law is symmetric about zero. The even moments come
// from the series of the moment generating function pi t / sin(pi t).
// The series of sin(pi t) / (pi t) gives the recurrence of the
// coefficients of that function.
fn logistic_moment(order: i64) -> f64 {
    if order % 2 == 1 { return 0.0; }
    let m = order / 2;
    // The coefficients of the series of sin(pi t) / (pi t). The
    // coefficient of the power t^(2j) is (-1)^j pi^(2j) / (2j+1)!.
    var cs: f64[>..] = [];
    var fac: f64 = 1.0;
    for j in 0..m + 1 {
        if j > 0 { fac = fac * (2.0 * j as f64) * (2.0 * j as f64 + 1.0); }
        var c = pow(PI, 2.0 * j as f64) / fac;
        if j % 2 == 1 { c = -c; }
        cs.push(c);
    }
    // The coefficients of the series of pi t / sin(pi t). The
    // coefficient of the power t^(2j) follows from the coefficients of
    // the other series.
    var a: f64[>..] = [];
    for i in 0..2 * m + 1 { a.push(0.0); }
    a[0] = 1.0;
    for j in 1..m + 1 {
        var s: f64 = 0.0;
        for t in 1..j + 1 { s = s + cs[t] * a[2 * (j - t)]; }
        a[2 * j] = -s;
    }
    // The moment is the coefficient times the factorial of the order.
    var fact: f64 = 1.0;
    for i in 1..order + 1 { fact = fact * (i as f64); }
    return fact * a[order];
}

fn dist_moment(d: Dist, order: i64) -> f64, i64 {
    if order < 0 || order > DS_MAX_ORDER { return nan(), ERR_DOMAIN; }
    if order == 0 { return 1.0, OK; }
    if d.kind == DS_NORMAL {
        // m(k) = mu m(k-1) + (k-1) sigma^2 m(k-2)
        var m0: f64 = 1.0;
        var m1: f64 = d.a;
        for k in 2..order + 1 {
            let m2 = d.a * m1 + (k - 1) as f64 * d.b * d.b * m0;
            m0 = m1;
            m1 = m2;
        }
        return m1, OK;
    }
    let k1 = (order + 1) as f64;
    if d.kind == DS_UNIFORM {
        return (pow(d.b, k1) - pow(d.a, k1)) / (k1 * (d.b - d.a)), OK;
    }
    if d.kind == DS_LOGISTIC {
        return logistic_moment(order), OK;
    }
    if d.kind == DS_BINOMIAL {
        var s: f64 = 0.0;
        var k: f64 = 0.0;
        while k <= d.a {
            let m, c = binom_pmf(k, d.a, d.b);
            if c != OK { return nan(), c; }
            s = s + m * pow(k, order as f64);
            k = k + 1.0;
        }
        return s, OK;
    }
    return nan(), ERR_DOMAIN;
}

// ---- the mixture ----------------------------------------------------------

// The check of the weights of a mixture. The weights are at least 0 and
// sum to 1. A count of components that differs from the count of weights
// gives the length fault.
fn mixture_check(comps: const Dist[:], w: const f64[:]) -> i64 {
    if comps.len == 0 { return ERR_EMPTY; }
    if w.len != comps.len { return ERR_LENGTH; }
    var s: f64 = 0.0;
    for v in w {
        if is_nan(v) { return ERR_NAN; }
        if is_inf(v) { return ERR_INFINITY; }
        if v < 0.0 { return ERR_DOMAIN; }
        s = s + v;
    }
    if abs(s - 1.0) > 1e-9 { return ERR_DOMAIN; }
    for c in comps {
        if c.kind == DS_BINOMIAL || c.kind == 0 { return ERR_DOMAIN; }
    }
    return OK;
}

// The mixture takes continuous components only, as in SciPy.
fn mixture_pdf(comps: const Dist[:], w: const f64[:], x: f64) -> f64, i64 {
    let k = mixture_check(comps, w);
    if k != OK { return nan(), k; }
    var s: f64 = 0.0;
    for i in 0..comps.len {
        let v, c = dist_pdf(comps[i], x);
        if c != OK { return nan(), c; }
        s = s + w[i] * v;
    }
    return s, OK;
}

fn mixture_cdf(comps: const Dist[:], w: const f64[:], x: f64) -> f64, i64 {
    let k = mixture_check(comps, w);
    if k != OK { return nan(), k; }
    var s: f64 = 0.0;
    for i in 0..comps.len {
        let v, c = dist_cdf(comps[i], x);
        if c != OK { return nan(), c; }
        s = s + w[i] * v;
    }
    return s, OK;
}

fn mixture_ccdf(comps: const Dist[:], w: const f64[:], x: f64) -> f64, i64 {
    let v, c = mixture_cdf(comps, w, x);
    if c != OK { return nan(), c; }
    return 1.0 - v, OK;
}

fn mixture_icdf(comps: const Dist[:], w: const f64[:], p: f64) -> f64, i64 {
    let k = mixture_check(comps, w);
    if k != OK { return nan(), k; }
    let cp = check_prob(p);
    if cp != OK { return nan(), cp; }
    if p == 0.0 { return -inf(), OK; }
    if p == 1.0 { return inf(), OK; }
    // The root lies between the lowest and the highest quantile of the
    // components at the same probability.
    var lo: f64 = inf();
    var hi: f64 = -inf();
    for i in 0..comps.len {
        let q, c = dist_icdf(comps[i], p);
        if c != OK { return nan(), c; }
        lo = min(lo, q);
        hi = max(hi, q);
    }
    if lo == hi { return lo, OK; }
    let x, c = brent_root(lo, hi, 1e-14, 1e-15) { v => mixture_cdf_value(comps, w, v) - p };
    if c != OK { return nan(), c; }
    return x, OK;
}

fn mixture_cdf_value(comps: const Dist[:], w: const f64[:], x: f64) -> f64 {
    let v, c = mixture_cdf(comps, w, x);
    if c != OK { return nan(); }
    return v;
}

fn mixture_iccdf(comps: const Dist[:], w: const f64[:], p: f64) -> f64, i64 {
    let c = check_prob(p);
    if c != OK { return nan(), c; }
    return mixture_icdf(comps, w, 1.0 - p);
}

// The raw moments of a mixture, up to the order n. The buffer holds the
// moments of order 1 to n. The check and the moments of every component
// are read once.
fn mixture_moments(out: f64[>..]&, comps: const Dist[:], w: const f64[:], n: i64) -> i64 {
    out.resize(0);
    let k = mixture_check(comps, w);
    if k != OK { return k; }
    for i in 0..n { out.push(0.0); }
    for i in 0..comps.len {
        for o in 1..n + 1 {
            let v, c = dist_moment(comps[i], o);
            if c != OK { return c; }
            out[o - 1] = out[o - 1] + w[i] * v;
        }
    }
    return OK;
}

fn mixture_moment(comps: const Dist[:], w: const f64[:], order: i64) -> f64, i64 {
    var m: f64[>..] = [];
    let k = mixture_moments(m, comps, w, order);
    if k != OK { return nan(), k; }
    return m[order - 1], OK;
}

fn mixture_mean(comps: const Dist[:], w: const f64[:]) -> f64, i64 {
    return mixture_moment(comps, w, 1);
}

fn mixture_variance(comps: const Dist[:], w: const f64[:]) -> f64, i64 {
    var m: f64[>..] = [];
    let k = mixture_moments(m, comps, w, 2);
    if k != OK { return nan(), k; }
    return m[1] - m[0] * m[0], OK;
}

fn mixture_standard_deviation(comps: const Dist[:], w: const f64[:]) -> f64, i64 {
    let v, c = mixture_variance(comps, w);
    if c != OK { return nan(), c; }
    return sqrt(v), OK;
}

fn mixture_skewness(comps: const Dist[:], w: const f64[:]) -> f64, i64 {
    var m: f64[>..] = [];
    let k = mixture_moments(m, comps, w, 3);
    if k != OK { return nan(), k; }
    let m1 = m[0];
    let v = m[1] - m1 * m1;
    return (m[2] - 3.0 * m1 * v - m1 * m1 * m1) / pow(v, 1.5), OK;
}

fn mixture_kurtosis(comps: const Dist[:], w: const f64[:]) -> f64, i64 {
    var m: f64[>..] = [];
    let k = mixture_moments(m, comps, w, 4);
    if k != OK { return nan(), k; }
    let m1 = m[0];
    let v = m[1] - m1 * m1;
    let m4c = m[3] - 4.0 * m1 * m[2] + 6.0 * m1 * m1 * m[1] - 3.0 * m1 * m1 * m1 * m1;
    return m4c / (v * v), OK;
}

fn mixture_support(comps: const Dist[:], w: const f64[:]) -> f64, f64, i64 {
    let k = mixture_check(comps, w);
    if k != OK { return nan(), nan(), k; }
    var lo: f64 = inf();
    var hi: f64 = -inf();
    for c in comps {
        let a, b, cs = dist_support(c);
        if cs != OK { return nan(), nan(), cs; }
        lo = min(lo, a);
        hi = max(hi, b);
    }
    return lo, hi, OK;
}
