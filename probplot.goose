// The probability plot and the probability plot correlation.
//
// The plot pairs the sorted sample with the quantiles of a family at the
// uniform order statistic medians of Filliben. The line through the
// pairs gives the slope, the intercept and the correlation.
//
// ppcc_plot gives the correlation for a grid of shapes. ppcc_max
// searches the shape with the best correlation by the bracket and the
// Brent minimizer of SciPy. The search starts from the bracket of the
// caller, which is from 0 to 1 by default. The tolerance is 1.48e-8.
//
// A sample above MAX_PPCC_POINTS values gives the too long fault in the
// two calls that search. A constant sample has no correlation and gives
// the undefined fault.
namespace statistics;

import std;
import math;

const PLOT_TOL: f64 = 1.48e-8;    // relative tolerance of the Brent search

// The medians of the order statistics of the uniform law, n of them.
fn plot_medians(n: i64, out: f64[>..]&) {
    out.resize(n, 0.0);
    let nf = n as f64;
    let last = pow(0.5, 1.0 / nf);
    out[n - 1] = last;
    out[0] = 1.0 - last;
    for k in 1..n - 1 {
        out[k] = ((k + 1) - 0.3175) / (nf + 0.365);
    }
}

// The quantile of the Tukey lambda family. Lambda 0 is the logistic law.
fn tukey_ppf(p: f64, lam: f64) -> f64 {
    if lam == 0.0 { return log(p) - log(1.0 - p); }
    return (pow(p, lam) - pow(1.0 - p, lam)) / lam;
}

// The quantile of a family at one probability inside the open unit range.
fn plot_ppf(family: i64, shape: f64, p: f64) -> f64, i64 {
    if family == PLOT_NORMAL { return norm_ppf(p); }
    if family == PLOT_EXPON { return expon_ppf(p); }
    if family == PLOT_UNIFORM { return unif_ppf(p); }
    if family == PLOT_TUKEY { return tukey_ppf(p, shape), OK; }
    if family == PLOT_LOGISTIC { return logi_ppf(p, 0.0, 1.0); }
    if family == PLOT_LAPLACE { return laplace_ppf(p); }
    if family == PLOT_CAUCHY { return cauchy_ppf(p); }
    if family == PLOT_GUMBEL_R { return gumbel_r_ppf(p); }
    if family == PLOT_GUMBEL_L { return gumbel_l_ppf(p); }
    if family == PLOT_RAYLEIGH { return rayleigh_ppf(p); }
    if family == PLOT_WEIBULL { return weibull_ppf(p, shape); }
    return nan(), ERR_DOMAIN;
}

// The code of a family and its shape.
fn plot_check(family: i64, shape: f64) -> i64 {
    if family < PLOT_NORMAL || family > PLOT_WEIBULL { return ERR_DOMAIN; }
    if family != PLOT_TUKEY && family != PLOT_WEIBULL { return OK; }
    if is_nan(shape) { return ERR_NAN; }
    if is_inf(shape) { return ERR_INFINITY; }
    if family == PLOT_WEIBULL && shape <= 0.0 { return ERR_DOMAIN; }
    return OK;
}

// The theoretical quantiles at the medians of Filliben, in out.
fn plot_quantiles(out: f64[>..]&, med: const f64[:], family: i64, shape: f64) -> i64 {
    out.resize(0);
    for p in med {
        let q, c = plot_ppf(family, shape, p);
        if c != OK { out.resize(0); return c; }
        out.push(q);
    }
    return OK;
}

// The sorted copy of a sample.
fn plot_sorted(out: f64[>..]&, xs: const f64[:]) {
    out.resize(0);
    for x in xs { out.push(x); }
    sort(out) { a, b => a < b };
}

// The probability plot of a sample. The lists osm and osr receive the
// theoretical quantiles and the sorted sample. The line is the least
// squares line of osr on osm.
fn probplot(osm: f64[>..]&, osr: f64[>..]&, xs: const f64[:], family: i64 = PLOT_NORMAL, shape: f64 = 0.0) -> ProbabilityPlot, i64 {
    osm.resize(0);
    osr.resize(0);
    let none = ProbabilityPlot { slope: nan(), intercept: nan(), r: nan() };
    let c = check_sample(xs, 2);
    if c != OK { return none, c; }
    let cf = plot_check(family, shape);
    if cf != OK { return none, cf; }
    var med: f64[>..] = [];
    plot_medians(xs.len, med);
    let cq = plot_quantiles(osm, med[..], family, shape);
    if cq != OK { return none, cq; }
    plot_sorted(osr, xs);
    let line, cl = linear_regression(osm[..], osr[..], false);
    if cl != OK { osm.resize(0); osr.resize(0); return none, cl; }
    return ProbabilityPlot { slope: line.slope, intercept: line.intercept, r: line.r }, OK;
}

// The correlation of the plot at one shape.
fn plot_corr(osr: const f64[:], med: const f64[:], family: i64, shape: f64) -> f64 {
    var q: f64[>..] = [];
    let c = plot_quantiles(q, med, family, shape);
    if c != OK { return nan(); }
    let line, cl = linear_regression(q[..], osr, false);
    if cl != OK { return nan(); }
    return line.r;
}

// The correlation of the plot for n shapes from a to b. The family is the
// Tukey lambda family or the Weibull family.
fn ppcc_plot(out: f64[>..]&, xs: const f64[:], a: f64, b: f64, n: i64, family: i64 = PLOT_TUKEY) -> i64 {
    out.resize(0);
    let c = check_sample(xs, 2);
    if c != OK { return c; }
    if xs.len > MAX_PPCC_POINTS || n > MAX_PPCC_GRID { return ERR_TOO_LONG; }
    if is_nan(a) || is_nan(b) { return ERR_NAN; }
    if is_inf(a) || is_inf(b) { return ERR_INFINITY; }
    if n < 1 || !(b > a) { return ERR_DOMAIN; }
    if family != PLOT_TUKEY && family != PLOT_WEIBULL { return ERR_DOMAIN; }
    var med: f64[>..] = [];
    plot_medians(xs.len, med);
    var osr: f64[>..] = [];
    plot_sorted(osr, xs);
    for k in 0..n {
        var s = a;
        if n > 1 { s = a + (b - a) * k / (n - 1); }
        if k == n - 1 && n > 1 { s = b; }
        if plot_check(family, s) != OK { out.resize(0); return ERR_DOMAIN; }
        let r = plot_corr(osr[..], med[..], family, s);
        if is_nan(r) { out.resize(0); return ERR_UNDEFINED; }
        out.push(r);
    }
    return OK;
}

// The shape with the best correlation. The search starts from the two
// points lo and hi and widens them to a bracket.
fn ppcc_max(xs: const f64[:], family: i64 = PLOT_TUKEY, lo: f64 = 0.0, hi: f64 = 1.0) -> f64, i64 {
    let c = check_sample(xs, 2);
    if c != OK { return nan(), c; }
    if xs.len > MAX_PPCC_POINTS { return nan(), ERR_TOO_LONG; }
    if is_nan(lo) || is_nan(hi) { return nan(), ERR_NAN; }
    if is_inf(lo) || is_inf(hi) { return nan(), ERR_INFINITY; }
    if family != PLOT_TUKEY && family != PLOT_WEIBULL { return nan(), ERR_DOMAIN; }
    var med: f64[>..] = [];
    plot_medians(xs.len, med);
    var osr: f64[>..] = [];
    plot_sorted(osr, xs);
    if is_nan(plot_corr(osr[..], med[..], family, lo)) { return nan(), ERR_UNDEFINED; }
    let s, cs = minimize_scalar(lo, hi, PLOT_TOL) { v => ppcc_loss(osr[..], med[..], family, v) };
    if cs != OK { return nan(), ERR_FIT; }
    return s, OK;
}

// One minus the correlation. A shape that gives no correlation counts as
// the worst fit.
fn ppcc_loss(osr: const f64[:], med: const f64[:], family: i64, s: f64) -> f64 {
    if plot_check(family, s) != OK { return 2.0; }
    let r = plot_corr(osr, med, family, s);
    if is_nan(r) { return 2.0; }
    return 1.0 - r;
}
