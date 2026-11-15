// Eight continuous families of the second tail group: crystalball,
// kstwobign, vonmises_line, ncx2, ncf, nct, ksone and kstwo. The names
// and the parameters follow SciPy 1.18.1, with a location of 0 and a
// scale of 1.
//
// Every family has three calls: pdf, cdf and ppf. Every call returns
// the result first and a fault code last.
//
// The non-central laws of chi-square and F are Poisson mixtures of the
// central laws. The non-central t law is an integral over the scale of a
// chi law. The law ksone has the sum of Birnbaum and Tingey. The law
// kstwo has the matrix method of Marsaglia, Tsang and Wang.
namespace statistics;

import std;
import math;

const NC_MAX: f64 = 100000.0;       // the largest non-centrality of a mixture

// The two sided bracket of a quantile of a law of the whole line. The
// block gives the cumulative function.
fn tl_ppf_line<F>(q: f64, center: f64) -> f64, i64 {
    var lo = center - 1.0;
    var hi = center + 1.0;
    var steps: i64 = 0;
    while steps < MAX_FAMILY_SEARCH && F(lo) > q {
        lo = center - 2.0 * (center - lo);
        steps += 1;
    }
    while steps < MAX_FAMILY_SEARCH && F(hi) < q {
        hi = center + 2.0 * (hi - center);
        steps += 1;
    }
    if steps >= MAX_FAMILY_SEARCH { return nan(), ERR_FIT; }
    let x, cx = brent_root(lo, hi, 1e-13, 4e-16) { v => F(v) - q };
    return x, cx;
}

// ---- crystalball ----------------------------------------------------------

fn crystalball_check(beta: f64, m: f64) -> i64 {
    let c1 = check_shape(beta);
    if c1 != OK { return c1; }
    let c2 = check_shape(m);
    if c2 != OK { return c2; }
    if m <= 1.0 { return ERR_DOMAIN; }
    return OK;
}

// The reciprocal of the constant of the density.
fn crystalball_norm(beta: f64, m: f64) -> f64 {
    return m / beta / (m - 1.0) * exp(-0.5 * beta * beta) + SQRT_2PI * tl_ncdf(beta);
}

fn crystalball_pdf(x: f64, beta: f64, m: f64) -> f64, i64 {
    let k = crystalball_check(beta, m);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    let n = 1.0 / crystalball_norm(beta, m);
    if x > -beta { return tl_out(n * exp(-0.5 * x * x)); }
    let l = m * log(m / beta) - 0.5 * beta * beta - m * log(m / beta - beta - x);
    return tl_out(n * exp(l));
}

fn crystalball_cdf(x: f64, beta: f64, m: f64) -> f64, i64 {
    let k = crystalball_check(beta, m);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    let n = 1.0 / crystalball_norm(beta, m);
    if x > -beta {
        let tail = m / beta * exp(-0.5 * beta * beta) / (m - 1.0);
        return tl_out(min(1.0, n * (tail + SQRT_2PI * (tl_ncdf(x) - tl_ncdf(-beta)))));
    }
    let l = m * log(m / beta) - 0.5 * beta * beta + (1.0 - m) * log(m / beta - beta - x) - log(m - 1.0);
    return tl_out(n * exp(l));
}

fn crystalball_ppf(p: f64, beta: f64, m: f64) -> f64, i64 {
    let k = crystalball_check(beta, m);
    if k != OK { return nan(), k; }
    let pc = check_prob(p);
    if pc != OK { return nan(), pc; }
    if p == 0.0 { return -inf(), OK; }
    if p == 1.0 { return inf(), OK; }
    let eb2 = exp(-0.5 * beta * beta);
    let c = m / beta * eb2 / (m - 1.0);
    let n = 1.0 / (c + SQRT_2PI * tl_ncdf(beta));
    if p < n * c {
        let inner = (m - 1.0) * pow(m / beta, -m) / eb2 * p / n;
        return tl_out(m / beta - beta - pow(inner, 1.0 / (1.0 - m)));
    }
    return tl_out(tl_nppf(tl_ncdf(-beta) + (p / n - c) / SQRT_2PI));
}

// ---- kstwobign: the limit of the Kolmogorov-Smirnov statistic -------------

// The cumulative function and the density at x above 0. Below 1 the form
// of the theta function converges fast. From 1 the alternating form does.
fn kstwobign_both(x: f64) -> f64, f64 {
    var cdf: f64 = 0.0;
    var pdf: f64 = 0.0;
    if x < 1.0 {
        for j in 1..9 {
            let o = (2 * j - 1) as f64;
            let a = o * o * PI * PI / 8.0;
            let t = SQRT_2PI / x * exp(-a / (x * x));
            cdf += t;
            pdf += t * (-1.0 / x + 2.0 * a / (x * x * x));
        }
        return cdf, pdf;
    }
    for j in 1..9 {
        let kk = j as f64;
        var sign = 1.0;
        if j % 2 == 0 { sign = -1.0; }
        let e = exp(-2.0 * kk * kk * x * x);
        cdf += sign * e;
        pdf += sign * kk * kk * e;
    }
    return 1.0 - 2.0 * cdf, 8.0 * x * pdf;
}

fn kstwobign_pdf(x: f64) -> f64, i64 {
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    let c, p = kstwobign_both(x);
    return tl_out(max(p, 0.0));
}

fn kstwobign_cdf(x: f64) -> f64, i64 {
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x <= 0.0 { return 0.0, OK; }
    let c, p = kstwobign_both(x);
    return tl_out(min(max(c, 0.0), 1.0));
}

fn kstwobign_cdf_value(x: f64) -> f64 {
    let v, c = kstwobign_cdf(x);
    if c != OK { return nan(); }
    return v;
}

fn kstwobign_ppf(q: f64) -> f64, i64 {
    let k = check_prob(q);
    if k != OK { return nan(), k; }
    if q == 0.0 { return 0.0, OK; }
    if q == 1.0 { return inf(), OK; }
    return fam_search_ppf(q, 0.0, 1.0) { v => kstwobign_cdf_value(v) };
}

// ---- vonmises_line: the von Mises law on the line from -pi to pi ----------

fn vonmises_line_pdf(x: f64, kappa: f64) -> f64, i64 {
    let k = check_shape_nonneg(kappa);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    if x < -PI || x > PI { return 0.0, OK; }
    return tl_out(exp(kappa * (cos(x) - 1.0)) / (2.0 * PI * vonmises_i0e(kappa)));
}

fn vonmises_line_cdf_value(x: f64, kappa: f64) -> f64 {
    if x <= -PI { return 0.0; }
    if x >= PI { return 1.0; }
    let c = 1.0 / (2.0 * PI * vonmises_i0e(kappa));
    // The law is symmetric, so the integral from zero gives the half.
    let v = fam_quad(0.0, x) { t => exp(kappa * (cos(t) - 1.0)) * c };
    return min(max(0.5 + v, 0.0), 1.0);
}

fn vonmises_line_cdf(x: f64, kappa: f64) -> f64, i64 {
    let k = check_shape_nonneg(kappa);
    if k != OK { return nan(), k; }
    let s = fam_x(x);
    if s != OK { return nan(), s; }
    return tl_out(vonmises_line_cdf_value(x, kappa));
}

fn vonmises_line_ppf(q: f64, kappa: f64) -> f64, i64 {
    let k = check_shape_nonneg(kappa);
    if k != OK { return nan(), k; }
    let pc = check_prob(q);
    if pc != OK { return nan(), pc; }
    if q == 0.0 { return -PI, OK; }
    if q == 1.0 { return PI, OK; }
    if q == 0.5 { return 0.0, OK; }
    let x, cx = brent_root(-PI, PI, 1e-14, 4e-16) { v => vonmises_line_cdf_value(v, kappa) - q };
    return x, cx;
}

// ---- the Poisson weights of the mixtures ----------------------------------

// The range of the Poisson counts that hold the mass of a mixture, for a
// mean of lam. The range spans twelve standard deviations.
fn nc_range(lam: f64) -> f64, f64 {
    let w = 12.0 * sqrt(lam) + 40.0;
    return max(0.0, floor(lam - w)), ceil(lam + w);
}

// The weight of the count j of the Poisson law of the mean lam.
fn nc_weight(j: f64, lam: f64) -> f64 {
    if lam == 0.0 {
        if j == 0.0 { return 1.0; }
        return 0.0;
    }
    return exp(-lam + j * log(lam) - fam_lngamma(j + 1.0));
}

// The weights of every count of the mixture, from lo to hi. A quantile
// search builds them once and reuses them at every step.
fn nc_weights(lam: f64, lo: i64, hi: i64, ws: f64[>..]&) {
    ws.resize(0);
    for j in lo..hi + 1 {
        ws.push(nc_weight(j, lam));
    }
}

// The weights of the mixture of the mean lam = nc / 2. The call returns
// the lowest count of the range, and it fills ws with the weights.
fn nc_prep(nc: f64, ws: f64[>..]&) -> i64 {
    let lam = 0.5 * nc;
    let lo, hi = nc_range(lam);
    nc_weights(lam, lo as i64, hi as i64, ws);
    return lo as i64;
}

// ---- ncx2: the non-central chi-square law ---------------------------------

fn ncx2_check(df: f64, nc: f64) -> i64 {
    let c1 = check_shape(df);
    if c1 != OK { return c1; }
    let c2 = check_shape_nonneg(nc);
    if c2 != OK { return c2; }
    if nc > NC_MAX { return ERR_DOMAIN; }
    return OK;
}

fn ncx2_pdf(x: f64, df: f64, nc: f64) -> f64, i64 {
    let k = tl_px(x, ncx2_check(df, nc));
    if k != OK { return nan(), k; }
    if x <= 0.0 && df >= 2.0 && !(x == 0.0 && df == 2.0) { return 0.0, OK; }
    if x < 0.0 { return 0.0, OK; }
    let lam = 0.5 * nc;
    let lo, hi = nc_range(lam);
    var s: f64 = 0.0;
    var j = lo;
    while j <= hi {
        let p, cp = chi2_pdf(x, df + 2.0 * j);
        if cp == OK { s += nc_weight(j, lam) * p; }
        j += 1.0;
    }
    return tl_out(s);
}

// The cumulative value from prepared weights. The cumulative values of
// the chi-square law of one count follow from the value of the count
// before: the shape grows by one, so one recurrence step costs one
// multiply and one subtract.
fn ncx2_cdf_core(x: f64, df: f64, lo: i64, ws: const f64[:]) -> f64 {
    if x <= 0.0 { return 0.0; }
    let z = 0.5 * x;
    let lz = log(z);
    var a = 0.5 * df + lo;
    let p0, c0 = gammainc(a, z);
    if c0 != OK { return 0.0; }
    let lg, sg, cg = lgamma(a + 1.0);
    if cg != OK { return 0.0; }
    // d holds z^a exp(-z) / gamma(a + 1), the step of the recurrence.
    var d = exp(a * lz - z - lg);
    var p = p0;
    var s = ws[0] * p0;
    for i in 1..ws.len {
        p = p - d;
        s += ws[i] * p;
        d = d * z / (a + 1.0);
        a = a + 1.0;
    }
    return min(s, 1.0);
}

fn ncx2_cdf_value(x: f64, df: f64, nc: f64) -> f64 {
    let lam = 0.5 * nc;
    let lo, hi = nc_range(lam);
    var ws: f64[>..] = [];
    nc_weights(lam, lo as i64, hi as i64, ws);
    return ncx2_cdf_core(x, df, lo as i64, ws[..]);
}

fn ncx2_cdf(x: f64, df: f64, nc: f64) -> f64, i64 {
    let k = tl_px(x, ncx2_check(df, nc));
    if k != OK { return nan(), k; }
    return tl_out(ncx2_cdf_value(x, df, nc));
}

fn ncx2_ppf(q: f64, df: f64, nc: f64) -> f64, i64 {
    let k = tl_pq(q, ncx2_check(df, nc));
    if k != OK { return nan(), k; }
    if q == 0.0 { return 0.0, OK; }
    if q == 1.0 { return inf(), OK; }
    // The weights serve every step of the search.
    var ws: f64[>..] = [];
    let lo = nc_prep(nc, ws);
    return fam_search_ppf(q, 0.0, df + nc + 1.0) { v => ncx2_cdf_core(v, df, lo, ws[..]) };
}

// ---- ncf: the non-central F law -------------------------------------------

fn ncf_check(dfn: f64, dfd: f64, nc: f64) -> i64 {
    let c1 = check_shape(dfn);
    if c1 != OK { return c1; }
    let c2 = check_shape(dfd);
    if c2 != OK { return c2; }
    let c3 = check_shape_nonneg(nc);
    if c3 != OK { return c3; }
    if nc > NC_MAX { return ERR_DOMAIN; }
    return OK;
}

// The density from prepared weights. The law of one count is a scaled
// F law, and every scaled F law of the mixture carries the beta density
// of the point u = dfn x / (dfn x + dfd) with the same u and the same
// second shape. The beta densities follow one from the next by a
// recurrence, so the loop holds no special call.
fn ncf_pdf_core(x: f64, dfn: f64, dfd: f64, lo: i64, ws: const f64[:]) -> f64 {
    let u = dfn * x / (dfn * x + dfd);
    let b = 0.5 * dfd;
    var a = 0.5 * dfn + lo;
    // q holds the beta density of u for the shape a.
    let lb, clb = lbeta(a, b);
    if clb != OK { return 0.0; }
    var q = exp((a - 1.0) * log(u) + (b - 1.0) * log(1.0 - u) - lb);
    if is_nan(q) { q = 0.0; }
    var s = ws[0] * q;
    for i in 1..ws.len {
        q = q * u * (a + b) / a;
        s += ws[i] * q;
        a = a + 1.0;
    }
    return s * dfn * dfd / ((dfn * x + dfd) * (dfn * x + dfd));
}

fn ncf_pdf(x: f64, dfn: f64, dfd: f64, nc: f64) -> f64, i64 {
    let k = tl_px(x, ncf_check(dfn, dfd, nc));
    if k != OK { return nan(), k; }
    if x < 0.0 { return 0.0, OK; }
    var ws: f64[>..] = [];
    let lo = nc_prep(nc, ws);
    return tl_out(ncf_pdf_core(x, dfn, dfd, lo, ws[..]));
}

// The cumulative value from prepared weights. The cumulative values of
// the beta law of one count follow from the value of the count before:
// the first shape grows by one, so one recurrence step costs one
// multiply and one subtract.
fn ncf_cdf_core(x: f64, dfn: f64, dfd: f64, lo: i64, ws: const f64[:]) -> f64 {
    if x <= 0.0 { return 0.0; }
    let u = dfn * x / (dfn * x + dfd);
    let b = 0.5 * dfd;
    var a = 0.5 * dfn + lo;
    let p0, c0 = betainc(a, b, u);
    if c0 != OK { return 0.0; }
    // d holds u^a (1 - u)^b / (a B(a, b)), the step of the recurrence.
    // The three logarithms of the gamma come direct: the core sits in
    // the hot loop of the quantile search, where the pair path of
    // lbeta costs more.
    let d0 = exp(a * log(u) + b * log(1.0 - u) - fam_lngamma(a) - fam_lngamma(b) + fam_lngamma(a + b)) / a;
    var d: f64 = d0;
    if is_nan(d) { d = 0.0; }
    var p = p0;
    var s = ws[0] * p0;
    for i in 1..ws.len {
        p = p - d;
        s += ws[i] * p;
        d = d * u * (a + b) / (a + 1.0);
        a = a + 1.0;
    }
    return min(s, 1.0);
}

fn ncf_cdf_value(x: f64, dfn: f64, dfd: f64, nc: f64) -> f64 {
    var ws: f64[>..] = [];
    let lo = nc_prep(nc, ws);
    return ncf_cdf_core(x, dfn, dfd, lo, ws[..]);
}

fn ncf_cdf(x: f64, dfn: f64, dfd: f64, nc: f64) -> f64, i64 {
    let k = tl_px(x, ncf_check(dfn, dfd, nc));
    if k != OK { return nan(), k; }
    return tl_out(ncf_cdf_value(x, dfn, dfd, nc));
}

fn ncf_ppf(q: f64, dfn: f64, dfd: f64, nc: f64) -> f64, i64 {
    let k = tl_pq(q, ncf_check(dfn, dfd, nc));
    if k != OK { return nan(), k; }
    if q == 0.0 { return 0.0, OK; }
    if q == 1.0 { return inf(), OK; }
    // The weights serve every step of the search.
    var ws: f64[>..] = [];
    let lo = nc_prep(nc, ws);
    return fam_search_ppf(q, 0.0, 2.0) { v => ncf_cdf_core(v, dfn, dfd, lo, ws[..]) };
}

// ---- nct: the non-central t law -------------------------------------------

fn nct_check(df: f64, nc: f64) -> i64 {
    let c1 = check_shape(df);
    if c1 != OK { return c1; }
    return check_scalar(nc);
}

// The edges of the panels of the integral over the scale s. The scale has
// a mean near 1 and a width near 1 / sqrt(2 df).
fn nct_edges(df: f64, edges: f64[>..]&) {
    edges.resize(0);
    let w = 1.0 / sqrt(2.0 * df);
    let zs: f64[15] = [-40.0, -20.0, -12.0, -8.0, -5.0, -3.0, -1.5, 0.0, 1.5, 3.0, 5.0, 8.0, 12.0, 20.0, 40.0];
    edges.push(0.0);
    for i in 0..zs.len {
        let e = 1.0 + w * zs[i];
        if e > edges[edges.len - 1] { edges.push(e); }
    }
}

// The log density of the scale s = sqrt(chi2 / df).
fn nct_lscale(s: f64, df: f64) -> f64 {
    return log(2.0) + 0.5 * df * log(0.5 * df) - fam_lngamma(0.5 * df) + (df - 1.0) * log(s) - 0.5 * df * s * s;
}

// The nodes of the fixed rule over the scale, with the density of the
// scale folded into the weight of every node. The quantile search builds
// them once and reuses them at every step. The integrand of the
// cumulative function is analytic on every panel, so one fixed rule of
// the panel is enough.
fn nct_prep(df: f64, sp: f64[>..]&, sw: f64[>..]&) {
    var edges: f64[>..] = [];
    nct_edges(df, edges);
    var gx: f64[>..] = [];
    var gw: f64[>..] = [];
    gauss_legendre(gx, gw, FAM_QUAD_NODES);
    let c = log(2.0) + 0.5 * df * log(0.5 * df) - fam_lngamma(0.5 * df);
    sp.resize(0);
    sw.resize(0);
    for i in 0..edges.len - 1 {
        let lo = edges[i];
        let hi = edges[i + 1];
        let w = (hi - lo) / 2.0;
        for j in 0..gx.len {
            let s = lo + w * (gx[j] + 1.0);
            sp.push(s);
            sw.push(gw[j] * w * exp(c + (df - 1.0) * log(s) - 0.5 * df * s * s));
        }
    }
}

// The cumulative value from prepared nodes. The weight of every node
// holds the density of the scale, so the loop evaluates one normal
// cumulative value per node.
fn nct_cdf_core(t: f64, nc: f64, sp: const f64[:], sw: const f64[:]) -> f64 {
    var total: f64 = 0.0;
    for i in 0..sp.len {
        total += sw[i] * tl_ncdf(t * sp[i] - nc);
    }
    return min(max(total, 0.0), 1.0);
}

// The cumulative function at t, as an integral over the scale.
fn nct_cdf_value(t: f64, df: f64, nc: f64) -> f64 {
    var sp: f64[>..] = [];
    var sw: f64[>..] = [];
    nct_prep(df, sp, sw);
    return nct_cdf_core(t, nc, sp[..], sw[..]);
}

fn nct_pdf(t: f64, df: f64, nc: f64) -> f64, i64 {
    let k = tl_px(t, nct_check(df, nc));
    if k != OK { return nan(), k; }
    var edges: f64[>..] = [];
    nct_edges(df, edges);
    var total: f64 = 0.0;
    for i in 0..edges.len - 1 {
        total += fam_quad(edges[i], edges[i + 1]) { s => s * exp(-0.5 * (t * s - nc) * (t * s - nc) + nct_lscale(s, df)) / SQRT_2PI };
    }
    return tl_out(total);
}

fn nct_cdf(t: f64, df: f64, nc: f64) -> f64, i64 {
    let k = tl_px(t, nct_check(df, nc));
    if k != OK { return nan(), k; }
    return tl_out(nct_cdf_value(t, df, nc));
}

fn nct_ppf(q: f64, df: f64, nc: f64) -> f64, i64 {
    let k = tl_pq(q, nct_check(df, nc));
    if k != OK { return nan(), k; }
    if q == 0.0 { return -inf(), OK; }
    if q == 1.0 { return inf(), OK; }
    // The nodes serve every step of the search.
    var sp: f64[>..] = [];
    var sw: f64[>..] = [];
    nct_prep(df, sp, sw);
    return tl_ppf_line(q, nc) { v => nct_cdf_core(v, nc, sp[..], sw[..]) };
}
