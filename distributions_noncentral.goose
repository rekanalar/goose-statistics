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
const KS_MAX_N: f64 = 1000.0;       // the largest sample size of ksone and kstwo

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

fn ncx2_cdf_value(x: f64, df: f64, nc: f64) -> f64 {
    if x <= 0.0 { return 0.0; }
    let lam = 0.5 * nc;
    let lo, hi = nc_range(lam);
    var s: f64 = 0.0;
    var j = lo;
    while j <= hi {
        let p, cp = chi2_cdf(x, df + 2.0 * j);
        if cp == OK { s += nc_weight(j, lam) * p; }
        j += 1.0;
    }
    return min(s, 1.0);
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
    return fam_search_ppf(q, 0.0, df + nc + 1.0) { v => ncx2_cdf_value(v, df, nc) };
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

fn ncf_pdf(x: f64, dfn: f64, dfd: f64, nc: f64) -> f64, i64 {
    let k = tl_px(x, ncf_check(dfn, dfd, nc));
    if k != OK { return nan(), k; }
    if x < 0.0 { return 0.0, OK; }
    let lam = 0.5 * nc;
    let lo, hi = nc_range(lam);
    var s: f64 = 0.0;
    var j = lo;
    while j <= hi {
        let r = dfn / (dfn + 2.0 * j);
        let p, cp = fdist_pdf(x * r, dfn + 2.0 * j, dfd);
        if cp == OK { s += nc_weight(j, lam) * r * p; }
        j += 1.0;
    }
    return tl_out(s);
}

fn ncf_cdf_value(x: f64, dfn: f64, dfd: f64, nc: f64) -> f64 {
    if x <= 0.0 { return 0.0; }
    let lam = 0.5 * nc;
    let lo, hi = nc_range(lam);
    let u = dfn * x / (dfn * x + dfd);
    var s: f64 = 0.0;
    var j = lo;
    while j <= hi {
        let p, cp = betainc(0.5 * dfn + j, 0.5 * dfd, u);
        if cp == OK { s += nc_weight(j, lam) * p; }
        j += 1.0;
    }
    return min(s, 1.0);
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
    return fam_search_ppf(q, 0.0, 2.0) { v => ncf_cdf_value(v, dfn, dfd, nc) };
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
    for i in 0..15 {
        let e = 1.0 + w * zs[i];
        if e > edges[edges.len - 1] { edges.push(e); }
    }
}

// The log density of the scale s = sqrt(chi2 / df).
fn nct_lscale(s: f64, df: f64) -> f64 {
    return log(2.0) + 0.5 * df * log(0.5 * df) - fam_lngamma(0.5 * df) + (df - 1.0) * log(s) - 0.5 * df * s * s;
}

// The cumulative function at t, as an integral over the scale.
fn nct_cdf_value(t: f64, df: f64, nc: f64) -> f64 {
    var edges: f64[>..] = [];
    nct_edges(df, edges);
    var total: f64 = 0.0;
    for i in 0..edges.len - 1 {
        total += fam_quad(edges[i], edges[i + 1]) { s => tl_ncdf(t * s - nc) * exp(nct_lscale(s, df)) };
    }
    return min(max(total, 0.0), 1.0);
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
    return tl_ppf_line(q, nc) { v => nct_cdf_value(v, df, nc) };
}

// ---- ksone: the one sided Kolmogorov-Smirnov statistic --------------------

fn ks_check(n: f64) -> i64 {
    let c = check_scalar(n);
    if c != OK { return c; }
    if n < 1.0 || n != floor(n) || n > KS_MAX_N { return ERR_DOMAIN; }
    return OK;
}

// The survival function of the sum of Birnbaum and Tingey, and the
// density as its negative derivative. The point lies from 0 to 1.
fn ksone_both(x: f64, n: f64) -> f64, f64 {
    var sf: f64 = 0.0;
    var dd: f64 = 0.0;
    let jmax = floor(n * (1.0 - x));
    var j: f64 = 0.0;
    while j <= jmax {
        let lc, cc = lchoose(n, j);
        if cc != OK { return nan(), nan(); }
        let a = x + j / n;
        let b = 1.0 - x - j / n;
        if j == 0.0 {
            sf += pow(1.0 - x, n);
            dd += -n * pow(1.0 - x, n - 1.0);
        } else if b > 0.0 {
            let l = exp(lc + (j - 1.0) * log(a) + (n - j) * log(b));
            sf += x * l;
            dd += l * (1.0 + x * (j - 1.0) / a - x * (n - j) / b);
        } else if n - j == 0.0 {
            let l = exp(lc + (j - 1.0) * log(a));
            sf += x * l;
            dd += l * (1.0 + x * (j - 1.0) / a);
        }
        j += 1.0;
    }
    return sf, -dd;
}

fn ksone_pdf(x: f64, n: f64) -> f64, i64 {
    let k = tl_px(x, ks_check(n));
    if k != OK { return nan(), k; }
    if x <= 0.0 || x >= 1.0 { return 0.0, OK; }
    let sf, d = ksone_both(x, n);
    return tl_out(max(d, 0.0));
}

fn ksone_cdf_value(x: f64, n: f64) -> f64 {
    if x <= 0.0 { return 0.0; }
    if x >= 1.0 { return 1.0; }
    let sf, d = ksone_both(x, n);
    return min(max(1.0 - sf, 0.0), 1.0);
}

fn ksone_cdf(x: f64, n: f64) -> f64, i64 {
    let k = tl_px(x, ks_check(n));
    if k != OK { return nan(), k; }
    return tl_out(ksone_cdf_value(x, n));
}

fn ksone_ppf(q: f64, n: f64) -> f64, i64 {
    let k = tl_pq(q, ks_check(n));
    if k != OK { return nan(), k; }
    if q == 0.0 { return 0.0, OK; }
    if q == 1.0 { return 1.0, OK; }
    let x, cx = brent_root(0.0, 1.0, 1e-14, 4e-16) { v => ksone_cdf_value(v, n) - q };
    return x, cx;
}

// ---- kstwo: the two sided Kolmogorov-Smirnov statistic --------------------

// The product c = a * b of two square matrices of the side m.
fn ks_mul(a: const f64[:], b: const f64[:], c: f64[>..]&, m: i64) {
    for i in 0..m {
        for j in 0..m {
            var s: f64 = 0.0;
            for l in 0..m { s += a[i * m + l] * b[l * m + j]; }
            c[i * m + j] = s;
        }
    }
}

// Scale a matrix by 1e-140 when its largest entry grows too large. The
// call returns the added power of ten.
fn ks_rescale(v: f64[>..]&) -> i64 {
    var top: f64 = 0.0;
    for i in 0..v.len { top = max(top, abs(v[i])); }
    if top < 1e140 { return 0; }
    for i in 0..v.len { v[i] = v[i] * 1e-140; }
    return 140;
}

// The cumulative function of the statistic, by the matrix method.
fn kstwo_cdf_value(n: f64, d: f64) -> f64 {
    if d <= 0.5 / n { return 0.0; }
    if d >= 1.0 { return 1.0; }
    let s = d * d * n;
    if s > 7.24 || (s > 3.76 && n > 99.0) {
        return 1.0 - 2.0 * exp(-(2.000071 + 0.331 / sqrt(n) + 1.409 / n) * s);
    }
    let k = floor(n * d) as i64 + 1;
    let m = 2 * k - 1;
    let h = k as f64 - n * d;
    var hm: f64[>..] = [];
    for i in 0..m {
        for j in 0..m {
            if i - j + 1 < 0 { hm.push(0.0); } else { hm.push(1.0); }
        }
    }
    for i in 0..m {
        hm[i * m] = hm[i * m] - pow(h, (i + 1) as f64);
        hm[(m - 1) * m + i] = hm[(m - 1) * m + i] - pow(h, (m - i) as f64);
    }
    if 2.0 * h - 1.0 > 0.0 { hm[(m - 1) * m] = hm[(m - 1) * m] + pow(2.0 * h - 1.0, m as f64); }
    for i in 0..m {
        for j in 0..m {
            if i - j + 1 > 0 {
                for g in 1..(i - j + 2) { hm[i * m + j] = hm[i * m + j] / g as f64; }
            }
        }
    }
    // The power of the matrix by repeated squaring, with a power of ten
    // kept apart.
    var r: f64[>..] = [];
    var p: f64[>..] = [];
    var t: f64[>..] = [];
    for i in 0..m * m {
        r.push(0.0);
        t.push(0.0);
        p.push(hm[i]);
    }
    for i in 0..m { r[i * m + i] = 1.0; }
    var er: i64 = 0;
    var ep: i64 = 0;
    var e = n as i64;
    while e > 0 {
        if e % 2 == 1 {
            ks_mul(r[..], p[..], t, m);
            for i in 0..m * m { r[i] = t[i]; }
            er += ep;
            er += ks_rescale(r);
        }
        e = e / 2;
        if e > 0 {
            ks_mul(p[..], p[..], t, m);
            for i in 0..m * m { p[i] = t[i]; }
            ep = ep * 2;
            ep += ks_rescale(p);
        }
    }
    var v = r[(k - 1) * m + k - 1];
    let ni = n as i64;
    for i in 1..(ni + 1) {
        v = v * i as f64 / n;
        if v < 1e-140 {
            v = v * 1e140;
            er -= 140;
        }
    }
    return min(max(v * pow(10.0, er as f64), 0.0), 1.0);
}

fn kstwo_cdf(x: f64, n: f64) -> f64, i64 {
    let k = tl_px(x, ks_check(n));
    if k != OK { return nan(), k; }
    return tl_out(kstwo_cdf_value(n, x));
}

// The density is the derivative of the cumulative function. The call
// takes a central difference.
fn kstwo_pdf(x: f64, n: f64) -> f64, i64 {
    let k = tl_px(x, ks_check(n));
    if k != OK { return nan(), k; }
    if x <= 0.0 || x >= 1.0 { return 0.0, OK; }
    let h = 1e-6;
    let hi = min(x + h, 1.0);
    let lo = max(x - h, 0.0);
    return tl_out(max((kstwo_cdf_value(n, hi) - kstwo_cdf_value(n, lo)) / (hi - lo), 0.0));
}

fn kstwo_ppf(q: f64, n: f64) -> f64, i64 {
    let k = tl_pq(q, ks_check(n));
    if k != OK { return nan(), k; }
    if q == 0.0 { return 0.5 / n, OK; }
    if q == 1.0 { return 1.0, OK; }
    let x, cx = brent_root(0.5 / n, 1.0, 1e-14, 4e-16) { v => kstwo_cdf_value(n, v) - q };
    return x, cx;
}
