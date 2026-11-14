// The histogram law of SciPy 1.18.1: rv_histogram.
//
// The law holds a histogram. The call takes the content of n bins and the
// n+1 bin edges. A flag says whether the content is a density or a count.
// A count becomes a density by the bin width. The density then scales so
// that the whole law has a mass of 1. The density is constant inside a bin
// and 0 outside the edges. The cumulative function is linear inside a bin.
//
// SciPy treats a missing flag as a density and warns when the bin widths
// vary. This package has no warning, so the flag has no default.
//
// Every call returns the result first and a fault code last.
namespace statistics;

import std;
import math;

// The check of the bins. The edges must rise. The content must be at
// least 0 and finite, and its total mass must be above 0. The check
// gives the total mass of the content.
fn hist_check(content: const f64[:], edges: const f64[:], density: bool, total: f64&) -> i64 {
    let n = content.len;
    if n == 0 { return ERR_EMPTY; }
    if n > MAX_BINS { return ERR_TOO_LONG; }
    if edges.len != n + 1 { return ERR_LENGTH; }
    for v in edges {
        if is_nan(v) { return ERR_NAN; }
        if is_inf(v) { return ERR_INFINITY; }
    }
    for i in 0..n {
        if edges[i + 1] <= edges[i] { return ERR_DOMAIN; }
        if is_nan(content[i]) { return ERR_NAN; }
        if is_inf(content[i]) { return ERR_INFINITY; }
        if content[i] < 0.0 { return ERR_DOMAIN; }
    }
    var sum: f64 = 0.0;
    for i in 0..n {
        var h = content[i];
        if !density { h = h / (edges[i + 1] - edges[i]); }
        sum = sum + h * (edges[i + 1] - edges[i]);
    }
    if !(sum > 0.0) || is_inf(sum) { return ERR_DOMAIN; }
    total = sum;
    return OK;
}

// The height of every bin after the scale to a mass of 1.
fn hist_heights(out: f64[>..]&, content: const f64[:], edges: const f64[:], density: bool, total: f64) {
    out.resize(0);
    for i in 0..content.len {
        var h = content[i];
        if !density { h = h / (edges[i + 1] - edges[i]); }
        out.push(h / total);
    }
}

// The bin that holds x, for x inside the edges. The last bin holds its
// upper edge, as SciPy does through its right-hand search.
fn hist_bin(edges: const f64[:], x: f64) -> i64 {
    var lo: i64 = 0;
    var hi: i64 = edges.len - 1;
    while hi - lo > 1 {
        let m = (lo + hi) / 2;
        if edges[m] <= x { lo = m; } else { hi = m; }
    }
    return lo;
}

fn histogram_pdf(content: const f64[:], edges: const f64[:], density: bool, x: f64) -> f64, i64 {
    var total: f64 = 0.0;
    let c = hist_check(content, edges, density, total);
    if c != OK { return nan(), c; }
    if is_nan(x) { return nan(), ERR_NAN; }
    if x < edges[0] || x >= edges[edges.len - 1] { return 0.0, OK; }
    var h: f64[>..] = [];
    hist_heights(h, content, edges, density, total);
    return h[hist_bin(edges, x)], OK;
}

fn histogram_cdf(content: const f64[:], edges: const f64[:], density: bool, x: f64) -> f64, i64 {
    var total: f64 = 0.0;
    let c = hist_check(content, edges, density, total);
    if c != OK { return nan(), c; }
    if is_nan(x) { return nan(), ERR_NAN; }
    if x <= edges[0] { return 0.0, OK; }
    if x >= edges[edges.len - 1] { return 1.0, OK; }
    var h: f64[>..] = [];
    hist_heights(h, content, edges, density, total);
    let k = hist_bin(edges, x);
    var s: f64 = 0.0;
    for i in 0..k { s = s + h[i] * (edges[i + 1] - edges[i]); }
    s = s + h[k] * (x - edges[k]);
    return min(s, 1.0), OK;
}

fn histogram_sf(content: const f64[:], edges: const f64[:], density: bool, x: f64) -> f64, i64 {
    let v, c = histogram_cdf(content, edges, density, x);
    if c != OK { return nan(), c; }
    return 1.0 - v, OK;
}

// The quantile. The cumulative function is linear in a bin, so the
// quantile solves in closed form inside the bin that holds the
// probability. A flat bin of height 0 is skipped.
fn histogram_ppf(content: const f64[:], edges: const f64[:], density: bool, p: f64) -> f64, i64 {
    var total: f64 = 0.0;
    let c = hist_check(content, edges, density, total);
    if c != OK { return nan(), c; }
    let cp = check_prob(p);
    if cp != OK { return nan(), cp; }
    let n = content.len;
    if p == 0.0 { return edges[0], OK; }
    if p == 1.0 { return edges[n], OK; }
    var h: f64[>..] = [];
    hist_heights(h, content, edges, density, total);
    var s: f64 = 0.0;
    for i in 0..n {
        let m = h[i] * (edges[i + 1] - edges[i]);
        if h[i] > 0.0 && s + m >= p {
            return edges[i] + (p - s) / h[i], OK;
        }
        s = s + m;
    }
    return edges[n], OK;
}

// The raw moment of a whole order, by the exact integral of each bin.
fn histogram_moment(content: const f64[:], edges: const f64[:], density: bool, order: i64) -> f64, i64 {
    var total: f64 = 0.0;
    let c = hist_check(content, edges, density, total);
    if c != OK { return nan(), c; }
    if order < 0 || order > 64 { return nan(), ERR_DOMAIN; }
    var h: f64[>..] = [];
    hist_heights(h, content, edges, density, total);
    let k = (order + 1) as f64;
    var s: f64 = 0.0;
    for i in 0..content.len {
        s = s + h[i] * (pow(edges[i + 1], k) - pow(edges[i], k)) / k;
    }
    return s, OK;
}

fn histogram_mean(content: const f64[:], edges: const f64[:], density: bool) -> f64, i64 {
    return histogram_moment(content, edges, density, 1);
}

fn histogram_var(content: const f64[:], edges: const f64[:], density: bool) -> f64, i64 {
    let m1, c1 = histogram_moment(content, edges, density, 1);
    if c1 != OK { return nan(), c1; }
    let m2, c2 = histogram_moment(content, edges, density, 2);
    if c2 != OK { return nan(), c2; }
    return m2 - m1 * m1, OK;
}

// The entropy of the density.
fn histogram_entropy(content: const f64[:], edges: const f64[:], density: bool) -> f64, i64 {
    var total: f64 = 0.0;
    let c = hist_check(content, edges, density, total);
    if c != OK { return nan(), c; }
    var h: f64[>..] = [];
    hist_heights(h, content, edges, density, total);
    var s: f64 = 0.0;
    for i in 0..content.len {
        if h[i] > 0.0 { s = s - h[i] * log(h[i]) * (edges[i + 1] - edges[i]); }
    }
    return s, OK;
}

// The draws of the law. A draw picks a bin by its mass and then a point
// inside the bin. The buffer holds size values.
fn histogram_rvs(out: f64[>..]&, content: const f64[:], edges: const f64[:], density: bool, size: i64, r: rng&) -> i64 {
    out.resize(0);
    var total: f64 = 0.0;
    let c = hist_check(content, edges, density, total);
    if c != OK { return c; }
    if size < 0 || size > 100000000 { return ERR_DOMAIN; }
    var h: f64[>..] = [];
    hist_heights(h, content, edges, density, total);
    for s in 0..size {
        let u = rs_unit(r);
        var acc: f64 = 0.0;
        var x: f64 = edges[edges.len - 1];
        for i in 0..h.len {
            let m = h[i] * (edges[i + 1] - edges[i]);
            if acc + m >= u {
                x = edges[i] + (u - acc) / h[i];
                break;
            }
            acc = acc + m;
        }
        out.push(x);
    }
    return OK;
}
