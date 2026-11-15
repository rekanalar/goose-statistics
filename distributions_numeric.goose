// The moments of the families whose law has no closed form.
//
// The rule is a Gauss-Legendre rule on panels. A law on 0 to inf is
// integrated in the logarithm of the point. The rule finds the peak of the
// weighted density first. It cuts the range where the weighted density
// falls below a small part of that peak. A law with a long tail keeps the
// accuracy of the rule.
//
// The scan for the peak widens when the peak lies at its edge. The peak
// then gets a golden section search and the borders of the range get a
// bisection. These steps let the rule work on a law with a narrow peak and
// on a law with a peak far from the origin.
//
// A law on the whole line is integrated in the point itself. The range
// reaches 30 on each side of the peak.
//
// The names and the parameters of the families follow SciPy.
namespace statistics;

import std;
import math;

const FAM_NODES: i64 = 16;    // the nodes of the rule on one panel
const FAM_PANELS: i64 = 100;   // the panels on the range of the integral
const FAM_SCAN: i64 = 400;     // the points of the scan for the peak
const FAM_TAIL: f64 = 1e-15;   // the part of the range that is cut
const FAM_ULO: f64 = -30.0;    // the lowest point of the scan, in logs
const FAM_UHI: f64 = 60.0;     // the highest point of the scan, in logs
const FAM_REACH: f64 = 30.0;   // the reach of a law on the whole line
const FAM_TRIES: i64 = 8;      // the widenings of the scan for the peak
const FAM_RAW_PANELS: i64 = 50; // the panels of the grid of the raw moments
const FAM_LINE_SCAN: i64 = 64;  // the points of the scan for the range of the line
const FAM_LINE_PANELS: i64 = 60; // the panels of the grid of the line

// The logarithm of one plus the exponential of a value. The form with
// the maximum of the two parts keeps the digits on each side. The direct
// form overflows for a large value.
fn fam_log1pexp(t: f64) -> f64, i64 {
    if t > 0.0 {
        let v, cv = log1p(exp(-t));
        if cv != OK { return nan(), cv; }
        return t + v, OK;
    }
    let v, cv = log1p(exp(t));
    if cv != OK { return nan(), cv; }
    return v, OK;
}

// The check of a result that the package computed. A result that is a
// NaN or an infinity gives the overflow fault. The code never returns
// such a value with the code OK.
fn fam_finite(v: f64) -> f64, i64 {
    if is_nan(v) { return nan(), ERR_OVERFLOW; }
    if is_inf(v) { return nan(), ERR_OVERFLOW; }
    return v, OK;
}

// The variance of a law from its first two raw moments.
fn fam_var(m1: f64, m2: f64) -> f64 {
    return m2 - m1 * m1;
}

// The skewness of a law from its first three raw moments.
fn fam_skew(m1: f64, m2: f64, m3: f64) -> f64 {
    let v = fam_var(m1, m2);
    return (m3 - 3.0 * m1 * m2 + 2.0 * m1 * m1 * m1) / (v * sqrt(v));
}

// The excess kurtosis of a law from its four raw moments.
fn fam_kurt(m1: f64, m2: f64, m3: f64, m4: f64) -> f64 {
    let v = fam_var(m1, m2);
    let m1s = m1 * m1;
    return (m4 - 4.0 * m1 * m3 + 6.0 * m1s * m2 - 3.0 * m1s * m1s) / (v * v) - 3.0;
}

// The logarithm of the weighted density of a law at the point u in the
// logarithm of the point. The scale gives the point of the law from u. The
// weight is the point to the power n. A point outside the range of a
// double gives the logarithm of zero.
fn fam_logw<F>(n: i64, scale: f64, u: f64) -> f64 {
    let x = scale * u;
    if x >= 700.0 { return -inf(); }
    let d = F(exp(x));
    if is_nan(d) { return -inf(); }
    if d <= 0.0 { return -inf(); }
    let l = log(scale) + (n + 1) as f64 * x + log(d);
    if is_nan(l) { return -inf(); }
    return l;
}

// The raw moments of a law on 0 to inf. The density of the law is the
// trailing block of the call. The rule works in the logarithm of the
// point. The scale gives the point of the law from the logarithm. A
// scale of one gives the point itself. The count gives the highest
// order that the rule computes.
//
// The weighted density of the order n is the weighted density of the
// order 1 plus (n - 1) times the point, so one scan and one grid serve
// every order. The scan holds the values of the first order, and the
// value of another order follows from it with one addition. The peak of
// an order comes from the scan, and the border of an order is the point
// where its weighted density falls under the cut. The grid runs over the
// range of the first order to the range of the last order, and the
// weighted density of the middle orders stays under the cut outside
// that range.
fn fam_raws<F>(count: i64, scale: f64) -> f64, f64, f64, f64, i64 {
    var gx: f64[>..] = [];
    var gw: f64[>..] = [];
    gauss_legendre(gx, gw, FAM_NODES);
    var m1: f64 = 0.0;
    var m2: f64 = 0.0;
    var m3: f64 = 0.0;
    var m4: f64 = 0.0;
    var ulo: f64 = FAM_ULO;
    var uhi: f64 = FAM_UHI;
    var done: i64 = 0;
    var tries: i64 = 0;
    var lo: f64 = 0.0;
    var hi: f64 = 0.0;
    var vals: f64[>..] = [];
    var stars: f64[>..] = [];
    while tries < FAM_TRIES && done == 0 {
        tries += 1;
        let step = (uhi - ulo) / FAM_SCAN as f64;
        vals.resize(0);
        // The scan of the weighted density of the first order.
        var pk1: f64 = -inf();
        var i1: i64 = 0;
        for k in 0..FAM_SCAN + 1 {
            let u = ulo + step * k as f64;
            let l1 = fam_logw(1, scale, u) { x => F(x) };
            vals.push(l1);
            if l1 > pk1 { pk1 = l1; i1 = k; }
        }
        // The scan found no density. The law lies outside the reach of
        // the rule. The call gives the fit fault.
        if is_inf(pk1) || is_nan(pk1) { break; }
        // The peak of every order, from the shared values.
        stars.resize(0);
        var edge_lo: bool = i1 == 0;
        var edge_hi: bool = i1 == FAM_SCAN;
        var ilast: i64 = i1;
        for n in 1..count + 1 {
            let w = (n - 1) as f64 * scale;
            var pk: f64 = -inf();
            var ik: i64 = 0;
            for k in 0..FAM_SCAN + 1 {
                let l = vals[k] + w * (ulo + step * k as f64);
                if l > pk { pk = l; ik = k; }
            }
            stars.push(pk);
            if ik == 0 { edge_lo = true; }
            if ik == FAM_SCAN { edge_hi = true; }
            ilast = ik;
        }
        // A peak lies at the edge of the scan. The scan widens.
        if edge_lo {
            if ulo > -690.0 { ulo = ulo * 2.0; continue; }
            if edge_hi && uhi < 690.0 { uhi = uhi * 2.0; continue; }
            break;
        }
        if edge_hi {
            if uhi < 690.0 { uhi = uhi * 2.0; continue; }
            break;
        }
        // The border of the first order on the left and the border of the
        // last order on the right, from the shared values. The walks run
        // past the ends of the scan, so a border outside the scan widens
        // it.
        let cut1 = stars[0] + log(FAM_TAIL);
        var li: i64 = i1;
        while li >= 0 {
            if vals[li] <= cut1 { break; }
            li -= 1;
        }
        let wc = (count - 1) as f64 * scale;
        let cutc = stars[count - 1] + log(FAM_TAIL);
        var ri: i64 = ilast;
        while ri <= FAM_SCAN {
            if vals[ri] + wc * (ulo + step * ri as f64) <= cutc { break; }
            ri += 1;
        }
        // A border lies outside the scan. The scan widens.
        if li < 0 || ri > FAM_SCAN {
            if li < 0 && ulo > -690.0 { ulo = ulo * 2.0; continue; }
            if uhi < 690.0 { uhi = uhi * 2.0; continue; }
            break;
        }
        lo = ulo + step * li as f64;
        hi = ulo + step * ri as f64;
        done = 1;
    }
    // The scan found no range. The law lies outside the reach of the
    // rule. The call gives the fit fault.
    if done == 0 {
        return m1, m2, m3, m4, ERR_FIT;
    }
    // The grid of the union range. The value of the first order comes
    // once per node, and the value of the order n follows with one
    // addition.
    let width = (hi - lo) / FAM_RAW_PANELS as f64;
    var tot: f64[>..] = [];
    for n in 1..count + 1 { tot.push(0.0); }
    for k in 0..FAM_RAW_PANELS {
        let left = lo + width * k as f64;
        for j in 0..FAM_NODES {
            let u = left + width * (gx[j] + 1.0) / 2.0;
            let l1 = fam_logw(1, scale, u) { x => F(x) };
            if !is_nan(l1) && !is_inf(l1) {
                let wj = gw[j] * width / 2.0;
                for n in 1..count + 1 {
                    let l = l1 + (n - 1) as f64 * scale * u - stars[n - 1];
                    if l > -700.0 {
                        tot[n - 1] += exp(l) * wj;
                    }
                }
            }
        }
    }
    for n in 1..count + 1 {
        let v = exp(stars[n - 1]) * tot[n - 1];
        if n == 1 { m1 = v; }
        else if n == 2 { m2 = v; }
        else if n == 3 { m3 = v; }
        else { m4 = v; }
    }
    return m1, m2, m3, m4, OK;
}

// The raw moments of a law on the whole line. The density of the law is
// the trailing block of the call. The rule works in the point itself. The
// count gives the highest order that the rule computes.
//
// A scan finds the range where the density carries weight, and one pass
// over the grid of that range serves every order: the value of the
// density at a node comes once, and the four powers of the node come
// from multiplications.
fn fam_raws_line<F>(count: i64) -> f64, f64, f64, f64, i64 {
    var gx: f64[>..] = [];
    var gw: f64[>..] = [];
    gauss_legendre(gx, gw, FAM_NODES);
    var m1: f64 = 0.0;
    var m2: f64 = 0.0;
    var m3: f64 = 0.0;
    var m4: f64 = 0.0;
    // The scan of the density over the reach. The values at the ends
    // serve the check of the reach.
    var vals: f64[>..] = [];
    var top: f64 = 0.0;
    let wscan = 2.0 * FAM_REACH / FAM_LINE_SCAN as f64;
    for k in 0..FAM_LINE_SCAN + 1 {
        let d = F(-FAM_REACH + wscan * k as f64);
        vals.push(d);
        if !is_nan(d) && !is_inf(d) && d > top { top = d; }
    }
    let e1 = vals[FAM_LINE_SCAN];
    let e2 = vals[0];
    // The density at the edge of the range decides if the law extends
    // past the reach of the rule. The call gives the fit fault then.
    if is_nan(e1) || is_nan(e2) || is_inf(e1) || is_inf(e2) { return m1, m2, m3, m4, ERR_FIT; }
    if e1 > top * 1e-9 || e2 > top * 1e-9 { return m1, m2, m3, m4, ERR_FIT; }
    // The range where the density is above the cut, with one step of
    // margin on both sides. A density below the cut everywhere keeps the
    // whole reach.
    let cut = top * 1e-13;
    var i0: i64 = 0;
    while i0 < FAM_LINE_SCAN && vals[i0] < cut { i0 += 1; }
    var i1: i64 = FAM_LINE_SCAN;
    while i1 > 0 && vals[i1] < cut { i1 -= 1; }
    if i0 > 0 { i0 -= 1; }
    if i1 < FAM_LINE_SCAN { i1 += 1; }
    var lo2: f64 = -FAM_REACH + wscan * i0 as f64;
    var hi2: f64 = -FAM_REACH + wscan * i1 as f64;
    if i0 >= i1 { lo2 = -FAM_REACH; hi2 = FAM_REACH; }
    // The single pass over the grid. One value of the density serves
    // every order.
    let width = (hi2 - lo2) / FAM_LINE_PANELS as f64;
    for k in 0..FAM_LINE_PANELS {
        let base = lo2 + width * k as f64;
        for j in 0..FAM_NODES {
            let x = base + width * (gx[j] + 1.0) / 2.0;
            let d = F(x);
            if !is_nan(d) && !is_inf(d) {
                let wj = d * gw[j] * width / 2.0;
                let x2 = x * x;
                if count >= 1 { m1 += wj * x; }
                if count >= 2 { m2 += wj * x2; }
                if count >= 3 { m3 += wj * x2 * x; }
                if count >= 4 { m4 += wj * x2 * x2; }
            }
        }
    }
    return m1, m2, m3, m4, OK;
}
