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
const FAM_GOLD: i64 = 40;      // the steps of the search for the peak
const FAM_BISECT: i64 = 40;    // the steps of the bisection of a border

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
// point. The scale gives the point of the law from the logarithm. A scale
// of one gives the point itself. The count gives the highest order that
// the rule computes.
fn fam_raws<F>(count: i64, scale: f64) -> f64, f64, f64, f64, i64 {
    var gx: f64[>..] = [];
    var gw: f64[>..] = [];
    gauss_legendre(gx, gw, FAM_NODES);
    var m1: f64 = 0.0;
    var m2: f64 = 0.0;
    var m3: f64 = 0.0;
    var m4: f64 = 0.0;
    var code: i64 = OK;
    for n in 1..count + 1 {
        // The scan of the logarithm of the weighted density. The scan
        // widens when the peak or a border lies at its edge.
        var ulo: f64 = FAM_ULO;
        var uhi: f64 = FAM_UHI;
        var lo: f64 = 0.0;
        var hi: f64 = 0.0;
        var lstar: f64 = 0.0;
        var done: i64 = 0;
        var tries: i64 = 0;
        while tries < FAM_TRIES && done == 0 {
            tries += 1;
            // The scan of the logarithm of the weighted density. The point
            // of the peak gives the scale of the integral.
            var peak: f64 = -inf();
            var peak_index: i64 = 0;
            for k in 0..FAM_SCAN + 1 {
                let u = ulo + (uhi - ulo) * k as f64 / FAM_SCAN as f64;
                let l = fam_logw(n, scale, u) { x => F(x) };
                if l > peak { peak = l; peak_index = k; }
            }
            // The scan found no density. The law lies outside the reach of
            // the rule. The call gives the fit fault.
            if is_inf(peak) || is_nan(peak) { code = ERR_FIT; break; }
            // The peak lies at the edge of the scan. The scan widens.
            if peak_index == 0 || peak_index == FAM_SCAN {
                if peak_index == 0 && ulo > -690.0 { ulo = ulo * 2.0; }
                else if uhi < 690.0 { uhi = uhi * 2.0; }
                continue;
            }
            let step = (uhi - ulo) / FAM_SCAN as f64;
            // The golden section search for the peak. The bracket holds the
            // two points of the scan beside the peak.
            var ga = ulo + step * (peak_index - 1) as f64;
            var gb = ulo + step * (peak_index + 1) as f64;
            let phi: f64 = 0.6180339887498949;
            var x1 = gb - phi * (gb - ga);
            var x2 = ga + phi * (gb - ga);
            var l1 = fam_logw(n, scale, x1) { x => F(x) };
            var l2 = fam_logw(n, scale, x2) { x => F(x) };
            for it in 0..FAM_GOLD {
                if l1 < l2 {
                    ga = x1;
                    x1 = x2;
                    l1 = l2;
                    x2 = ga + phi * (gb - ga);
                    l2 = fam_logw(n, scale, x2) { x => F(x) };
                } else {
                    gb = x2;
                    x2 = x1;
                    l2 = l1;
                    x1 = gb - phi * (gb - ga);
                    l1 = fam_logw(n, scale, x1) { x => F(x) };
                }
            }
            let ustar = (ga + gb) / 2.0;
            lstar = fam_logw(n, scale, ustar) { x => F(x) };
            // The range where the weighted density is above the cut.
            let thr = lstar + log(FAM_TAIL);
            // The border on the left. The walk gives the point of the scan
            // below the cut.
            var li: i64 = peak_index;
            while li >= 0 {
                let l = fam_logw(n, scale, ulo + step * li as f64) { x => F(x) };
                if l <= thr { break; }
                li -= 1;
            }
            // The border on the right.
            var ri: i64 = peak_index;
            while ri <= FAM_SCAN {
                let l = fam_logw(n, scale, ulo + step * ri as f64) { x => F(x) };
                if l <= thr { break; }
                ri += 1;
            }
            // A border lies outside the scan. The law extends past the
            // reach of the rule. The scan widens.
            if li < 0 || ri > FAM_SCAN {
                if li < 0 && ulo > -690.0 { ulo = ulo * 2.0; }
                else if uhi < 690.0 { uhi = uhi * 2.0; }
                continue;
            }
            // The bisection of the border on the left.
            var lref = if li + 1 > peak_index { ustar } else { ulo + step * (li + 1) as f64 };
            var la = ulo + step * li as f64;
            var lb = lref;
            for it in 0..FAM_BISECT {
                let mid = (la + lb) / 2.0;
                let lm = fam_logw(n, scale, mid) { x => F(x) };
                if lm > thr { lb = mid; } else { la = mid; }
            }
            lo = (la + lb) / 2.0;
            // The bisection of the border on the right.
            var rref = if ri - 1 < peak_index { ustar } else { ulo + step * (ri - 1) as f64 };
            var ra = ulo + step * ri as f64;
            var rb = rref;
            for it in 0..FAM_BISECT {
                let mid = (ra + rb) / 2.0;
                let lm = fam_logw(n, scale, mid) { x => F(x) };
                if lm > thr { rb = mid; } else { ra = mid; }
            }
            hi = (ra + rb) / 2.0;
            done = 1;
        }
        // The scan found no range. The law lies outside the reach of the
        // rule. The call gives the fit fault.
        if done == 0 { code = ERR_FIT; continue; }
        // The integral of the weighted density, scaled by the peak.
        let width = (hi - lo) / FAM_PANELS as f64;
        var total: f64 = 0.0;
        for k in 0..FAM_PANELS {
            let left = lo + width * k as f64;
            for j in 0..FAM_NODES {
                let u = left + width * (gx[j] + 1.0) / 2.0;
                let l = fam_logw(n, scale, u) { x => F(x) };
                if !is_nan(l) && !is_inf(l) {
                    total += exp(l - lstar) * gw[j] * width / 2.0;
                }
            }
        }
        let v = exp(lstar) * total;
        if n == 1 { m1 = v; }
        else if n == 2 { m2 = v; }
        else if n == 3 { m3 = v; }
        else { m4 = v; }
    }
    return m1, m2, m3, m4, code;
}

// The raw moments of a law on the whole line. The density of the law is
// the trailing block of the call. The rule works in the point itself. The
// count gives the highest order that the rule computes.
fn fam_raws_line<F>(count: i64) -> f64, f64, f64, f64, i64 {
    var gx: f64[>..] = [];
    var gw: f64[>..] = [];
    gauss_legendre(gx, gw, FAM_NODES);
    var m1: f64 = 0.0;
    var m2: f64 = 0.0;
    var m3: f64 = 0.0;
    var m4: f64 = 0.0;
    var code: i64 = OK;
    for n in 1..count + 1 {
        // The integral of the weighted density on each side of the peak.
        var total: f64 = 0.0;
        // The largest density on the range. The density at the edge of the
        // range is compared with it.
        var top: f64 = 0.0;
        for side in 0..2 {
            let left = if side == 0 { -FAM_REACH } else { 0.0 };
            let right = if side == 0 { 0.0 } else { FAM_REACH };
            let width = (right - left) / FAM_PANELS as f64;
            for k in 0..FAM_PANELS {
                let base = left + width * k as f64;
                for j in 0..FAM_NODES {
                    let x = base + width * (gx[j] + 1.0) / 2.0;
                    let d = F(x);
                    if !is_nan(d) && !is_inf(d) {
                        total += d * pow(x, n as f64) * gw[j] * width / 2.0;
                        if d > top { top = d; }
                    }
                }
            }
        }
        // The density at the edge of the range decides if the law extends
        // past the reach of the rule. The call gives the fit fault then.
        let e1 = F(FAM_REACH);
        let e2 = F(-FAM_REACH);
        if is_nan(e1) || is_nan(e2) || is_inf(e1) || is_inf(e2) { code = ERR_FIT; continue; }
        if e1 > top * 1e-9 || e2 > top * 1e-9 { code = ERR_FIT; continue; }
        if n == 1 { m1 = total; }
        else if n == 2 { m2 = total; }
        else if n == 3 { m3 = total; }
        else { m4 = total; }
    }
    return m1, m2, m3, m4, code;
}
