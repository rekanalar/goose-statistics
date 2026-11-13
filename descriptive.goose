// The descriptive helpers of a sample: the percentiles, the ranks of a
// score, the robust scale, the power means, the geometric scale, the
// expectiles, the L-moments and the z-scores.
//
// Every call returns the result first and a fault code last. A call
// sorts a copy of the sample, so the data of the caller never changes.
namespace statistics;

import std;
import math;

// The scale of the median absolute deviation that gives the deviation of
// a normal law. It is the quantile of the standard normal at 0.75.
const MAD_NORMAL: f64 = 0.6744897501960817;

// The smallest power that the linearized form of the power mean covers.
const PMEAN_EPS: f64 = 2e-6;

// The tolerance of the search of an expectile.
const EXPECTILE_TOL: f64 = 1e-14;

// Signature: scoreatpercentile(xs, per) -> value, code
// The percentile `per` of a sample, from 0 to 100. The value is the
// linear interpolation between the two order statistics that hold the
// percentile. A percentile outside 0 to 100 gives ERR_DOMAIN. An empty
// sample gives ERR_EMPTY.
fn scoreatpercentile(xs: const f64[:], per: f64) -> f64, i64 {
    return scoreatpercentile_method(xs, per, INTERP_FRACTION);
}

// Signature: scoreatpercentile_method(xs, per, method) -> value, code
// The percentile `per` of a sample, from 0 to 100. The method is
// INTERP_FRACTION for the linear interpolation, INTERP_LOWER for the
// lower value and INTERP_HIGHER for the higher value. A percentile
// outside 0 to 100 gives ERR_DOMAIN. An unknown method gives
// ERR_DOMAIN.
fn scoreatpercentile_method(xs: const f64[:], per: f64, method: i64) -> f64, i64 {
    let c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    let sc = check_scalar(per);
    if sc != OK { return nan(), sc; }
    if per < 0.0 || per > 100.0 { return nan(), ERR_DOMAIN; }
    if method < INTERP_FRACTION || method > INTERP_HIGHER { return nan(), ERR_DOMAIN; }
    let s = sort_copy(xs);
    return ds_pct_of_sorted(s[..], per, method), OK;
}

// Signature: scoreatpercentile_limit(xs, per, low, high) -> value, code
// The percentile `per` of the values of a sample that lie between low
// and high. The bounds are inclusive. A value outside the bounds leaves
// no value, which gives ERR_UNDEFINED. A low bound above the high bound
// gives ERR_DOMAIN.
fn scoreatpercentile_limit(xs: const f64[:], per: f64, low: f64, high: f64) -> f64, i64 {
    return scoreatpercentile_limit_method(xs, per, low, high, INTERP_FRACTION);
}

// Signature: scoreatpercentile_limit_method(xs, per, low, high, method) -> value, code
// The percentile `per` of the values of a sample that lie between low
// and high. The method is one of the three interpolation methods.
fn scoreatpercentile_limit_method(xs: const f64[:], per: f64, low: f64, high: f64, method: i64) -> f64, i64 {
    let c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    let sc = check_scalar(per);
    if sc != OK { return nan(), sc; }
    if per < 0.0 || per > 100.0 { return nan(), ERR_DOMAIN; }
    let lc = check_scalar(low);
    if lc != OK { return nan(), lc; }
    let hc = check_scalar(high);
    if hc != OK { return nan(), hc; }
    if low > high { return nan(), ERR_DOMAIN; }
    if method < INTERP_FRACTION || method > INTERP_HIGHER { return nan(), ERR_DOMAIN; }
    var keep: f64[>..] = [];
    for x in xs {
        if x >= low && x <= high { keep.push(x); }
    }
    if keep.len == 0 { return nan(), ERR_UNDEFINED; }
    let s = sort_copy(keep[..]);
    return ds_pct_of_sorted(s[..], per, method), OK;
}

// The percentile of a sorted sample. The arguments are checked.
fn ds_pct_of_sorted(s: const f64[:], per: f64, method: i64) -> f64 {
    let n = s.len;
    let v = (per / 100.0) * (n - 1);
    if v <= 0.0 { return s[0]; }
    if v >= (n - 1) { return s[n - 1]; }
    let p = floor(v);
    let j = p as i64;
    if method == INTERP_LOWER { return s[j]; }
    if method == INTERP_HIGHER { return s[j + 1]; }
    return lerp_exact(s[j], s[j + 1], v - p);
}

// Signature: percentileofscore(xs, score, kind) -> value, code
// The share of the values of a sample that lie at or below a score, as
// a percentage. The kind is PCT_RANK for the mean rank of the tied
// values, PCT_WEAK for the share at or below, PCT_STRICT for the share
// below and PCT_MEAN for the mean of the weak and the strict share. A
// NaN score gives ERR_NAN.
fn percentileofscore(xs: const f64[:], score: f64, kind: i64) -> f64, i64 {
    let c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    let sc = check_scalar(score);
    if sc != OK { return nan(), sc; }
    if kind < PCT_RANK || kind > PCT_MEAN { return nan(), ERR_DOMAIN; }
    var left: i64 = 0;
    var right: i64 = 0;
    for x in xs {
        if x < score { left += 1; }
        if x <= score { right += 1; }
    }
    let n = xs.len as f64;
    if kind == PCT_STRICT { return (left) * 100.0 / n, OK; }
    if kind == PCT_WEAK { return (right) * 100.0 / n, OK; }
    if kind == PCT_MEAN { return ((left + right)) * 50.0 / n, OK; }
    var plus: f64 = 0.0;
    if left < right { plus = 1.0; }
    return ((left + right) + plus) * 50.0 / n, OK;
}

// Signature: median_abs_deviation(xs, scale, use_mean) -> value, code
// The median of the absolute distance of the values to their center.
// The center is the median, or the mean when use_mean is true. The
// scale divides the result. A scale of 0 or below gives ERR_DOMAIN.
fn median_abs_deviation(xs: const f64[:], scale: f64, use_mean: bool) -> f64, i64 {
    let c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    let sc = check_scalar(scale);
    if sc != OK { return nan(), sc; }
    if scale <= 0.0 { return nan(), ERR_DOMAIN; }
    var center: f64 = 0.0;
    if use_mean {
        let m, c2 = exact_mean(xs);
        if c2 != 0 { return nan(), c2; }
        center = m;
    } else {
        let m, c2 = median(xs);
        if c2 != 0 { return nan(), c2; }
        center = m;
    }
    var dev: f64[>..] = [];
    for x in xs { dev.push(abs(x - center)); }
    let s = sort_copy(dev[..]);
    let mad = median_sorted(s[..]);
    return mad / scale, OK;
}

// Signature: median_abs_deviation_normal(xs) -> value, code
// The median absolute deviation with the scale of a normal law. The
// result is the deviation of the sample for a normal law.
fn median_abs_deviation_normal(xs: const f64[:]) -> f64, i64 {
    return median_abs_deviation(xs, MAD_NORMAL, false);
}

// Signature: pmean(xs, p) -> value, code
// The power mean of a sample. The exponent 0 gives the geometric mean.
// A value below 0 gives ERR_DOMAIN, as in SciPy, which gives NaN. The
// exponent must be finite.
fn pmean(xs: const f64[:], p: f64) -> f64, i64 {
    let c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    if is_nan(p) { return nan(), ERR_NAN; }
    if is_inf(p) { return nan(), ERR_DOMAIN; }
    for x in xs {
        if x < 0.0 { return nan(), ERR_DOMAIN; }
    }
    if p == 0.0 {
        let g, c2 = geometric_mean(xs);
        if c2 != 0 { return nan(), c2; }
        return g, OK;
    }
    if abs(p) <= PMEAN_EPS {
        let m0, c2 = geometric_mean(xs);
        if c2 != 0 { return nan(), c2; }
        let ln_avg, c3 = exact_mean_log(xs);
        if c3 != 0 { return nan(), c3; }
        let ln2_avg, c4 = exact_mean_log2(xs);
        if c4 != 0 { return nan(), c4; }
        return m0 * (1.0 + 0.5 * p * (ln2_avg - ln_avg * ln_avg)), OK;
    }
    var s: f64[>..] = [];
    for x in xs { s.push(pow(x, p)); }
    let m, c2 = exact_mean(s[..]);
    if c2 != 0 { return nan(), c2; }
    if m < 0.0 { return nan(), ERR_UNDEFINED; }
    let r = pow(m, 1.0 / p);
    if is_bad(r) { return nan(), ERR_OVERFLOW; }
    return r, OK;
}

// Signature: pmean_weighted(xs, ws, p) -> value, code
// The power mean of a sample with weights. A weight below 0 gives
// ERR_DOMAIN. The weights must not all be 0.
fn pmean_weighted(xs: const f64[:], ws: const f64[:], p: f64) -> f64, i64 {
    let c = check_weighted(xs, ws, 1);
    if c != 0 { return nan(), c; }
    if is_nan(p) { return nan(), ERR_NAN; }
    if is_inf(p) { return nan(), ERR_DOMAIN; }
    for x in xs {
        if x < 0.0 { return nan(), ERR_DOMAIN; }
    }
    // The terms of a value with weight 0 stay out of the mean, so a zero
    // value or an overflow there cannot reach the sum.
    var s: f64[>..] = [];
    var kept: f64[>..] = [];
    for i in 0..xs.len {
        if ws[i] == 0.0 { s.push(0.0); kept.push(0.0); continue; }
        var t: f64 = 0.0;
        if p == 0.0 {
            if xs[i] == 0.0 { return nan(), ERR_UNDEFINED; }
            t = log(xs[i]);
        } else {
            t = pow(xs[i], p);
        }
        if is_bad(t) {
            if xs[i] == 0.0 { return nan(), ERR_UNDEFINED; }
            return nan(), ERR_OVERFLOW;
        }
        s.push(t);
        kept.push(ws[i]);
    }
    let m, c2 = exact_mean_weighted(s[..], kept[..]);
    if c2 != 0 { return nan(), c2; }
    var r: f64 = 0.0;
    if p == 0.0 {
        r = exp(m);
    } else {
        if m < 0.0 { return nan(), ERR_UNDEFINED; }
        r = pow(m, 1.0 / p);
    }
    if is_bad(r) { return nan(), ERR_OVERFLOW; }
    return r, OK;
}

// The mean of the logarithms of a sample. A zero value gives no finite
// logarithm, so the mean has no defined value.
fn exact_mean_log(xs: const f64[:]) -> f64, i64 {
    var s: f64[>..] = [];
    for x in xs {
        let l = log(x);
        if is_inf(l) { return nan(), ERR_UNDEFINED; }
        s.push(l);
    }
    return exact_mean(s[..]);
}

// The mean of the squares of the logarithms of a sample.
fn exact_mean_log2(xs: const f64[:]) -> f64, i64 {
    var s: f64[>..] = [];
    for x in xs { s.push(log(x) * log(x)); }
    return exact_mean(s[..]);
}

// The mean of a sample with weights.
fn exact_mean_weighted(xs: const f64[:], ws: const f64[:]) -> f64, i64 {
    let top = top_weight(ws);
    if top <= 0.0 { return nan(), ERR_UNDEFINED; }
    var s1 = acc_new();
    var s2 = acc_new();
    for i in 0..xs.len {
        let w = ws[i] / top;
        acc_add(s1, w * xs[i]);
        acc_add(s2, w);
    }
    if s1.fault || s2.fault { return nan(), ERR_UNDEFINED; }
    let a = acc_total(s1);
    let b = acc_total(s2);
    if b == 0.0 { return nan(), ERR_UNDEFINED; }
    return a / b, OK;
}

// Signature: gstd(xs, ddof) -> value, code
// The geometric standard deviation of a sample: the exponential of the
// deviation of the logarithms. The divisor is n - ddof, as in SciPy, so
// the default ddof is 1. A value at or below 0 gives ERR_DOMAIN, as in
// SciPy, which gives NaN.
fn gstd(xs: const f64[:], ddof: i64) -> f64, i64 {
    let c = check_sample(xs, 2);
    if c != 0 { return nan(), c; }
    if ddof < 0 { return nan(), ERR_DOMAIN; }
    for x in xs {
        if x <= 0.0 { return nan(), ERR_DOMAIN; }
    }
    var s: f64[>..] = [];
    for x in xs { s.push(log(x)); }
    let sd, c2 = stdev_core(s[..], ddof);
    if c2 != 0 { return nan(), c2; }
    if sd == 0.0 { return nan(), ERR_UNDEFINED; }
    return exp(sd), OK;
}

// Signature: gzscore(out, xs, ddof) -> length, code
// The z-score of the logarithms of a sample. A value at or below 0
// gives ERR_DOMAIN. A constant sample gives ERR_UNDEFINED.
fn gzscore(out: f64[>..]&, xs: const f64[:], ddof: i64) -> i64, i64 {
    out.resize(0);
    let c = check_sample(xs, 2);
    if c != 0 { return 0, c; }
    if ddof < 0 { return 0, ERR_DOMAIN; }
    for x in xs {
        if x <= 0.0 { return 0, ERR_DOMAIN; }
    }
    var s: f64[>..] = [];
    for x in xs { s.push(log(x)); }
    let n, c2 = zscore(out, s[..], ddof);
    if c2 != 0 { out.resize(0); return 0, c2; }
    return n, OK;
}

// Signature: zmap(out, scores, compare, ddof) -> length, code
// The z-score of each score against a reference sample. The reference
// gives the mean and the deviation. A constant reference gives
// ERR_UNDEFINED.
fn zmap(out: f64[>..]&, scores: const f64[:], compare: const f64[:], ddof: i64) -> i64, i64 {
    out.resize(0);
    let c = check_sample(scores, 1);
    if c != 0 { return 0, c; }
    let cc = check_sample(compare, 2);
    if cc != 0 { return 0, cc; }
    if ddof < 0 { return 0, ERR_DOMAIN; }
    let sd, c2 = stdev_core(compare, ddof);
    if c2 != 0 { return 0, c2; }
    if sd == 0.0 { return 0, ERR_UNDEFINED; }
    let m, c3 = exact_mean(compare);
    if c3 != 0 { return 0, c3; }
    for x in scores { out.push((x - m) / sd); }
    return out.len, OK;
}

// Signature: describe(xs, ddof, bias) -> Summary, code
// The summary of a sample: the count, the minimum, the maximum, the
// mean, the variance, the skewness and the kurtosis. The divisor of
// the variance is n - ddof, so the default ddof is 1. The bias tells
// the skewness and the kurtosis use the biased form. The kurtosis is
// the Fisher form, as in SciPy. A constant sample gives ERR_UNDEFINED
// for the skewness and the kurtosis, as in SciPy, which gives NaN.
fn describe(xs: const f64[:], ddof: i64, bias: bool) -> Summary, i64 {
    let c = check_sample(xs, 1);
    if c != 0 { return Summary { .. }, c; }
    if ddof < 0 { return Summary { .. }, ERR_DOMAIN; }
    let n = xs.len;
    if n - ddof < 1 { return Summary { .. }, ERR_TOO_SHORT; }
    let lo, c2 = minimum(xs);
    if c2 != 0 { return Summary { .. }, c2; }
    let hi, c3 = maximum(xs);
    if c3 != 0 { return Summary { .. }, c3; }
    let m, c4 = exact_mean(xs);
    if c4 != 0 { return Summary { .. }, c4; }
    let v, c5 = variance_core(xs, ddof);
    if c5 != 0 { return Summary { .. }, c5; }
    let sk, c6 = skew(xs, bias);
    if c6 != 0 { return Summary { .. }, c6; }
    let ku, c7 = kurtosis(xs, true, bias);
    if c7 != 0 { return Summary { .. }, c7; }
    var s = Summary { .. };
    s.n = n;
    s.min = lo;
    s.max = hi;
    s.mean = m;
    s.variance = v;
    s.skewness = sk;
    s.kurtosis = ku;
    return s, OK;
}

// Signature: lmoment(out, xs, max_order, standardize) -> count, code
// The L-moments of a sample, up to the order max_order. The orders 3
// and up divide by the order 2 value when standardize is true, as in
// SciPy. An order above MAX_LORDER gives ERR_DOMAIN. An order above the
// number of values gives ERR_TOO_SHORT. A constant sample gives
// ERR_UNDEFINED for a standardized order of 3 or more.
fn lmoment(out: f64[>..]&, xs: const f64[:], max_order: i64, standardize: bool) -> i64, i64 {
    out.resize(0);
    let c = check_sample(xs, 1);
    if c != 0 { return 0, c; }
    if max_order < 1 { return 0, ERR_DOMAIN; }
    if max_order > MAX_LORDER { return 0, ERR_DOMAIN; }
    let n = xs.len;
    if max_order > n { return 0, ERR_TOO_SHORT; }
    let s = sort_copy(xs);
    var b: f64[>..] = [];
    for r in 0..max_order {
        let br, c2 = ds_lm_b(s[..], r);
        if c2 != 0 { return 0, c2; }
        b.push(br);
    }
    // The L-moment of order r + 1 is the sum of (-1)^(r-k) C(r,k) C(r+k,k) b_k.
    for r in 0..max_order {
        var v: f64 = 0.0;
        for k in 0..(r + 1) { v += ds_lm_prk(r, k) * b[k]; }
        out.push(v);
    }
    for v in out {
        if is_bad(v) { out.resize(0); return 0, ERR_OVERFLOW; }
    }
    if standardize && max_order >= 3 {
        let l2 = out[1];
        if l2 == 0.0 { out.resize(0); return 0, ERR_UNDEFINED; }
        for k in 2..max_order { out[k] = out[k] / l2; }
    }
    return out.len, OK;
}

// The b coefficient of the L-moments of a sorted sample: the sum of the
// binomial weights of the order statistics.
fn ds_lm_b(s: const f64[:], r: i64) -> f64, i64 {
    let n = s.len;
    var acc = acc_new();
    for j in 0..n {
        let c = ex_choose(j, r);
        if c < 0.0 { return nan(), ERR_UNDEFINED; }
        acc_add(acc, c * s[j]);
    }
    if acc.fault { return nan(), ERR_UNDEFINED; }
    let denom = ex_choose(n - 1, r) * n;
    if denom == 0.0 { return nan(), ERR_UNDEFINED; }
    return acc_total(acc) / denom, OK;
}

// The coefficient of b_k in the L-moment of order r + 1, from Hosking:
// (-1)^(r-k) C(r,k) C(r+k,k).
fn ds_lm_prk(r: i64, k: i64) -> f64 {
    let v = ex_choose(r, k) * ex_choose(r + k, k);
    if (r - k) % 2 != 0 { return -v; }
    return v;
}

// Signature: expectile(xs, alpha) -> value, code
// The expectile of a sample at the level alpha, from 0 to 1. The level
// 0.5 gives the mean. The value solves
// mean(abs(1[a <= t] - alpha) * (t - a)) = 0. The search runs on the
// range of the sample and stops after MAX_EXPECTILE_STEPS steps. A
// search that does not end gives ERR_FIT.
fn expectile(xs: const f64[:], alpha: f64) -> f64, i64 {
    return expectile_unweighted(xs, alpha);
}

// Signature: expectile_weighted(xs, ws, alpha) -> value, code
// The expectile of a sample with weights. A weight below 0 gives
// ERR_DOMAIN. Weights that sum to 0 give ERR_UNDEFINED.
fn expectile_weighted(xs: const f64[:], ws: const f64[:], alpha: f64) -> f64, i64 {
    return expectile_core(xs, ws, alpha, false);
}

// Signature: expectile_unweighted(xs, alpha) -> value, code
// The expectile of a sample. The weights are all 1.
fn expectile_unweighted(xs: const f64[:], alpha: f64) -> f64, i64 {
    return expectile_core(xs, xs, alpha, true);
}

// The search of the expectile. The flag equal tells the weights are all
// 1, which saves the work of the weight sum.
fn expectile_core(xs: const f64[:], ws: const f64[:], alpha: f64, equal: bool) -> f64, i64 {
    let c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    let sc = check_scalar(alpha);
    if sc != OK { return nan(), sc; }
    if alpha < 0.0 || alpha > 1.0 { return nan(), ERR_DOMAIN; }
    if !equal {
        let wc = check_weights(ws, xs.len, 1);
        if wc != 0 { return nan(), wc; }
    }
    var total: f64 = 0.0;
    if equal {
        total = xs.len;
    } else {
        for w in ws { total += w; }
    }
    if total <= 0.0 { return nan(), ERR_UNDEFINED; }
    let lo, c2 = minimum(xs);
    if c2 != 0 { return nan(), c2; }
    let hi, c3 = maximum(xs);
    if c3 != 0 { return nan(), c3; }
    var a = lo;
    var b = hi;
    if a == b { return a, OK; }
    var iter: i64 = 0;
    while iter < MAX_EXPECTILE_STEPS {
        let m = 0.5 * a + 0.5 * b;
        let f = ds_expectile_obj(xs, ws, m, alpha, equal);
        if f == 0.0 { return m, OK; }
        if f > 0.0 { b = m; } else { a = m; }
        if abs(b - a) <= EXPECTILE_TOL * (1.0 + abs(m)) { return 0.5 * a + 0.5 * b, OK; }
        iter += 1;
    }
    return nan(), ERR_FIT;
}

// The objective of the expectile at the point t. The value has the sign
// of the difference of the two sides of the equation.
fn ds_expectile_obj(xs: const f64[:], ws: const f64[:], t: f64, alpha: f64, equal: bool) -> f64 {
    var acc = acc_new();
    for i in 0..xs.len {
        let w = if equal { 1.0 } else { ws[i] };
        let side = if xs[i] <= t { 1.0 } else { 0.0 };
        acc_add(acc, w * abs(side - alpha) * (t - xs[i]));
    }
    if acc.fault { return 0.0; }
    return acc_total(acc);
}
