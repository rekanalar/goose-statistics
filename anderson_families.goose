// The Anderson-Darling test against the exponential, logistic, Gumbel
// and Weibull families, with the parameters fitted by maximum likelihood.
//
// Each call writes the critical values, the significance levels in
// percent and the fitted parameters into lists. The p-value is read from
// the critical values by linear interpolation, and it is clamped to the
// first and the last level. The statistic is
//     A2 = -N - sum((2 i - 1) / N * (log F(y_i) + log (1 - F(y_(N+1-i)))).
//
// The exponential fit keeps the location at 0 and takes the mean as the
// scale, so the sample must have a positive mean and no negative value.
// The logistic fit solves the two likelihood equations by Newton steps.
// SciPy stops its own solver at the tolerance 1e-5, so the statistics
// agree with SciPy to about that tolerance. The Gumbel fits solve the
// scale equation by the method of Brent.
//
// The Weibull fit has three parameters. The fit solves the likelihood
// equations of the shape and the location, and it keeps the root with
// the largest likelihood that has a shape of at least 1 and a location
// below the minimum. A sample with no such root gives the fit fault.
namespace statistics;

import std;
import math;

const ANDERSON_EXPON_A: f64[5] = [0.916, 1.062, 1.321, 1.591, 1.959];
const ANDERSON_GUMBEL_A: f64[5] = [0.474, 0.637, 0.757, 0.877, 1.038];
const ANDERSON_LOGISTIC_A: f64[6] = [0.426, 0.563, 0.660, 0.769, 0.906, 1.010];
const ANDERSON_WEIBULL_A: f64[88] = [
    0.292, 0.395, 0.467, 0.522, 0.617, 0.711, 0.836, 0.931,
    0.295, 0.399, 0.471, 0.527, 0.623, 0.719, 0.845, 0.941,
    0.298, 0.403, 0.476, 0.534, 0.631, 0.728, 0.856, 0.954,
    0.301, 0.408, 0.483, 0.541, 0.640, 0.738, 0.869, 0.969,
    0.305, 0.414, 0.490, 0.549, 0.650, 0.751, 0.885, 0.986,
    0.309, 0.421, 0.498, 0.559, 0.662, 0.765, 0.902, 1.007,
    0.314, 0.429, 0.508, 0.570, 0.676, 0.782, 0.923, 1.030,
    0.320, 0.438, 0.519, 0.583, 0.692, 0.802, 0.947, 1.057,
    0.327, 0.448, 0.532, 0.598, 0.711, 0.824, 0.974, 1.089,
    0.334, 0.469, 0.547, 0.615, 0.732, 0.850, 1.006, 1.125,
    0.342, 0.472, 0.563, 0.636, 0.757, 0.879, 1.043, 1.167];

// A value rounded to 3 decimals. A tie goes to the even digit, as in NumPy.
fn round3(v: f64) -> f64 {
    let y = v * 1000.0;
    let f = floor(y);
    let d = y - f;
    var r = f;
    if d > 0.5 { r = f + 1.0; }
    else if d == 0.5 && (f - 2.0 * floor(f / 2.0)) != 0.0 { r = f + 1.0; }
    return r / 1000.0;
}

// The logarithm of one minus exp(-a) for a > 0.
fn log_one_minus_exp_neg(a: f64) -> f64 {
    if a > 0.6931471805599453 { return log(1.0 - exp(-a)); }
    return log(-(exp(-a) - 1.0));
}

// The log of the logistic cdf and of its survival function.
fn logistic_log_cdf(z: f64) -> f64 {
    if z >= 0.0 { return -log(1.0 + exp(-z)); }
    return z - log(1.0 + exp(z));
}

// The two logarithms of a standard form at z for one family.
fn anderson_logs(family: i64, z: f64, lo: f64&, hi: f64&) {
    if family == ANDERSON_EXPON {
        if z <= 0.0 { lo = -inf(); } else { lo = log_one_minus_exp_neg(z); }
        hi = -z;
    } else if family == ANDERSON_LOGISTIC {
        lo = logistic_log_cdf(z);
        hi = logistic_log_cdf(-z);
    } else if family == ANDERSON_GUMBEL_R {
        let e = exp(-z);
        lo = -e;
        hi = log_one_minus_exp_neg(e);
    } else if family == ANDERSON_GUMBEL_L {
        let e = exp(z);
        lo = log_one_minus_exp_neg(e);
        hi = -e;
    } else {
        // the Weibull family reads the already powered score
        if z <= 0.0 { lo = -inf(); } else { lo = log_one_minus_exp_neg(z); }
        hi = -z;
    }
}

// The statistic from the sorted sample scaled to standard form.
fn anderson_stat(family: i64, zs: const f64[:]) -> f64 {
    let n = zs.len;
    let nf = n as f64;
    var lo: f64[>..] = [];
    var hi: f64[>..] = [];
    for k in 0..n {
        var l: f64 = 0.0;
        var h: f64 = 0.0;
        anderson_logs(family, zs[k], l, h);
        lo.push(l);
        hi.push(h);
    }
    return -nf - ad_sum(lo[..], hi[..]);
}

// The p-value from the critical values. The levels are fractions that
// fall as the critical values rise. The value is clamped at both ends.
fn anderson_pvalue(a2: f64, crit: const f64[:], level: const f64[:]) -> f64 {
    if a2 <= crit[0] { return level[0]; }
    for i in 1..crit.len {
        if a2 <= crit[i] {
            let t = (a2 - crit[i - 1]) / (crit[i] - crit[i - 1]);
            return level[i - 1] + t * (level[i] - level[i - 1]);
        }
    }
    return level[crit.len - 1];
}

// ---- the fits -------------------------------------------------------

// The fit of the scale of a Gumbel law for the maximum. The scale solves
// s = mean(x) - sum(x w) / sum(w) with w = exp(-x / s). The location is
// -s log(mean(w)). The values come shifted by their minimum.
fn gumbel_fit(loc: f64&, scale: f64&, xs: const f64[:]) -> i64 {
    let n = xs.len as f64;
    var mn = xs[0];
    var mx = xs[0];
    var sum: f64 = 0.0;
    for x in xs { if x < mn { mn = x; } if x > mx { mx = x; } sum += x; }
    let mean = sum / n;
    if mx == mn { return ERR_UNDEFINED; }
    let top = 4.0 * (mx - mn);
    let s, c = brent_root(1e-9 * (mx - mn), top, 1e-14, 4e-16) { s => gumbel_gap(xs, mean, mn, s) };
    if c != OK { return ERR_FIT; }
    var w: f64 = 0.0;
    for x in xs { w += exp(-(x - mn) / s); }
    scale = s;
    loc = mn - s * log(w / n);
    return OK;
}

// The residual of the scale equation of the Gumbel fit at the scale s.
fn gumbel_gap(xs: const f64[:], mean: f64, mn: f64, s: f64) -> f64 {
    var num: f64 = 0.0;
    var den: f64 = 0.0;
    for x in xs {
        let w = exp(-(x - mn) / s);
        num += x * w;
        den += w;
    }
    return s - mean + num / den;
}

// tanh(t / 2) without a tanh call: (1 - e) / (1 + e) with e = exp(-|t|).
fn half_tanh(t: f64) -> f64 {
    let e = exp(-abs(t));
    let v = (1.0 - e) / (1.0 + e);
    if t < 0.0 { return -v; }
    return v;
}

// One Newton solution of the logistic likelihood equations from the mean
// and the unbiased deviation. The step halves while the scale would not
// stay positive.
fn logistic_fit(loc: f64&, scale: f64&, xs: const f64[:]) -> i64 {
    let n = xs.len as f64;
    let mu, cm = exact_mean(xs);
    if cm != OK { return cm; }
    let v, cv = variance_core(xs, 1);
    if cv != OK { return cv; }
    if v <= 0.0 { return ERR_UNDEFINED; }
    var a = mu;
    var b = sqrt(v);
    for it in 0..200 {
        var f1: f64 = 0.0;
        var f2: f64 = n;
        var j11: f64 = 0.0;
        var j12: f64 = 0.0;
        var j21: f64 = 0.0;
        var j22: f64 = 0.0;
        for x in xs {
            let t = (x - a) / b;
            let th = half_tanh(t);
            let s = 0.5 * (1.0 - th);
            let g = -t * th;
            let ds = -s * (1.0 - s);
            let dg = -th - 0.5 * t * (1.0 - th * th);
            f1 += s;
            f2 += g;
            j11 += ds * (-1.0 / b);
            j12 += ds * (-t / b);
            j21 += dg * (-1.0 / b);
            j22 += dg * (-t / b);
        }
        f1 -= n / 2.0;
        let det = j11 * j22 - j12 * j21;
        if det == 0.0 || is_nan(det) { return ERR_FIT; }
        var da = (f1 * j22 - f2 * j12) / det;
        var db = (j11 * f2 - j21 * f1) / det;
        var h: f64 = 1.0;
        while b - h * db <= 0.0 && h > 1e-12 { h /= 2.0; }
        da *= h;
        db *= h;
        a -= da;
        b -= db;
        if abs(da) <= 1e-14 * (abs(a) + b) && abs(db) <= 1e-14 * b { loc = a; scale = b; return OK; }
    }
    return ERR_FIT;
}

// The shape that solves the first likelihood equation of the Weibull fit
// for the location u. The equation is 1/m - sum(d^m ln d) / sum(d^m) +
// mean(ln d) with d = x - u. A value that has no root gives 0.
fn weibull_shape(logs: const f64[:], lmax: f64, mean_log: f64) -> f64 {
    var hi: f64 = 2.0;
    var found = false;
    for k in 0..40 {
        if weibull_dm(logs, lmax, mean_log, hi) < 0.0 { found = true; break; }
        hi *= 2.0;
    }
    if !found { return 0.0; }
    let m, c = brent_root(1e-6, hi, 1e-15, 4e-16) { m => weibull_dm(logs, lmax, mean_log, m) };
    if c != OK { return 0.0; }
    return m;
}

fn weibull_dm(logs: const f64[:], lmax: f64, mean_log: f64, m: f64) -> f64 {
    var num: f64 = 0.0;
    var den: f64 = 0.0;
    for l in logs {
        let w = exp(m * (l - lmax));
        num += w * l;
        den += w;
    }
    return 1.0 / m - num / den + mean_log;
}

// The residual of the location equation and the profile log likelihood
// at the location u. The shape comes from the first equation.
fn weibull_profile(xs: const f64[:], u: f64, m: f64&, ll: f64&) -> f64 {
    let n = xs.len as f64;
    var logs: f64[>..] = [];
    var lmax = -inf();
    var mean_log: f64 = 0.0;
    var inv: f64 = 0.0;
    for x in xs {
        let l = log(x - u);
        logs.push(l);
        if l > lmax { lmax = l; }
        mean_log += l / n;
        inv += 1.0 / (x - u);
    }
    m = weibull_shape(logs[..], lmax, mean_log);
    if m == 0.0 { ll = -inf(); return nan(); }
    var sp: f64 = 0.0;
    var spm: f64 = 0.0;
    for l in logs {
        sp += exp(m * (l - lmax));
        spm += exp((m - 1.0) * (l - lmax));
    }
    // the scale s satisfies s^m = sum(d^m) / n
    let lns = lmax + log(sp / n) / m;
    ll = n * log(m) - n * m * lns + (m - 1.0) * mean_log * n - n;
    return (m - 1.0) / m * inv - n * exp(-lmax) * spm / sp;
}

// The three Weibull parameters: shape, location and scale.
fn weibull_fit(m: f64&, loc: f64&, scale: f64&, xs: const f64[:]) -> i64 {
    let n = xs.len as f64;
    var mn = xs[0];
    var mx = xs[0];
    for x in xs { if x < mn { mn = x; } if x > mx { mx = x; } }
    if mx == mn { return ERR_UNDEFINED; }
    let range = mx - mn;
    var best = -inf();
    var best_u: f64 = 0.0;
    var best_m: f64 = 0.0;
    var prev_u: f64 = 0.0;
    var prev_g: f64 = 0.0;
    var have = false;
    var d = 1e-7 * range;
    var step: i64 = 0;
    while step < 160 {
        let u = mn - d;
        var mm: f64 = 0.0;
        var ll: f64 = 0.0;
        let g = weibull_profile(xs, u, mm, ll);
        if !is_nan(g) {
            if have && prev_g * g < 0.0 {
                let r, c = brent_root(prev_u, u, 1e-15, 4e-16) { v => weibull_root_value(xs, v) };
                if c == OK {
                    var m2: f64 = 0.0;
                    var ll2: f64 = 0.0;
                    let g2 = weibull_profile(xs, r, m2, ll2);
                    if !is_nan(g2) && m2 >= 1.0 && r < mn && ll2 > best { best = ll2; best_u = r; best_m = m2; }
                }
            }
            have = true;
            prev_u = u;
            prev_g = g;
        } else { have = false; }
        d *= 1.2;
        step += 1;
    }
    if best == -inf() { return ERR_FIT; }
    var logs: f64[>..] = [];
    var sp: f64 = 0.0;
    for x in xs { sp += pow(x - best_u, best_m); }
    m = best_m;
    loc = best_u;
    scale = pow(sp / n, 1.0 / best_m);
    return OK;
}

fn weibull_root_value(xs: const f64[:], u: f64) -> f64 {
    var mm: f64 = 0.0;
    var ll: f64 = 0.0;
    return weibull_profile(xs, u, mm, ll);
}

// ---- the test -------------------------------------------------------

// The Weibull critical values for the shape m, by linear interpolation
// over c = 1/m in steps of 0.05. A larger c uses the last row.
fn weibull_critical(crit: f64[>..]&, m: f64) {
    crit.resize(0);
    let c = 1.0 / m;
    var pos = c / 0.05;
    if pos > 10.0 { pos = 10.0; }
    var i = floor(pos) as i64;
    if i > 9 { i = 9; }
    let t = pos - i;
    for j in 0..8 {
        let v = ANDERSON_WEIBULL_A[i * 8 + j] + t * (ANDERSON_WEIBULL_A[(i + 1) * 8 + j] - ANDERSON_WEIBULL_A[i * 8 + j]);
        crit.push(round3(v + 0.0005));
    }
}

// The test of a sample against one family. The lists receive the
// critical values, the levels in percent and the fitted parameters.
fn anderson_family(crit: f64[>..]&, level: f64[>..]&, params: f64[>..]&, xs: const f64[:], family: i64) -> TestResult, i64 {
    crit.resize(0);
    level.resize(0);
    params.resize(0);
    let c = check_sample(xs, 2);
    if c != OK { return fault_result(c); }
    var f = family;
    if f == ANDERSON_GUMBEL || f == ANDERSON_EXTREME1 { f = ANDERSON_GUMBEL_L; }
    if f != ANDERSON_EXPON && f != ANDERSON_LOGISTIC && f != ANDERSON_GUMBEL_L && f != ANDERSON_GUMBEL_R && f != ANDERSON_WEIBULL {
        return fault_result(ERR_DOMAIN);
    }
    let n = xs.len as f64;
    var sorted = sort_copy(xs);
    var zs: f64[>..] = [];
    var p1: f64 = 0.0;
    var p2: f64 = 0.0;
    var p3: f64 = 0.0;
    var cf: i64 = OK;
    if f == ANDERSON_EXPON {
        let mean, cm = exact_mean(xs);
        if cm != OK { return fault_result(cm); }
        if mean <= 0.0 || sorted[0] < 0.0 { return fault_result(ERR_DOMAIN); }
        p2 = mean;
        for y in sorted { zs.push(y / mean); }
        for v in ANDERSON_EXPON_A { crit.push(round3(v / (1.0 + 0.6 / n))); }
    } else if f == ANDERSON_LOGISTIC {
        cf = logistic_fit(p1, p2, xs);
        if cf != OK { return fault_result(cf); }
        for y in sorted { zs.push((y - p1) / p2); }
        for v in ANDERSON_LOGISTIC_A { crit.push(round3(v / (1.0 + 0.25 / n))); }
    } else if f == ANDERSON_GUMBEL_R || f == ANDERSON_GUMBEL_L {
        var data: f64[>..] = [];
        for y in xs { if f == ANDERSON_GUMBEL_R { data.push(y); } else { data.push(-y); } }
        cf = gumbel_fit(p1, p2, data[..]);
        if cf != OK { return fault_result(cf); }
        if f == ANDERSON_GUMBEL_L { p1 = -p1; }
        for y in sorted { zs.push((y - p1) / p2); }
        for v in ANDERSON_GUMBEL_A { crit.push(round3(v / (1.0 + 0.2 / sqrt(n)))); }
    } else {
        cf = weibull_fit(p1, p2, p3, xs);
        if cf != OK { return fault_result(cf); }
        for y in sorted { zs.push(pow((y - p2) / p3, p1)); }
        weibull_critical(crit, p1);
    }
    let a2 = anderson_stat(f, zs[..]);
    if is_nan(a2) || is_inf(a2) {
        crit.resize(0);
        if is_nan(a2) { return fault_result(ERR_UNDEFINED); }
        return fault_result(ERR_OVERFLOW);
    }
    var fracs: f64[>..] = [];
    if f == ANDERSON_EXPON { level.push(15.0); level.push(10.0); level.push(5.0); level.push(2.5); level.push(1.0); }
    else if f == ANDERSON_LOGISTIC { level.push(25.0); level.push(10.0); level.push(5.0); level.push(2.5); level.push(1.0); level.push(0.5); }
    else if f == ANDERSON_WEIBULL { level.push(0.5); level.push(0.75); level.push(0.85); level.push(0.9); level.push(0.95); level.push(0.975); level.push(0.99); level.push(0.995); }
    else { level.push(25.0); level.push(10.0); level.push(5.0); level.push(2.5); level.push(1.0); }
    for v in level {
        if f == ANDERSON_WEIBULL { fracs.push(1.0 - v); } else { fracs.push(v / 100.0); }
    }
    if f == ANDERSON_EXPON { params.push(0.0); params.push(p2); }
    else if f == ANDERSON_WEIBULL { params.push(p1); params.push(p2); params.push(p3); }
    else { params.push(p1); params.push(p2); }
    return good_result(a2, anderson_pvalue(a2, crit[..], fracs[..]), nan(), nan());
}

fn anderson_exponential(crit: f64[>..]&, level: f64[>..]&, params: f64[>..]&, xs: const f64[:]) -> TestResult, i64 {
    return anderson_family(crit, level, params, xs, ANDERSON_EXPON);
}

fn anderson_logistic(crit: f64[>..]&, level: f64[>..]&, params: f64[>..]&, xs: const f64[:]) -> TestResult, i64 {
    return anderson_family(crit, level, params, xs, ANDERSON_LOGISTIC);
}

fn anderson_gumbel_l(crit: f64[>..]&, level: f64[>..]&, params: f64[>..]&, xs: const f64[:]) -> TestResult, i64 {
    return anderson_family(crit, level, params, xs, ANDERSON_GUMBEL_L);
}

fn anderson_gumbel_r(crit: f64[>..]&, level: f64[>..]&, params: f64[>..]&, xs: const f64[:]) -> TestResult, i64 {
    return anderson_family(crit, level, params, xs, ANDERSON_GUMBEL_R);
}

fn anderson_weibull(crit: f64[>..]&, level: f64[>..]&, params: f64[>..]&, xs: const f64[:]) -> TestResult, i64 {
    return anderson_family(crit, level, params, xs, ANDERSON_WEIBULL);
}
