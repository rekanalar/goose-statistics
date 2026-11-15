// The Filliben statistic of a sample against a fitted family, and the
// fits of the chi-square, the student and the F families.
//
// The Filliben statistic of Filliben (1975) compares the sorted sample
// with the quantiles of the fitted family at the medians of the order
// statistics of the uniform law. The median of the k-th order
// statistic of n uniform draws is the median of the beta law with the
// shapes k and n + 1 - k. The statistic is the correlation of the two
// lists.
//
// The three fits follow the method of maximum likelihood. SciPy fits
// the location and the scale of these families as well, and its
// likelihood of the chi-square family has no unique optimum there. The
// package fixes the location at 0 for the chi-square and the F family,
// as it does for the gamma and the Weibull families. The student
// family keeps its location, because the law is symmetric about it.
// The searches use the method of Nelder and Mead. SciPy stops its
// search at another path through the same method, so the parameters
// agree to about 1e-4 and no better.
namespace statistics;

import std;
import math;

const FIT_NM_STEPS: i64 = 200;   // steps of the simplex search
// The largest degrees of freedom of a fit. A degree of freedom far above
// 1e6 gives the same statistic, because the law is flat there. The bound
// keeps the search out of the flat region, where no step improves the
// answer.
const FIT_DF_CAP: f64 = 1e6;

// ---- the Filliben statistic -----------------------------------------------

// The medians of the order statistics of the uniform law, n of them.
// The median of the k-th order statistic of n uniform draws is the
// median of the beta law with the shapes k + 1 and n - k.
fn filliben_medians(n: i64, out: f64[>..]&) {
    out.resize(0);
    for k in 0..n {
        // The median of the mirrored beta law is one minus the median,
        // so the second half reads the first half and no root find runs.
        if k * 2 >= n { out.push(1.0 - out[n - 1 - k]); continue; }
        let a = (k + 1) as f64;
        let b = (n - k) as f64;
        let v, c = betaincinv(a, b, 0.5);
        if c != OK { out.push(((k + 1) - 0.3175) / (n + 0.365)); }
        else { out.push(v); }
    }
}

// The correlation of two lists of the same length.
fn filliben_corr(xs: const f64[:], ys: const f64[:]) -> f64, i64 {
    let n = xs.len;
    if n < 2 || ys.len != n { return nan(), ERR_LENGTH; }
    var mx: f64 = 0.0;
    var my: f64 = 0.0;
    for i in 0..n { mx = mx + xs[i]; my = my + ys[i]; }
    mx = mx / n;
    my = my / n;
    var num: f64 = 0.0;
    var dx: f64 = 0.0;
    var dy: f64 = 0.0;
    for i in 0..n {
        let a = xs[i] - mx;
        let b = ys[i] - my;
        num = num + a * b;
        dx = dx + a * a;
        dy = dy + b * b;
    }
    if !(dx > 0.0) || !(dy > 0.0) { return nan(), ERR_UNDEFINED; }
    return num / sqrt(dx * dy), OK;
}

// The Filliben statistic of a sorted sample against a fitted family.
// The quantiles of the family at the medians of the order statistics
// give the second list.
// The medians come from the caller, so a Monte Carlo run computes them
// once and not once per resample.
fn filliben_statistic(s: const f64[:], family: i64, a: f64, b: f64, loc: f64, meds: const f64[:]) -> f64, i64 {
    let n = s.len;
    if n < 2 { return nan(), ERR_TOO_SHORT; }
    var q: f64[>..] = [];
    for i in 0..n {
        let v, c = ft_ppf(family, meds[i], a, b, loc);
        if c != OK { return nan(), c; }
        q.push(v);
    }
    return filliben_corr(s[..], q[..]);
}

// ---- the fits of the three families ---------------------------------------

// The negative log of the likelihood of a sample under the chi-square
// law with the degrees of freedom df and the scale.
fn fit_chi2_nll(s: const f64[:], df: f64, scale: f64) -> f64 {
    if !(df > 0.0) || !(scale > 0.0) { return inf(); }
    let df = if df > FIT_DF_CAP { FIT_DF_CAP } else { df };
    // The log density of the chi-square law has a closed form. The
    // constants of the degrees of freedom are fixed over the sample, so
    // they leave the loop.
    let g, sg, cg = lgamma(0.5 * df);
    if cg != OK { return inf(); }
    let a = 0.5 * df - 1.0;
    let half_df_log2 = 0.5 * df * log(2.0);
    var sum: f64 = 0.0;
    for v in s {
        if !(v > 0.0) { return inf(); }
        let z = v / scale;
        sum = sum + a * log(z) - 0.5 * z - half_df_log2 - g - log(scale);
    }
    return -sum;
}

// The negative log of the likelihood of a sample under the student law
// with the degrees of freedom df, the location and the scale.
fn fit_t_nll(s: const f64[:], p: const f64[:]) -> f64 {
    var df = p[0];
    let loc = p[1];
    let scale = p[2];
    if !(df > 0.0) || !(scale > 0.0) { return inf(); }
    if df > FIT_DF_CAP { df = FIT_DF_CAP; }
    // The log density of the student law has a closed form with the
    // gamma function. The constants of the degrees of freedom leave the
    // loop, as for the chi-square law.
    let g1, s1, c1 = lgamma(0.5 * (df + 1.0));
    let g2, s2, c2 = lgamma(0.5 * df);
    if c1 != OK || c2 != OK { return inf(); }
    let a = 0.5 * (df + 1.0);
    let head = g1 - g2 - 0.5 * log(df * PI) - log(scale);
    var sum: f64 = 0.0;
    for v in s {
        let z = (v - loc) / scale;
        sum = sum + head - a * log1p(z * z / df);
    }
    return -sum;
}

// The negative log of the likelihood of a sample under the F law with
// the degrees of freedom d1 and d2 and the scale.
fn fit_f_nll(s: const f64[:], p: const f64[:]) -> f64 {
    var d1 = p[0];
    var d2 = p[1];
    let scale = p[2];
    if !(d1 > 0.0) || !(d2 > 0.0) || !(scale > 0.0) { return inf(); }
    if d1 > FIT_DF_CAP { d1 = FIT_DF_CAP; }
    if d2 > FIT_DF_CAP { d2 = FIT_DF_CAP; }
    // The log density of the F law has a closed form with the gamma
    // function, as the two above.
    let g1, s1, c1 = lgamma(0.5 * d1);
    let g2, s2, c2 = lgamma(0.5 * d2);
    let g3, s3, c3 = lgamma(0.5 * (d1 + d2));
    if c1 != OK || c2 != OK || c3 != OK { return inf(); }
    let a = 0.5 * d1 - 1.0;
    let b = 0.5 * (d1 + d2);
    let head = g3 - g1 - g2 + 0.5 * d1 * log(d1 / d2) - log(scale);
    var sum: f64 = 0.0;
    for v in s {
        if !(v > 0.0) { return inf(); }
        let z = v / scale;
        sum = sum + head + a * log(z) - b * log1p(z * d1 / d2);
    }
    return -sum;
}

// ---- the simplex search ---------------------------------------------------

// The step of one parameter of the simplex. A parameter that is not
// positive moves by one, and a positive parameter moves by a tenth of
// its value.
fn fit_step_size(v: f64) -> f64 {
    if !(v > 0.0) { return 1.0; }
    return 0.1 * v;
}

// The order of the points of the simplex, from the best to the worst.
// The order is stable, so the search keeps its path.
fn fit_order(ord: i64[>..]&, vals: const f64[:]) {
    let n = ord.len;
    var i: i64 = 1;
    while i < n {
        let v = vals[ord[i]];
        var j = i;
        while j > 0 && vals[ord[j - 1]] > v {
            let t = ord[j];
            ord[j] = ord[j - 1];
            ord[j - 1] = t;
            j -= 1;
        }
        i += 1;
    }
}

// The search of the method of Nelder and Mead over the parameters. The
// buffer `x` holds the (n + 1) points of the simplex one after the
// other, each of n parameters. The buffer `v` holds the value of the
// objective at every point. The search ends with the best point in
// x[0..n].
fn fit_nelder<F>(x: f64[>..]&, v: f64[>..]&, n: i64, steps: i64) {
    // The order of the points, from the best to the worst.
    var ord: i64[>..] = [];
    for i in 0..n + 1 { ord.push(i); }
    var it: i64 = 0;
    while it < steps {
        fit_order(ord, v[..]);
        // The simplex is small enough when its points are close against
        // the scale of the parameters.
        var spread: f64 = 0.0;
        for k in 0..n {
            let d = abs(x[ord[0] * n + k] - x[ord[n] * n + k]) / (1.0 + abs(x[ord[0] * n + k]));
            if d > spread { spread = d; }
        }
        // The values of the objective shrink with the simplex. A spread
        // that is not positive means the objective is flat or not finite,
        // so no step can improve it.
        let vspread = v[ord[n]] - v[ord[0]];
        if spread < 1e-8 || !(vspread > 0.0) { break; }
        // The search stops when the values of the objective are close
        // against the size of the best value. This is the second rule of
        // the method of Nelder and Mead.
        if vspread < 1e-10 * (1.0 + abs(v[ord[0]])) { break; }
        let worst = ord[n];
        var centre: f64[>..] = [];
        for k in 0..n {
            var s: f64 = 0.0;
            for i in 0..n { s = s + x[ord[i] * n + k]; }
            centre.push(s / n as f64);
        }
        var refl: f64[>..] = [];
        for k in 0..n { refl.push(centre[k] + (centre[k] - x[worst * n + k])); }
        let fr = F(refl[..]);
        let best = ord[0];
        let second = ord[n - 1];
        if fr < v[best] {
            // The reflection beats the best point. Try to expand it.
            var exp: f64[>..] = [];
            for k in 0..n { exp.push(centre[k] + 2.0 * (centre[k] - x[worst * n + k])); }
            let fe = F(exp[..]);
            if fe < fr {
                for k in 0..n { x[worst * n + k] = exp[k]; }
                v[worst] = fe;
            } else {
                for k in 0..n { x[worst * n + k] = refl[k]; }
                v[worst] = fr;
            }
        } else if fr < v[second] {
            // The reflection beats the second worst point.
            for k in 0..n { x[worst * n + k] = refl[k]; }
            v[worst] = fr;
        } else {
            // The reflection is worse. Pull the worst point to the
            // centre, and shrink the whole simplex when that fails.
            var contr: f64[>..] = [];
            for k in 0..n { contr.push(centre[k] + 0.5 * (x[worst * n + k] - centre[k])); }
            let fc = F(contr[..]);
            if fc < v[worst] {
                for k in 0..n { x[worst * n + k] = contr[k]; }
                v[worst] = fc;
            } else {
                for i in 1..n + 1 {
                    let j = ord[i];
                    for k in 0..n {
                        x[j * n + k] = x[best * n + k] + 0.5 * (x[j * n + k] - x[best * n + k]);
                    }
                    v[j] = F(x[j * n..j * n + n]);
                }
            }
        }
        it += 1;
    }
    fit_order(ord, v[..]);
    // The best point moves to the front.
    let best = ord[0];
    if best != 0 {
        for k in 0..n {
            let t = x[k];
            x[k] = x[best * n + k];
            x[best * n + k] = t;
        }
        let tv = v[0];
        v[0] = v[best];
        v[best] = tv;
    }
}

// The maximum likelihood fit of the chi-square family. The location
// stays at 0. The parameters are the degrees of freedom and the scale.
fn fit_chi2(s: const f64[:]) -> f64, f64, i64 {
    var sum: f64 = 0.0;
    for v in s { sum += v; }
    let mean = sum / s.len;
    if !(mean > 0.0) { return nan(), nan(), ERR_FIT; }
    // The start follows the method of moments: the mean of the law is
    // the degrees of freedom times the scale.
    var full: f64[>..] = [];
    var vals: f64[>..] = [];
    let df0 = mean / 2.0;
    let sc0 = 2.0;
    full.push(df0);
    full.push(sc0);
    vals.push(fit_chi2_nll(s, df0, sc0));
    let df1 = df0 + fit_step_size(df0);
    full.push(df1);
    full.push(sc0);
    vals.push(fit_chi2_nll(s, df1, sc0));
    let sc1 = sc0 + fit_step_size(sc0);
    full.push(df0);
    full.push(sc1);
    vals.push(fit_chi2_nll(s, df0, sc1));
    fit_nelder(full, vals, 2, FIT_NM_STEPS) { p => fit_chi2_nll(s, p[0], p[1]) };
    let df = full[0];
    let scale = full[1];
    if !(df > 0.0) || !(scale > 0.0) || is_inf(df) || is_inf(scale) { return nan(), nan(), ERR_FIT; }
    return df, scale, OK;
}

// The maximum likelihood fit of the student family. The parameters are
// the degrees of freedom, the location and the scale.
fn fit_t(s: const f64[:]) -> f64, f64, f64, i64 {
    var sum: f64 = 0.0;
    for v in s { sum += v; }
    let m = sum / s.len;
    var ss: f64 = 0.0;
    for v in s { let d = v - m; ss += d * d; }
    let sd = sqrt(ss / (s.len - 1.0));
    if !(sd > 0.0) || is_inf(sd) { return nan(), nan(), nan(), ERR_FIT; }
    var full: f64[>..] = [];
    var vals: f64[>..] = [];
    // The start point: the degrees of freedom 5, the mean and the
    // deviation of the sample.
    var p0: f64[>..] = [];
    p0.push(5.0);
    p0.push(m);
    p0.push(sd);
    for p in p0 { full.push(p); }
    vals.push(fit_t_nll(s, p0[..]));
    for j in 0..3 {
        var pj: f64[>..] = [];
        for p in p0 { pj.push(p); }
        pj[j] = pj[j] + fit_step_size(pj[j]);
        for p in pj { full.push(p); }
        vals.push(fit_t_nll(s, pj[..]));
    }
    fit_nelder(full, vals, 3, FIT_NM_STEPS) { p => fit_t_nll(s, p) };
    let df = full[0];
    let loc = full[1];
    let scale = full[2];
    if !(df > 0.0) || !(scale > 0.0) || is_inf(df) || is_inf(scale) || is_inf(loc) {
        return nan(), nan(), nan(), ERR_FIT;
    }
    return df, scale, loc, OK;
}

// The maximum likelihood fit of the F family. The location stays at 0.
// The parameters are the two degrees of freedom and the scale.
fn fit_f(s: const f64[:]) -> f64, f64, f64, i64 {
    var sum: f64 = 0.0;
    for v in s { sum += v; }
    let mean = sum / s.len;
    if !(mean > 0.0) { return nan(), nan(), nan(), ERR_FIT; }
    // The start follows the method of moments: the mean of the law is
    // d2 / (d2 - 2) for d2 above 2, and the scale follows from it.
    var d2 = 4.0;
    if mean > 2.0 { d2 = 4.0 * mean / (2.0 * mean - 2.0); }
    if d2 <= 2.5 { d2 = 2.5; }
    var full: f64[>..] = [];
    var vals: f64[>..] = [];
    var p0: f64[>..] = [];
    p0.push(4.0);
    p0.push(d2);
    p0.push(mean * (d2 - 2.0) / d2);
    for p in p0 { full.push(p); }
    vals.push(fit_f_nll(s, p0[..]));
    for j in 0..3 {
        var pj: f64[>..] = [];
        for p in p0 { pj.push(p); }
        pj[j] = pj[j] + fit_step_size(pj[j]);
        for p in pj { full.push(p); }
        vals.push(fit_f_nll(s, pj[..]));
    }
    fit_nelder(full, vals, 3, FIT_NM_STEPS) { p => fit_f_nll(s, p) };
    let d1 = full[0];
    let d2f = full[1];
    let scale = full[2];
    if !(d1 > 0.0) || !(d2f > 0.0) || !(scale > 0.0) || is_inf(d1) || is_inf(d2f) || is_inf(scale) {
        return nan(), nan(), nan(), ERR_FIT;
    }
    return d1, d2f, scale, OK;
}
