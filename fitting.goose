// The fit of a family, the goodness-of-fit test and the survival data.
//
// SciPy raises FitError when a fit does not end. Goose returns faults as
// codes and has no exception, so this file gives the fault code ERR_FIT.
// The fits of the four families are closed forms. A resample that a
// closed form cannot fit gives ERR_FIT, and so does a fit that would
// need more than MAX_FIT_STEPS steps.
namespace statistics;

import std;
import math;

// Survival times in four parts. The struct holds views of the lists of the
// caller, so it never copies them and no call writes into them. The
// `interval` part holds pairs, low then high, one pair after the other.
// A `uncensored` time is an observed event. A `right` time is the last
// time a subject was seen with no event yet. Every value must be finite.
struct CensoredData {
    uncensored: const f64[:],
    left: const f64[:],
    right: const f64[:],
    interval: const f64[:],
}

// ---- the checks of CensoredData --------------------------------------

// The count of the observations of one struct.
fn ft_count(c: const CensoredData) -> i64 {
    return c.uncensored.len + c.left.len + c.right.len + c.interval.len / 2;
}

// Stage 1: a part, or the whole struct, longer than MAX_INPUT.
fn ft_cd_bound(c: const CensoredData) -> i64 {
    if c.uncensored.len > MAX_INPUT || c.left.len > MAX_INPUT { return ERR_TOO_LONG; }
    if c.right.len > MAX_INPUT || c.interval.len > MAX_INPUT { return ERR_TOO_LONG; }
    let total = c.uncensored.len + c.left.len + c.right.len + c.interval.len;
    if total > MAX_INPUT { return ERR_TOO_LONG; }
    return OK;
}

// Stage 3 and 4: NaN first, then infinity, over all four parts.
fn ft_cd_nan(c: const CensoredData) -> i64 {
    for x in c.uncensored { if is_nan(x) { return ERR_NAN; } }
    for x in c.left { if is_nan(x) { return ERR_NAN; } }
    for x in c.right { if is_nan(x) { return ERR_NAN; } }
    for x in c.interval { if is_nan(x) { return ERR_NAN; } }
    return OK;
}

fn ft_cd_inf(c: const CensoredData) -> i64 {
    for x in c.uncensored { if is_inf(x) { return ERR_INFINITY; } }
    for x in c.left { if is_inf(x) { return ERR_INFINITY; } }
    for x in c.right { if is_inf(x) { return ERR_INFINITY; } }
    for x in c.interval { if is_inf(x) { return ERR_INFINITY; } }
    return OK;
}

// Stage 5 and 6: the interval part holds pairs, and each pair has low <= high.
fn ft_cd_shape(c: const CensoredData) -> i64 {
    if c.interval.len % 2 != 0 { return ERR_LENGTH; }
    var i: i64 = 0;
    while i < c.interval.len {
        if c.interval[i] > c.interval[i + 1] { return ERR_DOMAIN; }
        i += 2;
    }
    return OK;
}

// The checks of two structs in the fixed order of the package: the length
// bound, the empty struct, NaN, infinity, the pair shape and the domain.
fn ft_check_censored(x: const CensoredData, y: const CensoredData) -> i64 {
    var c = ft_cd_bound(x);
    if c != OK { return c; }
    c = ft_cd_bound(y);
    if c != OK { return c; }
    if ft_count(x) == 0 || ft_count(y) == 0 { return ERR_EMPTY; }
    c = ft_cd_nan(x);
    if c != OK { return c; }
    c = ft_cd_nan(y);
    if c != OK { return c; }
    c = ft_cd_inf(x);
    if c != OK { return c; }
    c = ft_cd_inf(y);
    if c != OK { return c; }
    c = ft_cd_shape(x);
    if c != OK { return c; }
    return ft_cd_shape(y);
}

// ---- the logrank test --------------------------------------------------

// The count of the values in a sorted list that are at least t.
fn ft_count_ge(s: const f64[:], t: f64) -> i64 {
    var lo: i64 = 0;
    var hi = s.len;
    while lo < hi {
        let m = lo + (hi - lo) / 2;
        if s[m] < t { lo = m + 1; } else { hi = m; }
    }
    return s.len - lo;
}

// The sorted times of the uncensored and right-censored parts of one struct.
fn ft_times(c: const CensoredData) -> f64[>..] {
    var t: f64[>..] = [];
    for v in c.uncensored { t.push(v); }
    for v in c.right { t.push(v); }
    sort(t) { a, b => a < b };
    return t;
}

// The logrank test of two samples of right-censored times. The statistic
// is the signed square root of the usual chi-square statistic, so the
// sign tells which group has more deaths than expected. The at-risk
// counts follow the pooled sample: a subject is at risk at the time t of an
// event when its own time is at least t. The call takes the uncensored and
// right parts only. A left or an interval part gives ERR_DOMAIN, because the
// logrank test has no defined form for those parts. A sample with no
// event, or with no variance, gives ERR_UNDEFINED where SciPy gives NaN.
fn logrank(x: const CensoredData, y: const CensoredData, alt: i64 = ALT_TWO_SIDED) -> TestResult, i64 {
    let c = ft_check_censored(x, y);
    if c != OK { return fault_result(c); }
    if x.left.len > 0 || y.left.len > 0 || x.interval.len > 0 || y.interval.len > 0 {
        return fault_result(ERR_DOMAIN);
    }
    if check_alt(alt) != OK { return fault_result(ERR_DOMAIN); }

    let tx = ft_times(x);
    let ty = ft_times(y);
    var ev: f64[>..] = [];
    for v in x.uncensored { ev.push(v); }
    for v in y.uncensored { ev.push(v); }
    sort(ev) { a, b => a < b };

    var sum_exp: f64 = 0.0;
    var sum_var: f64 = 0.0;
    var i: i64 = 0;
    while i < ev.len {
        let t = ev[i];
        var j = i;
        while j < ev.len && ev[j] == t { j += 1; }
        let d = (j - i) as f64;
        let nx = ft_count_ge(tx[..], t) as f64;
        let ny = ft_count_ge(ty[..], t) as f64;
        let n = nx + ny;
        sum_exp += nx * (d / n);
        if n > 1.0 { sum_var += nx * ny * d * (n - d) / (n * n * (n - 1.0)); }
        i = j;
    }
    if !(sum_var > 0.0) { return fault_result(ERR_UNDEFINED); }
    let z = (x.uncensored.len - sum_exp) / sqrt(sum_var);
    return vr_finish_z(z, alt);
}

// ---- the goodness-of-fit test -------------------------------------------

// The cumulative function, the survival function and the quantile of a
// fitted family. The three numbers a, b and c are the parameters of
// the fit. The exponential family has the rate in a. The student
// family has the degrees of freedom, the scale and the location. The
// chi-square family has the degrees of freedom and the scale. The F
// family has the two degrees of freedom and the scale.
fn ft_cdf(family: i64, x: f64, a: f64, b: f64, c: f64) -> f64, i64 {
    if family == FAMILY_NORMAL { return norm_cdf(x, a, b); }
    if family == FAMILY_UNIFORM { return unif_cdf(x, a, b); }
    if family == FAMILY_EXPONENTIAL { return expon_cdf(x, a); }
    if family == FAMILY_CHI2 { return chi2_cdf(x / b, a); }
    if family == FAMILY_T { return tdist_cdf((x - c) / b, a); }
    if family == FAMILY_F { return fdist_cdf(x / c, a, b); }
    if family >= FAMILY_WEIBULL { return mle_cdf(family, x, a, b); }
    return lognorm_cdf(x, a, b);
}

fn ft_sf(family: i64, x: f64, a: f64, b: f64, c: f64) -> f64, i64 {
    if family == FAMILY_NORMAL { return norm_sf(x, a, b); }
    if family == FAMILY_UNIFORM { return unif_sf(x, a, b); }
    if family == FAMILY_EXPONENTIAL { return expon_sf(x, a); }
    if family == FAMILY_CHI2 { return chi2_sf(x / b, a); }
    if family == FAMILY_T { return tdist_sf((x - c) / b, a); }
    if family == FAMILY_F { return fdist_sf(x / c, a, b); }
    if family >= FAMILY_WEIBULL { return mle_sf(family, x, a, b); }
    return lognorm_sf(x, a, b);
}

fn ft_ppf(family: i64, p: f64, a: f64, b: f64, c: f64) -> f64, i64 {
    if family == FAMILY_NORMAL { return norm_ppf(p, a, b); }
    if family == FAMILY_UNIFORM { return unif_ppf(p, a, b); }
    if family == FAMILY_EXPONENTIAL { return expon_ppf(p, a); }
    if family == FAMILY_CHI2 {
        let v, k = chi2_ppf(p, a);
        if k != OK { return nan(), k; }
        return v * b, OK;
    }
    if family == FAMILY_T {
        let v, k = tdist_ppf(p, a);
        if k != OK { return nan(), k; }
        return v * b + c, OK;
    }
    if family == FAMILY_F {
        let v, k = fdist_ppf(p, a, b);
        if k != OK { return nan(), k; }
        return v * c, OK;
    }
    if family == FAMILY_GAMMA { return gammadist_ppf(p, a, b); }
    if family == FAMILY_BETA { return betadist_ppf(p, a, b); }
    if family == FAMILY_WEIBULL { return weibull_ppf(p, a, b); }
    return lognorm_ppf(p, a, b);
}

// The count of the distinct values of a sorted list. The fit of a family
// with n parameters needs n distinct values, or the search has no end.
fn ft_distinct(s: const f64[:]) -> i64 {
    var n: i64 = 0;
    var i: i64 = 0;
    while i < s.len {
        if i == 0 || s[i] != s[i - 1] { n += 1; }
        i += 1;
    }
    return n;
}

// The fit of a family to a sorted sample. The three numbers are the
// parameters of the family. The normal family takes the mean and the
// standard deviation with n - 1, as SciPy does. The exponential and
// the lognormal families fix the location at 0, as SciPy does with the
// fixed location 0; the lognormal takes the maximum likelihood scale.
// The uniform family takes the minimum and the maximum. The chi-square
// family takes the degrees of freedom and the scale. The student family
// takes the degrees of freedom, the scale and the location. The F family
// takes the two degrees of freedom and the scale. A fit that gives a
// scale that is not positive and finite gives ERR_FIT.
fn ft_fit(s: const f64[:], family: i64) -> f64, f64, f64, i64 {
    let n = s.len as f64;
    if family == FAMILY_WEIBULL {
        let a, b, c = mle_weibull(s);
        return a, b, 0.0, c;
    }
    if family == FAMILY_GAMMA {
        let a, b, c = mle_gamma(s);
        return a, b, 0.0, c;
    }
    if family == FAMILY_BETA {
        let a, b, c = mle_beta(s);
        return a, b, 0.0, c;
    }
    if family == FAMILY_CHI2 {
        // A value that is not positive has no density under this law,
        // so no fit exists. The check keeps the search short.
        for v in s { if !(v > 0.0) { return nan(), nan(), nan(), ERR_FIT; } }
        // The fit takes two parameters, so it needs two distinct values.
        if ft_distinct(s) < 2 { return nan(), nan(), nan(), ERR_FIT; }
        let a, b, c = fit_chi2(s);
        return a, b, 0.0, c;
    }
    if family == FAMILY_T {
        // The fit takes three parameters, so it needs three distinct
        // values. Without them the search has no end.
        if ft_distinct(s) < 3 { return nan(), nan(), nan(), ERR_FIT; }
        let a, b, c, k = fit_t(s);
        return a, b, c, k;
    }
    if family == FAMILY_F {
        // A value that is not positive has no density under this law,
        // so no fit exists. The check keeps the search short.
        for v in s { if !(v > 0.0) { return nan(), nan(), nan(), ERR_FIT; } }
        // The fit takes three parameters, so it needs three distinct
        // values. Without them the search has no end.
        if ft_distinct(s) < 3 { return nan(), nan(), nan(), ERR_FIT; }
        let a, b, c, k = fit_f(s);
        return a, b, c, k;
    }
    if family == FAMILY_UNIFORM {
        let lo = s[0];
        let hi = s[s.len - 1];
        if !(hi > lo) || is_inf(hi - lo) { return nan(), nan(), nan(), ERR_FIT; }
        return lo, hi, 0.0, OK;
    }
    var sum: f64 = 0.0;
    for v in s {
        if family == FAMILY_LOGNORMAL { sum += log(v); } else { sum += v; }
    }
    let m = sum / n;
    if is_bad(m) { return nan(), nan(), nan(), ERR_FIT; }
    if family == FAMILY_EXPONENTIAL {
        if !(m > 0.0) { return nan(), nan(), nan(), ERR_FIT; }
        let rate = 1.0 / m;
        if is_inf(rate) { return nan(), nan(), nan(), ERR_FIT; }
        return rate, nan(), 0.0, OK;
    }
    var ss: f64 = 0.0;
    for v in s {
        var d = v - m;
        if family == FAMILY_LOGNORMAL { d = log(v) - m; }
        ss += d * d;
    }
    var denom = n - 1.0;
    if family == FAMILY_LOGNORMAL { denom = n; }
    let sigma = sqrt(ss / denom);
    if !(sigma > 0.0) || is_inf(sigma) { return nan(), nan(), nan(), ERR_FIT; }
    return m, sigma, 0.0, OK;
}

// The Anderson sum of a sorted sample. The lower tail holds the log of
// the cumulative function of every value and the upper tail the log of
// the survival function of every value.
fn ad_sum(lo: const f64[:], hi: const f64[:]) -> f64 {
    let n = lo.len;
    let nf = n as f64;
    var sum: f64 = 0.0;
    for i in 0..n {
        sum += (2.0 * (i + 1) - 1.0) / nf * (lo[i] + hi[n - 1 - i]);
    }
    return sum;
}

// The statistic of a sorted sample against a fitted family.
fn ft_statistic(s: const f64[:], family: i64, a: f64, b: f64, c: f64, stat: i64, meds: const f64[:]) -> f64, i64 {
    let n = s.len;
    if stat == GOF_FILLIBEN { return filliben_statistic(s, family, a, b, c, meds); }
    let nf = n as f64;
    var acc: f64 = 0.0;
    var dplus: f64 = 0.0;
    var dminus: f64 = 0.0;
    if stat == GOF_ANDERSON {
        var lo: f64[>..] = [];
        var hi: f64[>..] = [];
        for k in 0..n {
            let lc, c1 = ft_cdf(family, s[k], a, b, c);
            let ls, c2 = ft_sf(family, s[k], a, b, c);
            if c1 != OK { return nan(), c1; }
            if c2 != OK { return nan(), c2; }
            // A zero probability makes the statistic infinite.
            if !(lc > 0.0) || !(ls > 0.0) { return nan(), ERR_UNDEFINED; }
            lo.push(log(lc));
            hi.push(log(ls));
        }
        acc = ad_sum(lo[..], hi[..]);
    }
    for i in 0..n {
        let fi = i as f64;
        if stat != GOF_ANDERSON {
            let f, c1 = ft_cdf(family, s[i], a, b, c);
            if c1 != OK { return nan(), c1; }
            if stat == GOF_KOLMOGOROV {
                dplus = max(dplus, (fi + 1.0) / nf - f);
                dminus = max(dminus, f - fi / nf);
            } else {
                let u = (2.0 * fi + 1.0) / (2.0 * nf) - f;
                acc += u * u;
            }
        }
    }
    var r: f64 = 0.0;
    if stat == GOF_ANDERSON { r = -nf - acc; }
    else if stat == GOF_KOLMOGOROV { r = max(dplus, dminus); }
    else { r = 1.0 / (12.0 * nf) + acc; }
    if is_bad(r) { return nan(), ERR_UNDEFINED; }
    return r, OK;
}


// One draw from a fitted family. The draw of the chi-square, student
// and F families inverts the cumulative function at a uniform draw.
fn ft_draw(family: i64, a: f64, b: f64, c: f64, r: rng&) -> f64 {
    if family == FAMILY_NORMAL { return a + b * rs_normal(r); }
    if family == FAMILY_UNIFORM { return a + (b - a) * rs_unit(r); }
    if family == FAMILY_EXPONENTIAL { return -log(1.0 - rs_unit(r)) / a; }
    if family >= FAMILY_CHI2 {
        let u = rs_unit(r);
        let v, k = ft_ppf(family, u, a, b, c);
        if k != OK { return nan(); }
        return v;
    }
    if family >= FAMILY_WEIBULL { return mle_draw(family, a, b, rs_unit(r)); }
    return exp(a + b * rs_normal(r));
}

// The goodness-of-fit test of a sample against the normal, uniform,
// exponential, lognormal, Weibull, gamma, beta, chi-square, student or
// F family, with the parameters fitted to the sample (the same fits as
// SciPy, except that the chi-square and F fits fix the location at 0).
// `statistic` is GOF_ANDERSON, GOF_KOLMOGOROV, GOF_CRAMER or
// GOF_FILLIBEN. The Filliben statistic is a correlation, so a small value
// is extreme, and the p-value counts the null values at most as large as
// the observed one. The p-value is the Monte Carlo rule of
// SciPy: each of `resamples` samples comes from the fitted family, is
// fitted again and gives one statistic. The p-value is (the count of
// null values at least as large as the observed one, plus 1) over
// (resamples plus 1). The fitted parameters go into `out`: the mean and
// the deviation, the lower and the upper bound, the rate, or the mean and
// the deviation of the logarithms. The draws come from `r`, resample by
// resample and in the order of the values. A family with no fit here gives
// ERR_DOMAIN. The exponential and lognormal families need positive values.
// The Anderson statistic of the uniform family is infinite and gives
// ERR_UNDEFINED. The draws of SciPy come from another stream. Each resample
// sorts and fits again, so the product of the resamples and the sample size
// stops at MAX_FIT_WORK. A larger job gives ERR_TOO_LONG.
fn goodness_of_fit(out: f64[>..]&, xs: const f64[:], family: i64, statistic: i64, resamples: i64, r: rng&) -> TestResult, i64 {
    out.resize(0);
    let c = check_sample(xs, 2);
    if c != OK { return fault_result(c); }
    if family < FAMILY_NORMAL || family > FAMILY_F { return fault_result(ERR_DOMAIN); }
    if statistic < GOF_ANDERSON || statistic > GOF_FILLIBEN { return fault_result(ERR_DOMAIN); }
    if resamples < 1 { return fault_result(ERR_DOMAIN); }
    if resamples > MAX_GOF_REPS || resamples > MAX_FIT_WORK / xs.len { return fault_result(ERR_TOO_LONG); }
    if family >= FAMILY_EXPONENTIAL && family != FAMILY_T {
        for x in xs { if !(x > 0.0) { return fault_result(ERR_DOMAIN); } }
    }
    if family == FAMILY_BETA {
        for x in xs { if !(x < 1.0) { return fault_result(ERR_DOMAIN); } }
    }
    let s = sort_copy(xs);
    if s[0] == s[s.len - 1] { return fault_result(ERR_UNDEFINED); }

    let a, b, c3, cf = ft_fit(s[..], family);
    if cf != OK { return fault_result(cf); }
    // The medians of the order statistics of the Filliben statistic. They
    // depend on the count alone, so one pass serves every resample.
    var meds: f64[>..] = [];
    if statistic == GOF_FILLIBEN { filliben_medians(s.len, meds); }
    let obs, co = ft_statistic(s[..], family, a, b, c3, statistic, meds[..]);
    if co != OK { return fault_result(co); }

    // The null values, one resample at a time. The count of the values that
    // reach the observed one, up to a relative bound of 100 machine epsilons.
    let gamma = abs(100.0 * EPS * obs);
    let less = statistic == GOF_FILLIBEN;
    var buf: f64[>..] = [];
    var count: i64 = 0;
    for k in 0..resamples {
        buf.resize(0);
        for i in 0..s.len {
            let v = ft_draw(family, a, b, c3, r);
            if is_bad(v) { return fault_result(ERR_FIT); }
            if family >= FAMILY_EXPONENTIAL && family != FAMILY_T && !(v > 0.0) { return fault_result(ERR_FIT); }
            buf.push(v);
        }
        sort(buf) { x, y => x < y };
        let ra, rb, rc, cr = ft_fit(buf[..], family);
        if cr != OK { return fault_result(cr); }
        let v, cs = ft_statistic(buf[..], family, ra, rb, rc, statistic, meds[..]);
        if cs != OK { return fault_result(cs); }
        if less { if v <= obs + gamma { count += 1; } }
        else if v >= obs - gamma { count += 1; }
    }

    // The parameters in the order of SciPy: the normal family gives the
    // location and the scale, the uniform family the location and the
    // width, the exponential family the location 0 and the scale, and
    // the lognormal family the shape, the location 0 and the scale.
    if family == FAMILY_NORMAL { out.push(a); out.push(b); }
    else if family == FAMILY_UNIFORM { out.push(a); out.push(b - a); }
    else if family == FAMILY_EXPONENTIAL { out.push(0.0); out.push(1.0 / a); }
    else if family == FAMILY_WEIBULL || family == FAMILY_GAMMA { out.push(a); out.push(0.0); out.push(b); }
    else if family == FAMILY_BETA { out.push(a); out.push(b); out.push(0.0); out.push(1.0); }
    else if family == FAMILY_CHI2 { out.push(a); out.push(0.0); out.push(b); }
    else if family == FAMILY_T { out.push(a); out.push(c3); out.push(b); }
    else if family == FAMILY_F { out.push(a); out.push(b); out.push(0.0); out.push(c3); }
    else { out.push(b); out.push(0.0); out.push(exp(a)); }
    let p = (count + 1) as f64 / (resamples + 1);
    return good_result(obs, p, nan(), nan());
}
