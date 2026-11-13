// The frequencies of a sample: the empirical distribution function, the
// estimated distribution function, the cumulative and the relative
// frequencies, and the tie correction of a rank test.
namespace statistics;

import std;
import math;

// The pairs of the a and b constants of the continuous methods of the
// estimated distribution function. The methods are CDF_INTERPOLATED,
// CDF_HAZEN, CDF_WEIBULL, CDF_LINEAR, CDF_MEDIAN_UNBIASED and
// CDF_NORMAL_UNBIASED.
fn fr_cdf_ab(method: i64) -> f64, f64 {
    if method == CDF_INTERPOLATED { return 0.0, 1.0; }
    if method == CDF_HAZEN { return 0.5, 0.5; }
    if method == CDF_WEIBULL { return 0.0, 0.0; }
    if method == CDF_LINEAR { return 1.0, 1.0; }
    if method == CDF_MEDIAN_UNBIASED { return 1.0 / 3.0, 1.0 / 3.0; }
    return 0.375, 0.375;
}

// The d constant of the discontinuous methods of the estimated
// distribution function. The methods are CDF_INVERTED,
// CDF_AVERAGED_INVERTED and CDF_CLOSEST.
fn fr_cdf_d(method: i64) -> f64 {
    if method == CDF_CLOSEST { return 0.5; }
    return 0.0;
}

// Whether a method of the estimated distribution function is
// discontinuous.
fn fr_cdf_discontinuous(method: i64) -> bool {
    return method == CDF_INVERTED || method == CDF_AVERAGED_INVERTED || method == CDF_CLOSEST;
}

// Signature: tiecorrect(ranks) -> value, code
// The tie correction of the Mann-Whitney and the Kruskal-Wallis tests:
// 1 - sum(t^3 - t) / (n^3 - n), where t are the counts of the tied
// ranks. A list of 0 or 1 values gives 1.
fn tiecorrect(ranks: const f64[:]) -> f64, i64 {
    let c = check_sample(ranks, 1);
    if c != 0 { return nan(), c; }
    let n = ranks.len;
    if n < 2 { return 1.0, OK; }
    let s = sort_copy(ranks);
    var acc = acc_new();
    var i: i64 = 0;
    while i < n {
        var j = i + 1;
        while j < n && s[j] == s[i] { j += 1; }
        let t = (j - i) as f64;
        acc_add(acc, t * t * t - t);
        i = j;
    }
    if acc.fault { return nan(), ERR_UNDEFINED; }
    let d = acc_total(acc);
    let nf = n as f64;
    return 1.0 - d / (nf * nf * nf - nf), OK;
}

// Signature: ecdf(quantiles, probs, xs) -> count, code
// The empirical distribution function of a sample. The output lists
// hold the sorted unique values and the share of the values at or below
// them. Ties count once.
fn ecdf(quantiles: f64[>..]&, probs: f64[>..]&, xs: const f64[:]) -> i64, i64 {
    quantiles.resize(0);
    probs.resize(0);
    let c = check_sample(xs, 1);
    if c != 0 { return 0, c; }
    let s = sort_copy(xs);
    let n = xs.len;
    var i: i64 = 0;
    while i < n {
        var j = i + 1;
        while j < n && s[j] == s[i] { j += 1; }
        quantiles.push(s[i]);
        probs.push(j / (n as f64));
        i = j;
    }
    return quantiles.len, OK;
}

// Signature: ecdf_value(xs, x) -> value, code
// The share of the values of a sample at or below the point x.
fn ecdf_value(xs: const f64[:], x: f64) -> f64, i64 {
    let c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    let sc = check_scalar(x);
    if sc != OK { return nan(), sc; }
    var count: i64 = 0;
    for v in xs {
        if v <= x { count += 1; }
    }
    return count / (xs.len as f64), OK;
}

// The times of a censored sample in one list. The second list holds 1
// for an event and 0 for a censored value. The list order holds the
// index of each time, sorted by the time.
fn fr_censored_times(times: f64[>..]&, events: f64[>..]&, order: i64[>..]&, uncensored: const f64[:], right: const f64[:]) {
    times.resize(0);
    events.resize(0);
    order.resize(0);
    for t in uncensored { times.push(t); events.push(1.0); }
    for t in right { times.push(t); events.push(0.0); }
    for i in 0..times.len { order.push(i); }
    sort(order) { a, b => times[a] < times[b] };
}

// Signature: ecdf_censored(quantiles, survival, uncensored, right) -> count, code
// The Kaplan-Meier estimate of the survival function of a censored
// sample. The list right holds the values that are known to be at or
// above the time. The output lists hold the times and the survival
// values. The last survival value is 0 only when the last time is an
// event, as in SciPy.
fn ecdf_censored(quantiles: f64[>..]&, survival: f64[>..]&, uncensored: const f64[:], right: const f64[:]) -> i64, i64 {
    quantiles.resize(0);
    survival.resize(0);
    let c = check_sample(uncensored, 1);
    if c != 0 { return 0, c; }
    let rc = check_sample(right, 0);
    if rc != 0 { return 0, rc; }
    let total = uncensored.len + right.len;
    if total > MAX_INPUT { return 0, ERR_TOO_LONG; }
    var times: f64[>..] = [];
    var events: f64[>..] = [];
    var order: i64[>..] = [];
    fr_censored_times(times, events, order, uncensored, right);
    // The risk set and the number of events at each time.
    var n = total as f64;
    var at_risk = n;
    var s = 1.0;
    var i: i64 = 0;
    while i < total {
        let t = times[order[i]];
        var d: i64 = 0;
        var c: i64 = 0;
        var j = i;
        while j < total && times[order[j]] == t {
            if events[order[j]] == 1.0 { d += 1; } else { c += 1; }
            j += 1;
        }
        if d > 0 {
            s = s * (1.0 - (d) / at_risk);
        }
        quantiles.push(t);
        survival.push(s);
        at_risk -= (d + c);
        i = j;
    }
    return quantiles.len, OK;
}

// Signature: ecdf_interval(low, high, xs, level) -> count, code
// The Greenwood confidence band of the empirical distribution function
// of a sample, at the level level. The output lists hold the lower and
// the upper bound at each unique value, clipped to 0 and 1. The formula
// is undefined at the largest value, where the bound is NaN, as in SciPy.
fn ecdf_interval(low: f64[>..]&, high: f64[>..]&, xs: const f64[:], level: f64) -> i64, i64 {
    low.resize(0);
    high.resize(0);
    let c = check_sample(xs, 1);
    if c != 0 { return 0, c; }
    let lc = check_level(level);
    if lc != 0 { return 0, lc; }
    let z, zc = norm_ppf(0.5 + level / 2.0);
    if zc != 0 { return 0, zc; }
    let s = sort_copy(xs);
    let n = xs.len;
    var cum: f64 = 0.0;
    var i: i64 = 0;
    while i < n {
        var j = i + 1;
        while j < n && s[j] == s[i] { j += 1; }
        let at_risk = (n - i) as f64;
        let d = (j - i) as f64;
        let p = j / (n as f64);
        let sf = 1.0 - p;
        if at_risk - d == 0.0 {
            low.push(nan());
            high.push(nan());
        } else {
            cum += d / (at_risk * (at_risk - d));
            let zse = z * sf * sqrt(cum);
            low.push(max(p - zse, 0.0));
            high.push(min(p + zse, 1.0));
        }
        i = j;
    }
    return low.len, OK;
}

// Signature: ecdf_censored_interval(low, high, uncensored, right, level) -> count, code
// The Greenwood confidence band of the Kaplan-Meier survival function of
// a censored sample, at the level level. The output lists hold the lower
// and the upper bound at each time, clipped to 0 and 1. The formula is
// undefined at the last time, where the bound is NaN, as in SciPy.
fn ecdf_censored_interval(low: f64[>..]&, high: f64[>..]&, uncensored: const f64[:], right: const f64[:], level: f64) -> i64, i64 {
    low.resize(0);
    high.resize(0);
    let c = check_sample(uncensored, 1);
    if c != 0 { return 0, c; }
    let rc = check_sample(right, 0);
    if rc != 0 { return 0, rc; }
    let lc = check_level(level);
    if lc != 0 { return 0, lc; }
    let z, zc = norm_ppf(0.5 + level / 2.0);
    if zc != 0 { return 0, zc; }
    var times: f64[>..] = [];
    var events: f64[>..] = [];
    var order: i64[>..] = [];
    fr_censored_times(times, events, order, uncensored, right);
    let total = times.len;
    var at_risk = total as f64;
    var s = 1.0;
    var cum: f64 = 0.0;
    var i: i64 = 0;
    while i < total {
        let t = times[order[i]];
        var d: i64 = 0;
        var cen: i64 = 0;
        var j = i;
        while j < total && times[order[j]] == t {
            if events[order[j]] == 1.0 { d += 1; } else { cen += 1; }
            j += 1;
        }
        if d > 0 { s = s * (1.0 - (d) / at_risk); }
        let df = d as f64;
        if at_risk - df == 0.0 {
            low.push(nan());
            high.push(nan());
        } else {
            cum += df / (at_risk * (at_risk - df));
            let zse = z * s * sqrt(cum);
            low.push(max(s - zse, 0.0));
            high.push(min(s + zse, 1.0));
        }
        at_risk -= (d + cen);
        i = j;
    }
    return low.len, OK;
}

// Signature: estimated_cdf(xs, y, method) -> value, code
// The estimated distribution function of a sample at the point y. The
// method is one of the nine methods of SciPy. A point below the
// smallest value gives 0 and a point above the largest gives 1.
fn estimated_cdf(xs: const f64[:], y: f64, method: i64) -> f64, i64 {
    let c = check_sample(xs, 1);
    if c != 0 { return nan(), c; }
    let sc = check_scalar(y);
    if sc != OK { return nan(), sc; }
    if method < CDF_INTERPOLATED || method > CDF_CLOSEST { return nan(), ERR_DOMAIN; }
    let s = sort_copy(xs);
    let n = s.len;
    if y < s[0] { return 0.0, OK; }
    if y > s[n - 1] { return 1.0, OK; }
    // The count of the values at or below y.
    var jp1: i64 = 0;
    while jp1 < n && s[jp1] <= y { jp1 += 1; }
    if fr_cdf_discontinuous(method) {
        let d = fr_cdf_d(method);
        return ((jp1) + d) / (n), OK;
    }
    // The continuous methods use the two neighbours of y.
    var j = jp1 - 1;
    if j < 0 { j = 0; }
    var delta: f64 = 1.0;
    if jp1 < n {
        let xj = s[j];
        let xjp1 = s[jp1];
        if xjp1 > xj { delta = (y - xj) / (xjp1 - xj); }
    }
    let a, b = fr_cdf_ab(method);
    let p = ((jp1) + delta - a) / ((n) + 1.0 - a - b);
    if p < 0.0 { return 0.0, OK; }
    if p > 1.0 { return 1.0, OK; }
    return p, OK;
}

// The counts of a histogram of a sample. The output list holds the count
// of each bin. The return holds the lower limit, the bin size and the
// count of the values outside the limits. Without a limit, the range of
// the sample grows by s on each side, where s is the range divided by
// twice the count of the bins minus one, as in SciPy. The limits are
// inclusive, and the last bin holds the values on the last edge.
fn fr_histogram(counts: f64[>..]&, xs: const f64[:], ws: const f64[:], nbins: i64, low: f64, high: f64, use_low: bool, use_high: bool) -> f64, f64, i64, i64 {
    counts.resize(0);
    let c = check_sample(xs, 1);
    if c != 0 { return nan(), nan(), 0, c; }
    let weighted = ws.len > 0;
    if weighted {
        let wc = check_weights(ws, xs.len, 1);
        if wc != 0 { return nan(), nan(), 0, wc; }
    }
    let bc = check_bins(nbins);
    if bc != 0 { return nan(), nan(), 0, bc; }
    let a, c2 = minimum(xs);
    if c2 != 0 { return nan(), nan(), 0, c2; }
    let b, c3 = maximum(xs);
    if c3 != 0 { return nan(), nan(), 0, c3; }
    var lo = a;
    var hi = b;
    if nbins > 1 {
        let s = (b - a) / (2.0 * (nbins - 1));
        lo = a - s;
        hi = b + s;
    }
    if use_low {
        let lc = check_scalar(low);
        if lc != OK { return nan(), nan(), 0, lc; }
        lo = low;
    }
    if use_high {
        let hc = check_scalar(high);
        if hc != OK { return nan(), nan(), 0, hc; }
        hi = high;
    }
    if hi <= lo { return nan(), nan(), 0, ERR_DOMAIN; }
    let binsize = (hi - lo) / nbins;
    // A range that overflows or a bin width that underflows has no valid edges.
    if is_bad(lo) || is_bad(hi) || is_bad(binsize) { return nan(), nan(), 0, ERR_OVERFLOW; }
    if binsize <= 0.0 { return nan(), nan(), 0, ERR_DOMAIN; }
    for i in 0..nbins { counts.push(0.0); }
    var extra: i64 = 0;
    for i in 0..xs.len {
        let x = xs[i];
        if x < lo || x > hi { extra += 1; continue; }
        // The first edge above x gives the bin. The last edge is inclusive.
        // Start from the guess and correct it against the same edge
        // expression, so the answer equals the linear scan at O(1) a point.
        var index = ((x - lo) / binsize) as i64;
        if index < 0 { index = 0; }
        if index > nbins - 1 { index = nbins - 1; }
        while index > 0 && x < lo + index * binsize { index -= 1; }
        while index < nbins - 1 && x >= lo + (index + 1) * binsize { index += 1; }
        counts[index] += if weighted { ws[i] } else { 1.0 };
    }
    return lo, binsize, extra, OK;
}

// Signature: cumfreq(counts, xs, numbins) -> Frequencies, code
// The cumulative frequencies of a sample. The bin count is numbins, and
// the range is the range of the sample. The output list holds the
// running sum of the counts of the bins.
fn cumfreq(counts: f64[>..]&, xs: const f64[:], numbins: i64) -> Frequencies, i64 {
    return cumfreq_limits(counts, xs, numbins, 0.0, 0.0, false, false);
}

// Signature: cumfreq_limits(counts, xs, numbins, low, high) -> Frequencies, code
// The cumulative frequencies of a sample over the range from low to
// high. The bounds are inclusive.
fn cumfreq_limits(counts: f64[>..]&, xs: const f64[:], numbins: i64, low: f64, high: f64) -> Frequencies, i64 {
    return cumfreq_limits(counts, xs, numbins, low, high, true, true);
}

// The cumulative frequencies with the bounds on or off.
fn cumfreq_limits(counts: f64[>..]&, xs: const f64[:], numbins: i64, low: f64, high: f64, use_low: bool, use_high: bool) -> Frequencies, i64 {
    counts.resize(0);
    var f = Frequencies { .. };
    var raw: f64[>..] = [];
    var empty: f64[>..] = [];
    let l, b, extra, c = fr_histogram(raw, xs, empty, numbins, low, high, use_low, use_high);
    if c != 0 { return f, c; }
    var run: f64 = 0.0;
    for v in raw {
        run += v;
        counts.push(run);
    }
    f.lowerlimit = l;
    f.binsize = b;
    f.extrapoints = extra;
    return f, OK;
}

// Signature: relfreq(freqs, xs, numbins) -> Frequencies, code
// The relative frequencies of a sample. The bin count is numbins, and
// the range is the range of the sample. The output list holds the share
// of each bin.
fn relfreq(freqs: f64[>..]&, xs: const f64[:], numbins: i64) -> Frequencies, i64 {
    return relfreq_limits(freqs, xs, numbins, 0.0, 0.0, false, false);
}

// Signature: relfreq_limits(freqs, xs, numbins, low, high) -> Frequencies, code
// The relative frequencies of a sample over the range from low to
// high. The bounds are inclusive.
fn relfreq_limits(freqs: f64[>..]&, xs: const f64[:], numbins: i64, low: f64, high: f64) -> Frequencies, i64 {
    return relfreq_limits(freqs, xs, numbins, low, high, true, true);
}

// The relative frequencies with the bounds on or off.
fn relfreq_limits(freqs: f64[>..]&, xs: const f64[:], numbins: i64, low: f64, high: f64, use_low: bool, use_high: bool) -> Frequencies, i64 {
    freqs.resize(0);
    var f = Frequencies { .. };
    var raw: f64[>..] = [];
    var empty: f64[>..] = [];
    let l, b, extra, c = fr_histogram(raw, xs, empty, numbins, low, high, use_low, use_high);
    if c != 0 { return f, c; }
    let total = xs.len as f64;
    for v in raw { freqs.push(v / total); }
    f.lowerlimit = l;
    f.binsize = b;
    f.extrapoints = extra;
    return f, OK;
}

// Signature: cumfreq_weighted(counts, xs, ws, numbins) -> Frequencies, code
// The cumulative frequencies of a sample with weights. The bin count is
// numbins, and the range is the range of the sample.
fn cumfreq_weighted(counts: f64[>..]&, xs: const f64[:], ws: const f64[:], numbins: i64) -> Frequencies, i64 {
    return cumfreq_weighted_core(counts, xs, ws, numbins, 0.0, 0.0, false, false);
}

// The cumulative frequencies with weights and the bounds on or off.
fn cumfreq_weighted_core(counts: f64[>..]&, xs: const f64[:], ws: const f64[:], numbins: i64, low: f64, high: f64, use_low: bool, use_high: bool) -> Frequencies, i64 {
    counts.resize(0);
    var f = Frequencies { .. };
    var raw: f64[>..] = [];
    let l, b, extra, c = fr_histogram(raw, xs, ws, numbins, low, high, use_low, use_high);
    if c != 0 { return f, c; }
    var run: f64 = 0.0;
    for v in raw {
        run += v;
        counts.push(run);
    }
    f.lowerlimit = l;
    f.binsize = b;
    f.extrapoints = extra;
    return f, OK;
}

// Signature: relfreq_weighted(freqs, xs, ws, numbins) -> Frequencies, code
// The relative frequencies of a sample with weights. The bin count is
// numbins, and the range is the range of the sample.
fn relfreq_weighted(freqs: f64[>..]&, xs: const f64[:], ws: const f64[:], numbins: i64) -> Frequencies, i64 {
    freqs.resize(0);
    var f = Frequencies { .. };
    var raw: f64[>..] = [];
    let l, b, extra, c = fr_histogram(raw, xs, ws, numbins, 0.0, 0.0, false, false);
    if c != 0 { return f, c; }
    var total: f64 = 0.0;
    for w in ws { total += w; }
    if total <= 0.0 { return f, ERR_UNDEFINED; }
    for v in raw { freqs.push(v / total); }
    f.lowerlimit = l;
    f.binsize = b;
    f.extrapoints = extra;
    return f, OK;
}
