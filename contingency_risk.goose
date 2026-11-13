// The odds ratio and the relative risk of a two by two table of counts,
// with their confidence intervals.
//
// The odds ratio has two forms. The sample form is the ratio of the
// products of the two diagonals of the table. The conditional form is the
// maximum likelihood estimate of that ratio under the conditional
// distribution of the table. The package finds the conditional form with
// the tilted hypergeometric distribution. The tilt w raises the weight of
// each count of the first type by the factor w, and the mean of the
// distribution grows with the tilt. The ratio is the tilt whose mean is
// the observed count.
//
// The interval of the sample form is the Wald interval. It uses the
// logarithm of the ratio and the standard error of that logarithm.
//
// The relative risk is the ratio of two risks. Its interval is the Katz
// interval. It uses the logarithm of the ratio too.
//
// A table with a zero row sum or a zero column sum has no ratio and
// gives the undefined fault.
namespace statistics;

import std;
import math;

// The checks of a table of four whole counts.
fn risk_table(table: const f64[:]) -> i64 {
    if table.len != 4 { return ERR_LENGTH; }
    for v in table {
        if is_nan(v) { return ERR_NAN; }
        if is_inf(v) { return ERR_INFINITY; }
        if v < 0.0 { return ERR_DOMAIN; }
        if v != floor(v) { return ERR_DOMAIN; }
    }
    return OK;
}

// The first and the last count of the support of the tilted
// hypergeometric.
fn nchyp_support(draws: f64, type1: f64, pop: f64) -> f64, f64 {
    let lo = max(0.0, draws - (pop - type1));
    let hi = min(draws, type1);
    return lo, hi;
}

// The log masses of the hypergeometric over the support. The tilt does
// not change them, so one array serves every step of a search.
fn nchyp_logpmf(out: f64[>..]&, draws: f64, type1: f64, pop: f64) {
    let lo, hi = nchyp_support(draws, type1, pop);
    let n = (hi - lo) as i64 + 1;
    out.resize(0);
    for k in 0..n {
        let d = hyper_logpmf(lo + k, pop, type1, draws);
        out.push(d.hi + d.lo);
    }
}

// The mean of the tilted hypergeometric. The weight of a count x is its
// mass times the tilt t to the power x. The weights are scaled by their
// largest value, so a large tilt does not overflow.
fn nchyp_mean(logs: const f64[:], lo: f64, t: f64) -> f64 {
    var m: f64 = -1.0e300;
    for j in 0..logs.len {
        let v = logs[j] + (lo + j) * t;
        if v > m { m = v; }
    }
    var num: f64 = 0.0;
    var den: f64 = 0.0;
    for j in 0..logs.len {
        let e = exp(logs[j] + (lo + j) * t - m);
        num += (lo + j) * e;
        den += e;
    }
    return num / den;
}

// The share of the tilted hypergeometric at or below the count k.
fn nchyp_cdf(logs: const f64[:], lo: f64, k: f64, t: f64) -> f64 {
    var m: f64 = -1.0e300;
    for j in 0..logs.len {
        let v = logs[j] + (lo + j) * t;
        if v > m { m = v; }
    }
    var num: f64 = 0.0;
    var den: f64 = 0.0;
    for j in 0..logs.len {
        let e = exp(logs[j] + (lo + j) * t - m);
        den += e;
        if lo + j <= k { num += e; }
    }
    return num / den;
}

// The tilt whose tilted hypergeometric has the mean x. The mean grows
// with the tilt.
fn nchyp_tilt_mean(x: f64, logs: const f64[:], lo: f64) -> f64, i64 {
    var tlo: f64 = 0.0;
    var thi: f64 = 0.0;
    var i: i64 = 0;
    if nchyp_mean(logs, lo, 0.0) < x {
        tlo = 0.0;
        thi = 1.0;
        while i < 60 && nchyp_mean(logs, lo, thi) < x {
            thi = thi * 2.0;
            i += 1;
        }
        if i >= 60 { return nan(), ERR_FIT; }
    } else {
        thi = 0.0;
        tlo = -1.0;
        while i < 60 && nchyp_mean(logs, lo, tlo) > x {
            tlo = tlo * 2.0;
            i += 1;
        }
        if i >= 60 { return nan(), ERR_FIT; }
    }
    let t, c = brent_root(tlo, thi, 1e-13, BRENT_RTOL) { v => nchyp_mean(logs, lo, v) - x };
    return t, c;
}

// The share of the tilted hypergeometric at the count x. The tail form
// gives the share above the count.
fn nchyp_share(logs: const f64[:], lo: f64, x: f64, tilt: f64, tail: bool) -> f64 {
    let c = nchyp_cdf(logs, lo, x, tilt);
    if tail { return 1.0 - c; }
    return c;
}

// The tilt whose tilted hypergeometric has the share alpha at or below the
// count x, or the share alpha above the count x. The tail form uses the
// share above the count. The share falls with the tilt.
fn nchyp_tilt(x: f64, alpha: f64, logs: const f64[:], lo: f64, tail: bool) -> f64, i64 {
    var tlo: f64 = -1.0;
    var thi: f64 = 1.0;
    var i: i64 = 0;
    var j: i64 = 0;
    var found = false;
    while i < 60 {
        let s = nchyp_share(logs, lo, x, thi, tail);
        if tail {
            if s >= alpha { found = true; break; }
        } else {
            if s <= alpha { found = true; break; }
        }
        thi = thi * 2.0;
        i += 1;
    }
    if !found {
        if tail { return -1.0 / 0.0, OK; }
        return 1.0 / 0.0, OK;
    }
    while j < 60 {
        let s = nchyp_share(logs, lo, x, tlo, tail);
        if tail {
            if s <= alpha { break; }
        } else {
            if s >= alpha { break; }
        }
        tlo = tlo * 2.0;
        j += 1;
    }
    if j >= 60 { return -1.0 / 0.0, OK; }
    let t, c = brent_root(tlo, thi, 1e-13, BRENT_RTOL) { v => nchyp_share(logs, lo, x, v, tail) - alpha };
    return t, c;
}

// The odds ratio of a table of four counts, with its interval. The kind
// selects the sample form or the conditional form. The alternative
// selects the two-sided interval or a one-sided interval. The level is
// the share of the interval.
fn odds_ratio(table: const f64[:], kind: i64 = ODDS_SAMPLE, alt: i64 = ALT_TWO_SIDED, level: f64 = 0.95) -> IntervalResult, i64 {
    let c = risk_table(table);
    if c != OK { return interval_fault(c), c; }
    let ca = check_alt(alt);
    if ca != OK { return interval_fault(ca), ca; }
    if kind != ODDS_SAMPLE && kind != ODDS_CONDITIONAL { return interval_fault(ERR_DOMAIN), ERR_DOMAIN; }
    if !(level > 0.0 && level < 1.0) { return interval_fault(ERR_DOMAIN), ERR_DOMAIN; }
    let a = table[0];
    let b = table[1];
    let cc = table[2];
    let d = table[3];
    if a + b == 0.0 || cc + d == 0.0 || a + cc == 0.0 || b + d == 0.0 {
        return interval_fault(ERR_UNDEFINED), ERR_UNDEFINED;
    }
    if kind == ODDS_SAMPLE { return odds_ratio_sample(a, b, cc, d, alt, level); }
    return odds_ratio_conditional(a, b, cc, d, alt, level);
}

// The sample odds ratio with the Wald interval. The standard error of the
// logarithm of the ratio is the root of the sum of the reciprocals of the
// four counts.
fn odds_ratio_sample(a: f64, b: f64, cc: f64, d: f64, alt: i64, level: f64) -> IntervalResult, i64 {
    if b == 0.0 || cc == 0.0 {
        if a == 0.0 || d == 0.0 { return interval_fault(ERR_UNDEFINED), ERR_UNDEFINED; }
        return interval_fault(ERR_OVERFLOW), ERR_OVERFLOW;
    }
    let or = a * d / (b * cc);
    if is_nan(or) { return interval_fault(ERR_UNDEFINED), ERR_UNDEFINED; }
    if is_inf(or) { return interval_fault(ERR_OVERFLOW), ERR_OVERFLOW; }
    let se = sqrt(1.0 / a + 1.0 / b + 1.0 / cc + 1.0 / d);
    let z = odds_ratio_z(alt, level);
    let lo = odds_ratio_low(or, se, z, alt);
    let hi = odds_ratio_high(or, se, z, alt);
    if is_nan(lo) || is_nan(hi) { return interval_fault(ERR_UNDEFINED), ERR_UNDEFINED; }
    return IntervalResult { statistic: or, low: lo, high: hi }, OK;
}

// The z value of the interval of a ratio.
fn odds_ratio_z(alt: i64, level: f64) -> f64 {
    if alt == ALT_LESS || alt == ALT_GREATER {
        let z, c = norm_ppf(level);
        return z;
    }
    let z, c = norm_ppf(0.5 * level + 0.5);
    return z;
}

// The lower end of the Wald interval of a ratio.
fn odds_ratio_low(or: f64, se: f64, z: f64, alt: i64) -> f64 {
    if alt == ALT_LESS { return 0.0; }
    let l = log(or) - z * se;
    return exp(l);
}

// The upper end of the Wald interval of a ratio.
fn odds_ratio_high(or: f64, se: f64, z: f64, alt: i64) -> f64 {
    if alt == ALT_GREATER { return 1.0 / 0.0; }
    let h = log(or) + z * se;
    return exp(h);
}

// The conditional odds ratio with its exact interval. The ratio is the
// tilt of the tilted hypergeometric whose mean is the observed count.
fn odds_ratio_conditional(a: f64, b: f64, cc: f64, d: f64, alt: i64, level: f64) -> IntervalResult, i64 {
    let x = a;
    let pop = a + b + cc + d;
    let draws = a + b;
    let type1 = a + cc;
    let lo, hi = nchyp_support(draws, type1, pop);
    if hi - lo + 1.0 > MAX_TABLE_CELLS { return interval_fault(ERR_TOO_LONG), ERR_TOO_LONG; }
    var logs: f64[>..] = [];
    nchyp_logpmf(logs, draws, type1, pop);
    if x == lo { return odds_ratio_conditional_at(x, 0.0, draws, type1, pop, alt, level, logs[..]); }
    if x == hi { return interval_fault(ERR_OVERFLOW), ERR_OVERFLOW; }
    let t, c = nchyp_tilt_mean(x, logs[..], lo);
    if c != OK { return interval_fault(c), c; }
    let or = exp(t);
    if is_inf(or) { return interval_fault(ERR_OVERFLOW), ERR_OVERFLOW; }
    return odds_ratio_conditional_at(x, or, draws, type1, pop, alt, level, logs[..]);
}

// The interval of the conditional odds ratio. The upper end is the tilt
// whose share at the observed count is the share of the tail. The lower
// end is the tilt whose share above the observed count is that share.
fn odds_ratio_conditional_at(x: f64, or: f64, draws: f64, type1: f64, pop: f64, alt: i64, level: f64, logs: const f64[:]) -> IntervalResult, i64 {
    var alpha: f64 = 0.5 * (1.0 - level);
    if alt == ALT_LESS || alt == ALT_GREATER { alpha = 1.0 - level; }
    var low: f64 = 0.0;
    var high: f64 = 1.0 / 0.0;
    let lo, hi = nchyp_support(draws, type1, pop);
    if alt != ALT_LESS {
        let t, c = nchyp_tilt(x - 1.0, alpha, logs[..], lo, true);
        if c != OK { return interval_fault(c), c; }
        low = exp(t);
    }
    if alt != ALT_GREATER {
        let t, c = nchyp_tilt(x, alpha, logs[..], lo, false);
        if c != OK { return interval_fault(c), c; }
        high = exp(t);
    }
    if is_nan(low) || is_nan(high) { return interval_fault(ERR_UNDEFINED), ERR_UNDEFINED; }
    return IntervalResult { statistic: or, low: low, high: high }, OK;
}

// The relative risk of two groups with their Katz interval. The exposed
// group holds the cases and the total, the control group holds the same.
// The interval uses the logarithm of the ratio.
fn relative_risk(exposed_cases: f64, exposed_total: f64, control_cases: f64, control_total: f64, level: f64 = 0.95) -> IntervalResult, i64 {
    let c = risk_count(exposed_cases);
    if c != OK { return interval_fault(c), c; }
    let c2 = risk_count(exposed_total);
    if c2 != OK { return interval_fault(c2), c2; }
    let c3 = risk_count(control_cases);
    if c3 != OK { return interval_fault(c3), c3; }
    let c4 = risk_count(control_total);
    if c4 != OK { return interval_fault(c4), c4; }
    if exposed_total < 1.0 || control_total < 1.0 { return interval_fault(ERR_DOMAIN), ERR_DOMAIN; }
    if exposed_cases > exposed_total || control_cases > control_total {
        return interval_fault(ERR_DOMAIN), ERR_DOMAIN;
    }
    if !(level > 0.0 && level < 1.0) { return interval_fault(ERR_DOMAIN), ERR_DOMAIN; }
    if exposed_cases == 0.0 && control_cases == 0.0 {
        return interval_fault(ERR_UNDEFINED), ERR_UNDEFINED;
    }
    if control_cases == 0.0 { return interval_fault(ERR_OVERFLOW), ERR_OVERFLOW; }
    let rr = (exposed_cases / exposed_total) / (control_cases / control_total);
    if is_nan(rr) { return interval_fault(ERR_UNDEFINED), ERR_UNDEFINED; }
    if is_inf(rr) { return interval_fault(ERR_OVERFLOW), ERR_OVERFLOW; }
    if exposed_cases == 0.0 { return interval_fault(ERR_UNDEFINED), ERR_UNDEFINED; }
    let se = sqrt(1.0 / exposed_cases - 1.0 / exposed_total + 1.0 / control_cases - 1.0 / control_total);
    let z = odds_ratio_z(ALT_TWO_SIDED, level);
    let low = rr * exp(-z * se);
    let high = rr * exp(z * se);
    if is_nan(low) || is_nan(high) { return interval_fault(ERR_UNDEFINED), ERR_UNDEFINED; }
    return IntervalResult { statistic: rr, low: low, high: high }, OK;
}

// The check of a count of a table.
fn risk_count(v: f64) -> i64 {
    if is_nan(v) { return ERR_NAN; }
    if is_inf(v) { return ERR_INFINITY; }
    if v < 0.0 { return ERR_DOMAIN; }
    if v != floor(v) { return ERR_DOMAIN; }
    return OK;
}
