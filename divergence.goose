// The Cressie-Read power divergence test of the package. The statistic
// holds the Pearson, the log-likelihood, the Freeman-Tukey, the modified
// log-likelihood and the Neyman statistics as special cases of one power.
//
// Each call returns a TestResult and then a fault code. A fault gives NaN
// in every field of the result.
namespace statistics;

import std;
import math;

// The named powers of the divergence.
const LAMBDA_PEARSON: f64 = 1.0;
const LAMBDA_LOG_LIKELIHOOD: f64 = 0.0;
const LAMBDA_FREEMAN_TUKEY: f64 = -0.5;
const LAMBDA_MOD_LOG_LIKELIHOOD: f64 = -1.0;
const LAMBDA_NEYMAN: f64 = -2.0;
const LAMBDA_CRESSIE_READ: f64 = 0.6666666666666666;

// x times log(y), with 0 for x = 0 whatever y is, as the reference does.
fn xlogy_term(x: f64, y: f64) -> f64 {
    if x == 0.0 && !is_nan(y) { return 0.0; }
    return x * log(y);
}

// One term of the statistic for the observed count o and the expected
// count e. A term that is not finite means the statistic has no value.
fn power_term(o: f64, e: f64, lam: f64) -> f64 {
    if lam == 0.0 { return 2.0 * xlogy_term(o, o / e); }
    if lam == -1.0 { return 2.0 * xlogy_term(e, e / o); }
    return 2.0 * o * (pow(o / e, lam) - 1.0) / (lam * (lam + 1.0));
}

// The statistic of the lists, or a fault. The expected counts have the
// same length as the observed counts.
fn power_stat(obs: const f64[:], expected: const f64[:], lam: f64) -> f64, i64 {
    var acc = acc_new();
    let n = obs.len;
    for i in 0..n {
        let t = power_term(obs[i], expected[i], lam);
        if is_nan(t) || is_inf(t) { return nan(), ERR_UNDEFINED; }
        acc_add(acc, t);
    }
    if acc.fault { return nan(), ERR_OVERFLOW; }
    let s = acc_total(acc);
    if is_inf(s) { return nan(), ERR_OVERFLOW; }
    return s, OK;
}

// The power divergence test of the observed counts against the expected
// counts. An empty expected list means equal expected counts. The p-value
// uses k - 1 - ddof degrees of freedom.
fn power_divergence(obs: const f64[:], expected: const f64[:], lam: f64 = 1.0, ddof: i64 = 0) -> TestResult, i64 {
    let c = check_sample(obs, 2);
    if c != OK { return fault_result(c); }
    if is_nan(lam) { return fault_result(ERR_NAN); }
    if is_inf(lam) { return fault_result(ERR_INFINITY); }
    if expected.len > 0 {
        let ce = check_sample(expected, 0);
        if ce != OK { return fault_result(ce); }
        if expected.len != obs.len { return fault_result(ERR_LENGTH); }
    }
    let k = obs.len;
    if ddof < 0 || ddof > k - 2 { return fault_result(ERR_DOMAIN); }
    var total = acc_new();
    for o in obs {
        if o < 0.0 || o > MAX_COUNT { return fault_result(ERR_DOMAIN); }
        acc_add(total, o);
    }
    let sum_obs = acc_total(total);
    var sum_exp = sum_obs;
    var flat: f64[>..] = [];
    if expected.len > 0 {
        var te = acc_new();
        for e in expected {
            if e < 0.0 || e > MAX_COUNT { return fault_result(ERR_DOMAIN); }
            acc_add(te, e);
        }
        sum_exp = acc_total(te);
        if abs(sum_obs - sum_exp) > 1e-8 * max(sum_obs, sum_exp) { return fault_result(ERR_DOMAIN); }
    } else {
        for i in 0..k { flat.push(sum_obs / k); }
    }
    if sum_obs == 0.0 { return fault_result(ERR_UNDEFINED); }
    var s: f64 = 0.0;
    var cs: i64 = OK;
    if expected.len > 0 {
        let r, c2 = power_stat(obs, expected, lam);
        s = r;
        cs = c2;
    } else {
        let r, c2 = power_stat(obs, flat[..], lam);
        s = r;
        cs = c2;
    }
    if cs != OK { return fault_result(cs); }
    return chi_result(s, (k - 1 - ddof));
}

// The power divergence test against equal expected counts.
fn power_divergence_uniform(obs: const f64[:], lam: f64 = 1.0, ddof: i64 = 0) -> TestResult, i64 {
    let none: f64[0] = [];
    return power_divergence(obs, none[..], lam, ddof);
}
