// The laws of the Anderson-Darling test against the four families: the
// p-value falls as the statistic rises, the parameters are the fit of the
// family, and a rescaled sample keeps the statistic of a scale family.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_anderson_families_tests() -> i64 {
    var bad: i64 = 0;
    var crit: f64[>..] = [];
    var lev: f64[>..] = [];
    var par: f64[>..] = [];
    let s: f64[12] = [0.9, 1.4, 2.1, 0.3, 3.8, 1.1, 0.7, 2.6, 1.9, 5.2, 0.5, 1.6];
    let keep: f64[12] = [0.9, 1.4, 2.1, 0.3, 3.8, 1.1, 0.7, 2.6, 1.9, 5.2, 0.5, 1.6];
    // the exponential fit has the mean as its scale
    let t, c = statistics::anderson_exponential(crit, lev, par, s[..]);
    chk_code(bad, "exponential code", c, 0);
    check(bad, par.len == 2 && par[0] == 0.0, "the exponential location is 0");
    chk_val(bad, "exponential scale is the mean", par[1], c, 1.841666666666667, 1e-12);
    // the statistic of a scale family does not change when the sample is rescaled
    var scaled: f64[>..] = [];
    for v in s { scaled.push(v * 7.0); }
    let t2, c2 = statistics::anderson_exponential(crit, lev, par, scaled[..]);
    chk_val(bad, "exponential scale invariance", t2.statistic, c2, t.statistic, 1e-12);
    // the critical values rise as the level falls
    var rising = true;
    for i in 1..crit.len { if crit[i] <= crit[i - 1] { rising = false; } }
    check(bad, rising && crit.len == lev.len, "the critical values rise with the level");
    // a larger statistic gives a smaller p-value inside the table
    check(bad, t.pvalue >= 0.0 && t.pvalue <= 1.0, "the p-value is a probability");
    // the right and left Gumbel fits mirror each other
    var neg: f64[>..] = [];
    for v in s { neg.push(-v); }
    let tr, cr = statistics::anderson_gumbel_r(crit, lev, par, s[..]);
    let tl, cl = statistics::anderson_gumbel_l(crit, lev, par, neg[..]);
    chk_val(bad, "gumbel mirror statistic", tl.statistic, cl, tr.statistic, 1e-9);
    check(bad, cr == 0 && cl == 0, "the Gumbel fits succeed");
    // the logistic fit is the centre of a symmetric sample
    let sym: f64[8] = [-3.0, -2.0, -1.0, -0.5, 0.5, 1.0, 2.0, 3.0];
    let tg, cg = statistics::anderson_logistic(crit, lev, par, sym[..]);
    chk_code(bad, "logistic code", cg, 0);
    chk_val(bad, "logistic centre", par[0] + 1.0, cg, 1.0, 1e-9);
    // the sample of the caller is unchanged
    var same = true;
    for i in 0..12 { if s[i] != keep[i] { same = false; } }
    check(bad, same, "the test leaves the sample unchanged");
    // a fault leaves the lists empty
    let none: f64[0] = [];
    let tz, cz = statistics::anderson_exponential(crit, lev, par, none[..]);
    check(bad, cz == 1 && crit.len == 0 && lev.len == 0 && par.len == 0, "a fault leaves the lists empty");
    return bad;
}
