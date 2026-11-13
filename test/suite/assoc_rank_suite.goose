// The laws of the rank association calls: the correlation lies in -1 to 1,
// the p-value lies in 0 to 1, the weighted tau lies in -1 to 1, and the
// Spearman rho of a sample against itself is 1.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_assoc_rank_tests() -> i64 {
    var bad: i64 = 0;
    let x: f64[6] = [2.1, 3.4, 1.9, 5.6, 4.4, 3.3];
    let y: f64[6] = [1.2, 3.9, 2.2, 5.1, 4.8, 2.9];
    let keep: f64[6] = [2.1, 3.4, 1.9, 5.6, 4.4, 3.3];
    // a sample against itself
    let a1, c1 = statistics::spearmanrho(x[..], x[..]);
    chk_val(bad, "spearmanrho of a sample against itself", a1.statistic, c1, 1.0, 1e-12);
    // the correlation lies in -1 to 1 and the p-value in 0 to 1
    for alt in 0..3 {
        let r, c = statistics::spearmanrho(x[..], y[..], alt);
        chk_code(bad, "spearmanrho code", c, 0);
        if r.statistic < -1.0 || r.statistic > 1.0 { bad += 1; print("FAIL spearmanrho range\n"); }
        if r.pvalue < 0.0 || r.pvalue > 1.0 { bad += 1; print("FAIL spearmanrho p-value range\n"); }
    }
    // the alternative splits the p-value: the two-sided value is twice
    // the one-sided value of the same tail
    let two, ct2 = statistics::spearmanrho(x[..], y[..], statistics::ALT_TWO_SIDED);
    let less, cl = statistics::spearmanrho(x[..], y[..], statistics::ALT_LESS);
    let gre, cg = statistics::spearmanrho(x[..], y[..], statistics::ALT_GREATER);
    var tail = less.pvalue;
    if gre.pvalue < tail { tail = gre.pvalue; }
    chk_val(bad, "the two-sided p-value is twice the smaller tail", two.pvalue, ct2 + cl + cg, 2.0 * tail, 1e-12);
    // the xi correlation is the same with and without the continuous law
    let d1, cd1 = statistics::chatterjeexi(x[..], y[..], false);
    let d2, cd2 = statistics::chatterjeexi(x[..], y[..], true);
    chk_val(bad, "the xi statistic is the same", d1.statistic, cd1 + cd2, d2.statistic, 1e-12);
    // the weighted tau of a sample against itself is 1
    let w, cw = statistics::weightedtau(x[..], x[..]);
    chk_val(bad, "weightedtau of a sample against itself", w, cw, 1.0, 1e-12);
    let w2, cw2 = statistics::weightedtau(x[..], x[..], false);
    chk_val(bad, "weightedtau not additive against itself", w2, cw2, 1.0, 1e-12);
    // the average of the two directions makes the call symmetric
    let r1, cr1 = statistics::weightedtau(x[..], y[..]);
    let r2, cr2 = statistics::weightedtau(y[..], x[..]);
    chk_val(bad, "the weighted tau is symmetric", r1 - r2, cr1 + cr2, 0.0, 1e-12);
    // the samples of the caller are unchanged
    var same = true;
    for i in 0..6 { if x[i] != keep[i] { same = false; } }
    check(bad, same, "the rank calls leave the samples unchanged");
    // the point biserial correlation of a list with two values
    let flags: f64[6] = [0.0, 0.0, 1.0, 1.0, 0.0, 1.0];
    let b1, cb1 = statistics::pointbiserialr(flags[..], y[..]);
    chk_code(bad, "pointbiserialr code", cb1, 0);
    if b1.statistic < -1.0 || b1.statistic > 1.0 { bad += 1; print("FAIL pointbiserialr range\n"); }
    return bad;
}
