// The tests of the fitting calls. The vectors live in fitting_vectors.goose.
namespace suite;

import std;
import math;
import statistics from "../..";

fn cens(u: const f64[:], l: const f64[:], r: const f64[:], iv: const f64[:]) -> statistics::CensoredData {
    return statistics::CensoredData { uncensored: u, left: l, right: r, interval: iv };
}

fn run_fitting_tests() -> i64 {
    var bad: i64 = 0;
    var none: f64[>..] = [];
    let ux: f64[4] = [2.0, 3.0, 5.0, 8.0];
    let rx: f64[3] = [4.0, 6.0, 9.0];
    let uy: f64[4] = [1.0, 2.0, 2.0, 4.0];
    let ry: f64[3] = [3.0, 5.0, 7.0];
    let x = cens(ux[..], none[..], rx[..], none[..]);
    let y = cens(uy[..], none[..], ry[..], none[..]);

    // logrank: the measured pair, and the faults in the fixed order.
    let r1, c1 = statistics::logrank(x, y);
    check(bad, c1 == 0 && near_tol(r1.statistic, -0.8042527331835678, 1e-12)
        && near_tol(r1.pvalue, 0.42125102949556004, 1e-12), "logrank, the measured pair");
    check(bad, is_nan(r1.df) && is_nan(r1.df2), "logrank, df stays NaN");
    let rl, cl = statistics::logrank(x, y, statistics::ALT_LESS);
    let rg, cg = statistics::logrank(x, y, statistics::ALT_GREATER);
    check(bad, cl == 0 && cg == 0 && near_tol(rl.pvalue + rg.pvalue, 1.0, 1e-12)
        && near_tol(r1.pvalue, 2.0 * min(rl.pvalue, rg.pvalue), 1e-12), "logrank, the one-sided p-values");
    let ra, ca = statistics::logrank(y, x);
    check(bad, ca == 0 && near_tol(ra.statistic, -r1.statistic, 1e-12), "logrank, swapped groups flip the sign");
    let _e1, ce1 = statistics::logrank(x, y, 3);
    check(bad, ce1 == 5, "logrank, an unknown alternative");
    let _e1b, ce1b = statistics::logrank(x, y, -1);
    check(bad, ce1b == 5, "logrank, a negative alternative");

    let empty = cens(none[..], none[..], none[..], none[..]);
    let e2, ce2 = statistics::logrank(empty, y);
    check(bad, ce2 == 1 && is_nan(e2.statistic) && is_nan(e2.pvalue), "logrank, an empty group");
    let _e3, ce3 = statistics::logrank(x, empty);
    check(bad, ce3 == 1, "logrank, an empty second group");

    let nanu: f64[2] = [1.0, statistics::nan()];
    let infr: f64[2] = [2.0, 1.0 / 0.0];
    let _e4, ce4 = statistics::logrank(cens(nanu[..], none[..], rx[..], none[..]), y);
    check(bad, ce4 == 3, "logrank, NaN in the uncensored part");
    let _e5, ce5 = statistics::logrank(x, cens(uy[..], none[..], infr[..], none[..]));
    check(bad, ce5 == 4, "logrank, infinity in the right part");
    let _e5b, ce5b = statistics::logrank(cens(infr[..], none[..], rx[..], none[..]), cens(nanu[..], none[..], ry[..], none[..]));
    check(bad, ce5b == 3, "logrank, NaN outranks infinity");

    let left: f64[1] = [1.5];
    let _e6, ce6 = statistics::logrank(cens(ux[..], left[..], rx[..], none[..]), y);
    check(bad, ce6 == 5, "logrank, a left part in x");
    let _e6b, ce6b = statistics::logrank(x, cens(uy[..], left[..], ry[..], none[..]));
    check(bad, ce6b == 5, "logrank, a left part in y");
    let pair: f64[2] = [1.0, 2.0];
    let _e7, ce7 = statistics::logrank(x, cens(uy[..], none[..], ry[..], pair[..]));
    check(bad, ce7 == 5, "logrank, an interval part");
    let back: f64[2] = [3.0, 2.0];
    let _e8, ce8 = statistics::logrank(x, cens(uy[..], none[..], ry[..], back[..]));
    check(bad, ce8 == 5, "logrank, an interval pair with low above high");
    let odd: f64[3] = [1.0, 2.0, 3.0];
    let _e9, ce9 = statistics::logrank(x, cens(uy[..], none[..], ry[..], odd[..]));
    check(bad, ce9 == 6, "logrank, an interval list of odd length");
    let nint: f64[2] = [1.0, statistics::nan()];
    let _e9b, ce9b = statistics::logrank(x, cens(uy[..], none[..], ry[..], nint[..]));
    check(bad, ce9b == 3, "logrank, NaN in the interval part");

    // No event at all, and no variance, give the undefined fault.
    let _e10, ce10 = statistics::logrank(cens(none[..], none[..], rx[..], none[..]), cens(none[..], none[..], ry[..], none[..]));
    check(bad, ce10 == 7, "logrank, no event");
    let one: f64[1] = [2.0];
    let early: f64[1] = [1.0];
    let _e11, ce11 = statistics::logrank(cens(one[..], none[..], none[..], none[..]), cens(none[..], none[..], early[..], none[..]));
    check(bad, ce11 == 7, "logrank, one event and no variance");
    let e12, ce12a = statistics::logrank(cens(one[..], none[..], none[..], none[..]), cens(none[..], none[..], one[..], none[..]));
    check(bad, ce12a == 0 && near_tol(e12.statistic, 1.0, 1e-12) && near_tol(e12.pvalue, 0.31731050786291415, 1e-12), "logrank, one event and one tied censored time");

    // The length bound comes before the other checks.
    var big: f64[>..] = [];
    for i in 0..4194305 { big.push(i); }
    let _e12, cb12 = statistics::logrank(cens(big[..], none[..], none[..], none[..]), y);
    check(bad, cb12 == 9, "logrank, a part above the bound");
    let _e13, ce13 = statistics::logrank(x, cens(none[..], none[..], none[..], big[..]));
    check(bad, ce13 == 9, "logrank, an interval part above the bound");
    big.resize(0);

    // No call writes into the lists of the caller.
    let ux2 = ux;
    let ry2 = ry;
    let _n1, _n2 = statistics::logrank(x, y);
    var same = true;
    for i in 0..4 { if ux[i] != ux2[i] { same = false; } }
    for i in 0..3 { if ry[i] != ry2[i] { same = false; } }
    check(bad, same, "logrank, no write into the lists");

    bad += run_gof_tests();
    return bad;
}

fn run_gof_tests() -> i64 {
    var bad: i64 = 0;
    var out: f64[>..] = [];
    let normal: f64[10] = [9.1, 10.4, 11.2, 8.7, 10.9, 9.6, 12.1, 10.2, 9.9, 11.0];
    var g = rng { s: 7 };

    // The normal family: the fit goes into out, and the p-value stays above 0.05.
    let r1, c1 = statistics::goodness_of_fit(out, normal[..], statistics::FAMILY_NORMAL, statistics::GOF_ANDERSON, 999, g);
    let m, _cm = statistics::mean(normal[..]);
    let sd, _cs = statistics::stdev(normal[..]);
    check(bad, c1 == 0 && out.len == 2 && near_tol(out[0], m, 1e-12) && near_tol(out[1], sd, 1e-12),
        "goodness_of_fit, the fit of the normal family");
    check(bad, c1 == 0 && r1.pvalue > 0.05 && r1.pvalue <= 1.0 && r1.statistic > 0.0
        && is_nan(r1.df) && is_nan(r1.df2), "goodness_of_fit, normal data keep a large p-value");

    // The same seed gives the same answer. Another seed gives another stream.
    var g1 = rng { s: 99 };
    var g2 = rng { s: 99 };
    var o2: f64[>..] = [];
    let a1, ca1 = statistics::goodness_of_fit(out, normal[..], statistics::FAMILY_NORMAL, statistics::GOF_CRAMER, 500, g1);
    let a2, ca2 = statistics::goodness_of_fit(o2, normal[..], statistics::FAMILY_NORMAL, statistics::GOF_CRAMER, 500, g2);
    check(bad, ca1 == 0 && ca2 == 0 && a1.pvalue == a2.pvalue && a1.statistic == a2.statistic, "goodness_of_fit, one seed one answer");
    // A seed of 0 is not special.
    var g0 = rng { s: 0 };
    let z1, cz1 = statistics::goodness_of_fit(out, normal[..], statistics::FAMILY_NORMAL, statistics::GOF_KOLMOGOROV, 200, g0);
    check(bad, cz1 == 0 && z1.pvalue > 0.0, "goodness_of_fit, a seed of 0");

    // The fits of the other families, and the size of out.
    let pos: f64[6] = [0.5, 1.5, 0.8, 2.9, 1.1, 0.4];
    let r3, c3 = statistics::goodness_of_fit(out, pos[..], statistics::FAMILY_EXPONENTIAL, statistics::GOF_CRAMER, 200, g);
    let pm, _pc = statistics::mean(pos[..]);
    check(bad, c3 == 0 && out.len == 2 && near_tol(out[0], 0.0, 1e-12) && near_tol(out[1], pm, 1e-12), "goodness_of_fit, the scale of the exponential family");
    let r4, c4 = statistics::goodness_of_fit(out, pos[..], statistics::FAMILY_LOGNORMAL, statistics::GOF_ANDERSON, 200, g);
    var lm: f64 = 0.0;
    var ls: f64 = 0.0;
    for v in pos { lm += log(v) / 6.0; }
    for v in pos { ls += (log(v) - lm) * (log(v) - lm) / 6.0; }
    check(bad, c4 == 0 && out.len == 3 && near_tol(out[0], sqrt(ls), 1e-12) && out[1] == 0.0 && near_tol(out[2], exp(lm), 1e-12), "goodness_of_fit, the lognormal fit");
    let r5, c5 = statistics::goodness_of_fit(out, pos[..], statistics::FAMILY_UNIFORM, statistics::GOF_KOLMOGOROV, 200, g);
    check(bad, c5 == 0 && out.len == 2 && out[0] == 0.4 && near_tol(out[1], 2.5, 1e-12), "goodness_of_fit, the location and the width of the uniform family");
    check(bad, r3.pvalue > 0.0 && r4.pvalue > 0.0 && r5.pvalue > 0.0, "goodness_of_fit, p-values stay positive");

    // A sample far from the family gives a small p-value.
    let skew: f64[12] = [0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.5, 9.0, 12.0, 20.0, 31.0];
    let r6, c6 = statistics::goodness_of_fit(out, skew[..], statistics::FAMILY_NORMAL, statistics::GOF_ANDERSON, 999, g);
    check(bad, c6 == 0 && r6.pvalue < 0.05, "goodness_of_fit, skewed data against the normal family");

    // The faults. The out buffer is empty after each of them.
    let c7: i64 = 0;
    var junk: f64[>..] = [1.0, 2.0];
    let neg: f64[3] = [1.0, -1.0, 2.0];
    let _f1, cf1 = statistics::goodness_of_fit(junk, neg[..], statistics::FAMILY_GAMMA, statistics::GOF_ANDERSON, 99, g);
    check(bad, cf1 == 5 && junk.len == 0, "goodness_of_fit, the gamma family needs positive values");
    for fam in [7, 8, 9] {
        let _f, cf = statistics::goodness_of_fit(out, normal[..], fam, statistics::GOF_ANDERSON, 99, g);
        check(bad, cf == 5, "goodness_of_fit, a family with no fit");
    }
    let _f2, cf2 = statistics::goodness_of_fit(out, normal[..], 99, statistics::GOF_ANDERSON, 99, g);
    check(bad, cf2 == 5, "goodness_of_fit, an unknown family");
    let _f2b, cf2b = statistics::goodness_of_fit(out, normal[..], -1, statistics::GOF_ANDERSON, 99, g);
    check(bad, cf2b == 5, "goodness_of_fit, a negative family");
    let _f3, cf3 = statistics::goodness_of_fit(out, normal[..], statistics::FAMILY_NORMAL, 3, 99, g);
    check(bad, cf3 == 5, "goodness_of_fit, an unknown statistic");
    let _f3b, cf3b = statistics::goodness_of_fit(out, normal[..], statistics::FAMILY_NORMAL, -1, 99, g);
    check(bad, cf3b == 5, "goodness_of_fit, a negative statistic");
    let _f4, cf4 = statistics::goodness_of_fit(out, normal[..], statistics::FAMILY_NORMAL, statistics::GOF_ANDERSON, 0, g);
    check(bad, cf4 == 5, "goodness_of_fit, no resample");
    let _f5, cf5 = statistics::goodness_of_fit(out, normal[..], statistics::FAMILY_NORMAL, statistics::GOF_ANDERSON, 1000001, g);
    check(bad, cf5 == 9, "goodness_of_fit, too many resamples");
    let _f5b, cf5b = statistics::goodness_of_fit(out, normal[..], statistics::FAMILY_NORMAL, statistics::GOF_ANDERSON, 9223372036854775807, g);
    check(bad, cf5b == 9, "goodness_of_fit, the largest count of resamples");
    // The work bound of the fit: the resamples times the sample size.
    var big: f64[>..] = [];
    for i in 0..1000 { big.push((i * 37) % 101 + 0.25 * (i % 7)); }
    let _f5c, cf5c = statistics::goodness_of_fit(out, big[..], statistics::FAMILY_NORMAL, statistics::GOF_KOLMOGOROV, 2001, g);
    check(bad, cf5c == 9, "goodness_of_fit, the work bound of the fit");
    check(bad, out.len == 0, "goodness_of_fit, the work bound leaves out empty");
    var none: f64[>..] = [];
    let _f6, cf6 = statistics::goodness_of_fit(out, none[..], statistics::FAMILY_NORMAL, statistics::GOF_ANDERSON, 99, g);
    check(bad, cf6 == 1, "goodness_of_fit, an empty sample");
    let nanv: f64[3] = [1.0, statistics::nan(), 2.0];
    let _f7, cf7 = statistics::goodness_of_fit(out, nanv[..], statistics::FAMILY_NORMAL, statistics::GOF_ANDERSON, 99, g);
    check(bad, cf7 == 3, "goodness_of_fit, NaN");
    let infv: f64[3] = [1.0, 1.0 / 0.0, 2.0];
    let _f8, cf8 = statistics::goodness_of_fit(out, infv[..], statistics::FAMILY_NORMAL, statistics::GOF_ANDERSON, 99, g);
    check(bad, cf8 == 4, "goodness_of_fit, infinity");
    let one: f64[1] = [1.0];
    let _f9, cf9 = statistics::goodness_of_fit(out, one[..], statistics::FAMILY_NORMAL, statistics::GOF_ANDERSON, 99, g);
    check(bad, cf9 == 2, "goodness_of_fit, one value");
    let constant: f64[4] = [2.0, 2.0, 2.0, 2.0];
    let _f10, cf10 = statistics::goodness_of_fit(out, constant[..], statistics::FAMILY_NORMAL, statistics::GOF_ANDERSON, 99, g);
    check(bad, cf10 == 7 && out.len == 0, "goodness_of_fit, a constant sample");
    let neg: f64[3] = [1.0, -2.0, 3.0];
    let _f11, cf11 = statistics::goodness_of_fit(out, neg[..], statistics::FAMILY_EXPONENTIAL, statistics::GOF_ANDERSON, 99, g);
    check(bad, cf11 == 5, "goodness_of_fit, a negative value for the exponential family");
    let _f12, cf12 = statistics::goodness_of_fit(out, neg[..], statistics::FAMILY_LOGNORMAL, statistics::GOF_ANDERSON, 99, g);
    check(bad, cf12 == 5, "goodness_of_fit, a negative value for the lognormal family");
    let zero: f64[3] = [0.0, 2.0, 3.0];
    let _f13, cf13 = statistics::goodness_of_fit(out, zero[..], statistics::FAMILY_EXPONENTIAL, statistics::GOF_ANDERSON, 99, g);
    check(bad, cf13 == 5, "goodness_of_fit, a zero for the exponential family");
    let _f14, cf14 = statistics::goodness_of_fit(out, pos[..], statistics::FAMILY_UNIFORM, statistics::GOF_ANDERSON, 99, g);
    check(bad, cf14 == 7, "goodness_of_fit, the Anderson statistic of the uniform family");

    // A fit that cannot end gives the fit fault: the sum overflows.
    let huge: f64[3] = [1.7e308, 1.6e308, 1.7e308];
    let _f15, cf15 = statistics::goodness_of_fit(out, huge[..], statistics::FAMILY_EXPONENTIAL, statistics::GOF_CRAMER, 99, g);
    check(bad, cf15 == 10 && out.len == 0, "goodness_of_fit, an exponential fit that overflows");
    let _f16, cf16 = statistics::goodness_of_fit(out, huge[..], statistics::FAMILY_NORMAL, statistics::GOF_CRAMER, 99, g);
    check(bad, cf16 == 10, "goodness_of_fit, a normal fit that overflows");
    let wide: f64[2] = [-1.7e308, 1.7e308];
    let _f18, cf18 = statistics::goodness_of_fit(out, wide[..], statistics::FAMILY_UNIFORM, statistics::GOF_KOLMOGOROV, 99, g);
    check(bad, cf18 == 10, "goodness_of_fit, a uniform width that overflows");

    // The sample is not changed.
    let keep: f64[10] = normal;
    let _k, _kc = statistics::goodness_of_fit(out, normal[..], statistics::FAMILY_NORMAL, statistics::GOF_KOLMOGOROV, 50, g);
    var same = true;
    for i in 0..10 { if keep[i] != normal[i] { same = false; } }
    check(bad, same && c7 == 0, "goodness_of_fit, no write into the sample");

    // The length bound comes first, before the work.
    var big: f64[>..] = [];
    for i in 0..4194305 { big.push(i); }
    let _f19, cf19 = statistics::goodness_of_fit(out, big[..], statistics::FAMILY_NORMAL, statistics::GOF_ANDERSON, 99, g);
    check(bad, cf19 == 9, "goodness_of_fit, a sample above the bound");
    return bad;
}
