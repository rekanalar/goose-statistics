// Scenario tests and hostile input tests of the distributions. Every call
// must answer a fault with NaN and a code, answer a value with code 0,
// hold no NaN in a value, give a whole number for a count, and finish in a
// short time on a crafted input.
namespace suite;

import std;
import math;
import os;
import statistics from "../..";

// The count of seeded inputs for each call: 100000, or the number in the
// variable STATISTICS_HOSTILE_N. The mutation probes use a smaller count.
fn dist_count() -> i64 {
    var buf: u8[>..] = [];
    if !env("STATISTICS_HOSTILE_N", buf) { return 100000; }
    var n: i64 = 0;
    for ch in buf {
        if ch < 48 || ch > 57 { return 100000; }
        n = n * 10 + (ch - 48);
    }
    if n < 1 { return 100000; }
    return n;
}

// The count of seeded inputs for the calls that cost a continued fraction,
// a search or a sum. Each one costs a fraction of a millisecond, and a
// full count of them takes longer than the whole test program may run.
fn dist_slow_count() -> i64 {
    var n = dist_count() / 50;
    if n < 200 { return 200; }
    return n;
}

// True when the call costs a continued fraction, a search or a sum. The
// cumulative function and the quantile of every family cost one, and the
// moments of the hypergeometric cost a sum of its masses.
fn dist_is_slow(fam: i64, op: i64) -> bool {
    if op == 1 || op == 2 || op == 3 { return true; }
    if fam == 14 && (op == 6 || op == 7) { return true; }
    return false;
}

// A seeded random number in the unit range.
fn dist_unit(st: u64&) -> f64 {
    return (rnd(st) % 1000001) / 1000000.0;
}

// A seeded random number that is often out of range.
fn dist_any(st: u64&) -> f64 {
    let r = rnd(st) % 8;
    let u = (rnd(st) % 2000001 - 1000000) as f64 / 1000000.0;
    if r == 0 { return u; }
    if r == 1 { return abs(u); }
    if r == 2 { return abs(u) * 1000.0; }
    if r == 3 { return abs(u) * 1e-6; }
    if r == 4 { return u * 1e6; }
    if r == 5 { return 0.0; }
    if r == 6 { return floor(u * 100.0); }
    return (rnd(st) % 1000);
}

// The arguments of one family, in the range that the family needs. The
// first is the value of a density or a mass, the second the probability of
// a quantile, and the last three the parameters.
fn dist_args(fam: i64, op: i64, st: u64&) -> f64, f64, f64, f64, f64 {
    var x = 0.0;
    var q = dist_unit(st);
    var a = 0.0;
    var b = 0.0;
    var c = 0.0;
    if fam == 0 {
        a = dist_any(st) * 3.0;
        b = abs(dist_any(st)) + 0.5;
        x = dist_any(st) * 4.0;
    } else if fam == 1 {
        a = dist_any(st);
        b = a + abs(dist_any(st)) + 0.5;
        x = dist_any(st) * 3.0;
    } else if fam == 2 {
        a = abs(dist_any(st)) + 0.01;
        x = dist_any(st) * 5.0;
    } else if fam == 3 {
        a = dist_any(st);
        b = abs(dist_any(st)) + 0.1;
        x = exp(dist_any(st) * 5.0);
    } else if fam == 4 {
        a = abs(dist_any(st)) + 0.05;
        b = abs(dist_any(st)) + 0.05;
        x = dist_any(st) * 5.0;
    } else if fam == 5 {
        a = abs(dist_any(st)) + 0.05;
        b = abs(dist_any(st)) + 0.05;
        x = dist_any(st) * 5.0;
    } else if fam == 6 {
        a = abs(dist_any(st)) + 0.05;
        b = abs(dist_any(st)) + 0.05;
        x = dist_unit(st);
    } else if fam == 7 {
        a = abs(dist_any(st)) + 0.5;
        x = dist_any(st) * 10.0;
    } else if fam == 8 {
        a = abs(dist_any(st)) + 0.5;
        x = dist_any(st) * 10.0;
    } else if fam == 9 {
        a = abs(dist_any(st)) + 0.5;
        b = abs(dist_any(st)) + 0.5;
        x = dist_any(st) * 10.0;
    } else if fam == 10 {
        a = floor(abs(dist_any(st)) * 50.0);
        b = dist_unit(st);
        x = floor(abs(dist_any(st)) * 50.0);
    } else if fam == 11 {
        a = abs(dist_any(st)) * 20.0;
        x = floor(abs(dist_any(st)) * 50.0);
    } else if fam == 12 {
        a = dist_unit(st) * 0.999 + 0.001;
        x = floor(abs(dist_any(st)) * 50.0) + 1.0;
    } else if fam == 13 {
        a = abs(dist_any(st)) * 20.0 + 0.5;
        b = dist_unit(st) * 0.999 + 0.001;
        x = floor(abs(dist_any(st)) * 50.0);
    } else {
        a = floor(abs(dist_any(st)) * 60.0) + 1.0;
        b = floor(dist_unit(st) * a);
        c = floor(dist_unit(st) * a);
        x = floor(abs(dist_any(st)) * 60.0);
    }
    if op == 3 { q = dist_unit(st); }
    return x, q, a, b, c;
}

// The seeded inputs of every call. A fault must hold NaN. A value must
// hold no NaN. A count of a discrete family must be a whole number.
fn hostile_distributions() -> i64 {
    var bad: i64 = 0;
    for fam in 0..15 {
        for op in 0..8 {
            print("hostile ", dist_family_name(fam), " ", dist_op_name(op), "\n");
            var st: u64 = 20261112 + (fam as u64) * 7919 + (op as u64) * 104729;
            var faults: i64 = 0;
            var values: i64 = 0;
            var count = dist_count();
            if dist_is_slow(fam, op) { count = dist_slow_count(); }
            for i in 0..count {
                let x, q, a, b, c = dist_args(fam, op, st);
                let v, code = dist_call(fam, op, x, q, a, b, c);
                var ok = code >= 0 && code <= 9;
                if code != 0 {
                    faults += 1;
                    ok = ok && is_nan(v);
                } else {
                    values += 1;
                    ok = ok && !is_nan(v);
                    // only the quantile of a discrete family is a count
                    if fam >= 10 && op == 3 {
                        ok = ok && floor(v) == v;
                    }
                }
                if !ok {
                    bad += 1;
                    if bad <= 20 {
                        print("FAIL ", dist_family_name(fam), " ", dist_op_name(op), " x ", x, " q ", q, " a ", a, " b ", b, " c ", c, " gave ", v, " code ", code, "\n");
                    }
                }
            }
            if values == 0 {
                bad += 1;
                print("FAIL ", dist_family_name(fam), " ", dist_op_name(op), " gave no value\n");
            }
        }
    }
    return bad;
}

// The crafted inputs: the shapes at the bounds, the end points of the
// probability range and the long supports. Each call must take at most
// fifty times the time of a random call.
fn dist_time() -> i64 {
    var bad: i64 = 0;
    let t0 = clock();
    for i in 0..50 {
        let v, c = statistics::norm_cdf(0.5 + i * 0.01);
    }
    let base = (clock() - t0) / 50.0;
    let shapes: const f64[8] = [1.0, 2.0, 3.0, 10.0, 100.0, 1e4, 1e8, 1e10];
    for k in 0..8 {
        let s = shapes[k];
        let t1 = clock();
        let _a, c1 = statistics::gammadist_cdf(s, s);
        let _b, c2 = statistics::betadist_cdf(0.5, s, 0.3);
        let _c, c3 = statistics::tdist_cdf(s, s);
        let _d, c4 = statistics::fdist_cdf(s, s, s);
        let _e, c5 = statistics::hypergeom_cdf(3.0, 200.0, 100.0, 60.0);
        let _f, c6 = statistics::poisson_ppf(0.5, min(s, 1e6));
        let spent = clock() - t1;
        if spent > 0.5 + 50.0 * base {
            bad += 1;
            print("FAIL time at the shape ", s, ": ", spent, " seconds\n");
        }
        if c1 != 0 || c2 != 0 || c3 != 0 || c4 != 0 || c5 != 0 || c6 != 0 {
            bad += 1;
            print("FAIL codes at the shape ", s, ": ", c1, c2, c3, c4, c5, c6, "\n");
        }
    }
    return bad;
}

// The fault codes of the calls that must reject an argument.
fn dist_cases() -> i64 {
    var bad: i64 = 0;
    let nan = 0.0 / 0.0;
    let inf = 1.0 / 0.0;
    // a NaN gives the NaN code
    let v1, c1 = statistics::norm_pdf(nan);
    if c1 != 3 || !is_nan(v1) { bad += 1; print("FAIL norm_pdf of a NaN: ", v1, " ", c1, "\n"); }
    let v2, c2 = statistics::binom_pmf(nan, 10.0, 0.5);
    if c2 != 3 || !is_nan(v2) { bad += 1; print("FAIL binom_pmf of a NaN: ", v2, " ", c2, "\n"); }
    // an infinity gives the infinity code
    let v3, c3 = statistics::gammadist_pdf(inf, 1.0, 1.0);
    if c3 != 4 || !is_nan(v3) { bad += 1; print("FAIL gammadist_pdf of an infinity: ", v3, " ", c3, "\n"); }
    let v4, c4 = statistics::expon_pdf(1.0, inf);
    if c4 != 4 || !is_nan(v4) { bad += 1; print("FAIL expon_pdf of an infinite rate: ", v4, " ", c4, "\n"); }
    // a shape at or below 0 gives the domain code
    let v5, c5 = statistics::gammadist_pdf(1.0, 0.0, 1.0);
    if c5 != 5 || !is_nan(v5) { bad += 1; print("FAIL gammadist_pdf of a zero shape: ", v5, " ", c5, "\n"); }
    let v6, c6 = statistics::betadist_cdf(0.5, -1.0, 1.0);
    if c6 != 5 || !is_nan(v6) { bad += 1; print("FAIL betadist_cdf of a negative shape: ", v6, " ", c6, "\n"); }
    // a scale of 0 gives the domain code
    let v7, c7 = statistics::norm_pdf(1.0, 0.0, 0.0);
    if c7 != 5 || !is_nan(v7) { bad += 1; print("FAIL norm_pdf of a zero scale: ", v7, " ", c7, "\n"); }
    // degrees of freedom below 1 give the domain code
    let v8, c8 = statistics::tdist_cdf(1.0, 0.5);
    if c8 != 5 || !is_nan(v8) { bad += 1; print("FAIL tdist_cdf of a small df: ", v8, " ", c8, "\n"); }
    // a probability outside the unit range gives the domain code
    let v9, c9 = statistics::norm_ppf(1.5);
    if c9 != 5 || !is_nan(v9) { bad += 1; print("FAIL norm_ppf of 1.5: ", v9, " ", c9, "\n"); }
    let v10, c10 = statistics::binom_ppf(-0.1, 10.0, 0.5);
    if c10 != 5 || !is_nan(v10) { bad += 1; print("FAIL binom_ppf of a negative probability: ", v10, " ", c10, "\n"); }
    // a count with a fraction gives the domain code
    let v11, c11 = statistics::hypergeom_pmf(1.5, 10.0, 5.0, 3.0);
    if c11 != 5 || !is_nan(v11) { bad += 1; print("FAIL hypergeom_pmf of a fraction: ", v11, " ", c11, "\n"); }
    // a draw count above the population gives the domain code
    let v12, c12 = statistics::hypergeom_pmf(1.0, 10.0, 5.0, 20.0);
    if c12 != 5 || !is_nan(v12) { bad += 1; print("FAIL hypergeom_pmf above the population: ", v12, " ", c12, "\n"); }
    // a shape with no moment gives the undefined code
    let v13, c13 = statistics::tdist_var(2.0);
    if c13 != 7 || !is_nan(v13) { bad += 1; print("FAIL tdist_var at df 2: ", v13, " ", c13, "\n"); }
    let v14, c14 = statistics::fdist_kurt(8.0, 8.0);
    if c14 != 7 || !is_nan(v14) { bad += 1; print("FAIL fdist_kurt at df2 8: ", v14, " ", c14, "\n"); }
    // the moment poles sit at df 4 for the F variance and the t kurtosis
    let v21, c21 = statistics::fdist_var(3.0, 4.0);
    if c21 != 7 || !is_nan(v21) { bad += 1; print("FAIL fdist_var at df2 4: ", v21, " ", c21, "\n"); }
    let v22, c22 = statistics::tdist_kurt(4.0);
    if c22 != 7 || !is_nan(v22) { bad += 1; print("FAIL tdist_kurt at df 4: ", v22, " ", c22, "\n"); }
    // a Weibull mean with a large gamma argument still holds: 100 factorial
    let v23, c23 = statistics::weibull_mean(0.01, 1.0);
    if c23 != 0 || abs(v23 / 9.332621544394415e157 - 1.0) > 1e-13 { bad += 1; print("FAIL weibull_mean at shape 0.01: ", v23, " ", c23, "\n"); }
    // the probability 0 gives the smallest count of the support
    let v24, c24 = statistics::geom_ppf(0.0, 0.5);
    if c24 != 0 || v24 != 1.0 { bad += 1; print("FAIL geom_ppf at 0: ", v24, " ", c24, "\n"); }
    let v25, c25 = statistics::binom_ppf(0.0, 10.0, 0.5);
    if c25 != 0 || v25 != 0.0 { bad += 1; print("FAIL binom_ppf at 0: ", v25, " ", c25, "\n"); }
    // a support above the bound gives the too long code
    let v15, c15 = statistics::hypergeom_cdf(1.0, 400000.0, 200000.0, 200000.0);
    if c15 != 9 || !is_nan(v15) { bad += 1; print("FAIL hypergeom_cdf above the bound: ", v15, " ", c15, "\n"); }
    // the values at the ends of the support
    let v16, c16 = statistics::binom_pmf(0.0, 0.0, 0.5);
    if c16 != 0 || v16 != 1.0 { bad += 1; print("FAIL binom_pmf at n 0: ", v16, " ", c16, "\n"); }
    let v17, c17 = statistics::geom_pmf(0.0, 0.5);
    if c17 != 0 || v17 != 0.0 { bad += 1; print("FAIL geom_pmf at 0: ", v17, " ", c17, "\n"); }
    let v18, c18 = statistics::hypergeom_pmf(0.0, 20.0, 10.0, 4.0);
    if c18 != 0 { bad += 1; print("FAIL hypergeom_pmf at the mode: ", v18, " ", c18, "\n"); }
    // the quantile of the discrete families is a whole number
    let v19, c19 = statistics::poisson_ppf(0.5, 3.0);
    if c19 != 0 || floor(v19) != v19 { bad += 1; print("FAIL poisson_ppf is not a count: ", v19, " ", c19, "\n"); }
    let v20, c20 = statistics::hypergeom_ppf(0.5, 20.0, 10.0, 4.0);
    if c20 != 0 || floor(v20) != v20 { bad += 1; print("FAIL hypergeom_ppf is not a count: ", v20, " ", c20, "\n"); }
    return bad;
}

// A call that returns a fault code 0 must return a finite value. The grid
// holds the extreme arguments that once gave NaN or an infinity with no fault.
fn extreme_grid() -> i64 {
    var bad: i64 = 0;
    let g: f64[8] = [0.0, -1.0, 1e-300, 1e-5, 0.5, 3.0, 1e10, 1e300];
    for fam in 0..10 {
        for op in 0..8 {
            for i in 0..8 {
                for j in 0..8 {
                    for k in 0..8 {
                        let v, c = dist_call(fam, op, g[i], 0.5, g[j], g[k], 0.0);
                        if c == 0 && (is_nan(v) || is_inf(v)) {
                            bad += 1;
                            if bad < 20 { print("FAIL extreme ", dist_family_name(fam), " ", dist_op_name(op), " ", g[i], " ", g[j], " ", g[k], " gives ", v, " with no fault\n"); }
                        }
                    }
                }
            }
        }
    }
    return bad;
}

// The cases of the security review.
fn review_cases() -> i64 {
    var bad: i64 = 0;
    let inf = 1.0 / 0.0;
    // a pair product above 1e299 keeps a finite result
    let a1, b1 = statistics::norm_pdf(0.0, 1.0, 1e-300);
    if b1 != 0 || a1 != 0.0 { bad += 1; print("FAIL norm_pdf with a huge z: ", a1, " ", b1, "\n"); }
    let a2, b2 = statistics::tdist_pdf(1e300, 1.0);
    if b2 != 0 || a2 != 0.0 { bad += 1; print("FAIL tdist_pdf at 1e300: ", a2, " ", b2, "\n"); }
    let a3, b3 = statistics::lognorm_cdf(1.0, 1.0, 1e-300);
    if b3 != 0 || a3 != 0.0 { bad += 1; print("FAIL lognorm_cdf with a tiny sigma: ", a3, " ", b3, "\n"); }
    // a quantile that does not fit gives the overflow code
    let a4, b4 = statistics::tdist_ppf(1e-300, 1.0);
    if b4 != 8 || !is_nan(a4) { bad += 1; print("FAIL tdist_ppf of 1e-300: ", a4, " ", b4, "\n"); }
    // the moments of the beta and F families keep their digits
    let a5, b5 = statistics::betadist_var(1e-300, 1e-300);
    if b5 != 0 || abs(a5 - 0.25) > 1e-15 { bad += 1; print("FAIL betadist_var of tiny shapes: ", a5, " ", b5, "\n"); }
    let a6, b6 = statistics::fdist_var(1e300, 1e300);
    if b6 != 0 || abs(a6 / 4e-300 - 1.0) > 1e-14 { bad += 1; print("FAIL fdist_var of huge degrees: ", a6, " ", b6, "\n"); }
    let a7, b7 = statistics::betadist_skew(1e300, 1e300);
    if b7 != 0 || a7 != 0.0 { bad += 1; print("FAIL betadist_skew of equal shapes: ", a7, " ", b7, "\n"); }
    // a density with an infinite value gives the overflow code
    let a8, b8 = statistics::weibull_pdf(0.0, 0.5, 1.0);
    if b8 != 8 || !is_nan(a8) { bad += 1; print("FAIL weibull_pdf at the pole: ", a8, " ", b8, "\n"); }
    let a9, b9 = statistics::gammadist_pdf(0.0, 0.5, 1.0);
    if b9 != 8 || !is_nan(a9) { bad += 1; print("FAIL gammadist_pdf at the pole: ", a9, " ", b9, "\n"); }
    let a10, b10 = statistics::chi2_pdf(0.0, 1.0);
    if b10 != 8 || !is_nan(a10) { bad += 1; print("FAIL chi2_pdf at the pole: ", a10, " ", b10, "\n"); }
    let a11, b11 = statistics::betadist_pdf(0.0, 0.5, 1.0);
    if b11 != 8 || !is_nan(a11) { bad += 1; print("FAIL betadist_pdf at the pole: ", a11, " ", b11, "\n"); }
    let a12, b12 = statistics::fdist_pdf(0.0, 1.0, 1.0);
    if b12 != 8 || !is_nan(a12) { bad += 1; print("FAIL fdist_pdf at the pole: ", a12, " ", b12, "\n"); }
    // an infinite scale gives the infinity code in every call
    let a13, b13 = statistics::norm_pdf(1.0, 0.0, inf);
    if b13 != 4 || !is_nan(a13) { bad += 1; print("FAIL norm_pdf of an infinite sigma: ", a13, " ", b13, "\n"); }
    // an infinite x gives the infinity code in the cumulative calls
    let a14, b14 = statistics::fdist_cdf(inf, 1.0, 1.0);
    if b14 != 4 || !is_nan(a14) { bad += 1; print("FAIL fdist_cdf of an infinity: ", a14, " ", b14, "\n"); }
    let a15, b15 = statistics::betadist_sf(-inf, 1.0, 1.0);
    if b15 != 4 || !is_nan(a15) { bad += 1; print("FAIL betadist_sf of an infinity: ", a15, " ", b15, "\n"); }
    // the infinity comes before the domain
    let a16, b16 = statistics::gammadist_pdf(inf, 0.0, 1.0);
    if b16 != 4 { bad += 1; print("FAIL gammadist_pdf infinity before domain: ", b16, "\n"); }
    let a17, b17 = statistics::tdist_cdf(inf, 0.5);
    if b17 != 4 { bad += 1; print("FAIL tdist_cdf infinity before domain: ", b17, "\n"); }
    let a18, b18 = statistics::lognorm_pdf(inf, 0.0, 0.0);
    if b18 != 4 { bad += 1; print("FAIL lognorm_pdf infinity before domain: ", b18, "\n"); }
    // the NaN comes before the infinity
    let a19, b19 = statistics::weibull_pdf(inf, nan_value(), 1.0);
    if b19 != 3 { bad += 1; print("FAIL weibull_pdf NaN before infinity: ", b19, "\n"); }
    // the tails of the incomplete beta at the ends of the range
    let a20, b20 = statistics::betainc_upper(2.0, 3.0, 1.0);
    if b20 != 0 || a20 != 0.0 { bad += 1; print("FAIL betainc_upper at 1: ", a20, " ", b20, "\n"); }
    let a21, b21 = statistics::betainc_upper(2.0, 3.0, 0.0);
    if b21 != 0 || a21 != 1.0 { bad += 1; print("FAIL betainc_upper at 0: ", a21, " ", b21, "\n"); }
    return bad;
}

fn nan_value() -> f64 { return 0.0 / 0.0; }

fn run_distribution_extra() -> i64 {
    var bad: i64 = 0;
    bad += dist_cases();
    bad += review_cases();
    bad += extreme_grid();
    bad += dist_time();
    bad += hostile_distributions();
    return bad;
}
