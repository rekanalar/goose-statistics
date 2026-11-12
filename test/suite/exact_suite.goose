// The tests of the exact calls: the oracle rows, the tied rows, the faults
// and the no-write checks. The random rows live in exact_vectors.goose,
// which the generator test/make_exact_vectors.py writes.
namespace suite;

import std;
import math;
import statistics from "../..";

// One result against the measured row.
fn ex_res(bad: i64&, name: u8[:], r: statistics::TestResult, c: i64, stat: f64, p: f64, ts: f64, tp: f64) {
    check(bad, c == 0 && near_tol(r.statistic, stat, ts) && near_tol(r.pvalue, p, tp), name);
}

// The Poisson window of this size is 10 to 40 values wide, so no fault shows.
fn run_exact_tests() -> i64 {
    var bad: i64 = 0;
    let a: f64[8] = [2.1, 3.4, 1.9, 5.6, 4.4, 3.8, 2.7, 4.9];
    let b: f64[8] = [3.9, 5.2, 4.8, 6.1, 5.5, 4.1, 6.6, 5.0];
    let c3: f64[8] = [1.2, 2.5, 2.2, 3.1, 1.8, 2.9, 2.4, 3.3];
    let np: f64[2] = [3.5, 1.5];
    let nan: f64 = 0.0 / 0.0;
    let inf: f64 = 1.0 / 0.0;

    // The oracle rows of the 2 by 2 tables.
    {
        let r, c = statistics::barnard_exact(7.0, 3.0, 2.0, 8.0);
        ex_res(bad, "barnard two-sided", r, c, 2.247332874877474, 0.03345870971679689, 1e-12, 1e-9);
        let rg, cg = statistics::barnard_exact(7.0, 3.0, 2.0, 8.0, statistics::ALT_GREATER);
        ex_res(bad, "barnard greater", rg, cg, 2.247332874877474, 0.016786917914384802, 1e-12, 1e-9);
        let rl, cl = statistics::barnard_exact(7.0, 3.0, 2.0, 8.0, statistics::ALT_LESS);
        ex_res(bad, "barnard less", rl, cl, 2.247332874877474, 1.0, 1e-12, 1e-9);
        let rd, cd = statistics::barnard_exact(7.0, 12.0, 8.0, 3.0, statistics::ALT_LESS);
        ex_res(bad, "barnard documented row", rd, cd, -1.8943380760602064, 0.034076716365767874, 1e-12, 1e-9);
        let s, cs = statistics::boschloo_exact(7.0, 3.0, 2.0, 8.0);
        ex_res(bad, "boschloo two-sided", s, cs, 0.03488925934746368, 0.03196459926398076, 1e-12, 1e-9);
        let sg, csg = statistics::boschloo_exact(7.0, 3.0, 2.0, 8.0, statistics::ALT_GREATER);
        ex_res(bad, "boschloo greater", sg, csg, 0.03488925934746368, 0.01598229963199038, 1e-12, 1e-9);
        let sl, csl = statistics::boschloo_exact(7.0, 3.0, 2.0, 8.0, statistics::ALT_LESS);
        ex_res(bad, "boschloo less", sl, csl, 0.9972612526792093, 0.9883427185178026, 1e-12, 1e-9);
        let sd, csd = statistics::boschloo_exact(74.0, 31.0, 43.0, 32.0, statistics::ALT_GREATER);
        check(bad, csd == 0 && near_tol(sd.statistic, 0.0483, 2e-3) && near_tol(sd.pvalue, 0.0355, 2e-3), "boschloo documented row");
        // A column of zeros: Barnard gives a p-value of 1, Boschloo too.
        let z, cz = statistics::barnard_exact(0.0, 3.0, 0.0, 4.0);
        check(bad, cz == statistics::ERR_UNDEFINED || (cz == 0 && z.pvalue == 1.0), "barnard empty column");
        let z2, cz2 = statistics::boschloo_exact(0.0, 3.0, 0.0, 4.0);
        check(bad, cz2 == statistics::ERR_UNDEFINED || (cz2 == 0 && z2.pvalue == 1.0), "boschloo empty column");
        // A table that is too large, a grid that is too large, a bad count.
        let t1, ct1 = statistics::barnard_exact(300.0, 200.0, 200.0, 300.0);
        fault_is(bad, "barnard table above the bound", t1, ct1, statistics::ERR_TOO_LONG);
        let t2, ct2 = statistics::boschloo_exact(7.0, 3.0, 2.0, 8.0, statistics::ALT_TWO_SIDED, 2000);
        fault_is(bad, "boschloo grid above the bound", t2, ct2, statistics::ERR_TOO_LONG);
        let t3, ct3 = statistics::barnard_exact(-1.0, 3.0, 2.0, 8.0);
        fault_is(bad, "barnard negative count", t3, ct3, statistics::ERR_DOMAIN);
        let t4, ct4 = statistics::barnard_exact(1.5, 3.0, 2.0, 8.0);
        fault_is(bad, "barnard fractional count", t4, ct4, statistics::ERR_DOMAIN);
        let t5, ct5 = statistics::barnard_exact(nan, 3.0, 2.0, 8.0);
        fault_is(bad, "barnard NaN count", t5, ct5, statistics::ERR_NAN);
        let t6, ct6 = statistics::boschloo_exact(inf, 3.0, 2.0, 8.0);
        fault_is(bad, "boschloo infinite count", t6, ct6, statistics::ERR_INFINITY);
        let t7, ct7 = statistics::barnard_exact(7.0, 3.0, 2.0, 8.0, 9);
        fault_is(bad, "barnard unknown alternative", t7, ct7, statistics::ERR_DOMAIN);
        let t8, ct8 = statistics::boschloo_exact(7.0, 3.0, 2.0, 8.0, statistics::ALT_TWO_SIDED, 0);
        fault_is(bad, "boschloo zero grid", t8, ct8, statistics::ERR_DOMAIN);
    }

    // The one-sample Cramer-von Mises test.
    {
        let r, c = statistics::cramervonmises(a[..], statistics::FAMILY_NORMAL, np[..]);
        ex_res(bad, "cramervonmises oracle", r, c, 0.02081411830801615, 0.9983318698471384, 1e-12, 1e-9);
        let one: f64[1] = [1.0];
        let t1, c1 = statistics::cramervonmises(one[..], statistics::FAMILY_NORMAL, np[..]);
        fault_is(bad, "cramervonmises one value", t1, c1, statistics::ERR_TOO_SHORT);
        let e: f64[0] = [];
        let t2, c2 = statistics::cramervonmises(e[..], statistics::FAMILY_NORMAL, np[..]);
        fault_is(bad, "cramervonmises empty", t2, c2, statistics::ERR_EMPTY);
        let bn: f64[3] = [1.0, nan, 2.0];
        let t3, c3b = statistics::cramervonmises(bn[..], statistics::FAMILY_NORMAL, np[..]);
        fault_is(bad, "cramervonmises NaN", t3, c3b, statistics::ERR_NAN);
        let bi: f64[3] = [1.0, inf, 2.0];
        let t4, c4 = statistics::cramervonmises(bi[..], statistics::FAMILY_NORMAL, np[..]);
        fault_is(bad, "cramervonmises infinity", t4, c4, statistics::ERR_INFINITY);
        let t5, c5 = statistics::cramervonmises(a[..], 99, np[..]);
        fault_is(bad, "cramervonmises unknown family", t5, c5, statistics::ERR_DOMAIN);
        let t6, c6 = statistics::cramervonmises(a[..], -1, np[..]);
        fault_is(bad, "cramervonmises negative family", t6, c6, statistics::ERR_DOMAIN);
        let one_p: f64[1] = [3.5];
        let t7, c7 = statistics::cramervonmises(a[..], statistics::FAMILY_NORMAL, one_p[..]);
        fault_is(bad, "cramervonmises wrong parameter count", t7, c7, statistics::ERR_LENGTH);
        let bad_p: f64[2] = [3.5, -1.0];
        let t8, c8 = statistics::cramervonmises(a[..], statistics::FAMILY_NORMAL, bad_p[..]);
        check(bad, c8 != 0, "cramervonmises negative scale");
        let nan_p: f64[2] = [nan, 1.0];
        let t9, c9 = statistics::cramervonmises(a[..], statistics::FAMILY_NORMAL, nan_p[..]);
        check(bad, c9 != 0, "cramervonmises NaN parameter");
        // A constant sample is a valid input for a one-sample test.
        let same_v: f64[5] = [2.0, 2.0, 2.0, 2.0, 2.0];
        let t10, c10 = statistics::cramervonmises(same_v[..], statistics::FAMILY_NORMAL, np[..]);
        check(bad, c10 == 0 && t10.pvalue >= 0.0 && t10.pvalue <= 1.0, "cramervonmises constant sample");
    }

    // The two-sample Cramer-von Mises test.
    {
        let r, c = statistics::cramervonmises_2samp(a[..], b[..]);
        ex_res(bad, "cvm2 auto", r, c, 0.546875, 0.032012432012432014, 1e-12, 1e-12);
        let re, ce = statistics::cramervonmises_2samp(a[..], b[..], statistics::MW_EXACT);
        ex_res(bad, "cvm2 exact", re, ce, 0.546875, 0.032012432012432014, 1e-12, 1e-12);
        let ra, ca = statistics::cramervonmises_2samp(a[..], b[..], statistics::MW_ASYMPTOTIC);
        ex_res(bad, "cvm2 asymptotic", ra, ca, 0.546875, 0.029770401611030284, 1e-12, 1e-9);
        // Tied values.
        let x: f64[5] = [1.0, 2.0, 2.0, 3.0, 4.0];
        let y: f64[5] = [2.0, 3.0, 3.0, 5.0, 6.0];
        let rt, ct = statistics::cramervonmises_2samp(x[..], y[..]);
        check(bad, ct == 0 && rt.pvalue > 0.0 && rt.pvalue <= 1.0 && rt.statistic >= 0.0, "cvm2 tied values");
        let t1, c1 = statistics::cramervonmises_2samp(a[..], b[..], 7);
        fault_is(bad, "cvm2 unknown method", t1, c1, statistics::ERR_DOMAIN);
        let one: f64[1] = [1.0];
        let t2, c2 = statistics::cramervonmises_2samp(one[..], b[..]);
        fault_is(bad, "cvm2 one value", t2, c2, statistics::ERR_TOO_SHORT);
        let e: f64[0] = [];
        let t3, c3b = statistics::cramervonmises_2samp(a[..], e[..]);
        fault_is(bad, "cvm2 empty", t3, c3b, statistics::ERR_EMPTY);
        let bn: f64[3] = [1.0, nan, 2.0];
        let t4, c4 = statistics::cramervonmises_2samp(a[..], bn[..]);
        fault_is(bad, "cvm2 NaN", t4, c4, statistics::ERR_NAN);
        let bi: f64[3] = [1.0, inf, 2.0];
        let t5, c5 = statistics::cramervonmises_2samp(bi[..], b[..]);
        fault_is(bad, "cvm2 infinity", t5, c5, statistics::ERR_INFINITY);
        // The exact form above 100 values refuses. The asymptotic form works.
        var big: f64[>..] = [];
        for i in 0..150 { big.push(i * 0.37); }
        let t6, c6 = statistics::cramervonmises_2samp(big[..], b[..], statistics::MW_EXACT);
        fault_is(bad, "cvm2 exact above the bound", t6, c6, statistics::ERR_TOO_LONG);
        let t7, c7 = statistics::cramervonmises_2samp(big[..], b[..], statistics::MW_AUTO);
        check(bad, c7 == 0 && t7.pvalue >= 0.0 && t7.pvalue <= 1.0, "cvm2 auto above 20 values");
        // Two equal constant samples give an undefined or a certain result.
        let k1: f64[4] = [1.0, 1.0, 1.0, 1.0];
        let t8, c8 = statistics::cramervonmises_2samp(k1[..], k1[..]);
        check(bad, (c8 == 0 && t8.pvalue == 1.0) || c8 == statistics::ERR_UNDEFINED, "cvm2 constant samples");
    }

    // The Epps-Singleton test.
    {
        let r, c = statistics::epps_singleton_2samp(a[..], b[..]);
        ex_res(bad, "epps oracle", r, c, 5.013687598904544, 0.28589594331337015, 1e-9, 1e-9);
        check(bad, c == 0 && r.df == 4.0 && is_nan(r.df2), "epps df is the rank");
        let k: f64[8] = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0];
        let t1, c1 = statistics::epps_singleton_2samp(k[..], k[..]);
        fault_is(bad, "epps constant samples", t1, c1, statistics::ERR_UNDEFINED);
        let sh: f64[4] = [1.0, 2.0, 3.0, 4.0];
        let t2, c2 = statistics::epps_singleton_2samp(sh[..], b[..]);
        fault_is(bad, "epps short sample", t2, c2, statistics::ERR_TOO_SHORT);
        let e: f64[0] = [];
        let t3, c3b = statistics::epps_singleton_2samp(e[..], b[..]);
        fault_is(bad, "epps empty", t3, c3b, statistics::ERR_EMPTY);
        let bn: f64[8] = [1.0, 2.0, nan, 4.0, 5.0, 6.0, 7.0, 8.0];
        let t4, c4 = statistics::epps_singleton_2samp(bn[..], b[..]);
        fault_is(bad, "epps NaN", t4, c4, statistics::ERR_NAN);
        let bi: f64[8] = [1.0, 2.0, inf, 4.0, 5.0, 6.0, 7.0, 8.0];
        let t5, c5 = statistics::epps_singleton_2samp(a[..], bi[..]);
        fault_is(bad, "epps infinity", t5, c5, statistics::ERR_INFINITY);
        let neg: f64[2] = [0.4, -0.8];
        let t6, c6 = statistics::epps_singleton_2samp(a[..], b[..], neg[..]);
        fault_is(bad, "epps negative point", t6, c6, statistics::ERR_DOMAIN);
        let zer: f64[1] = [0.0];
        let t7, c7 = statistics::epps_singleton_2samp(a[..], b[..], zer[..]);
        fault_is(bad, "epps zero point", t7, c7, statistics::ERR_DOMAIN);
        var many: f64[>..] = [];
        for i in 0..40 { many.push(0.1 + i * 0.05); }
        let t8, c8 = statistics::epps_singleton_2samp(a[..], b[..], many[..]);
        fault_is(bad, "epps too many points", t8, c8, statistics::ERR_TOO_LONG);
        // A given default list gives the default answer.
        let dflt: f64[2] = [0.4, 0.8];
        let r2, c2b = statistics::epps_singleton_2samp(a[..], b[..], dflt[..]);
        check(bad, c2b == 0 && r2.statistic == r.statistic && r2.pvalue == r.pvalue, "epps default points");
    }

    // The Baumgartner-Weiss-Schindler test.
    {
        var g = rng { s: 1 };
        let r, c = statistics::bws_test(a[..], b[..], g);
        ex_res(bad, "bws exact", r, c, 2.978111049107143, 0.03154623154623155, 1e-12, 1e-12);
        let rg, cg = statistics::bws_test(a[..], b[..], g, statistics::ALT_GREATER);
        ex_res(bad, "bws greater", rg, cg, -2.978111049107143, 0.9846930846930847, 1e-12, 1e-12);
        let rl, cl = statistics::bws_test(a[..], b[..], g, statistics::ALT_LESS);
        ex_res(bad, "bws less", rl, cl, -2.978111049107143, 0.015384615384615385, 1e-12, 1e-12);
        // The generator stays untouched in the exact form.
        check(bad, g.s == 1, "bws exact leaves the generator alone");
        let dx: f64[7] = [1.0, 2.0, 3.0, 4.0, 6.0, 7.0, 8.0];
        let dy: f64[7] = [5.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0];
        let rd, cd = statistics::bws_test(dx[..], dy[..], g);
        ex_res(bad, "bws documented row", rd, cd, 5.132167152575315, 0.002913752913752914, 1e-12, 1e-12);
        // Tied values use the midranks.
        let tx: f64[4] = [1.0, 2.0, 2.0, 3.0];
        let ty: f64[4] = [2.0, 3.0, 4.0, 4.0];
        let rt, ct = statistics::bws_test(tx[..], ty[..], g);
        check(bad, ct == 0 && rt.pvalue > 0.0 && rt.pvalue <= 1.0, "bws tied values");
        // The seeded draw above the exact bound: the same seed, the same answer.
        var big_x: f64[>..] = [];
        var big_y: f64[>..] = [];
        for i in 0..30 { big_x.push(i * 1.1); big_y.push(i * 1.1 + 5.0); }
        var g1 = rng { s: 5 };
        var g2 = rng { s: 5 };
        var g3 = rng { s: 6 };
        let m1, cm1 = statistics::bws_test(big_x[..], big_y[..], g1, statistics::ALT_TWO_SIDED, 2000);
        let m2, cm2 = statistics::bws_test(big_x[..], big_y[..], g2, statistics::ALT_TWO_SIDED, 2000);
        let m3, cm3 = statistics::bws_test(big_x[..], big_y[..], g3, statistics::ALT_TWO_SIDED, 2000);
        check(bad, cm1 == 0 && cm2 == 0 && cm3 == 0 && m1.pvalue == m2.pvalue && m1.statistic == m2.statistic, "bws one seed one answer");
        check(bad, m1.pvalue >= 1.0 / 2001.0 && m1.pvalue <= 1.0 && m3.pvalue >= 1.0 / 2001.0, "bws drawn p-value range");
        check(bad, g1.s != 5, "bws drawn form advances the generator");
        let sm: f64[1] = [1.0];
        let t1, c1 = statistics::bws_test(sm[..], b[..], g);
        fault_is(bad, "bws one value", t1, c1, statistics::ERR_TOO_SHORT);
        let e: f64[0] = [];
        let t2, c2 = statistics::bws_test(a[..], e[..], g);
        fault_is(bad, "bws empty", t2, c2, statistics::ERR_EMPTY);
        let bn: f64[3] = [1.0, nan, 2.0];
        let t3, c3b = statistics::bws_test(bn[..], b[..], g);
        fault_is(bad, "bws NaN", t3, c3b, statistics::ERR_NAN);
        let bi: f64[3] = [1.0, inf, 2.0];
        let t4, c4 = statistics::bws_test(a[..], bi[..], g);
        fault_is(bad, "bws infinity", t4, c4, statistics::ERR_INFINITY);
        let t5, c5 = statistics::bws_test(a[..], b[..], g, 9);
        fault_is(bad, "bws unknown alternative", t5, c5, statistics::ERR_DOMAIN);
        let t6, c6 = statistics::bws_test(big_x[..], big_y[..], g, statistics::ALT_TWO_SIDED, 0);
        fault_is(bad, "bws zero resamples", t6, c6, statistics::ERR_DOMAIN);
        let t7, c7 = statistics::bws_test(big_x[..], big_y[..], g, statistics::ALT_TWO_SIDED, statistics::MAX_RESAMPLES + 1);
        fault_is(bad, "bws resamples above the bound", t7, c7, statistics::ERR_TOO_LONG);
    }

    // The k-sample Anderson-Darling test.
    {
        var all: f64[>..] = [];
        for v in a { all.push(v); }
        for v in b { all.push(v); }
        for v in c3 { all.push(v); }
        let sz: i64[3] = [8, 8, 8];
        var out: f64[>..] = [];
        let r, c = statistics::anderson_ksamp(out, all[..], sz[..]);
        ex_res(bad, "anderson midrank", r, c, 5.518090763382168, 0.0010495357151611717, 1e-11, 1e-9);
        let want: f64[7] = [0.44925883860929594, 1.3052767, 1.9434184, 2.57696569, 3.41634856, 4.07210043, 5.56419101];
        check(bad, out.len == 7 && near_tol(out[0], want[0], 1e-6) && near_tol(out[3], want[3], 1e-6) && near_tol(out[6], want[6], 1e-6), "anderson critical values");
        let rr, cr = statistics::anderson_ksamp(out, all[..], sz[..], statistics::KSAMP_RIGHT);
        ex_res(bad, "anderson right", rr, cr, 5.432977778355046, 0.0011501345058158967, 1e-11, 1e-9);
        let rc, cc = statistics::anderson_ksamp(out, all[..], sz[..], statistics::KSAMP_CONTINUOUS);
        ex_res(bad, "anderson continuous", rc, cc, 5.432977778355046, 0.0011501345058158967, 1e-11, 1e-9);
        // A p-value above the table is capped at 0.25. One below is floored.
        let lo: f64[8] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let lo2: f64[8] = [1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5];
        var pair: f64[>..] = [];
        for v in lo { pair.push(v); }
        for v in lo2 { pair.push(v); }
        let sz2: i64[2] = [8, 8];
        let rp, cp = statistics::anderson_ksamp(out, pair[..], sz2[..]);
        check(bad, cp == 0 && rp.pvalue == 0.25, "anderson p-value capped");
        // Faults.
        let t1, c1 = statistics::anderson_ksamp(out, all[..], sz[..], 9);
        fault_is(bad, "anderson unknown variant", t1, c1, statistics::ERR_DOMAIN);
        let sz_bad: i64[3] = [8, 8, 7];
        let t2, c2 = statistics::anderson_ksamp(out, all[..], sz_bad[..]);
        fault_is(bad, "anderson sizes do not match", t2, c2, statistics::ERR_LENGTH);
        let sz_one: i64[1] = [24];
        let t3, c3b = statistics::anderson_ksamp(out, all[..], sz_one[..]);
        fault_is(bad, "anderson one sample", t3, c3b, statistics::ERR_TOO_SHORT);
        let few: f64[3] = [1.0, 2.0, 3.0];
        let sz_few: i64[2] = [2, 1];
        let t5, c5 = statistics::anderson_ksamp(out, few[..], sz_few[..]);
        fault_is(bad, "anderson fewer than four values", t5, c5, statistics::ERR_TOO_SHORT);
        let sz_zero: i64[3] = [8, 0, 16];
        let t4, c4 = statistics::anderson_ksamp(out, all[..], sz_zero[..]);
        check(bad, c4 == statistics::ERR_EMPTY || c4 == statistics::ERR_DOMAIN || c4 == statistics::ERR_TOO_SHORT, "anderson empty group");
        let e: f64[0] = [];
        let sz_e: i64[0] = [];
        let t5, c5 = statistics::anderson_ksamp(out, e[..], sz_e[..]);
        fault_is(bad, "anderson empty", t5, c5, statistics::ERR_EMPTY);
        var bn: f64[>..] = [];
        for v in all { bn.push(v); }
        bn[3] = nan;
        let t6, c6 = statistics::anderson_ksamp(out, bn[..], sz[..]);
        fault_is(bad, "anderson NaN", t6, c6, statistics::ERR_NAN);
        bn[3] = inf;
        let t7, c7 = statistics::anderson_ksamp(out, bn[..], sz[..]);
        fault_is(bad, "anderson infinity", t7, c7, statistics::ERR_INFINITY);
        let k: f64[6] = [2.0, 2.0, 2.0, 2.0, 2.0, 2.0];
        let sz3: i64[2] = [3, 3];
        let t8, c8 = statistics::anderson_ksamp(out, k[..], sz3[..]);
        fault_is(bad, "anderson constant pool", t8, c8, statistics::ERR_UNDEFINED);
    }

    // The Kolmogorov-Smirnov dispatchers.
    {
        var out: f64[>..] = [];
        let r, c = statistics::kstest_2samp(out, a[..], b[..]);
        ex_res(bad, "kstest_2samp", r, c, 0.625, 0.08702408702408701, 1e-12, 1e-12);
        check(bad, out.len == 2 && out[0] == 3.8 && out[1] == 1.0 && is_nan(r.df) && is_nan(r.df2), "kstest_2samp location and sign");
        let r1, c1 = statistics::kstest(out, a[..], statistics::FAMILY_NORMAL, np[..]);
        ex_res(bad, "kstest one sample", r1, c1, 0.1430611921955091, 0.9880479270584139, 1e-12, 1e-12);
        check(bad, out.len == 2 && out[0] == 1.9 && out[1] == -1.0, "kstest location and sign");
        let t1, ct1 = statistics::kstest(out, a[..], 99, np[..]);
        fault_is(bad, "kstest unknown family", t1, ct1, statistics::ERR_DOMAIN);
        check(bad, out.len == 0, "kstest fault clears out");
        let one_p: f64[1] = [3.5];
        let t2, ct2 = statistics::kstest(out, a[..], statistics::FAMILY_NORMAL, one_p[..]);
        fault_is(bad, "kstest wrong parameter count", t2, ct2, statistics::ERR_LENGTH);
        let e: f64[0] = [];
        let t3, ct3 = statistics::kstest(out, e[..], statistics::FAMILY_NORMAL, np[..]);
        fault_is(bad, "kstest empty", t3, ct3, statistics::ERR_EMPTY);
        let bn: f64[3] = [1.0, nan, 2.0];
        let t4, ct4 = statistics::kstest_2samp(out, a[..], bn[..]);
        fault_is(bad, "kstest_2samp NaN", t4, ct4, statistics::ERR_NAN);
        let t5, ct5 = statistics::kstest_2samp(out, e[..], b[..]);
        fault_is(bad, "kstest_2samp empty", t5, ct5, statistics::ERR_EMPTY);
        let t6, ct6 = statistics::kstest_2samp(out, a[..], b[..], 9);
        fault_is(bad, "kstest_2samp unknown alternative", t6, ct6, statistics::ERR_DOMAIN);
        // Tied values: the location is the first point of the largest gap.
        let tx: f64[4] = [1.0, 2.0, 2.0, 3.0];
        let ty: f64[4] = [2.0, 3.0, 4.0, 4.0];
        let rt, ctt = statistics::kstest_2samp(out, tx[..], ty[..]);
        let r_ref, c_ref = statistics::ks_2samp(tx[..], ty[..]);
        check(bad, ctt == 0 && c_ref == 0 && rt.statistic == r_ref.statistic && rt.pvalue == r_ref.pvalue && out.len == 2, "kstest_2samp agrees with ks_2samp on ties");
    }

    // The Poisson means test.
    {
        let r, c = statistics::poisson_means_test(10.0, 2.0, 18.0, 2.5);
        ex_res(bad, "poisson_means oracle", r, c, -0.9484873498302991, 0.355227490108269, 1e-12, 1e-9);
        let r2, c2 = statistics::poisson_means_test(0.0, 100.0, 3.0, 100.0);
        ex_res(bad, "poisson_means documented row", r2, c2, -1.7320508075688772, 0.08837900929018157, 1e-12, 1e-7);
        let r3, c3b = statistics::poisson_means_test(3.0, 1.0, 0.0, 2.0, 5.0);
        check(bad, c3b == 0 && r3.statistic == 0.0 && r3.pvalue == 1.0, "poisson_means a rate below the difference");
        let t1, c1 = statistics::poisson_means_test(-1.0, 2.0, 18.0, 2.5);
        fault_is(bad, "poisson_means negative count", t1, c1, statistics::ERR_DOMAIN);
        let t2, c2b = statistics::poisson_means_test(1.5, 2.0, 18.0, 2.5);
        fault_is(bad, "poisson_means fractional count", t2, c2b, statistics::ERR_DOMAIN);
        let t3, c3 = statistics::poisson_means_test(10.0, 0.0, 18.0, 2.5);
        fault_is(bad, "poisson_means zero size", t3, c3, statistics::ERR_DOMAIN);
        let t4, c4 = statistics::poisson_means_test(10.0, 2.0, 18.0, 2.5, -1.0);
        fault_is(bad, "poisson_means negative difference", t4, c4, statistics::ERR_DOMAIN);
        let t5, c5 = statistics::poisson_means_test(nan, 2.0, 18.0, 2.5);
        fault_is(bad, "poisson_means NaN", t5, c5, statistics::ERR_NAN);
        let t6, c6 = statistics::poisson_means_test(10.0, inf, 18.0, 2.5);
        fault_is(bad, "poisson_means infinite size", t6, c6, statistics::ERR_INFINITY);
        let t7, c7 = statistics::poisson_means_test(10.0, 2.0, 18.0, 2.5, 0.0, 9);
        fault_is(bad, "poisson_means unknown alternative", t7, c7, statistics::ERR_DOMAIN);
        // Two zero counts: the variance is 0, so the statistic is undefined.
        let t8, c8 = statistics::poisson_means_test(0.0, 1.0, 0.0, 1.0);
        check(bad, c8 == statistics::ERR_UNDEFINED || (c8 == 0 && t8.pvalue == 1.0), "poisson_means two zero counts");
        // A large mean makes the window too wide.
        let t9, c9 = statistics::poisson_means_test(1.0e6, 10.0, 1.0e6, 10.0);
        check(bad, c9 == statistics::ERR_TOO_LONG, "poisson_means window above the bound");
    }

    // The Page trend test.
    {
        var data: f64[>..] = [];
        for i in 0..8 {
            data.push(a[i]);
            data.push(b[i]);
            data.push(c3[i]);
        }
        let r, c = statistics::page_trend_test(data[..], 3);
        ex_res(bad, "page auto", r, c, 90.0, 0.9479184527892085, 1e-12, 1e-11);
        check(bad, c == 0 && r.df == statistics::PAGE_EXACT, "page auto picks the exact form");
        let ra, ca = statistics::page_trend_test(data[..], 3, false, [], statistics::PAGE_ASYMPTOTIC);
        ex_res(bad, "page asymptotic", ra, ca, 90.0, 0.9331927987311419, 1e-12, 1e-11);
        check(bad, ca == 0 && ra.df == statistics::PAGE_ASYMPTOTIC, "page asymptotic marks its method");
        // The documented table: ratings of three methods by ten students.
        let tab: f64[30] = [3.0, 4.0, 3.0, 2.0, 2.0, 4.0, 3.0, 3.0, 5.0, 1.0, 3.0, 2.0, 2.0, 3.0, 2.0,
                            2.0, 4.0, 5.0, 1.0, 2.0, 4.0, 3.0, 4.0, 4.0, 2.0, 4.0, 5.0, 1.0, 3.0, 4.0];
        let rd, cd = statistics::page_trend_test(tab[..], 3);
        ex_res(bad, "page documented row", rd, cd, 133.5, 0.0018191161948127822, 1e-12, 1e-11);
        let rda, cda = statistics::page_trend_test(tab[..], 3, false, [], statistics::PAGE_ASYMPTOTIC);
        ex_res(bad, "page documented asymptotic", rda, cda, 133.5, 0.0012693433690751756, 1e-12, 1e-11);
        // A column count that does not divide the data.
        let t1, c1 = statistics::page_trend_test(data[..], 5);
        fault_is(bad, "page columns do not divide", t1, c1, statistics::ERR_LENGTH);
        let t2, c2 = statistics::page_trend_test(data[..], 2);
        fault_is(bad, "page two columns", t2, c2, statistics::ERR_TOO_SHORT);
        let one_row: f64[3] = [1.0, 2.0, 3.0];
        let t3, c3b = statistics::page_trend_test(one_row[..], 3);
        fault_is(bad, "page one row", t3, c3b, statistics::ERR_TOO_SHORT);
        let e: f64[0] = [];
        let t4, c4 = statistics::page_trend_test(e[..], 3);
        fault_is(bad, "page empty", t4, c4, statistics::ERR_EMPTY);
        var bn: f64[>..] = [];
        for v in data { bn.push(v); }
        bn[4] = nan;
        let t5, c5 = statistics::page_trend_test(bn[..], 3);
        fault_is(bad, "page NaN", t5, c5, statistics::ERR_NAN);
        bn[4] = inf;
        let t6, c6 = statistics::page_trend_test(bn[..], 3);
        fault_is(bad, "page infinity", t6, c6, statistics::ERR_INFINITY);
        let t7, c7 = statistics::page_trend_test(data[..], 3, false, [], 9);
        fault_is(bad, "page unknown method", t7, c7, statistics::ERR_DOMAIN);
        let pbad: f64[3] = [1.0, 1.0, 3.0];
        let t8, c8 = statistics::page_trend_test(data[..], 3, false, pbad[..]);
        fault_is(bad, "page predicted ranks are no permutation", t8, c8, statistics::ERR_DOMAIN);
        let pshort: f64[2] = [1.0, 2.0];
        let t9, c9 = statistics::page_trend_test(data[..], 3, false, pshort[..]);
        fault_is(bad, "page predicted ranks length", t9, c9, statistics::ERR_LENGTH);
        let rk: f64[6] = [1.0, 2.0, 4.0, 1.0, 2.0, 3.0];
        let t10, c10 = statistics::page_trend_test(rk[..], 3, true);
        fault_is(bad, "page ranks above the column count", t10, c10, statistics::ERR_DOMAIN);
        // The exact form above 12 columns refuses. The asymptotic form works.
        var wide: f64[>..] = [];
        for i in 0..28 { wide.push((i * 7 % 13) * 1.0); }
        let t11, c11 = statistics::page_trend_test(wide[..], 14, false, [], statistics::PAGE_EXACT);
        fault_is(bad, "page exact above the column bound", t11, c11, statistics::ERR_TOO_LONG);
        let t12, c12 = statistics::page_trend_test(wide[..], 14);
        check(bad, c12 == 0 && t12.df == statistics::PAGE_ASYMPTOTIC, "page auto above 8 columns");
    }

    // No call writes into its input lists.
    {
        var xa: f64[>..] = [];
        var xb: f64[>..] = [];
        for v in a { xa.push(v); }
        for v in b { xb.push(v); }
        var g = rng { s: 3 };
        var out: f64[>..] = [];
        let sz: i64[2] = [8, 8];
        var pool: f64[>..] = [];
        for v in a { pool.push(v); }
        for v in b { pool.push(v); }
        let _1, _c1 = statistics::cramervonmises(xa[..], statistics::FAMILY_NORMAL, np[..]);
        let _2, _c2 = statistics::cramervonmises_2samp(xa[..], xb[..]);
        let _3, _c3 = statistics::epps_singleton_2samp(xa[..], xb[..]);
        let _4, _c4 = statistics::bws_test(xa[..], xb[..], g);
        let _5, _c5 = statistics::anderson_ksamp(out, pool[..], sz[..]);
        let _6, _c6 = statistics::kstest_2samp(out, xa[..], xb[..]);
        let _7, _c7 = statistics::kstest(out, xa[..], statistics::FAMILY_NORMAL, np[..]);
        let _8, _c8 = statistics::page_trend_test(pool[..], 4);
        var same_a = true;
        var same_b = true;
        for i in 0..8 {
            if xa[i] != a[i] { same_a = false; }
            if xb[i] != b[i] { same_b = false; }
        }
        for i in 0..8 {
            if pool[i] != a[i] || pool[8 + i] != b[i] { same_a = false; }
        }
        check(bad, same_a && same_b && sz[0] == 8 && sz[1] == 8 && np[0] == 3.5 && np[1] == 1.5, "the exact calls leave their inputs alone");
    }
    return bad;
}
