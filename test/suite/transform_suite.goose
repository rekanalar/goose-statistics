// One test for each documented case of the transform calls. The measured
// rows live in transform_vectors.goose, which the generator
// test/make_transform_vectors.py writes.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_transform_tests() -> i64 {
    var bad: i64 = 0;
    let pos: f64[8] = [0.8, 1.5, 2.2, 3.1, 4.7, 6.3, 9.8, 12.5];
    let mixed: f64[8] = [-2.1, -0.4, 0.3, 1.1, 1.9, 2.6, 4.8, 7.5];
    var out: f64[>..] = [];

    // The oracle rows of SciPy for the two samples.
    let l1, c1 = statistics::boxcox_llf(0.3, pos[..]);
    check(bad, c1 == 0 && near_tol(l1, -9.301671184737357, 1e-12), "boxcox_llf, 0.3");
    let l2, c2 = statistics::boxcox_llf(0.0, pos[..]);
    check(bad, c2 == 0 && near_tol(l2, -9.298102800427483, 1e-12), "boxcox_llf, 0");
    let l3, c3 = statistics::yeojohnson_llf(0.4, mixed[..]);
    check(bad, c3 == 0 && near_tol(l3, -8.560033420240845, 1e-12), "yeojohnson_llf, mixed signs");
    let bv: f64[8] = [-0.21583850725792886, 0.4311564515228515, 0.8895156400367473, 1.3471128222404196,
                      1.969501527378815, 2.456697948733182, 3.2773532576754634, 3.7780116774108063];
    let r1, c4 = statistics::boxcox(out, pos[..], 0.3);
    check(bad, c4 == 0 && r1 == 0.3 && near_list(out[..], bv[..], 1e-12), "boxcox, lambda 0.3");
    let yv: f64[8] = [-2.9720754847428625, -0.4376682261785951, 0.280350850198276, 0.8982753492378878,
                      1.40587727318528, 1.7947331922020555, 2.816637831516918, 3.8309518948453003];
    let r2, c5 = statistics::yeojohnson(out, mixed[..], 0.5);
    check(bad, c5 == 0 && r2 == 0.5 && near_list(out[..], yv[..], 1e-12), "yeojohnson, lambda 0.5");
    let m1, d1 = statistics::boxcox_normmax(pos[..], statistics::NORMMAX_MLE);
    check(bad, d1 == 0 && near_tol(m1, 0.14670215286264066, 1e-6), "boxcox_normmax, mle");
    let m2, d2 = statistics::boxcox_normmax(pos[..]);
    check(bad, d2 == 0 && near_tol(m2, 0.165857035985584, 1e-6), "boxcox_normmax, the default is pearsonr");
    let m3, d3 = statistics::yeojohnson_normmax(mixed[..]);
    check(bad, d3 == 0 && near_tol(m3, 0.7205299104935587, 1e-6), "yeojohnson_normmax");
    let lo, hi, lm, d4 = statistics::boxcox_interval(pos[..], 0.05);
    check(bad, d4 == 0 && near_tol(lo, -0.7115331485527032, 1e-6) && near_tol(hi, 1.0680425539507548, 1e-6)
        && near_tol(lm, 0.14670215286264066, 1e-6), "boxcox_interval, alpha 0.05");

    // The lambda 0 and the lambdas near 0 agree: the transform is continuous.
    var z0: f64[>..] = [];
    var z1: f64[>..] = [];
    let _a, ca = statistics::boxcox(z0, pos[..], 0.0);
    let _b, cb = statistics::boxcox(z1, pos[..], 1e-12);
    var cont = ca == 0 && cb == 0 && near_list(z0[..], z1[..], 1e-9);
    for i in 0..8 { if !near_tol(z0[i], log(pos[i]), 1e-15) { cont = false; } }
    check(bad, cont, "boxcox, lambda 0 is the log");
    let _c, cc = statistics::boxcox(z1, pos[..], 1e-30);
    check(bad, cc == 0 && near_list(z0[..], z1[..], 1e-12), "boxcox, a tiny lambda");
    let ll0, cl0 = statistics::boxcox_llf(0.0, pos[..]);
    let ll1, cl1 = statistics::boxcox_llf(1e-9, pos[..]);
    check(bad, cl0 == 0 && cl1 == 0 && near_tol(ll0, ll1, 1e-6), "boxcox_llf, continuity at 0");
    let ly0, cy0 = statistics::yeojohnson_llf(0.0, pos[..]);
    let ly1, cy1 = statistics::yeojohnson_llf(1e-9, pos[..]);
    check(bad, cy0 == 0 && cy1 == 0 && near_tol(ly0, ly1, 1e-6), "yeojohnson_llf, continuity at 0");
    let _d, cd = statistics::yeojohnson(z0, mixed[..], 2.0);
    let _e, ce = statistics::yeojohnson(z1, mixed[..], 2.0 + 1e-12);
    check(bad, cd == 0 && ce == 0 && near_list(z0[..], z1[..], 1e-9), "yeojohnson, continuity at 2");

    // The fit modes use the lambda of the search.
    let f1, fc1 = statistics::boxcox_fit(out, pos[..]);
    check(bad, fc1 == 0 && near_tol(f1, 0.14670215286264066, 1e-6) && out.len == 8
        && near_tol(out[0], -0.21953070812164815, 1e-6), "boxcox_fit");
    let f2, fc2 = statistics::yeojohnson_fit(out, mixed[..]);
    check(bad, fc2 == 0 && near_tol(f2, 0.7205299104935587, 1e-6) && out.len == 8
        && near_tol(out[0], -2.5423608997559173, 1e-6), "yeojohnson_fit");

    // The normal plot tables.
    let n1, pc1 = statistics::boxcox_normplot(out, pos[..], -1.0, 1.0, 5);
    check(bad, pc1 == 0 && n1 == 5 && out.len == 10 && out[0] == -1.0 && out[8] == 1.0
        && near_tol(out[1], 0.9016204286380205, 1e-10), "boxcox_normplot, five rows");
    let n2, pc2 = statistics::yeojohnson_normplot(out, mixed[..], -1.0, 1.0, 5);
    check(bad, pc2 == 0 && n2 == 5 && near_tol(out[5], 0.8782878739874824, 1e-10), "yeojohnson_normplot, five rows");

    // A value at or below 0 in boxcox gives the domain fault.
    let nonpos: f64[3] = [1.0, 0.0, 2.0];
    let negv: f64[3] = [1.0, -1.0, 2.0];
    let _g1, g1c = statistics::boxcox(out, nonpos[..], 0.5);
    check(bad, g1c == 5 && out.len == 0, "boxcox, a zero");
    let _g2, g2c = statistics::boxcox(out, negv[..], 0.5);
    check(bad, g2c == 5, "boxcox, a negative value");
    let _g3, g3c = statistics::boxcox_llf(0.5, negv[..]);
    check(bad, g3c == 5, "boxcox_llf, a negative value");
    let _g4, g4c = statistics::boxcox_normmax(nonpos[..]);
    check(bad, g4c == 5, "boxcox_normmax, a zero");
    let _g5, g5c = statistics::boxcox_fit(out, negv[..]);
    check(bad, g5c == 5, "boxcox_fit, a negative value");
    let _g6, _g6b, _g6c, g6d = statistics::boxcox_interval(negv[..], 0.05);
    check(bad, g6d == 5, "boxcox_interval, a negative value");
    let _g7, g7c = statistics::boxcox_normplot(out, negv[..], -1.0, 1.0, 5);
    check(bad, g7c == 5, "boxcox_normplot, a negative value");
    let _g8, g8c = statistics::yeojohnson(out, negv[..], 0.5);
    check(bad, g8c == 0, "yeojohnson takes negative values");

    // A constant sample is undefined.
    let flat: f64[4] = [3.0, 3.0, 3.0, 3.0];
    let _h1, h1c = statistics::boxcox_llf(0.5, flat[..]);
    let _h2, h2c = statistics::yeojohnson_llf(0.5, flat[..]);
    let _h3, h3c = statistics::boxcox_normmax(flat[..]);
    let _h4, h4c = statistics::boxcox_normmax(flat[..], statistics::NORMMAX_MLE);
    let _h5, h5c = statistics::yeojohnson_normmax(flat[..]);
    let _h6, h6c = statistics::boxcox_fit(out, flat[..]);
    let _h7, h7c = statistics::yeojohnson_fit(out, flat[..]);
    let _h8, _h8b, _h8c, h8d = statistics::boxcox_interval(flat[..], 0.05);
    let _h9, h9c = statistics::boxcox_normplot(out, flat[..], -1.0, 1.0, 5);
    let _h10, h10c = statistics::yeojohnson_normplot(out, flat[..], -1.0, 1.0, 5);
    check(bad, h1c == 7 && h2c == 7 && h3c == 7 && h4c == 7 && h5c == 7 && h6c == 7 && h7c == 7
        && h8d == 7 && h9c == 7 && h10c == 7, "a constant sample is undefined");
    let zeros: f64[3] = [0.0, 0.0, 0.0];
    let _h11, h11c = statistics::yeojohnson_normmax(zeros[..]);
    check(bad, h11c == 7, "yeojohnson_normmax, zeros");

    // The input checks, in the fixed order.
    let none: f64[0] = [];
    let one: f64[1] = [2.0];
    let nan_x: f64[3] = [1.0, statistics::nan(), 2.0];
    let inf_x: f64[3] = [1.0, 1.0 / 0.0, 2.0];
    let neg_nan: f64[3] = [-1.0, statistics::nan(), 2.0];
    let _i1, i1c = statistics::boxcox_llf(0.5, none[..]);
    let _i2, i2c = statistics::boxcox_llf(0.5, nan_x[..]);
    let _i3, i3c = statistics::boxcox_llf(0.5, inf_x[..]);
    let _i4, i4c = statistics::boxcox_llf(0.5, one[..]);
    let _i5, i5c = statistics::boxcox_llf(0.5, neg_nan[..]);
    check(bad, i1c == 1 && i2c == 3 && i3c == 4 && i4c == 2 && i5c == 3, "boxcox_llf, the order of the checks");
    let _j1, j1c = statistics::yeojohnson_llf(0.5, none[..]);
    let _j2, j2c = statistics::yeojohnson_llf(0.5, nan_x[..]);
    let _j3, j3c = statistics::yeojohnson_llf(0.5, inf_x[..]);
    let _j4, j4c = statistics::yeojohnson_llf(0.5, one[..]);
    check(bad, j1c == 1 && j2c == 3 && j3c == 4 && j4c == 2, "yeojohnson_llf, the order of the checks");
    let _k1, k1c = statistics::boxcox(out, none[..], 0.5);
    let _k2, k2c = statistics::boxcox(out, nan_x[..], 0.5);
    let _k3, k3c = statistics::yeojohnson(out, inf_x[..], 0.5);
    let k4, k4c = statistics::boxcox(out, one[..], 0.5);
    check(bad, k1c == 1 && k2c == 3 && k3c == 4 && k4c == 0 && k4 == 0.5 && out.len == 1, "boxcox, the checks of a sample");
    let _k5, k5c = statistics::boxcox_normmax(none[..]);
    let _k6, k6c = statistics::yeojohnson_normmax(nan_x[..]);
    let _k7, k7c = statistics::yeojohnson_normmax(inf_x[..]);
    let _k8, k8c = statistics::yeojohnson_normmax(one[..]);
    check(bad, k5c == 1 && k6c == 3 && k7c == 4 && k8c == 2, "the searches check the sample");
    let nanl = statistics::nan();
    let infl = 1.0 / 0.0;
    let _m1, m1c = statistics::boxcox_llf(nanl, pos[..]);
    let _m2, m2c = statistics::boxcox(out, pos[..], infl);
    let _m3, m3c = statistics::yeojohnson(out, pos[..], nanl);
    let _m4, m4c = statistics::yeojohnson_llf(infl, pos[..]);
    check(bad, m1c == 3 && m2c == 4 && m3c == 3 && m4c == 4, "a NaN or infinite lambda");
    let _m5, m5c = statistics::boxcox_normmax(pos[..], 2);
    let _m6, m6c = statistics::boxcox_normmax(pos[..], -1);
    check(bad, m5c == 5 && m6c == 5, "boxcox_normmax, an unknown method");
    let _m7, _m7b, _m7c, m7d = statistics::boxcox_interval(pos[..], 0.0);
    let _m8, _m8b, _m8c, m8d = statistics::boxcox_interval(pos[..], 1.0);
    let _m9, _m9b, _m9c, m9d = statistics::boxcox_interval(pos[..], nanl);
    check(bad, m7d == 5 && m8d == 5 && m9d == 3, "boxcox_interval, the level");

    // The table arguments.
    let _p1, p1c = statistics::boxcox_normplot(out, pos[..], 1.0, -1.0, 5);
    let _p2, p2c = statistics::boxcox_normplot(out, pos[..], 1.0, 1.0, 5);
    let _p3, p3c = statistics::boxcox_normplot(out, pos[..], -1.0, 1.0, 1);
    let _p4, p4c = statistics::boxcox_normplot(out, pos[..], -1.0, 1.0, 0);
    let _p5, p5c = statistics::boxcox_normplot(out, pos[..], nanl, 1.0, 5);
    let _p6, p6c = statistics::boxcox_normplot(out, pos[..], -1.0, infl, 5);
    let _p7, p7c = statistics::yeojohnson_normplot(out, mixed[..], 1.0, -1.0, 5);
    let _p8, p8c = statistics::boxcox_normplot(out, pos[..], -1.0, 1.0, statistics::MAX_INPUT);
    let _p9, p9c = statistics::yeojohnson_normplot(out, mixed[..], -1.0, 1.0, statistics::MAX_INPUT);
    check(bad, p1c == 5 && p2c == 5 && p3c == 5 && p4c == 5 && p5c == 3 && p6c == 4 && p7c == 5
        && p8c == 9 && p9c == 9 && out.len == 0, "the plot tables, the arguments");
    let two, p10c = statistics::boxcox_normplot(out, pos[..], -1.0, 1.0, 2);
    check(bad, p10c == 0 && two == 2 && out.len == 4 && out[0] == -1.0 && out[2] == 1.0, "a plot table of two rows");

    // A sample above the input bound gives ERR_TOO_LONG before any work.
    var long: f64[>..] = [];
    for i in 0..statistics::MAX_INPUT + 1 { long.push(1.0 + (i % 7)); }
    let _q1, q1c = statistics::boxcox_llf(0.5, long[..]);
    let _q2, q2c = statistics::yeojohnson_llf(0.5, long[..]);
    let _q3, q3c = statistics::boxcox(out, long[..], 0.5);
    let _q4, q4c = statistics::yeojohnson(out, long[..], 0.5);
    let _q5, q5c = statistics::boxcox_normmax(long[..]);
    let _q6, q6c = statistics::yeojohnson_normmax(long[..]);
    let _q7, q7c = statistics::boxcox_fit(out, long[..]);
    let _q8, _q8b, _q8c, q8d = statistics::boxcox_interval(long[..], 0.05);
    let _q9, q9c = statistics::boxcox_normplot(out, long[..], -1.0, 1.0, 5);
    let _q10, q10c = statistics::yeojohnson_normplot(out, long[..], -1.0, 1.0, 5);
    check(bad, q1c == 9 && q2c == 9 && q3c == 9 && q4c == 9 && q5c == 9 && q6c == 9 && q7c == 9
        && q8d == 9 && q9c == 9 && q10c == 9, "the calls refuse a sample above the input bound");

    // The overflow cases. No result is NaN or infinite with the code OK.
    let wide: f64[4] = [1e-300, 1e300, 5.0, 7.0];
    let _r1, r1c = statistics::boxcox_normmax(wide[..], statistics::NORMMAX_MLE);
    let r2o, r2c = statistics::boxcox_llf(40.0, wide[..]);
    check(bad, (r1c == 0 || r1c == 8 || r1c == 10) && (r2c == 0 || r2c == 8) && (r2c != 0 || !is_nan(r2o)), "wide data, no silent NaN");
    let _r3, r3c = statistics::boxcox(out, wide[..], 5.0);
    check(bad, r3c == 8 && out.len == 0, "boxcox, an overflow");
    let _r4, r4c = statistics::yeojohnson(out, wide[..], 40.0);
    check(bad, r4c == 8, "yeojohnson, an overflow");
    let _r5, r5c = statistics::boxcox_normmax(wide[..]);
    check(bad, r5c == 0 || r5c == 8 || r5c == 10, "boxcox_normmax, wide data ends with a code");
    let _r6, r6c = statistics::yeojohnson_normmax(wide[..]);
    check(bad, r6c == 0 || r6c == 8 || r6c == 10, "yeojohnson_normmax, wide data ends with a code");
    var huge: f64[3] = [1e300, 2e300, -1e300];
    let _r7, r7c = statistics::yeojohnson_normmax(huge[..]);
    check(bad, r7c == 0 || r7c == 8 || r7c == 10, "yeojohnson_normmax, huge values");
    let _r8, r8c = statistics::yeojohnson(out, huge[..], 1.0);
    check(bad, r8c == 0 && is_finite_list(out[..]), "yeojohnson, huge values at lambda 1");

    // No call writes into the sample of the caller.
    let keep: f64[8] = [0.8, 1.5, 2.2, 3.1, 4.7, 6.3, 9.8, 12.5];
    let keep0: f64[8] = [0.8, 1.5, 2.2, 3.1, 4.7, 6.3, 9.8, 12.5];
    let _s1, _sc1 = statistics::boxcox_llf(0.3, keep[..]);
    let _s2, _sc2 = statistics::boxcox(out, keep[..], 0.3);
    let _s3, _sc3 = statistics::boxcox_normmax(keep[..]);
    let _s4, _sc4 = statistics::boxcox_normmax(keep[..], statistics::NORMMAX_MLE);
    let _s5, _sc5 = statistics::boxcox_fit(out, keep[..]);
    let _s6, _s6b, _s6c, _sc6 = statistics::boxcox_interval(keep[..], 0.05);
    let _s7, _sc7 = statistics::boxcox_normplot(out, keep[..], -1.0, 1.0, 9);
    let _s8, _sc8 = statistics::yeojohnson_llf(0.3, keep[..]);
    let _s9, _sc9 = statistics::yeojohnson(out, keep[..], 0.3);
    let _s10, _sc10 = statistics::yeojohnson_normmax(keep[..]);
    let _s11, _sc11 = statistics::yeojohnson_fit(out, keep[..]);
    let _s12, _sc12 = statistics::yeojohnson_normplot(out, keep[..], -1.0, 1.0, 9);
    check(bad, near_list(keep[..], keep0[..], 0.0), "the calls leave the sample unchanged");

    // The order of the values does not change the table or the lambda.
    let rev: f64[8] = [12.5, 9.8, 6.3, 4.7, 3.1, 2.2, 1.5, 0.8];
    let v1, vc1 = statistics::boxcox_normmax(pos[..]);
    let v2, vc2 = statistics::boxcox_normmax(rev[..]);
    check(bad, vc1 == 0 && vc2 == 0 && v1 == v2, "boxcox_normmax, the order of the values");

    // A lambda near zero keeps the series form, not the plain log.
    let big: f64[3] = [1e10, 2e10, 3e10];
    var tb: f64[>..] = [];
    let lt = 1e-11;
    let _lt, ltc = statistics::boxcox(tb, big[..], lt);
    let g0 = log(1e10);
    let w0 = g0 * (1.0 + lt * g0 / 2.0 + lt * lt * g0 * g0 / 6.0);
    check(bad, ltc == 0 && tb.len == 3 && near_tol(tb[0], w0, 1e-13), "boxcox, a lambda near zero");

    // The root finder needs a sign change, even when the product is small.
    let _rt, rtc = statistics::brent_root(0.5, 1.0, 1e-12, 1e-12) { x => x };
    check(bad, rtc == statistics::ERR_DOMAIN, "brent_root, no sign change");
    return bad;
}

fn is_finite_list(xs: const f64[:]) -> bool {
    for x in xs {
        if is_nan(x) || abs(x) > 1.7976931348623157e308 { return false; }
    }
    return true;
}
