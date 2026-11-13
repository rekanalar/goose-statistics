// The fault edges of the descriptive calls: the empty sample, the NaN,
// the infinity, the domain errors and the input that stays unchanged.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_descriptive_tests() -> i64 {
    var bad: i64 = 0;
    var out: f64[>..] = [];
    let empty: f64[0] = [];
    let one: f64[1] = [2.0];
    let nan_sample: f64[3] = [1.0, 0.0 / 0.0, 3.0];
    let inf_sample: f64[3] = [1.0, 1.0 / 0.0, 3.0];
    let const_sample: f64[4] = [2.0, 2.0, 2.0, 2.0];
    let zero: f64[2] = [0.0, 1.0];
    let neg: f64[2] = [-1.0, 2.0];

    // scoreatpercentile
    let _v1, c1 = statistics::scoreatpercentile(empty[..], 50.0);
    chk_code(bad, "scoreatpercentile empty", c1, statistics::ERR_EMPTY);
    let _v2, c2 = statistics::scoreatpercentile(one[..], 150.0);
    chk_code(bad, "scoreatpercentile above 100", c2, statistics::ERR_DOMAIN);
    let _v3, c3 = statistics::scoreatpercentile(one[..], -1.0);
    chk_code(bad, "scoreatpercentile below 0", c3, statistics::ERR_DOMAIN);
    let _v4, c4 = statistics::scoreatpercentile(nan_sample[..], 50.0);
    chk_code(bad, "scoreatpercentile NaN", c4, statistics::ERR_NAN);
    let _v5, c5 = statistics::scoreatpercentile(inf_sample[..], 50.0);
    chk_code(bad, "scoreatpercentile infinity", c5, statistics::ERR_INFINITY);
    let _v6, c6 = statistics::scoreatpercentile_method(one[..], 50.0, 9);
    chk_code(bad, "scoreatpercentile unknown method", c6, statistics::ERR_DOMAIN);

    // scoreatpercentile_limit and scoreatpercentile_limit_method
    let _v6b, c6b = statistics::scoreatpercentile_limit(empty[..], 50.0, 0.0, 10.0);
    chk_code(bad, "scoreatpercentile_limit empty", c6b, statistics::ERR_EMPTY);
    let _v6c, c6c = statistics::scoreatpercentile_limit(one[..], 50.0, 10.0, 0.0);
    chk_code(bad, "scoreatpercentile_limit bounds that do not ascend", c6c, statistics::ERR_DOMAIN);
    let _v6d, c6d = statistics::scoreatpercentile_limit(one[..], 50.0, 0.0 / 0.0, 10.0);
    chk_code(bad, "scoreatpercentile_limit NaN bound", c6d, statistics::ERR_NAN);
    let _v6e, c6e = statistics::scoreatpercentile_limit(one[..], 50.0, 0.0, 1.0 / 0.0);
    chk_code(bad, "scoreatpercentile_limit infinite bound", c6e, statistics::ERR_INFINITY);
    let _v6f, c6f = statistics::scoreatpercentile_limit(one[..], 150.0, 0.0, 10.0);
    chk_code(bad, "scoreatpercentile_limit above 100", c6f, statistics::ERR_DOMAIN);
    let _v6g, c6g = statistics::scoreatpercentile_limit(one[..], 50.0, 3.0, 4.0);
    chk_code(bad, "scoreatpercentile_limit no value in the bounds", c6g, statistics::ERR_UNDEFINED);
    let _v6h, c6h = statistics::scoreatpercentile_limit_method(one[..], 50.0, 0.0, 10.0, 9);
    chk_code(bad, "scoreatpercentile_limit_method unknown method", c6h, statistics::ERR_DOMAIN);
    let _v6i, c6i = statistics::scoreatpercentile_limit_method(one[..], 50.0, 10.0, 0.0, 0);
    chk_code(bad, "scoreatpercentile_limit_method bounds that do not ascend", c6i, statistics::ERR_DOMAIN);
    let _v6j, c6j = statistics::scoreatpercentile_limit_method(nan_sample[..], 50.0, 0.0, 10.0, 0);
    chk_code(bad, "scoreatpercentile_limit_method NaN", c6j, statistics::ERR_NAN);

    // percentileofscore
    let _v7, c7 = statistics::percentileofscore(empty[..], 1.0, statistics::PCT_RANK);
    chk_code(bad, "percentileofscore empty", c7, statistics::ERR_EMPTY);
    let _v8, c8 = statistics::percentileofscore(one[..], 0.0 / 0.0, statistics::PCT_RANK);
    chk_code(bad, "percentileofscore NaN score", c8, statistics::ERR_NAN);
    let _v9, c9 = statistics::percentileofscore(one[..], 1.0, 9);
    chk_code(bad, "percentileofscore unknown kind", c9, statistics::ERR_DOMAIN);

    // median_abs_deviation
    let _v10, c10 = statistics::median_abs_deviation(empty[..], 1.0, false);
    chk_code(bad, "median_abs_deviation empty", c10, statistics::ERR_EMPTY);
    let _v11, c11 = statistics::median_abs_deviation(one[..], 0.0, false);
    chk_code(bad, "median_abs_deviation zero scale", c11, statistics::ERR_DOMAIN);
    let _v12, c12 = statistics::median_abs_deviation(one[..], -1.0, false);
    chk_code(bad, "median_abs_deviation negative scale", c12, statistics::ERR_DOMAIN);
    let _v13, c13 = statistics::median_abs_deviation(nan_sample[..], 1.0, false);
    chk_code(bad, "median_abs_deviation NaN", c13, statistics::ERR_NAN);

    // A zero value or an overflow gives a fault, never a silent NaN.
    let zv: f64[3] = [0.0, 2.0, 3.0];
    let unit: f64[3] = [1.0, 1.0, 1.0];
    let zw: f64[3] = [0.0, 1.0, 1.0];
    let huge: f64[2] = [1e300, 1e300];
    let _p1, cp1 = statistics::pmean_weighted(zv[..], unit[..], -1.0);
    chk_code(bad, "pmean_weighted zero value, negative exponent", cp1, statistics::ERR_UNDEFINED);
    let _p2, cp2 = statistics::pmean_weighted(zv[..], unit[..], 0.0);
    chk_code(bad, "pmean_weighted zero value, exponent 0", cp2, statistics::ERR_UNDEFINED);
    let _p3, cp3 = statistics::pmean_weighted(huge[..], unit[0..2], 2.0);
    chk_code(bad, "pmean_weighted overflow", cp3, statistics::ERR_OVERFLOW);
    let p4, cp4 = statistics::pmean_weighted(zv[..], zw[..], -1.0);
    chk_val(bad, "pmean_weighted skips a value of weight 0", p4, cp4, 2.4, 1e-12);

    // A result that overflows gives a fault, never a silent infinity or NaN.
    let big2: f64[2] = [1e308, 1.7e308];
    let _e0, ce0 = statistics::expectile(big2[..], 0.5);
    check(bad, ce0 != 0 || !is_inf(_e0), "expectile midpoint does not overflow");
    let _e1, ce1 = statistics::pmean(one[..], 1e300);
    check(bad, ce1 != 0, "pmean with a huge exponent gives a fault");
    var lm_in: f64[>..] = [];
    for i in 0..100 { lm_in.push(1e300); }
    var lm_out: f64[>..] = [];
    let _l0, cl0 = statistics::lmoment(lm_out, lm_in[..], 16, false);
    check(bad, cl0 != 0 || (lm_out.len == 16 && !is_nan(lm_out[5])), "lmoment gives no silent NaN");

    // pmean
    let _v14, c14 = statistics::pmean(empty[..], 1.0);
    chk_code(bad, "pmean empty", c14, statistics::ERR_EMPTY);
    let _v15, c15 = statistics::pmean(neg[..], 1.0);
    chk_code(bad, "pmean negative value", c15, statistics::ERR_DOMAIN);
    let _v16, c16 = statistics::pmean(one[..], 1.0 / 0.0);
    chk_code(bad, "pmean infinite exponent", c16, statistics::ERR_DOMAIN);
    let _v17, c17 = statistics::pmean(one[..], 0.0 / 0.0);
    chk_code(bad, "pmean NaN exponent", c17, statistics::ERR_NAN);
    // A zero value is allowed, as in SciPy.
    let v18, c18 = statistics::pmean(zero[..], 1.0);
    chk_val(bad, "pmean zero value", v18, c18, 0.5, 1e-12);
    let _v19, c19 = statistics::pmean_weighted(one[..], const_sample[..], 1.0);
    chk_code(bad, "pmean_weighted length", c19, statistics::ERR_LENGTH);

    // gstd and gzscore
    let _v20, c20 = statistics::gstd(neg[..], 0);
    chk_code(bad, "gstd negative value", c20, statistics::ERR_DOMAIN);
    let _v21, c21 = statistics::gstd(zero[..], 0);
    chk_code(bad, "gstd zero value", c21, statistics::ERR_DOMAIN);
    let _v22, c22 = statistics::gstd(one[..], 1);
    chk_code(bad, "gstd one value ddof 1", c22, statistics::ERR_TOO_SHORT);
    let _v23, c23 = statistics::gzscore(out, const_sample[..], 0);
    chk_code(bad, "gzscore constant sample", c23, statistics::ERR_UNDEFINED);

    // zmap
    let _v24, c24 = statistics::zmap(out, empty[..], one[..], 0);
    chk_code(bad, "zmap empty scores", c24, statistics::ERR_EMPTY);
    let _v25, c25 = statistics::zmap(out, const_sample[..], const_sample[..], 0);
    chk_code(bad, "zmap constant reference", c25, statistics::ERR_UNDEFINED);
    let _v26, c26 = statistics::zmap(out, one[..], empty[..], 0);
    chk_code(bad, "zmap empty reference", c26, statistics::ERR_EMPTY);

    // describe
    let _v27, c27 = statistics::describe(one[..], 1, false);
    chk_code(bad, "describe one value ddof 1", c27, statistics::ERR_TOO_SHORT);
    let _v28, c28 = statistics::describe(const_sample[..], 0, false);
    chk_code(bad, "describe constant sample", c28, statistics::ERR_UNDEFINED);
    let _v29, c29 = statistics::describe(nan_sample[..], 0, false);
    chk_code(bad, "describe NaN", c29, statistics::ERR_NAN);

    // lmoment
    let _v30, c30 = statistics::lmoment(out, empty[..], 4, true);
    chk_code(bad, "lmoment empty", c30, statistics::ERR_EMPTY);
    let _v31, c31 = statistics::lmoment(out, one[..], 4, true);
    chk_code(bad, "lmoment order above n", c31, statistics::ERR_TOO_SHORT);
    let _v32, c32 = statistics::lmoment(out, one[..], 17, true);
    chk_code(bad, "lmoment order above the bound", c32, statistics::ERR_DOMAIN);
    let _v33, c33 = statistics::lmoment(out, const_sample[..], 4, true);
    chk_code(bad, "lmoment constant sample", c33, statistics::ERR_UNDEFINED);

    // expectile
    let _v34, c34 = statistics::expectile(empty[..], 0.5);
    chk_code(bad, "expectile empty", c34, statistics::ERR_EMPTY);
    let _v35, c35 = statistics::expectile(one[..], 1.5);
    chk_code(bad, "expectile alpha above 1", c35, statistics::ERR_DOMAIN);
    let _v36, c36 = statistics::expectile(one[..], -0.5);
    chk_code(bad, "expectile alpha below 0", c36, statistics::ERR_DOMAIN);
    let _v37, c37 = statistics::expectile_weighted(one[..], const_sample[..], 0.5);
    chk_code(bad, "expectile_weighted length", c37, statistics::ERR_LENGTH);

    // The calls leave the list of the caller unchanged.
    let keep: f64[4] = [2.5, -1.0, 4.0, 3.25];
    let before: f64[4] = [2.5, -1.0, 4.0, 3.25];
    let _w1, _c38 = statistics::describe(keep[..], 0, true);
    let _w2, _c39 = statistics::lmoment(out, keep[..], 4, true);
    let _w3, _c40 = statistics::scoreatpercentile(keep[..], 50.0);
    var same = true;
    for i in 0..keep.len { if keep[i] != before[i] { same = false; } }
    check(bad, same, "the descriptive calls leave the sample unchanged");
    return bad;
}
