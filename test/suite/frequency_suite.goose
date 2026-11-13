// The fault edges of the frequency calls: the empty sample, the NaN,
// the infinity, the domain errors and the input that stays unchanged.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_frequency_tests() -> i64 {
    var bad: i64 = 0;
    var out: f64[>..] = [];
    var out2: f64[>..] = [];
    let empty: f64[0] = [];
    let one: f64[1] = [2.0];
    let nan_sample: f64[3] = [1.0, 0.0 / 0.0, 3.0];
    let inf_sample: f64[3] = [1.0, 1.0 / 0.0, 3.0];
    let const_sample: f64[4] = [2.0, 2.0, 2.0, 2.0];
    let three: f64[3] = [1.0, 2.0, 3.0];

    // tiecorrect
    let _v1, c1 = statistics::tiecorrect(empty[..]);
    chk_code(bad, "tiecorrect empty", c1, statistics::ERR_EMPTY);
    let _v2, c2 = statistics::tiecorrect(nan_sample[..]);
    chk_code(bad, "tiecorrect NaN", c2, statistics::ERR_NAN);
    let _v3, c3 = statistics::tiecorrect(inf_sample[..]);
    chk_code(bad, "tiecorrect infinity", c3, statistics::ERR_INFINITY);
    let v4, c4 = statistics::tiecorrect(one[..]);
    chk_val(bad, "tiecorrect one value", v4, c4, 1.0, 1e-12);

    // ecdf and ecdf_value
    let _v5, c5 = statistics::ecdf(out, out2, empty[..]);
    chk_code(bad, "ecdf empty", c5, statistics::ERR_EMPTY);
    let _v6, c6 = statistics::ecdf(out, out2, nan_sample[..]);
    chk_code(bad, "ecdf NaN", c6, statistics::ERR_NAN);
    let _v7, c7 = statistics::ecdf_value(empty[..], 1.0);
    chk_code(bad, "ecdf_value empty", c7, statistics::ERR_EMPTY);
    let _v8, c8 = statistics::ecdf_value(one[..], 0.0 / 0.0);
    chk_code(bad, "ecdf_value NaN point", c8, statistics::ERR_NAN);
    let _v9, c9 = statistics::ecdf_value(one[..], 1.0 / 0.0);
    chk_code(bad, "ecdf_value infinite point", c9, statistics::ERR_INFINITY);

    // ecdf_censored and its band
    let _v10, c10 = statistics::ecdf_censored(out, out2, empty[..], three[..]);
    chk_code(bad, "ecdf_censored empty", c10, statistics::ERR_EMPTY);
    let _v11, c11 = statistics::ecdf_censored_interval(out, out2, empty[..], three[..], 0.95);
    chk_code(bad, "ecdf_censored_interval empty", c11, statistics::ERR_EMPTY);
    let _v12, c12 = statistics::ecdf_censored_interval(out, out2, three[..], three[..], 1.5);
    chk_code(bad, "ecdf_censored_interval level above 1", c12, statistics::ERR_DOMAIN);

    // ecdf_interval
    let _v13, c13 = statistics::ecdf_interval(out, out2, empty[..], 0.95);
    chk_code(bad, "ecdf_interval empty", c13, statistics::ERR_EMPTY);
    let _v14, c14 = statistics::ecdf_interval(out, out2, three[..], 0.0);
    chk_code(bad, "ecdf_interval level 0", c14, statistics::ERR_DOMAIN);
    let _v15, c15 = statistics::ecdf_interval(out, out2, three[..], 1.5);
    chk_code(bad, "ecdf_interval level above 1", c15, statistics::ERR_DOMAIN);

    // estimated_cdf
    let _v16, c16 = statistics::estimated_cdf(empty[..], 1.0, statistics::CDF_LINEAR);
    chk_code(bad, "estimated_cdf empty", c16, statistics::ERR_EMPTY);
    let _v17, c17 = statistics::estimated_cdf(one[..], 0.0 / 0.0, statistics::CDF_LINEAR);
    chk_code(bad, "estimated_cdf NaN point", c17, statistics::ERR_NAN);
    let _v18, c18 = statistics::estimated_cdf(one[..], 1.0, 99);
    chk_code(bad, "estimated_cdf unknown method", c18, statistics::ERR_DOMAIN);

    // cumfreq and relfreq
    let _v19, c19 = statistics::cumfreq(out, empty[..], 4);
    chk_code(bad, "cumfreq empty", c19, statistics::ERR_EMPTY);
    let _v20, c20 = statistics::cumfreq(out, const_sample[..], 0);
    chk_code(bad, "cumfreq zero bins", c20, statistics::ERR_DOMAIN);
    let _v21, c21 = statistics::cumfreq(out, const_sample[..], 1048577);
    chk_code(bad, "cumfreq bins above the bound", c21, statistics::ERR_TOO_LONG);
    let _v22, c22 = statistics::relfreq(out, empty[..], 4);
    chk_code(bad, "relfreq empty", c22, statistics::ERR_EMPTY);
    let _v23, c23 = statistics::relfreq_limits(out, three[..], 4, 5.0, 1.0);
    chk_code(bad, "relfreq_limits low above high", c23, statistics::ERR_DOMAIN);
    let _v24, c24 = statistics::cumfreq_weighted(out, three[..], one[..], 4);
    chk_code(bad, "cumfreq_weighted length", c24, statistics::ERR_LENGTH);
    let _v25, c25 = statistics::relfreq_weighted(out, three[..], one[..], 4);
    chk_code(bad, "relfreq_weighted length", c25, statistics::ERR_LENGTH);

    // The calls leave the list of the caller unchanged.
    let keep: f64[3] = [1.0, 2.0, 3.0];
    let before: f64[3] = [1.0, 2.0, 3.0];
    let _w1, _c26 = statistics::ecdf(out, out2, keep[..]);
    let _w2, _c27 = statistics::tiecorrect(keep[..]);
    let _w3, _c28 = statistics::estimated_cdf(keep[..], 2.0, statistics::CDF_LINEAR);
    var same = true;
    for i in 0..keep.len { if keep[i] != before[i] { same = false; } }
    check(bad, same, "the frequency calls leave the sample unchanged");

    // A range that overflows, or a width that underflows, gives a fault.
    let wide: f64[2] = [-1e308, 1e308];
    let _f1, cf1 = statistics::cumfreq(out, wide[..], 4);
    chk_code(bad, "cumfreq range that overflows", cf1, statistics::ERR_OVERFLOW);
    let tiny: f64[2] = [0.0, 4.9e-324];
    let _f2, cf2 = statistics::relfreq(out, tiny[..], 4);
    chk_code(bad, "relfreq width that underflows", cf2, statistics::ERR_DOMAIN);
    return bad;
}
