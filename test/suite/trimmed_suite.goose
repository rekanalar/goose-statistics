// The fault edges of the trimmed calls: the empty sample, the NaN, the
// infinity, the domain errors and the input that stays unchanged.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_trimmed_tests() -> i64 {
    var bad: i64 = 0;
    var out: f64[>..] = [];
    let empty: f64[0] = [];
    let one: f64[1] = [2.0];
    let three: f64[3] = [1.0, 2.0, 3.0];
    let nan_sample: f64[3] = [1.0, 0.0 / 0.0, 3.0];
    let inf_sample: f64[3] = [1.0, 1.0 / 0.0, 3.0];

    // tmin and tmax
    let _v1, c1 = statistics::tmin(empty[..], 0.0, true);
    chk_code(bad, "tmin empty", c1, statistics::ERR_EMPTY);
    let _v2, c2 = statistics::tmin(one[..], 0.0 / 0.0, true);
    chk_code(bad, "tmin NaN limit", c2, statistics::ERR_NAN);
    let _v3, c3 = statistics::tmin(one[..], 1.0 / 0.0, true);
    chk_code(bad, "tmin infinite limit", c3, statistics::ERR_INFINITY);
    let _v4, c4 = statistics::tmin(three[..], 10.0, true);
    chk_code(bad, "tmin no value above the limit", c4, statistics::ERR_UNDEFINED);
    let _v5, c5 = statistics::tmax(empty[..], 0.0, true);
    chk_code(bad, "tmax empty", c5, statistics::ERR_EMPTY);
    let _v6, c6 = statistics::tmax(one[..], 0.0 / 0.0, true);
    chk_code(bad, "tmax NaN limit", c6, statistics::ERR_NAN);
    let _v7, c7 = statistics::tmax(three[..], 0.0, true);
    chk_code(bad, "tmax no value below the limit", c7, statistics::ERR_UNDEFINED);

    // tvar, tstd and tsem
    let _v8, c8 = statistics::tvar(empty[..], 0.0, 0.0, false, false, true, true, 1);
    chk_code(bad, "tvar empty", c8, statistics::ERR_EMPTY);
    let _v9, c9 = statistics::tvar(one[..], 0.0, 0.0, false, false, true, true, 1);
    chk_code(bad, "tvar one value ddof 1", c9, statistics::ERR_TOO_SHORT);
    let _v10, c10 = statistics::tvar(three[..], 5.0, 6.0, true, true, true, true, 1);
    chk_code(bad, "tvar empty selection", c10, statistics::ERR_UNDEFINED);
    let _v11, c11 = statistics::tvar(three[..], 6.0, 5.0, true, true, true, true, 1);
    chk_code(bad, "tvar low above high", c11, statistics::ERR_DOMAIN);
    let _v12, c12 = statistics::tstd(three[..], 5.0, 6.0, true, true, true, true, 1);
    chk_code(bad, "tstd empty selection", c12, statistics::ERR_UNDEFINED);
    let _v13, c13 = statistics::tsem(three[..], 5.0, 6.0, true, true, true, true, 1);
    chk_code(bad, "tsem empty selection", c13, statistics::ERR_UNDEFINED);
    let _v14, c14 = statistics::tvar(nan_sample[..], 0.0, 0.0, false, false, true, true, 1);
    chk_code(bad, "tvar NaN", c14, statistics::ERR_NAN);
    let _v15, c15 = statistics::tstd(inf_sample[..], 0.0, 0.0, false, false, true, true, 1);
    chk_code(bad, "tstd infinity", c15, statistics::ERR_INFINITY);

    // trimboth and trim1
    let _v16, c16 = statistics::trimboth(out, empty[..], 0.1);
    chk_code(bad, "trimboth empty", c16, statistics::ERR_EMPTY);
    let _v17, c17 = statistics::trimboth(out, three[..], 0.5);
    chk_code(bad, "trimboth proportion of 0.5", c17, statistics::ERR_DOMAIN);
    let _v18, c18 = statistics::trimboth(out, three[..], -0.1);
    chk_code(bad, "trimboth negative proportion", c18, statistics::ERR_DOMAIN);
    let _v19, c19 = statistics::trimboth(out, three[..], 1.0);
    chk_code(bad, "trimboth proportion of 1", c19, statistics::ERR_DOMAIN);
    let _v20, c20 = statistics::trim1(out, three[..], 1.0, statistics::TAIL_RIGHT);
    chk_code(bad, "trim1 proportion of 1", c20, statistics::ERR_DOMAIN);
    let _v21, c21 = statistics::trim1(out, three[..], 0.25, 9);
    chk_code(bad, "trim1 unknown tail", c21, statistics::ERR_DOMAIN);
    // The left tail of a sample of three values keeps every value.
    let v22, c22 = statistics::trim1(out, three[..], 0.25, statistics::TAIL_LEFT);
    chk_code(bad, "trim1 left tail", c22, statistics::OK);
    let w22: f64[3] = [1.0, 2.0, 3.0];
    chk_list(bad, "trim1 left tail keeps every value", out[..], c22, w22[..], 1e-12);

    // The calls leave the list of the caller unchanged.
    let keep: f64[3] = [1.0, 2.0, 3.0];
    let before: f64[3] = [1.0, 2.0, 3.0];
    let _w1, _c23 = statistics::tvar(keep[..], 0.0, 0.0, false, false, true, true, 1);
    let _w2, _c24 = statistics::trimboth(out, keep[..], 0.1);
    let _w3, _c25 = statistics::trim1(out, keep[..], 0.1, statistics::TAIL_RIGHT);
    var same = true;
    for i in 0..keep.len { if keep[i] != before[i] { same = false; } }
    check(bad, same, "the trimmed calls leave the sample unchanged");
    return bad;
}
