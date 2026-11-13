// The fault edges of the entropy calls: the empty list, the NaN, the
// infinity, the negative probability and the input that stays unchanged.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_entropy_tests() -> i64 {
    var bad: i64 = 0;
    let empty: f64[0] = [];
    let half: f64[2] = [0.5, 0.5];
    let other: f64[2] = [0.9, 0.1];
    let neg: f64[2] = [-0.5, 1.5];
    let zeros: f64[2] = [0.0, 0.0];
    let zero_one: f64[2] = [0.0, 1.0];
    let nan_list: f64[2] = [0.5, 0.0 / 0.0];
    let inf_list: f64[2] = [0.5, 1.0 / 0.0];
    let one: f64[1] = [2.0];
    let const_sample: f64[4] = [2.0, 2.0, 2.0, 2.0];
    let ties: f64[9] = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
    let three_p: f64[3] = [0.5, 0.25, 0.25];

    // A share that underflows to 0 adds 0. A sample whose cube or fourth
    // power overflows gives a fault, never a silent NaN.
    let mixed: f64[3] = [1e300, 1e-300, 1.0];
    let wide: f64[4] = [1e200, -1e200, 0.0, 0.0];
    let e0, ce0 = statistics::entropy(mixed[..], 0.0);
    chk_val(bad, "entropy with a share that underflows", e0, ce0, 6.907755278982137e-298, 1e-9);
    let _k0, ck0 = statistics::kstat(wide[..], 4);
    chk_code(bad, "kstat fourth order overflow", ck0, statistics::ERR_OVERFLOW);

    // A window whose double wraps gives the domain fault.
    let ten: f64[10] = [0.1, 0.5, 0.9, 1.4, 2.0, 2.2, 3.1, 3.5, 4.2, 5.0];
    let _w0, cw0 = statistics::differential_entropy(ten[..], statistics::ENTROPY_VAN_ES, 4611686018427387904, 0.0);
    chk_code(bad, "differential_entropy window that wraps", cw0, statistics::ERR_DOMAIN);

    // entropy
    let _v1, c1 = statistics::entropy(empty[..], 0.0);
    chk_code(bad, "entropy empty", c1, statistics::ERR_EMPTY);
    let _v2, c2 = statistics::entropy(nan_list[..], 0.0);
    chk_code(bad, "entropy NaN", c2, statistics::ERR_NAN);
    let _v3, c3 = statistics::entropy(inf_list[..], 0.0);
    chk_code(bad, "entropy infinity", c3, statistics::ERR_INFINITY);
    let _v4, c4 = statistics::entropy(neg[..], 0.0);
    chk_code(bad, "entropy negative probability", c4, statistics::ERR_DOMAIN);
    let _v5, c5 = statistics::entropy(zeros[..], 0.0);
    chk_code(bad, "entropy zero sum", c5, statistics::ERR_UNDEFINED);
    let _v6, c6 = statistics::entropy(half[..], 1.0);
    chk_code(bad, "entropy base 1", c6, statistics::ERR_DOMAIN);
    let _v7, c7 = statistics::entropy(half[..], -2.0);
    chk_code(bad, "entropy negative base", c7, statistics::ERR_DOMAIN);

    // relative_entropy
    let _v8, c8 = statistics::relative_entropy(half[..], three_p[..], 0.0);
    chk_code(bad, "relative_entropy length", c8, statistics::ERR_LENGTH);
    let _v9, c9 = statistics::relative_entropy(half[..], zero_one[..], 0.0);
    chk_code(bad, "relative_entropy zero q", c9, statistics::ERR_UNDEFINED);
    let _v10, c10 = statistics::relative_entropy(neg[..], other[..], 0.0);
    chk_code(bad, "relative_entropy negative probability", c10, statistics::ERR_DOMAIN);

    // kstat and kstatvar
    let _v11, c11 = statistics::kstat(empty[..], 1);
    chk_code(bad, "kstat empty", c11, statistics::ERR_EMPTY);
    let _v12, c12 = statistics::kstat(one[..], 2);
    chk_code(bad, "kstat too few values", c12, statistics::ERR_TOO_SHORT);
    let _v13, c13 = statistics::kstat(one[..], 5);
    chk_code(bad, "kstat unknown order", c13, statistics::ERR_DOMAIN);
    let _v14, c14 = statistics::kstatvar(one[..], 3);
    chk_code(bad, "kstatvar unknown order", c14, statistics::ERR_DOMAIN);
    let _v15, c15 = statistics::kstatvar(one[..], 2);
    chk_code(bad, "kstatvar too few values", c15, statistics::ERR_TOO_SHORT);

    // differential_entropy
    let _v16, c16 = statistics::differential_entropy(empty[..], statistics::ENTROPY_AUTO, 0, 0.0);
    chk_code(bad, "differential_entropy empty", c16, statistics::ERR_EMPTY);
    let _v17, c17 = statistics::differential_entropy(one[..], statistics::ENTROPY_AUTO, 0, 0.0);
    chk_code(bad, "differential_entropy one value", c17, statistics::ERR_DOMAIN);
    let _v18, c18 = statistics::differential_entropy(ties[..], statistics::ENTROPY_VASICEK, 0, 0.0);
    chk_code(bad, "differential_entropy ties", c18, statistics::ERR_UNDEFINED);
    let _v18b, c18b = statistics::differential_entropy(const_sample[..], statistics::ENTROPY_VASICEK, 0, 0.0);
    chk_code(bad, "differential_entropy constant sample", c18b, statistics::ERR_DOMAIN);
    let _v19, c19 = statistics::differential_entropy(const_sample[..], 9, 0, 0.0);
    chk_code(bad, "differential_entropy unknown method", c19, statistics::ERR_DOMAIN);
    let _v20, c20 = statistics::differential_entropy(const_sample[..], statistics::ENTROPY_VASICEK, 2, 0.0);
    chk_code(bad, "differential_entropy window too wide", c20, statistics::ERR_DOMAIN);

    // The calls leave the list of the caller unchanged.
    let keep: f64[2] = [0.5, 0.5];
    let before: f64[2] = [0.5, 0.5];
    let _w1, _c21 = statistics::entropy(keep[..], 0.0);
    let _w2, _c22 = statistics::relative_entropy(keep[..], other[..], 0.0);
    var same = true;
    for i in 0..keep.len { if keep[i] != before[i] { same = false; } }
    check(bad, same, "the entropy calls leave the list unchanged");
    return bad;
}
