// The laws of the contingency calls: the margins add to the total, the
// expected counts keep the margins of the table, the measures of
// association lie in 0 to 1, and no call changes the table of the caller.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_contingency_tests() -> i64 {
    var bad: i64 = 0;
    var out: f64[>..] = [];
    var m0: f64[>..] = [];
    var m1: f64[>..] = [];
    let t: f64[6] = [15.0, 20.0, 10.0, 30.0, 25.0, 12.0];
    let keep: f64[6] = [15.0, 20.0, 10.0, 30.0, 25.0, 12.0];
    let sh: i64[2] = [2, 3];
    let c1 = statistics::margins(m0, t[..], sh[..], 0);
    let c2 = statistics::margins(m1, t[..], sh[..], 1);
    chk_code(bad, "margins rows", c1, 0);
    chk_code(bad, "margins columns", c2, 0);
    chk_val(bad, "margin rows add to the total", m0[0] + m0[1], 0, 112.0, 0.0);
    chk_val(bad, "margin columns add to the total", m1[0] + m1[1] + m1[2], 0, 112.0, 0.0);
    let ce = statistics::expected_freq(out, t[..], sh[..]);
    chk_code(bad, "expected_freq", ce, 0);
    let e: f64[6] = [out[0], out[1], out[2], out[3], out[4], out[5]];
    chk_val(bad, "expected row 0", e[0] + e[1] + e[2], 0, m0[0], 1e-12);
    chk_val(bad, "expected column 1", e[1] + e[4], 0, m1[1], 1e-12);
    // the three measures lie in 0 to 1 and a perfect table gives 1
    for m in 0..3 {
        let v, cv = statistics::association(t[..], 2, 3, m);
        chk_code(bad, "association code", cv, 0);
        if v < 0.0 || v > 1.0 { bad += 1; print("FAIL association range ", v, "\n"); }
    }
    let diag: f64[4] = [10.0, 0.0, 0.0, 10.0];
    let v1, cv1 = statistics::association(diag[..], 2, 2, statistics::ASSOC_CRAMER);
    chk_val(bad, "cramer of a perfect table", v1, cv1, 1.0, 1e-12);
    // the continuity correction changes a 2 by 2 table and no other
    let a, ca = statistics::association(t[..], 2, 3, 0, false);
    let b, cb = statistics::association(t[..], 2, 3, 0, true);
    chk_val(bad, "correction leaves a 2 by 3 table", a, ca + cb, b, 0.0);
    let sq: f64[4] = [10.0, 20.0, 30.0, 25.0];
    let a2, ca2 = statistics::association(sq[..], 2, 2, 0, false);
    let b2, cb2 = statistics::association(sq[..], 2, 2, 0, true);
    if !(b2 < a2) { bad += 1; print("FAIL correction lowers a 2 by 2 table\n"); }
    // the table of the caller is unchanged
    var same = true;
    for i in 0..6 { if t[i] != keep[i] { same = false; } }
    check(bad, same, "the contingency calls leave the table unchanged");
    // the check order: length bound, empty, NaN, infinity, length, domain
    let nanv: f64[4] = [1.0, 0.0 / 0.0, -1.0, 2.0];
    let _v3, c3 = statistics::association(nanv[..], 2, 2);
    chk_code(bad, "NaN before the negative count", c3, statistics::ERR_NAN);
    let shbig: i64[2] = [2000, 2];
    let _v4, c4 = statistics::association(t[..], 2000, 2);
    chk_code(bad, "too many rows", c4, statistics::ERR_TOO_LONG);
    let c5 = statistics::margins(out, t[..], shbig[..], 0);
    chk_code(bad, "margins axis too long", c5, statistics::ERR_TOO_LONG);
    return bad;
}
