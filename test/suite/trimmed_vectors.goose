// The measured rows of the trimmed calls: the answers of SciPy 1.18.1.
// The generator test/make_shape_vectors.py writes this file.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_trimmed_vectors() -> i64 {
    var bad: i64 = 0;
    var out: f64[>..] = [];
    var out2: f64[>..] = [];
    let s2: f64[10] = [2.5, -1.0, 4.0, 3.25, 10.5, 0.5, 7.75, 3.0, 8.0, 6.0];
    let v1, c1 = statistics::tmin(s2[..], -1.0e300, true);
    chk_val(bad, "tmin", v1, c1, -1.0, 1e-12);
    let v2, c2 = statistics::tmin(s2[..], 3.0, true);
    chk_val(bad, "tmin lower 3", v2, c2, 3.0, 1e-12);
    let v3, c3 = statistics::tmin(s2[..], 3.0, false);
    chk_val(bad, "tmin lower 3 exclusive", v3, c3, 3.25, 1e-12);
    let v4, c4 = statistics::tmax(s2[..], 1.0e300, true);
    chk_val(bad, "tmax", v4, c4, 10.5, 1e-12);
    let v5, c5 = statistics::tmax(s2[..], 4.0, true);
    chk_val(bad, "tmax upper 4", v5, c5, 4.0, 1e-12);
    let v6, c6 = statistics::tmax(s2[..], 4.0, false);
    chk_val(bad, "tmax upper 4 exclusive", v6, c6, 3.25, 1e-12);
    let v7, c7 = statistics::tvar(s2[..], 0.0, 0.0, false, false, true, true, 1);
    chk_val(bad, "tvar", v7, c7, 12.816666666666668, 1e-12);
    let v8, c8 = statistics::tstd(s2[..], 0.0, 0.0, false, false, true, true, 1);
    chk_val(bad, "tstd", v8, c8, 3.58003724375413, 1e-12);
    let v9, c9 = statistics::tstd(s2[..], 2.0, 6.0, true, true, true, true, 1);
    chk_val(bad, "tstd limits", v9, c9, 1.3693063937629153, 1e-12);
    let v10, c10 = statistics::tsem(s2[..], 0.0, 0.0, false, false, true, true, 1);
    chk_val(bad, "tsem", v10, c10, 1.1321071798494464, 1e-12);
    let v11, c11 = statistics::tsem(s2[..], 2.0, 6.0, true, true, true, true, 1);
    chk_val(bad, "tsem limits", v11, c11, 0.6123724356957945, 1e-12);
    let _r12, c12 = statistics::trimboth(out, s2[..], 0.1);
    let w12: f64[8] = [0.5, 2.5, 3.0, 3.25, 4.0, 6.0, 7.75, 8.0];
    chk_list(bad, "trimboth 0.1", out[..], c12, w12[..], 1e-12);
    let _r13, c13 = statistics::trimboth(out, s2[..], 0.25);
    let w13: f64[6] = [2.5, 3.0, 3.25, 4.0, 6.0, 7.75];
    chk_list(bad, "trimboth 0.25", out[..], c13, w13[..], 1e-12);
    let _r14, c14 = statistics::trim1(out, s2[..], 0.25, 1);
    let w14: f64[8] = [-1.0, 0.5, 2.5, 3.0, 3.25, 4.0, 6.0, 7.75];
    chk_list(bad, "trim1 right", out[..], c14, w14[..], 1e-12);
    let _r15, c15 = statistics::trim1(out, s2[..], 0.6, 1);
    let w15: f64[4] = [-1.0, 0.5, 2.5, 3.0];
    chk_list(bad, "trim1 right large", out[..], c15, w15[..], 1e-12);
    let _r16, c16 = statistics::trim1(out, s2[..], 0.25, 0);
    let w16: f64[8] = [2.5, 3.0, 3.25, 4.0, 6.0, 7.75, 8.0, 10.5];
    chk_list(bad, "trim1 left", out[..], c16, w16[..], 1e-12);
    return bad;
}
