// The measured rows of the circular calls: the answers of SciPy 1.18.1.
// The generator test/make_shape_vectors.py writes this file.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_circular_vectors() -> i64 {
    var bad: i64 = 0;
    var out: f64[>..] = [];
    var out2: f64[>..] = [];
    let ang: f64[7] = [0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0];
    let v1, c1 = statistics::circmean(ang[..], 0.0, 6.283185307179586);
    chk_val(bad, "circmean", v1, c1, 0.8704198904898331, 1e-12);
    let v2, c2 = statistics::circmean(ang[..], -3.141592653589793, 3.141592653589793);
    chk_val(bad, "circmean low -pi", v2, c2, 0.870419890489833, 1e-12);
    let v3, c3 = statistics::circvar(ang[..], 0.0, 6.283185307179586);
    chk_val(bad, "circvar", v3, c3, 0.8588578576450885, 1e-12);
    let v4, c4 = statistics::circstd(ang[..], 0.0, 6.283185307179586, false);
    chk_val(bad, "circstd", v4, c4, 1.9788824091337363, 1e-12);
    let v5, c5 = statistics::circstd(ang[..], 0.0, 6.283185307179586, true);
    chk_val(bad, "circstd normalized", v5, c5, 1.9788824091337363, 1e-12);
    let v6, c6 = statistics::circmean(ang[..], 0.0, 10.0);
    chk_val(bad, "circmean range 10", v6, c6, 2.072678936445207, 1e-12);
    let v7, c7 = statistics::circstd(ang[..], 0.0, 10.0, false);
    chk_val(bad, "circstd range 10", v7, c7, 1.843997257204314, 1e-12);
    let dra: f64[4] = [3.0, 4.0, 6.0, -8.0];
    let _r8, c8 = statistics::directional_stats(out, dra[..], 2, true);
    let w8: f64[2] = [1.0, 0.0];
    chk_list(bad, "directional a direction", out[..], c8, w8[..], 1e-12);
    let v9, c9 = statistics::directional_stats(out, dra[..], 2, true);
    chk_val(bad, "directional a resultant", v9, c9, 0.6, 1e-12);
    let drb: f64[6] = [0.8660254, 0.5, 0.0, 0.8660254, -0.5, 0.0];
    let _r10, c10 = statistics::directional_stats(out, drb[..], 3, true);
    let w10: f64[3] = [1.0, 0.0, 0.0];
    chk_list(bad, "directional b direction", out[..], c10, w10[..], 1e-12);
    let v11, c11 = statistics::directional_stats(out, drb[..], 3, true);
    chk_val(bad, "directional b resultant", v11, c11, 0.8660254028383291, 1e-12);
    let drc: f64[9] = [1.0, 2.0, 2.0, 2.0, 0.5, 1.0, -1.0, 3.0, 0.5];
    let _r12, c12 = statistics::directional_stats(out, drc[..], 3, true);
    let w12: f64[3] = [0.3742585803780969, 0.7628419064037253, 0.5272596522091386];
    chk_list(bad, "directional c direction", out[..], c12, w12[..], 1e-12);
    let v13, c13 = statistics::directional_stats(out, drc[..], 3, true);
    chk_val(bad, "directional c resultant", v13, c13, 0.7961138967939102, 1e-12);
    return bad;
}
