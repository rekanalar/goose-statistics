// The fault edges of the circular calls: the empty sample, the NaN, the
// infinity, the domain errors and the input that stays unchanged.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_circular_tests() -> i64 {
    var bad: i64 = 0;
    var dir: f64[>..] = [];
    let empty: f64[0] = [];
    let one: f64[1] = [2.0];
    let three: f64[3] = [1.0, 2.0, 3.0];
    let nan_sample: f64[3] = [1.0, 0.0 / 0.0, 3.0];
    let inf_sample: f64[3] = [1.0, 1.0 / 0.0, 3.0];
    // Two opposite angles, so the resultant is zero.
    let opposite: f64[2] = [0.0, 3.141592653589793];
    let zero_row: f64[2] = [0.0, 0.0];

    // circmean, circvar and circstd
    let _v1, c1 = statistics::circmean(empty[..], 0.0, 6.283185307179586);
    chk_code(bad, "circmean empty", c1, statistics::ERR_EMPTY);
    let _v2, c2 = statistics::circmean(three[..], 0.0, 0.0);
    chk_code(bad, "circmean empty range", c2, statistics::ERR_DOMAIN);
    let _v3, c3 = statistics::circmean(three[..], 6.283185307179586, 0.0);
    chk_code(bad, "circmean high below low", c3, statistics::ERR_DOMAIN);
    let _v4, c4 = statistics::circmean(nan_sample[..], 0.0, 6.283185307179586);
    chk_code(bad, "circmean NaN", c4, statistics::ERR_NAN);
    let _v5, c5 = statistics::circmean(inf_sample[..], 0.0, 6.283185307179586);
    chk_code(bad, "circmean infinity", c5, statistics::ERR_INFINITY);
    // The resultant of two opposite angles is not zero in floating point,
    // so the calls give a value, as SciPy does.
    let v6, c6 = statistics::circmean(opposite[..], 0.0, 6.283185307179586);
    chk_val(bad, "circmean opposite angles", v6, c6, 1.5707963267948966, 1e-9);
    let v7, c7 = statistics::circvar(opposite[..], 0.0, 6.283185307179586);
    chk_val(bad, "circvar opposite angles", v7, c7, 0.9999999999999999, 1e-9);
    let v8, c8 = statistics::circstd(opposite[..], 0.0, 6.283185307179586, false);
    chk_val(bad, "circstd opposite angles", v8, c8, 8.640816650440966, 1e-9);
    let v9, c9 = statistics::circstd(one[..], 0.0, 6.283185307179586, false);
    chk_val(bad, "circstd one value", v9, c9, 0.0, 1e-12);

    // directional_stats
    let _v10, c10 = statistics::directional_stats(dir, empty[..], 2, true);
    chk_code(bad, "directional_stats empty", c10, statistics::ERR_EMPTY);
    let _v11, c11 = statistics::directional_stats(dir, three[..], 2, true);
    chk_code(bad, "directional_stats odd rows", c11, statistics::ERR_LENGTH);
    let _v12, c12 = statistics::directional_stats(dir, three[..], 0, true);
    chk_code(bad, "directional_stats zero dimensions", c12, statistics::ERR_DOMAIN);
    let _v13, c13 = statistics::directional_stats(dir, zero_row[..], 2, true);
    chk_code(bad, "directional_stats zero row", c13, statistics::ERR_UNDEFINED);

    // The calls leave the list of the caller unchanged.
    let keep: f64[3] = [0.1, 0.5, 1.0];
    let before: f64[3] = [0.1, 0.5, 1.0];
    let _w1, _c14 = statistics::circmean(keep[..], 0.0, 6.283185307179586);
    let _w2, _c15 = statistics::circstd(keep[..], 0.0, 6.283185307179586, false);
    var same = true;
    for i in 0..keep.len { if keep[i] != before[i] { same = false; } }
    check(bad, same, "the circular calls leave the sample unchanged");

    // An angle that overflows when scaled gives a fault, never a silent NaN.
    let huge_angle: f64[1] = [1e308];
    let _h1, ch1 = statistics::circmean(huge_angle[..], 0.0, 1.0);
    chk_code(bad, "circmean angle that overflows", ch1, statistics::ERR_OVERFLOW);
    let _h2, ch2 = statistics::circvar(huge_angle[..], 0.0, 1.0);
    chk_code(bad, "circvar angle that overflows", ch2, statistics::ERR_OVERFLOW);
    return bad;
}
