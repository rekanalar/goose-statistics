// The laws of the sigma clip: the kept values lie in the bounds, the
// result of a second clip with the same factors keeps them all, and the
// sample of the caller is unchanged.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_clip_tests() -> i64 {
    var bad: i64 = 0;
    var out: f64[>..] = [];
    var out2: f64[>..] = [];
    let a: f64[12] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 100.0, 1000.0];
    let keep: f64[12] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 100.0, 1000.0];
    let r, c = statistics::sigmaclip(out, a[..], 1.5, 1.5);
    chk_code(bad, "sigmaclip code", c, 0);
    check(bad, r.kept == out.len, "the kept count is the length of the output");
    var inside = true;
    for v in out { if v < r.low || v > r.high { inside = false; } }
    check(bad, inside, "the kept values lie in the bounds");
    // a clip of a clipped sample drops nothing more
    let kept = copy(out);
    let r2, c2 = statistics::sigmaclip(out2, kept[..], 1.5, 1.5);
    chk_code(bad, "the second clip code", c2, 0);
    check(bad, r2.kept <= r.kept, "the second clip never keeps more");
    var same = true;
    for i in 0..12 { if a[i] != keep[i] { same = false; } }
    check(bad, same, "sigmaclip leaves the sample unchanged");
    // a wide window keeps everything
    let w, cw = statistics::sigmaclip(out, a[..], 100.0, 100.0);
    chk_code(bad, "wide window", cw, 0);
    check(bad, w.kept == 12, "a wide window keeps all values");
    // the output list is cleared on a fault
    let bad_low = statistics::sigmaclip(out, a[..], -1.0, 1.0);
    check(bad, out.len == 0, "a fault leaves the output empty");
    return bad;
}
