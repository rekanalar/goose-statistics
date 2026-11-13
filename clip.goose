// The sigma clip of the package: an iterative removal of the values that
// lie far from the mean.
//
// The call returns the result first and a fault code last. A fault gives
// NaN in the bounds and 0 kept values.
namespace statistics;

import std;
import math;

// Clips the sample. A pass keeps the values in [mean - low * std,
// mean + high * std], with the mean and the deviation of the values still
// kept. The passes stop when no value drops. The kept values go into out.
// Each pass that continues drops at least one value, so the loop ends
// within the length of the sample.
fn sigmaclip(out: f64[>..]&, xs: const f64[:], low: f64 = 4.0, high: f64 = 4.0) -> ClipResult, i64 {
    out.resize(0);
    let fault = ClipResult { low: nan(), high: nan(), kept: 0 };
    let c = check_sample(xs, 1);
    if c != OK { return fault, c; }
    if is_nan(low) || is_nan(high) { return fault, ERR_NAN; }
    if is_inf(low) || is_inf(high) { return fault, ERR_INFINITY; }
    if low < 0.0 || high < 0.0 { return fault, ERR_DOMAIN; }
    var cur: f64[>..] = [];
    for x in xs { cur.push(x); }
    var lo: f64 = 0.0;
    var hi: f64 = 0.0;
    // Each pass that continues drops at least one value, so the passes
    // are at most the length of the sample.
    while true {
        if cur.len == 0 { return fault, ERR_UNDEFINED; }
        let m, cm = exact_mean(cur[..]);
        let sd, cs = stdev_core(cur[..], 0);
        if cm != OK { return fault, cm; }
        if cs != OK { return fault, cs; }
        lo = m - sd * low;
        hi = m + sd * high;
        if is_inf(lo) || is_inf(hi) { return fault, ERR_OVERFLOW; }
        var kept: i64 = 0;
        let len = cur.len;
        for k in 0..len {
            let x = cur[k];
            if x >= lo && x <= hi { cur[kept] = x; kept += 1; }
        }
        let dropped = len - kept;
        cur.resize(kept);
        if dropped == 0 { break; }
    }
    for x in cur { out.push(x); }
    return ClipResult { low: lo, high: hi, kept: cur.len }, OK;
}
