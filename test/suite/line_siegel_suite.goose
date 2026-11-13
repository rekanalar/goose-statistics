// The laws of the Siegel estimate: the estimate follows a line that the
// points sit near, the two methods give a slope near each other, and the
// samples of the caller stay unchanged.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_line_siegel_tests() -> i64 {
    var bad: i64 = 0;
    let x: f64[8] = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
    let y: f64[8] = [1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0];
    let keep: f64[8] = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
    // a straight line gives its slope and its intercept
    let s, b, c = statistics::siegelslopes(y[..], x[..]);
    chk_code(bad, "siegelslopes code", c, 0);
    chk_val(bad, "the slope of a straight line", s, c, 2.0, 1e-12);
    chk_val(bad, "the intercept of a straight line", b, c, 1.0, 1e-12);
    // the separate method gives the same slope
    let s2, b2, c2 = statistics::siegelslopes(y[..], x[..], statistics::SIEGEL_SEPARATE);
    chk_code(bad, "siegelslopes separate code", c2, 0);
    chk_val(bad, "the separate slope", s2, c2, s, 1e-12);
    chk_val(bad, "the separate intercept", b2, c2, b, 1e-12);
    // an outlier moves the estimate less than the least squares line
    let out = [1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 500.0];
    let s3, b3, c3 = statistics::siegelslopes(out[..], x[..]);
    chk_code(bad, "siegelslopes with an outlier", c3, 0);
    check(bad, abs(s3 - 2.0) < 0.5, "an outlier moves the slope little");
    // the samples of the caller are unchanged
    var same = true;
    for i in 0..8 { if x[i] != keep[i] { same = false; } }
    check(bad, same, "the Siegel estimate leaves the samples unchanged");
    // a constant first sample has no slope
    let flat: f64[4] = [1.0, 1.0, 1.0, 1.0];
    let xs: f64[4] = [0.0, 1.0, 2.0, 3.0];
    let _s4, _b4, c4 = statistics::siegelslopes(xs[..], flat[..]);
    chk_code(bad, "a constant first sample", c4, 7);
    return bad;
}
