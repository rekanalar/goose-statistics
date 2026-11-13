// The laws of Somers D: the statistic lies from -1 to 1, a perfect order
// gives 1, a reversed order gives -1, and the two alternatives split the
// two-sided p-value.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_assoc_somers_tests() -> i64 {
    var bad: i64 = 0;
    let x: f64[6] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    let up: f64[6] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    let down: f64[6] = [6.0, 5.0, 4.0, 3.0, 2.0, 1.0];
    let keep: f64[6] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    let a, ca = statistics::somersd(x[..], up[..]);
    chk_code(bad, "perfect order code", ca, 0);
    chk_val(bad, "perfect order", a.statistic, ca, 1.0, 1e-12);
    let b, cb = statistics::somersd(x[..], down[..]);
    chk_val(bad, "reversed order", b.statistic, cb, -1.0, 1e-12);
    // the alternatives split the two-sided p-value
    let t: f64[8] = [3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0];
    let u: f64[8] = [2.0, 7.0, 1.0, 8.0, 2.0, 8.0, 1.0, 8.0];
    let two, c0 = statistics::somersd(t[..], u[..], statistics::ALT_TWO_SIDED);
    let lo, c1 = statistics::somersd(t[..], u[..], statistics::ALT_LESS);
    let hi, c2 = statistics::somersd(t[..], u[..], statistics::ALT_GREATER);
    chk_val(bad, "tails add to one", lo.pvalue + hi.pvalue, c1, 1.0, 1e-12);
    chk_val(bad, "two-sided is twice the smaller tail", two.pvalue, c0, 2.0 * min(lo.pvalue, hi.pvalue), 1e-12);
    chk_val(bad, "the alternatives share the statistic", lo.statistic, c2, hi.statistic, 0.0);
    check(bad, abs(two.statistic) <= 1.0, "the statistic lies from -1 to 1");
    // the list form agrees with the table form
    let xs: f64[7] = [1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0];
    let ys: f64[7] = [1.0, 1.0, 2.0, 1.0, 2.0, 2.0, 2.0];
    let tab: f64[4] = [2.0, 1.0, 1.0, 3.0];
    let ta, cta = statistics::somersd_table(tab[..], 2, 2);
    let tb, ctb = statistics::somersd(xs[..], ys[..]);
    chk_code(bad, "table code", cta, 0);
    chk_code(bad, "list code", ctb, 0);
    chk_val(bad, "list and table agree", tb.statistic, ctb, ta.statistic, 1e-12);
    chk_val(bad, "list and table agree on p", tb.pvalue, ctb, ta.pvalue, 1e-12);
    // a table of two rows and three columns reads the cells in row order
    let wide: f64[6] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    let wa, cwa = statistics::somersd_table(wide[..], 2, 3);
    chk_code(bad, "wide table code", cwa, 0);
    let wa2, cwa2 = statistics::somersd_table(wide[..], 3, 2);
    chk_code(bad, "transposed table code", cwa2, 0);
    check(bad, wa.statistic != wa2.statistic, "the rows and the columns of the table differ");
    // D depends on the order of the samples: the denominator holds the
    // squares of the row sums, and a transpose changes them
    let ax: f64[7] = [1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0];
    let ay: f64[7] = [1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0];
    let fa, cfa = statistics::somersd(ax[..], ay[..]);
    let fb, cfb = statistics::somersd(ay[..], ax[..]);
    chk_code(bad, "first order code", cfa, 0);
    chk_code(bad, "second order code", cfb, 0);
    check(bad, fa.statistic != fb.statistic, "D depends on the order of the samples");
    chk_val(bad, "the first order value", fa.statistic, cfa, 0.6666666666666666, 1e-12);
    chk_val(bad, "the second order value", fb.statistic, cfb, 0.8, 1e-12);
    // the samples of the caller are unchanged
    var same = true;
    for i in 0..6 { if x[i] != keep[i] { same = false; } }
    check(bad, same, "somersd leaves the sample unchanged");
    return bad;
}
