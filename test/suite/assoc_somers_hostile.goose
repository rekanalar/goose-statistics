// Hostile input tests of Somers D. Each answer is a statistic in the unit
// range and a p-value in the unit range with code 0, or a fault code.
namespace suite;

import std;
import statistics from "../..";

fn run_assoc_somers_hostile() -> i64 {
    var bad: i64 = 0;
    var st: u64 = 30977;
    let n = tests_count();
    for k in 0..n {
        let hostile = rnd(st) % 3 != 0;
        let xs = hostile_sample(st, hostile);
        let ys = hostile_sample(st, hostile);
        let alt = rnd(st) % 4;
        let r, c = statistics::somersd(xs[..], ys[..], alt);
        hostile_report(bad, sound_corr(r, c), "somersd", k, c);
        let rows = 1 + rnd(st) % 4;
        let cols = 1 + rnd(st) % 4;
        let cells = hostile_sample(st, hostile);
        let r2, c2 = statistics::somersd_table(cells[..], rows, cols, alt % 3);
        hostile_report(bad, sound_corr(r2, c2), "somersd_table", k, c2);
    }
    // a table of one row has no order and gives the undefined fault
    let one_row: f64[2] = [1.0, 2.0];
    let r1, c1 = statistics::somersd_table(one_row[..], 1, 2);
    chk_code(bad, "a table of one row", c1, 7);
    // a table of one column has no order and gives the undefined fault
    let r2, c2 = statistics::somersd_table(one_row[..], 2, 1);
    chk_code(bad, "a table of one column", c2, 7);
    // a table whose counts share no order gives the undefined fault
    let diag: f64[4] = [0.0, 1.0, 0.0, 1.0];
    let r3, c3 = statistics::somersd_table(diag[..], 2, 2);
    chk_code(bad, "a table with no order", c3, 7);
    // a table of the wrong shape gives the length fault
    let r4, c4 = statistics::somersd_table(one_row[..], 2, 2);
    chk_code(bad, "a table of the wrong shape", c4, 6);
    // a table of more cells than the bound gives the too long fault
    var wide_x: f64[>..] = [];
    var wide_y: f64[>..] = [];
    for k in 0..1025 { wide_x.push(k as f64); wide_y.push(k as f64); }
    let rw, cw = statistics::somersd(wide_x[..], wide_y[..]);
    chk_code(bad, "a table above the cell bound", cw, 9);
    return bad;
}
