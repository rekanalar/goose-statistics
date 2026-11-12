// Helpers of the test suites.
namespace suite;

import std;
import statistics from "../..";

// The next number of a seeded generator. The state wraps on purpose.
fn rnd(st: u64&) -> i64 {
    st = st * 6364136223846793005 + 1442695040888963407;
    return (st >> 33) as i64;
}

// True when a equals e within a relative 1e-9. A tiny e uses an absolute
// bound, so a p-value of 0 and a p-value of 1e-320 both pass.
fn near(a: f64, e: f64) -> bool {
    return near_tol(a, e, 1e-9);
}

// True when a equals e within the relative bound tol. A value of e near 0
// uses an absolute bound of 1e-300.
fn near_tol(a: f64, e: f64, tol: f64) -> bool {
    if a == e { return true; }
    return abs(a - e) <= tol * abs(e) + 1e-300;
}

// True when two lists hold the same values within the relative bound tol.
fn near_list(got: const f64[:], want: const f64[:], tol: f64) -> bool {
    if got.len != want.len { return false; }
    for i in 0..got.len {
        if !near_tol(got[i], want[i], tol) { return false; }
    }
    return true;
}

// The mean of a list, or NaN for an empty list. The resampling rows use it
// as their statistic.
fn mean_of(xs: const f64[:]) -> f64 {
    if xs.len == 0 { return 0.0 / 0.0; }
    var t: f64 = 0.0;
    for v in xs { t += v; }
    return t / xs.len;
}

// A faulted result has the code and NaN in every field.
fn fault_is(bad: i64&, name: u8[:], r: statistics::TestResult, c: i64, want: i64) {
    if c == want && is_nan(r.statistic) && is_nan(r.pvalue) && is_nan(r.df) && is_nan(r.df2) { return; }
    bad += 1;
    print("FAIL ", name, " code ", c, " want ", want, "\n");
}
