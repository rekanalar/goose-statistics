// Hostile input tests of the odds ratio and the relative risk. A good
// answer holds a finite estimate and an interval that starts at zero and
// ends above its start. A fault gives NaN in every field.
namespace suite;

import std;
import math;
import statistics from "../..";

fn risk_ok(r: statistics::IntervalResult, c: i64) -> bool {
    if c == 0 {
        if is_nan(r.statistic) || is_nan(r.low) || is_nan(r.high) { return false; }
        if r.statistic < 0.0 { return false; }
        if r.low < 0.0 { return false; }
        if r.low > r.high { return false; }
        return true;
    }
    return c >= 1 && c <= 10 && is_nan(r.statistic) && is_nan(r.low) && is_nan(r.high);
}

fn run_contingency_risk_hostile() -> i64 {
    var bad: i64 = 0;
    var st: u64 = 4271;
    let n = tests_count();
    for k in 0..n {
        var tab: f64[4] = [0.0, 0.0, 0.0, 0.0];
        for i in 0..4 { tab[i] = (rnd(st) % 21) as f64; }
        // the conditional form costs a search, so it runs once in fifty rounds
        var kind = statistics::ODDS_SAMPLE;
        if k % 50 == 0 { kind = statistics::ODDS_CONDITIONAL; }
        let alt = rnd(st) % 3;
        let level = 0.5 + 0.5 * (rnd(st) % 50) as f64 / 50.0;
        let r, c = statistics::odds_ratio(tab[..], kind, alt, level);
        hostile_report(bad, risk_ok(r, c), "odds_ratio", k, c);
        let ec = (rnd(st) % 21) as f64;
        let et = ec + (rnd(st) % 21) as f64;
        let cc = (rnd(st) % 21) as f64;
        let ct = cc + (rnd(st) % 21) as f64;
        let r2, c2 = statistics::relative_risk(ec, et, cc, ct, level);
        hostile_report(bad, risk_ok(r2, c2), "relative_risk", k, c2);
    }
    // a table of the wrong length
    let short: f64[3] = [1.0, 2.0, 3.0];
    let r3, c3 = statistics::odds_ratio(short[..]);
    chk_code(bad, "a table of three counts", c3, 6);
    // a count that is not a whole number
    let frac: f64[4] = [1.0, 2.0, 3.0, 4.5];
    let r4, c4 = statistics::odds_ratio(frac[..]);
    chk_code(bad, "a fractional count", c4, 5);
    // a negative count
    let neg: f64[4] = [-1.0, 2.0, 3.0, 4.0];
    let r5, c5 = statistics::odds_ratio(neg[..]);
    chk_code(bad, "a negative count", c5, 5);
    // a kind that the call does not hold
    let r6, c6 = statistics::odds_ratio(short[..], 0, 0, 0.95);
    chk_code(bad, "a table of three counts", c6, 6);
    let drug: f64[4] = [7.0, 15.0, 58.0, 472.0];
    let r7, c7 = statistics::odds_ratio(drug[..], 9);
    chk_code(bad, "an unknown kind", c7, 5);
    // a level outside the unit range
    let r8, c8 = statistics::odds_ratio(drug[..], 0, 0, 1.5);
    chk_code(bad, "a level above one", c8, 5);
    // a relative risk with more cases than people
    let r9, c9 = statistics::relative_risk(20.0, 10.0, 1.0, 10.0);
    chk_code(bad, "more cases than people", c9, 5);
    // a relative risk with a group of no people
    let r10, c10 = statistics::relative_risk(1.0, 0.0, 1.0, 10.0);
    chk_code(bad, "an empty group", c10, 5);
    // a conditional estimate over a support above the bound gives the too long fault
    let wide_t: f64[4] = [1.0e9, 1.0e9, 1.0e9, 1.0e9];
    let rw, cw = statistics::odds_ratio(wide_t[..], statistics::ODDS_CONDITIONAL);
    chk_code(bad, "conditional odds ratio over a wide support", cw, 9);
    // a control group with no case gives the overflow fault
    let rz, cz = statistics::relative_risk(5.0, 10.0, 0.0, 10.0);
    chk_code(bad, "relative risk with no control case", cz, 8);
    // both groups with no case give the undefined fault
    let rz2, cz2 = statistics::relative_risk(0.0, 10.0, 0.0, 10.0);
    chk_code(bad, "relative risk with no case at all", cz2, 7);
    return bad;
}
