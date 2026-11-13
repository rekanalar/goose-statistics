// The tests of the odds ratio and the relative risk. The sample odds
// ratio is the ratio of the products of the diagonals. The conditional
// odds ratio of a symmetric table is one. The interval of the sample form
// grows with the level.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_contingency_risk_tests() -> i64 {
    var bad: i64 = 0;
    let drug: f64[4] = [7.0, 15.0, 58.0, 472.0];
    // the sample odds ratio of a known table
    let s, cs = statistics::odds_ratio(drug[..], statistics::ODDS_SAMPLE);
    chk_code(bad, "sample code", cs, 0);
    chk_val(bad, "sample statistic", s.statistic, cs, 3.7977011494252872, 1e-12);
    chk_val(bad, "sample low", s.low, cs, 1.4869302660564159, 1e-12);
    chk_val(bad, "sample high", s.high, cs, 9.699536252360435, 1e-12);
    // the statistic is the ratio of the products of the diagonals
    chk_val(bad, "the diagonal product", s.statistic, cs, 7.0 * 472.0 / (15.0 * 58.0), 1e-12);
    // a larger level gives a wider interval
    let s9, c9 = statistics::odds_ratio(drug[..], statistics::ODDS_SAMPLE, statistics::ALT_TWO_SIDED, 0.99);
    chk_code(bad, "level code", c9, 0);
    check(bad, s9.low < s.low, "a larger level gives a lower bound");
    check(bad, s9.high > s.high, "a larger level gives a higher bound");
    // the one-sided intervals keep one bound at the edge
    let sl, cl = statistics::odds_ratio(drug[..], statistics::ODDS_SAMPLE, statistics::ALT_LESS);
    let sg, cg = statistics::odds_ratio(drug[..], statistics::ODDS_SAMPLE, statistics::ALT_GREATER);
    chk_code(bad, "less code", cl, 0);
    chk_code(bad, "greater code", cg, 0);
    chk_val(bad, "the lower bound of the less interval", sl.low, cl, 0.0, 1e-12);
    check(bad, is_inf(sg.high), "the upper bound of the greater interval is infinite");
    chk_val(bad, "the statistic of the less interval", sl.statistic, cl, s.statistic, 0.0);
    // the conditional odds ratio of a balanced table is one
    let bal: f64[4] = [5.0, 5.0, 5.0, 5.0];
    let cnd, cc = statistics::odds_ratio(bal[..], statistics::ODDS_CONDITIONAL);
    chk_code(bad, "conditional code", cc, 0);
    chk_val(bad, "conditional statistic", cnd.statistic, cc, 1.0, 1e-9);
    check(bad, cnd.low < 1.0 && cnd.high > 1.0, "the interval holds the value one");
    // the conditional odds ratio of a table with a zero cell is the ratio
    // of the conditional maximum likelihood
    let cc2, cc3 = statistics::odds_ratio(drug[..], statistics::ODDS_CONDITIONAL);
    chk_code(bad, "conditional drug code", cc3, 0);
    chk_val(bad, "conditional drug statistic", cc2.statistic, cc3, 3.7836687705553493, 1e-9);
    // the relative risk of two equal risks is one
    let rr, cr = statistics::relative_risk(2.0, 10.0, 4.0, 20.0);
    chk_code(bad, "relative risk code", cr, 0);
    chk_val(bad, "relative risk statistic", rr.statistic, cr, 1.0, 1e-12);
    check(bad, rr.low < 1.0 && rr.high > 1.0, "the interval of the relative risk holds the value one");
    // the relative risk is the ratio of the two risks
    let rr2, cr2 = statistics::relative_risk(27.0, 122.0, 44.0, 487.0);
    chk_code(bad, "relative risk drug code", cr2, 0);
    chk_val(bad, "relative risk drug statistic", rr2.statistic, cr2, 2.4495156482861398, 1e-12);
    chk_val(bad, "relative risk drug low", rr2.low, cr2, 1.5836990926700116, 1e-12);
    chk_val(bad, "relative risk drug high", rr2.high, cr2, 3.7886786315466354, 1e-12);
    // the sample of the caller is unchanged
    let keep: f64[4] = [7.0, 15.0, 58.0, 472.0];
    var same = true;
    for i in 0..4 { if drug[i] != keep[i] { same = false; } }
    check(bad, same, "the odds ratio leaves the table unchanged");
    return bad;
}
