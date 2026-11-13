// The laws of the power divergence: the powers 1, 0 and -1 are the
// Pearson, the log-likelihood and the modified log-likelihood statistics,
// and the call agrees with chisquare for the power 1.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_divergence_tests() -> i64 {
    var bad: i64 = 0;
    let obs: f64[6] = [16.0, 18.0, 16.0, 14.0, 12.0, 12.0];
    let exp: f64[6] = [16.0, 16.0, 16.0, 16.0, 16.0, 8.0];
    let none: f64[0] = [];
    let r1, c1 = statistics::power_divergence(obs[..], exp[..], statistics::LAMBDA_PEARSON);
    let r2, c2 = statistics::chisquare(obs[..], exp[..]);
    chk_val(bad, "the Pearson power equals chisquare", r1.statistic, c1 + c2, r2.statistic, 1e-13);
    chk_val(bad, "the same p-value", r1.pvalue, c1 + c2, r2.pvalue, 1e-13);
    chk_val(bad, "the degrees of freedom", r1.df, c1, 5.0, 0.0);
    // the limit of the power near 0 is the log-likelihood statistic
    let g0, cg0 = statistics::power_divergence(obs[..], exp[..], 0.0);
    let g1, cg1 = statistics::power_divergence(obs[..], exp[..], 1e-6);
    chk_val(bad, "the power near 0", g1.statistic, cg0 + cg1, g0.statistic, 1e-5);
    let m0, cm0 = statistics::power_divergence(obs[..], exp[..], -1.0);
    let m1, cm1 = statistics::power_divergence(obs[..], exp[..], -1.0 + 1e-6);
    chk_val(bad, "the power near -1", m1.statistic, cm0 + cm1, m0.statistic, 1e-5);
    // the wrapper with equal expected counts
    let u1, cu1 = statistics::power_divergence_uniform(obs[..], 0.0);
    let u2, cu2 = statistics::power_divergence(obs[..], none[..], 0.0);
    chk_val(bad, "uniform equals an empty expected list", u1.statistic, cu1 + cu2, u2.statistic, 0.0);
    // a bigger ddof gives fewer degrees of freedom and a smaller p-value
    let d1, cd1 = statistics::power_divergence(obs[..], exp[..], 1.0, 1);
    if !(d1.pvalue < r1.pvalue) { bad += 1; print("FAIL ddof lowers the p-value\n"); }
    chk_code(bad, "ddof code", cd1, 0);
    return bad;
}
