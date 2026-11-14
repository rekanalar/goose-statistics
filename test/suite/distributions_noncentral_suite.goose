// The rows of the eight laws of the second tail group that the vector file does not hold: the round
// trip of the quantile and the bounds of the shapes.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_noncentral_tests() -> i64 {
    var bad: i64 = 0;
    let q0, kq0 = statistics::crystalball_ppf(0.35, 1.5, 3.0);
    let r0, kr0 = statistics::crystalball_cdf(q0, 1.5, 3.0);
    chk_val(bad, "crystalball cdf of the quantile", r0, kq0 + kr0, 0.35, 1e-7);
    let q1, kq1 = statistics::kstwobign_ppf(0.35);
    let r1, kr1 = statistics::kstwobign_cdf(q1);
    chk_val(bad, "kstwobign cdf of the quantile", r1, kq1 + kr1, 0.35, 1e-7);
    let q2, kq2 = statistics::vonmises_line_ppf(0.35, 2.0);
    let r2, kr2 = statistics::vonmises_line_cdf(q2, 2.0);
    chk_val(bad, "vonmises_line cdf of the quantile", r2, kq2 + kr2, 0.35, 1e-7);
    let q3, kq3 = statistics::ncx2_ppf(0.35, 3.0, 2.0);
    let r3, kr3 = statistics::ncx2_cdf(q3, 3.0, 2.0);
    chk_val(bad, "ncx2 cdf of the quantile", r3, kq3 + kr3, 0.35, 1e-7);
    let q4, kq4 = statistics::ncf_ppf(0.35, 4.0, 7.0, 1.5);
    let r4, kr4 = statistics::ncf_cdf(q4, 4.0, 7.0, 1.5);
    chk_val(bad, "ncf cdf of the quantile", r4, kq4 + kr4, 0.35, 1e-7);
    let q5, kq5 = statistics::nct_ppf(0.35, 6.0, 1.0);
    let r5, kr5 = statistics::nct_cdf(q5, 6.0, 1.0);
    chk_val(bad, "nct cdf of the quantile", r5, kq5 + kr5, 0.35, 1e-7);
    let q6, kq6 = statistics::ksone_ppf(0.35, 10.0);
    let r6, kr6 = statistics::ksone_cdf(q6, 10.0);
    chk_val(bad, "ksone cdf of the quantile", r6, kq6 + kr6, 0.35, 1e-7);
    let q7, kq7 = statistics::kstwo_ppf(0.35, 10.0);
    let r7, kr7 = statistics::kstwo_cdf(q7, 10.0);
    chk_val(bad, "kstwo cdf of the quantile", r7, kq7 + kr7, 0.35, 1e-7);
    // A shape outside its range gives the domain fault.
    let _d1, cd1 = statistics::crystalball_pdf(1.0, 1.5, 1.0);
    chk_code(bad, "crystalball m 1", cd1, 5);
    let _d2, cd2 = statistics::ksone_pdf(0.2, 0.0);
    chk_code(bad, "ksone n 0", cd2, 5);
    let _d3, cd3 = statistics::kstwo_cdf(0.2, 2000.0);
    chk_code(bad, "kstwo n above the bound", cd3, 5);
    let _d4, cd4 = statistics::ncx2_pdf(1.0, 3.0, -1.0);
    chk_code(bad, "ncx2 negative non-centrality", cd4, 5);
    let s1, k1 = statistics::vonmises_line_ppf(0.0, 2.0);
    chk_val(bad, "vonmises_line ppf 0", s1, k1, -3.141592653589793, 1e-12);
    return bad;
}
