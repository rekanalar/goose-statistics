// Hostile input for the eight laws of the second tail group. A NaN gives the NaN fault, an infinity
// gives the infinity fault and a probability outside 0 to 1 gives the
// domain fault.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_noncentral_hostile() -> i64 {
    var bad: i64 = 0;
    let n0, kn0 = statistics::crystalball_pdf(0.0 / 0.0, 1.5, 3.0);
    chk_code(bad, "crystalball a NaN point", kn0, 3);
    let zi0, kzi0 = statistics::crystalball_cdf(1.0 / 0.0, 1.5, 3.0);
    chk_code(bad, "crystalball an infinite point", kzi0, 4);
    let p0, kp0 = statistics::crystalball_ppf(1.5, 1.5, 3.0);
    chk_code(bad, "crystalball a probability above 1", kp0, 5);
    let n1, kn1 = statistics::kstwobign_pdf(0.0 / 0.0);
    chk_code(bad, "kstwobign a NaN point", kn1, 3);
    let zi1, kzi1 = statistics::kstwobign_cdf(1.0 / 0.0);
    chk_code(bad, "kstwobign an infinite point", kzi1, 4);
    let p1, kp1 = statistics::kstwobign_ppf(1.5);
    chk_code(bad, "kstwobign a probability above 1", kp1, 5);
    let n2, kn2 = statistics::vonmises_line_pdf(0.0 / 0.0, 2.0);
    chk_code(bad, "vonmises_line a NaN point", kn2, 3);
    let zi2, kzi2 = statistics::vonmises_line_cdf(1.0 / 0.0, 2.0);
    chk_code(bad, "vonmises_line an infinite point", kzi2, 4);
    let p2, kp2 = statistics::vonmises_line_ppf(1.5, 2.0);
    chk_code(bad, "vonmises_line a probability above 1", kp2, 5);
    let n3, kn3 = statistics::ncx2_pdf(0.0 / 0.0, 3.0, 2.0);
    chk_code(bad, "ncx2 a NaN point", kn3, 3);
    let zi3, kzi3 = statistics::ncx2_cdf(1.0 / 0.0, 3.0, 2.0);
    chk_code(bad, "ncx2 an infinite point", kzi3, 4);
    let p3, kp3 = statistics::ncx2_ppf(1.5, 3.0, 2.0);
    chk_code(bad, "ncx2 a probability above 1", kp3, 5);
    let n4, kn4 = statistics::ncf_pdf(0.0 / 0.0, 4.0, 7.0, 1.5);
    chk_code(bad, "ncf a NaN point", kn4, 3);
    let zi4, kzi4 = statistics::ncf_cdf(1.0 / 0.0, 4.0, 7.0, 1.5);
    chk_code(bad, "ncf an infinite point", kzi4, 4);
    let p4, kp4 = statistics::ncf_ppf(1.5, 4.0, 7.0, 1.5);
    chk_code(bad, "ncf a probability above 1", kp4, 5);
    let n5, kn5 = statistics::nct_pdf(0.0 / 0.0, 6.0, 1.0);
    chk_code(bad, "nct a NaN point", kn5, 3);
    let zi5, kzi5 = statistics::nct_cdf(1.0 / 0.0, 6.0, 1.0);
    chk_code(bad, "nct an infinite point", kzi5, 4);
    let p5, kp5 = statistics::nct_ppf(1.5, 6.0, 1.0);
    chk_code(bad, "nct a probability above 1", kp5, 5);
    let n6, kn6 = statistics::ksone_pdf(0.0 / 0.0, 10.0);
    chk_code(bad, "ksone a NaN point", kn6, 3);
    let zi6, kzi6 = statistics::ksone_cdf(1.0 / 0.0, 10.0);
    chk_code(bad, "ksone an infinite point", kzi6, 4);
    let p6, kp6 = statistics::ksone_ppf(1.5, 10.0);
    chk_code(bad, "ksone a probability above 1", kp6, 5);
    let n7, kn7 = statistics::kstwo_pdf(0.0 / 0.0, 10.0);
    chk_code(bad, "kstwo a NaN point", kn7, 3);
    let zi7, kzi7 = statistics::kstwo_cdf(1.0 / 0.0, 10.0);
    chk_code(bad, "kstwo an infinite point", kzi7, 4);
    let p7, kp7 = statistics::kstwo_ppf(1.5, 10.0);
    chk_code(bad, "kstwo a probability above 1", kp7, 5);
    // A random sweep: every answer with the code OK is finite.
    var st: u64 = 4711;
    for k in 0..tests_count() / 20 {
        let x = (rnd(st) % 801) as f64 / 100.0 - 4.0;
        let v0, c0 = statistics::crystalball_cdf(x, 1.5, 3.0);
        hostile_report(bad, c0 == 0 && v0 >= 0.0 && v0 <= 1.0, "crystalball cdf", k, c0);
        let v1, c1 = statistics::kstwobign_cdf(x);
        hostile_report(bad, c1 == 0 && v1 >= 0.0 && v1 <= 1.0, "kstwobign cdf", k, c1);
        let v2, c2 = statistics::vonmises_line_cdf(x, 2.0);
        hostile_report(bad, c2 == 0 && v2 >= 0.0 && v2 <= 1.0, "vonmises_line cdf", k, c2);
        let v3, c3 = statistics::ncx2_cdf(x, 3.0, 2.0);
        hostile_report(bad, c3 == 0 && v3 >= 0.0 && v3 <= 1.0, "ncx2 cdf", k, c3);
        let v4, c4 = statistics::ncf_cdf(x, 4.0, 7.0, 1.5);
        hostile_report(bad, c4 == 0 && v4 >= 0.0 && v4 <= 1.0, "ncf cdf", k, c4);
        let v5, c5 = statistics::nct_cdf(x, 6.0, 1.0);
        hostile_report(bad, c5 == 0 && v5 >= 0.0 && v5 <= 1.0, "nct cdf", k, c5);
    }
    return bad;
}
