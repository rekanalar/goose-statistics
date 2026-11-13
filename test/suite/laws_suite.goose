// The laws of the truncated law, the order statistic and the Dirichlet
// multinomial: the cdf is a monotone function from 0 to 1, the quantile
// inverts the cdf, and the probabilities of all outcomes sum to 1.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_laws_tests() -> i64 {
    var bad: i64 = 0;
    // the quantile inverts the cdf of a truncated law
    for k in 1..10 {
        let p = k / 10.0;
        let x, cx = statistics::truncate_ppf(p, 0, 0.0, 1.0, -1.0, 2.0);
        let back, cb = statistics::truncate_cdf(x, 0, 0.0, 1.0, -1.0, 2.0);
        chk_code(bad, "truncate round trip", cb, 0);
        chk_val(bad, "truncate round trip", back, cx, p, 1e-12);
    }
    // the cdf at the bounds is 0 and 1, and the density is 0 outside
    let c0, k0 = statistics::truncate_cdf(-1.0, 0, 0.0, 1.0, -1.0, 2.0);
    chk_val(bad, "truncate cdf at the low bound", c0, k0, 0.0, 0.0);
    let c1, k1 = statistics::truncate_cdf(2.0, 0, 0.0, 1.0, -1.0, 2.0);
    chk_val(bad, "truncate cdf at the high bound", c1, k1, 1.0, 0.0);
    let d0, kd = statistics::truncate_pdf(3.0, 0, 0.0, 1.0, -1.0, 2.0);
    chk_val(bad, "truncate density outside", d0, kd, 0.0, 0.0);
    // an unbounded truncation is the family itself
    let u, ku = statistics::truncate_pdf(0.3, 0, 0.0, 1.0);
    let w, kw = statistics::norm_pdf(0.3);
    chk_val(bad, "truncate without bounds", u, ku, w, 1e-15);
    // the mean of the symmetric truncation is the centre
    let m, km = statistics::truncate_mean(0, 0.0, 1.0, -2.0, 2.0);
    chk_val(bad, "truncate symmetric mean", m + 1.0, km, 1.0, 1e-12);
    // the variance of a truncation is smaller than the variance of the law
    let v, kv = statistics::truncate_var(0, 0.0, 1.0, -1.0, 1.0);
    check(bad, kv == 0 && v < 1.0 && v > 0.0, "truncation shrinks the variance");
    // the order statistic of rank 1 and n 1 is the family
    let o, ko = statistics::order_statistic_cdf(0.3, 0, 0.0, 1.0, 1, 1);
    let f, kf = statistics::norm_cdf(0.3);
    chk_val(bad, "order statistic of one draw", o, ko, f, 1e-14);
    // the means of the order statistics of a symmetric law mirror
    let lo, kl = statistics::order_statistic_mean(0, 0.0, 1.0, 1, 5);
    let hi, kh = statistics::order_statistic_mean(0, 0.0, 1.0, 5, 5);
    chk_val(bad, "order statistic mirror", lo + hi + 1.0, kl, 1.0, 1e-12);
    // the quantile of the order statistic inverts its cdf
    let q, kq = statistics::order_statistic_ppf(0.4, 2, 0.0, 1.0, 3, 7);
    let qb, kb = statistics::order_statistic_cdf(q, 2, 0.0, 1.0, 3, 7);
    chk_val(bad, "order statistic round trip", qb, kq, 0.4, 1e-12);
    // the probabilities of all outcomes of four draws sum to 1
    let al: f64[3] = [1.5, 2.0, 0.5];
    var total: f64 = 0.0;
    for a in 0..5 {
        for b in 0..5 - a {
            let x: f64[3] = [a, b, 4 - a - b];
            let p, kp = statistics::dirichlet_multinomial_pmf(x[..], al[..], 4);
            chk_code(bad, "dirichlet multinomial outcome", kp, 0);
            total += p;
        }
    }
    chk_val(bad, "dirichlet multinomial total", total, 0, 1.0, 1e-12);
    // the mean sums to n and the covariance rows sum to 0
    var out: f64[>..] = [];
    let km2 = statistics::dirichlet_multinomial_mean(out, al[..], 4);
    chk_val(bad, "dirichlet multinomial mean sum", out[0] + out[1] + out[2], km2, 4.0, 1e-12);
    let kc = statistics::dirichlet_multinomial_cov(out, al[..], 4);
    chk_code(bad, "dirichlet multinomial covariance", kc, 0);
    check(bad, out.len == 9, "the covariance holds k squared values");
    chk_val(bad, "dirichlet multinomial covariance row", out[0] + out[1] + out[2] + 1.0, kc, 1.0, 1e-12);
    // a fault leaves the output empty
    let kz = statistics::dirichlet_multinomial_var(out, al[..], -1);
    check(bad, kz != 0 && out.len == 0, "a fault leaves the output empty");
    return bad;
}
