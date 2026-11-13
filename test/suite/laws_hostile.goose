// Hostile input tests of the truncated law, the order statistic and the
// Dirichlet multinomial. Each answer is a finite value with code 0, or a
// fault code with no value.
namespace suite;

import std;
import statistics from "../..";

fn run_laws_hostile() -> i64 {
    var bad: i64 = 0;
    var st: u64 = 664201;
    let n = tests_count();
    var out: f64[>..] = [];
    for k in 0..n {
        let hostile = rnd(st) % 3 != 0;
        let xs = hostile_sample(st, hostile);
        let ys = hostile_sample(st, hostile);
        let fam = rnd(st) % 9;
        var a: f64 = 0.0;
        var b: f64 = 0.0;
        var lb: f64 = 0.0;
        var ub: f64 = 0.0;
        if xs.len > 0 { a = xs[0]; lb = xs[xs.len - 1]; }
        if ys.len > 0 { b = ys[0]; ub = ys[ys.len - 1]; }
        let p = (rnd(st) % 12) / 10.0 - 0.1;
        let v1, c1 = statistics::truncate_pdf(a, fam, b, ub, lb, ub);
        hostile_report(bad, sound_finite(v1, c1), "truncate_pdf", k, c1);
        let v2, c2 = statistics::truncate_cdf(a, fam, 0.0, 1.0, lb, ub);
        hostile_report(bad, sound_finite(v2, c2), "truncate_cdf", k, c2);
        let v3, c3 = statistics::truncate_ppf(p, fam, 0.0, 1.0, lb, ub);
        hostile_report(bad, sound_finite(v3, c3), "truncate_ppf", k, c3);
        let v4, c4 = statistics::truncate_mean(fam, a, b, lb, ub);
        hostile_report(bad, sound_finite(v4, c4), "truncate_mean", k, c4);
        let r = rnd(st) % 8;
        let sz = rnd(st) % 8;
        let v5, c5 = statistics::order_statistic_pdf(a, fam, 0.0, 1.0, r, sz);
        hostile_report(bad, sound_finite(v5, c5), "order_statistic_pdf", k, c5);
        let v6, c6 = statistics::order_statistic_ppf(p, fam, 0.0, 1.0, r, sz);
        hostile_report(bad, sound_finite(v6, c6), "order_statistic_ppf", k, c6);
        let v7, c7 = statistics::dirichlet_multinomial_logpmf(xs[..], ys[..], rnd(st) % 20);
        hostile_report(bad, sound_finite(v7, c7), "dirichlet_multinomial_logpmf", k, c7);
        let c8 = statistics::dirichlet_multinomial_cov(out, ys[..], rnd(st) % 20);
        hostile_report(bad, c8 == 0 || (c8 >= 1 && c8 <= 10 && out.len == 0), "dirichlet_multinomial_cov", k, c8);
    }
    // more categories than the bound give the too long fault
    var wide: f64[>..] = [];
    for i in 0..statistics::MAX_DIRMULT_K + 1 { wide.push(1.0); }
    let kw = statistics::dirichlet_multinomial_var(out, wide[..], 2);
    chk_code(bad, "dirichlet multinomial too many categories", kw, 9);
    // parameters whose sum overflows give the overflow fault
    let ovb: f64[2] = [1e200, 1e200];
    chk_code(bad, "dirichlet multinomial variance that overflows", statistics::dirichlet_multinomial_var(out, ovb[..], 10), 8);
    chk_code(bad, "dirichlet multinomial covariance that overflows", statistics::dirichlet_multinomial_cov(out, ovb[..], 10), 8);
    let ovinf: f64[2] = [1.7e308, 1.7e308];
    chk_code(bad, "dirichlet multinomial mean that overflows", statistics::dirichlet_multinomial_mean(out, ovinf[..], 10), 8);
    // the logistic density far below its location is 0 and not a fault
    let far, fc = statistics::truncate_pdf(-1000.0, 2, 0.0, 1.0, -2000.0, 0.0);
    chk_code(bad, "logistic density far below the location", fc, 0);
    check(bad, far == 0.0, "the logistic density far below the location is 0");
    return bad;
}
