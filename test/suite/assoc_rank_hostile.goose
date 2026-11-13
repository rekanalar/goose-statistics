// Hostile input tests of the rank association calls. Each answer is a
// statistic in -1 to 1 with a p-value in 0 to 1 and code 0, or NaN with a
// known fault code.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_assoc_rank_hostile() -> i64 {
    var bad: i64 = 0;
    var st: u64 = 313131;
    let n = tests_count();
    for k in 0..n {
        let hostile = rnd(st) % 3 != 0;
        let xs = hostile_sample(st, hostile);
        var ys = hostile_sample(st, hostile);
        if rnd(st) % 2 == 0 {
            ys.resize(0);
            for v in xs { ys.push(if hostile { hostile_value(st) } else { v }); }
        }
        let alt = rnd(st) % 4 - 1;
        let cont = rnd(st) % 2 == 0;
        let a1, c1 = statistics::spearmanrho(xs[..], ys[..], alt);
        hostile_report(bad, sound_corr(a1, c1), "spearmanrho", k, c1);
        let a2, c2 = statistics::pointbiserialr(xs[..], ys[..], alt);
        hostile_report(bad, sound_corr(a2, c2), "pointbiserialr", k, c2);
        let a3, c3 = statistics::chatterjeexi(xs[..], ys[..], cont);
        hostile_report(bad, sound_corr(a3, c3), "chatterjeexi", k, c3);
        let w1, c4 = statistics::weightedtau(xs[..], ys[..]);
        hostile_report(bad, sound_value(w1, c4), "weightedtau", k, c4);
        let w2, c5 = statistics::weightedtau(xs[..], ys[..], rnd(st) % 2 == 0);
        hostile_report(bad, sound_value(w2, c5), "weightedtau not additive", k, c5);
        let w3, c6 = statistics::weightedtau_ranked(xs[..], ys[..], xs[..]);
        hostile_report(bad, sound_value(w3, c6), "weightedtau ranked", k, c6);
    }
    // ranks that hold a NaN or an infinity give a fault and no index error
    let rk_x: f64[3] = [1.0, 2.0, 3.0];
    let rk_y: f64[3] = [1.0, 3.0, 2.0];
    let rk_nan: f64[3] = [0.0 / 0.0, 1.0, 2.0];
    let rk_w, rk_c = statistics::weightedtau_ranked(rk_x[..], rk_y[..], rk_nan[..]);
    chk_code(bad, "weightedtau ranked with a NaN rank", rk_c, 3);
    let rk_inf: f64[3] = [1.0 / 0.0, 1.0, 2.0];
    let rk_w2, rk_c2 = statistics::weightedtau_ranked(rk_x[..], rk_y[..], rk_inf[..]);
    chk_code(bad, "weightedtau ranked with an infinite rank", rk_c2, 4);
    return bad;
}
