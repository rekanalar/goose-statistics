// Hostile input tests of the frequency calls. Each answer must be a value with
// code 0, or NaN with a known fault code.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_frequency_hostile() -> i64 {
    var bad: i64 = 0;
    var st: u64 = 246813585;
    let n = tests_count();
    var out: f64[>..] = [];
    var out2: f64[>..] = [];
    for k in 0..n {
        let hostile = rnd(st) % 3 != 0;
        let xs = hostile_sample(st, hostile);
        var ys = hostile_sample(st, hostile);
        if rnd(st) % 2 == 0 {
            ys.resize(0);
            for v in xs { ys.push(if hostile { hostile_value(st) } else { rnd(st) % 7 - 3 }); }
        }
        let per = (rnd(st) % 120) as f64;
        let alpha = (rnd(st) % 120) as f64 / 100.0;
        let level = (rnd(st) % 120) as f64 / 100.0;
        let method = rnd(st) % 12;
        let kind = rnd(st) % 8;
        let stat = rnd(st) % 10;
        let tail = rnd(st) % 4;
        let prop = (rnd(st) % 90) as f64 / 100.0;

        let v17, c17 = statistics::estimated_cdf(xs[..], alpha, method);
        hostile_report(bad, sound_value(v17, c17), "estimated_cdf", k, c17);
        let v18, c18 = statistics::tiecorrect(xs[..]);
        hostile_report(bad, sound_value(v18, c18), "tiecorrect", k, c18);
        let n19, c19 = statistics::ecdf(out, out2, xs[..]);
        hostile_report(bad, sound_count(n19, c19, out.len), "ecdf", k, c19);
        let n20, c20 = statistics::ecdf_interval(out, out2, xs[..], level);
        hostile_report(bad, sound_count(n20, c20, out.len), "ecdf_interval", k, c20);
        let n21, c21 = statistics::ecdf_censored(out, out2, xs[..], ys[..]);
        hostile_report(bad, sound_count(n21, c21, out.len), "ecdf_censored", k, c21);
        let n22, c22 = statistics::ecdf_censored_interval(out, out2, xs[..], ys[..], level);
        hostile_report(bad, sound_count(n22, c22, out.len), "ecdf_censored_interval", k, c22);
        let f23, c23 = statistics::cumfreq(out, xs[..], (rnd(st) % 6) as i64);
        hostile_report(bad, sound_freqs(f23, c23), "cumfreq", k, c23);
        let f24, c24 = statistics::relfreq(out, xs[..], (rnd(st) % 6) as i64);
        hostile_report(bad, sound_freqs(f24, c24), "relfreq", k, c24);
        let v40, c40 = statistics::ecdf_value(xs[..], alpha);
        hostile_report(bad, sound_value(v40, c40), "ecdf_value", k, c40);
        let f45, c45 = statistics::cumfreq_limits(out, xs[..], (rnd(st) % 6) as i64, alpha, level);
        hostile_report(bad, sound_freqs(f45, c45), "cumfreq_limits", k, c45);
        let f46, c46 = statistics::relfreq_limits(out, xs[..], (rnd(st) % 6) as i64, alpha, level);
        hostile_report(bad, sound_freqs(f46, c46), "relfreq_limits", k, c46);
        let ws = weights_of(st, xs.len);
        let f47, c47 = statistics::cumfreq_weighted(out, xs[..], ws[..], (rnd(st) % 6) as i64);
        hostile_report(bad, sound_freqs(f47, c47), "cumfreq_weighted", k, c47);
        let f48, c48 = statistics::relfreq_weighted(out, xs[..], ws[..], (rnd(st) % 6) as i64);
        hostile_report(bad, sound_freqs(f48, c48), "relfreq_weighted", k, c48);
    }
    return bad;
}
