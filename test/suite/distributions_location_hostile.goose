// Hostile input tests of the Cauchy, Laplace, Rayleigh and Gumbel calls.
// Each answer is a finite value with code 0, or NaN with a known fault code.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_location_hostile() -> i64 {
    var bad: i64 = 0;
    var st: u64 = 90210137;
    let n = tests_count();
    for k in 0..n {
        let x = hostile_value(st);
        let loc = hostile_value(st);
        let scale = hostile_value(st);
        let p = (rnd(st) % 130 - 10) / 100.0;
        let v1, c1 = statistics::cauchy_pdf(x, loc, scale);
        hostile_report(bad, sound_finite(v1, c1), "cauchy_pdf", k, c1);
        let v2, c2 = statistics::cauchy_cdf(x, loc, scale);
        hostile_report(bad, sound_finite(v2, c2), "cauchy_cdf", k, c2);
        let v3, c3 = statistics::cauchy_sf(x, loc, scale);
        hostile_report(bad, sound_finite(v3, c3), "cauchy_sf", k, c3);
        let v4, c4 = statistics::cauchy_ppf(p, loc, scale);
        hostile_report(bad, sound_finite(v4, c4), "cauchy_ppf", k, c4);
        let v5, c5 = statistics::cauchy_mean(loc, scale);
        hostile_report(bad, sound_finite(v5, c5), "cauchy_mean", k, c5);
        let v6, c6 = statistics::laplace_pdf(x, loc, scale);
        hostile_report(bad, sound_finite(v6, c6), "laplace_pdf", k, c6);
        let v7, c7 = statistics::laplace_cdf(x, loc, scale);
        hostile_report(bad, sound_finite(v7, c7), "laplace_cdf", k, c7);
        let v8, c8 = statistics::laplace_sf(x, loc, scale);
        hostile_report(bad, sound_finite(v8, c8), "laplace_sf", k, c8);
        let v9, c9 = statistics::laplace_ppf(p, loc, scale);
        hostile_report(bad, sound_finite(v9, c9), "laplace_ppf", k, c9);
        let v10, c10 = statistics::laplace_var(loc, scale);
        hostile_report(bad, sound_finite(v10, c10), "laplace_var", k, c10);
        let v11, c11 = statistics::rayleigh_pdf(x, loc, scale);
        hostile_report(bad, sound_finite(v11, c11), "rayleigh_pdf", k, c11);
        let v12, c12 = statistics::rayleigh_cdf(x, loc, scale);
        hostile_report(bad, sound_finite(v12, c12), "rayleigh_cdf", k, c12);
        let v13, c13 = statistics::rayleigh_sf(x, loc, scale);
        hostile_report(bad, sound_finite(v13, c13), "rayleigh_sf", k, c13);
        let v14, c14 = statistics::rayleigh_ppf(p, loc, scale);
        hostile_report(bad, sound_finite(v14, c14), "rayleigh_ppf", k, c14);
        let v15, c15 = statistics::rayleigh_mean(loc, scale);
        hostile_report(bad, sound_finite(v15, c15), "rayleigh_mean", k, c15);
        let v16, c16 = statistics::rayleigh_var(loc, scale);
        hostile_report(bad, sound_finite(v16, c16), "rayleigh_var", k, c16);
        let v17, c17 = statistics::gumbel_r_pdf(x, loc, scale);
        hostile_report(bad, sound_finite(v17, c17), "gumbel_r_pdf", k, c17);
        let v18, c18 = statistics::gumbel_r_cdf(x, loc, scale);
        hostile_report(bad, sound_finite(v18, c18), "gumbel_r_cdf", k, c18);
        let v19, c19 = statistics::gumbel_r_sf(x, loc, scale);
        hostile_report(bad, sound_finite(v19, c19), "gumbel_r_sf", k, c19);
        let v20, c20 = statistics::gumbel_r_ppf(p, loc, scale);
        hostile_report(bad, sound_finite(v20, c20), "gumbel_r_ppf", k, c20);
        let v21, c21 = statistics::gumbel_r_mean(loc, scale);
        hostile_report(bad, sound_finite(v21, c21), "gumbel_r_mean", k, c21);
        let v22, c22 = statistics::gumbel_r_var(loc, scale);
        hostile_report(bad, sound_finite(v22, c22), "gumbel_r_var", k, c22);
        let v23, c23 = statistics::gumbel_l_pdf(x, loc, scale);
        hostile_report(bad, sound_finite(v23, c23), "gumbel_l_pdf", k, c23);
        let v24, c24 = statistics::gumbel_l_cdf(x, loc, scale);
        hostile_report(bad, sound_finite(v24, c24), "gumbel_l_cdf", k, c24);
        let v25, c25 = statistics::gumbel_l_sf(x, loc, scale);
        hostile_report(bad, sound_finite(v25, c25), "gumbel_l_sf", k, c25);
        let v26, c26 = statistics::gumbel_l_ppf(p, loc, scale);
        hostile_report(bad, sound_finite(v26, c26), "gumbel_l_ppf", k, c26);
        let v27, c27 = statistics::gumbel_l_mean(loc, scale);
        hostile_report(bad, sound_finite(v27, c27), "gumbel_l_mean", k, c27);
        let v28, c28 = statistics::gumbel_l_skew(loc, scale);
        hostile_report(bad, sound_finite(v28, c28), "gumbel_l_skew", k, c28);
    }
    return bad;
}
