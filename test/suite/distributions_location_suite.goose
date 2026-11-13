// The laws of the Cauchy, Laplace, Rayleigh and Gumbel families: the
// cumulative function and the survival function add to 1, the quantile
// inverts the cumulative function, and a default equals the explicit form.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_location_tests() -> i64 {
    var bad: i64 = 0;
    let xs: f64[7] = [-4.0, -1.5, -0.1, 0.3, 1.0, 2.2, 6.0];
    for x in xs {
        let a1, ca1 = statistics::cauchy_cdf(x, 0.5, 2.0);
        let b1, cb1 = statistics::cauchy_sf(x, 0.5, 2.0);
        chk_val(bad, "cauchy cdf plus sf", a1 + b1, ca1 + cb1, 1.0, 1e-14);
        let q1, cq1 = statistics::cauchy_ppf(a1, 0.5, 2.0);
        chk_val(bad, "cauchy ppf of cdf", q1, cq1, x, 1e-9);
        let a2, ca2 = statistics::laplace_cdf(x, 0.5, 2.0);
        let b2, cb2 = statistics::laplace_sf(x, 0.5, 2.0);
        chk_val(bad, "laplace cdf plus sf", a2 + b2, ca2 + cb2, 1.0, 1e-14);
        let q2, cq2 = statistics::laplace_ppf(a2, 0.5, 2.0);
        chk_val(bad, "laplace ppf of cdf", q2, cq2, x, 1e-9);
        let a3, ca3 = statistics::gumbel_r_cdf(x, 0.5, 2.0);
        let b3, cb3 = statistics::gumbel_r_sf(x, 0.5, 2.0);
        chk_val(bad, "gumbel_r cdf plus sf", a3 + b3, ca3 + cb3, 1.0, 1e-14);
        let q3, cq3 = statistics::gumbel_r_ppf(a3, 0.5, 2.0);
        chk_val(bad, "gumbel_r ppf of cdf", q3, cq3, x, 1e-9);
        let a4, ca4 = statistics::gumbel_l_cdf(x, 0.5, 2.0);
        let b4, cb4 = statistics::gumbel_l_sf(x, 0.5, 2.0);
        chk_val(bad, "gumbel_l cdf plus sf", a4 + b4, ca4 + cb4, 1.0, 1e-14);
        let q4, cq4 = statistics::gumbel_l_ppf(a4, 0.5, 2.0);
        chk_val(bad, "gumbel_l ppf of cdf", q4, cq4, x, 1e-9);
        let a5, ca5 = statistics::rayleigh_cdf(x, -5.0, 2.0);
        let b5, cb5 = statistics::rayleigh_sf(x, -5.0, 2.0);
        chk_val(bad, "rayleigh cdf plus sf", a5 + b5, ca5 + cb5, 1.0, 1e-14);
        let q5, cq5 = statistics::rayleigh_ppf(a5, -5.0, 2.0);
        chk_val(bad, "rayleigh ppf of cdf", q5, cq5, x, 1e-8);
    }
    // the defaults are the standard form
    let d1, cd1 = statistics::laplace_pdf(0.3);
    let d2, cd2 = statistics::laplace_pdf(0.3, 0.0, 1.0);
    chk_val(bad, "laplace default", d1, cd1 + cd2, d2, 0.0);
    // the mirror of the Gumbel family
    let r1, cr1 = statistics::gumbel_r_pdf(0.8);
    let r2, cr2 = statistics::gumbel_l_pdf(-0.8);
    chk_val(bad, "gumbel mirror", r1, cr1 + cr2, r2, 1e-14);
    let s1, cs1 = statistics::gumbel_r_skew();
    let s2, cs2 = statistics::gumbel_l_skew();
    chk_val(bad, "gumbel skew sign", s1 + s2, cs1 + cs2, 0.0, 1e-15);
    // the order of the checks: NaN first, then infinity, then domain
    let _v1, c1 = statistics::laplace_pdf(0.0 / 0.0, 0.0, -1.0);
    chk_code(bad, "NaN before domain", c1, statistics::ERR_NAN);
    let _v2, c2 = statistics::rayleigh_pdf(1.0 / 0.0, 0.0, -1.0);
    chk_code(bad, "infinity before domain", c2, statistics::ERR_INFINITY);
    let _v3, c3 = statistics::gumbel_l_cdf(0.0, 1.0 / 0.0, 1.0);
    chk_code(bad, "infinite location", c3, statistics::ERR_INFINITY);
    let _v4, c4 = statistics::cauchy_mean(0.0, 0.0);
    chk_code(bad, "cauchy moment with a bad scale", c4, statistics::ERR_DOMAIN);
    return bad;
}
