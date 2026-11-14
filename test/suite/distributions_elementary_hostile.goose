// Hostile input for the elementary families. A NaN gives the NaN fault, an
// infinity gives the infinity fault and a probability outside 0 to 1 gives
// the domain fault. No answer with the code OK is NaN or infinite.
namespace suite;

import std;
import math;
import statistics from "../..";

// The code of the density, the cumulative function and the survival
// function of one value for the eight families.
fn fam_elem_value(bad: i64&, x: f64, want: i64, what: u8[:]) {
    let a1, b1 = statistics::anglit_pdf(x);
    let a2, b2 = statistics::arcsine_cdf(x);
    let a3, b3 = statistics::cosine_sf(x);
    let a4, b4 = statistics::gibrat_pdf(x);
    let a5, b5 = statistics::halfcauchy_cdf(x);
    let a6, b6 = statistics::halflogistic_sf(x);
    let a7, b7 = statistics::halfnorm_pdf(x);
    let a8, b8 = statistics::hypsecant_cdf(x);
    chk_code(bad, what, b1, want);
    chk_code(bad, what, b2, want);
    chk_code(bad, what, b3, want);
    chk_code(bad, what, b4, want);
    chk_code(bad, what, b5, want);
    chk_code(bad, what, b6, want);
    chk_code(bad, what, b7, want);
    chk_code(bad, what, b8, want);
}

// The code of the quantile of one probability for the eight families.
fn fam_elem_prob(bad: i64&, p: f64, want: i64, what: u8[:]) {
    let a1, b1 = statistics::anglit_ppf(p);
    let a2, b2 = statistics::arcsine_ppf(p);
    let a3, b3 = statistics::cosine_ppf(p);
    let a4, b4 = statistics::gibrat_ppf(p);
    let a5, b5 = statistics::halfcauchy_ppf(p);
    let a6, b6 = statistics::halflogistic_ppf(p);
    let a7, b7 = statistics::halfnorm_ppf(p);
    let a8, b8 = statistics::hypsecant_ppf(p);
    chk_code(bad, what, b1, want);
    chk_code(bad, what, b2, want);
    chk_code(bad, what, b3, want);
    chk_code(bad, what, b4, want);
    chk_code(bad, what, b5, want);
    chk_code(bad, what, b6, want);
    chk_code(bad, what, b7, want);
    chk_code(bad, what, b8, want);
}

fn run_distributions_elementary_hostile() -> i64 {
    var bad: i64 = 0;
    fam_elem_value(bad, (0.0 / 0.0), 3, "a NaN value");
    fam_elem_value(bad, (1.0 / 0.0), 4, "an infinite value");
    fam_elem_value(bad, -(1.0 / 0.0), 4, "a negative infinite value");
    fam_elem_prob(bad, (0.0 / 0.0), 3, "a NaN probability");
    fam_elem_prob(bad, -0.1, 5, "a probability below 0");
    fam_elem_prob(bad, 1.1, 5, "a probability above 1");
    fam_elem_prob(bad, (1.0 / 0.0), 5, "an infinite probability");
    // An infinite bound of the support does not fit in a result.
    let g, cg = statistics::gibrat_ppf(1.0);
    chk_code(bad, "the quantile of 1 of an unbounded law", cg, 8);
    let h, ch = statistics::hypsecant_ppf(0.0);
    chk_code(bad, "the quantile of 0 of a law on the whole line", ch, 8);
    // A value far in the tail gives a density of 0 and no fault.
    let t1, ct1 = statistics::hypsecant_pdf(1e5);
    chk_val(bad, "hypsecant far tail", t1, ct1, 0.0, 1e-300);
    let t2, ct2 = statistics::halfnorm_pdf(1e5);
    chk_val(bad, "halfnorm far tail", t2, ct2, 0.0, 1e-300);
    // The density at the pole of the arcsine law does not fit.
    let ap, cap = statistics::arcsine_pdf(0.0);
    chk_code(bad, "the arcsine density at the pole", cap, 8);
    // A random sweep: every answer with the code OK is finite.
    var st: u64 = 5531;
    for k in 0..tests_count() {
        let x = ((rnd(st) % 2001) as f64 - 1000.0) / 100.0;
        let p = (rnd(st) % 1001) as f64 / 1000.0;
        let v1, c1 = statistics::cosine_pdf(x);
        let v2, c2 = statistics::halflogistic_cdf(x);
        let v3, c3 = statistics::gibrat_ppf(p);
        let v4, c4 = statistics::cosine_ppf(p);
        hostile_report(bad, c1 == 0 && !is_nan(v1) && !is_inf(v1), "cosine pdf", k, c1);
        hostile_report(bad, c2 == 0 && !is_nan(v2) && !is_inf(v2), "halflogistic cdf", k, c2);
        hostile_report(bad, (c3 == 0 && !is_nan(v3) && !is_inf(v3)) || c3 == 8, "gibrat ppf", k, c3);
        hostile_report(bad, c4 == 0 && !is_nan(v4) && !is_inf(v4), "cosine ppf", k, c4);
    }
    return bad;
}
