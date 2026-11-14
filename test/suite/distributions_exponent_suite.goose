// The rows of the families of a power or an exponent that the vector file
// does not hold: the support of each law, the end points of the quantile
// and the shape bounds.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_exponent_tests() -> i64 {
    var bad: i64 = 0;
    // A value outside the support gives a density of 0 and not a fault.
    let v1, c1 = statistics::exponpow_pdf(-1.0, 1.5);
    chk_val(bad, "exponpow pdf below the support", v1, c1, 0.0, 1e-12);
    let v2, c2 = statistics::exponweib_pdf(-1.0, 1.5, 2.0);
    chk_val(bad, "exponweib pdf below the support", v2, c2, 0.0, 1e-12);
    let v5, c5 = statistics::powerlognorm_pdf(-1.0, 1.4, 2.0);
    chk_val(bad, "powerlognorm pdf below the support", v5, c5, 0.0, 1e-12);
    // The cumulative function of a value below the support is 0.
    let v6, c6 = statistics::exponpow_cdf(-1.0, 1.5);
    chk_val(bad, "exponpow cdf below the support", v6, c6, 0.0, 1e-12);
    let v7, c7 = statistics::exponweib_cdf(-1.0, 1.5, 2.0);
    chk_val(bad, "exponweib cdf below the support", v7, c7, 0.0, 1e-12);
    let v8, c8 = statistics::powerlognorm_cdf(-1.0, 1.4, 2.0);
    chk_val(bad, "powerlognorm cdf below the support", v8, c8, 0.0, 1e-12);
    // The cumulative function of a value above the support is 1.
    let v9, c9 = statistics::exponweib_cdf(1.0e10, 1.5, 2.0);
    chk_val(bad, "exponweib cdf above the support", v9, c9, 1.0, 1e-12);
    let v10, c10 = statistics::exponpow_cdf(1.0e10, 1.5);
    chk_val(bad, "exponpow cdf above the support", v10, c10, 1.0, 1e-12);
    // The quantile of the end points gives the bounds of the support.
    let e1, ce1 = statistics::exponpow_ppf(0.0, 1.5);
    chk_val(bad, "exponpow ppf 0", e1, ce1, 0.0, 1e-12);
    let e2, ce2 = statistics::powerlognorm_ppf(0.0, 1.4, 2.0);
    chk_val(bad, "powerlognorm ppf 0", e2, ce2, 0.0, 1e-12);
    let e3, ce3 = statistics::exponweib_ppf(0.0, 1.5, 2.0);
    chk_val(bad, "exponweib ppf 0", e3, ce3, 0.0, 1e-12);
    // The cumulative function of the quantile of a probability is that
    // probability.
    let q2, kq2 = statistics::powernorm_ppf(0.4, 1.6);
    chk_code(bad, "powernorm ppf code", kq2, 0);
    let r2, kr2 = statistics::powernorm_cdf(q2, 1.6);
    chk_val(bad, "powernorm cdf of the quantile", r2, kr2, 0.4, 1e-9);
    let q3, kq3 = statistics::powerlognorm_ppf(0.6, 1.4, 2.0);
    chk_code(bad, "powerlognorm ppf code", kq3, 0);
    let r3, kr3 = statistics::powerlognorm_cdf(q3, 1.4, 2.0);
    chk_val(bad, "powerlognorm cdf of the quantile", r3, kr3, 0.6, 1e-9);
    let q4, kq4 = statistics::exponweib_ppf(0.25, 1.5, 2.0);
    chk_code(bad, "exponweib ppf code", kq4, 0);
    let r4, kr4 = statistics::exponweib_cdf(q4, 1.5, 2.0);
    chk_val(bad, "exponweib cdf of the quantile", r4, kr4, 0.25, 1e-9);
    // A shape of 0 or below gives the domain fault.
    let _z3, cz3 = statistics::exponpow_pdf(1.0, 0.0);
    chk_code(bad, "exponpow shape 0", cz3, 5);
    let _z4, cz4 = statistics::exponweib_pdf(1.0, 0.0, 2.0);
    chk_code(bad, "exponweib shape 0", cz4, 5);
    let _z5, cz5 = statistics::powernorm_pdf(1.0, -1.0);
    chk_code(bad, "powernorm shape below 0", cz5, 5);
    let _z6, cz6 = statistics::powerlognorm_pdf(1.0, 1.4, 0.0);
    chk_code(bad, "powerlognorm scale 0", cz6, 5);
    // The law of the power of the normal law has its mean near 0.
    let m1, cm1 = statistics::powernorm_mean(1.6);
    chk_code(bad, "powernorm mean code", cm1, 0);
    check(bad, abs(m1) < 0.5, "the mean of the power of the normal law is near 0");
    // The law of the power of the log-normal law has a heavy tail. Its
    // fourth raw moment is much larger than the fourth moment of the
    // log-normal law itself.
    let m4, cm4 = statistics::powerlognorm_kurt(1.4, 2.0);
    chk_code(bad, "powerlognorm kurt code", cm4, 0);
    check(bad, m4 > 1000.0, "the kurtosis of the power of the log-normal law is large");
    return bad;
}
