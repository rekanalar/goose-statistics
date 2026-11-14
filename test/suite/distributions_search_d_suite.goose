// The rows of the stable law that the vector file does not hold: the
// bounds of the shapes and the round trip of the quantile.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_search_d_tests() -> i64 {
    var bad: i64 = 0;
    // The cumulative function of the quantile gives the probability back.
    let q1, kq1 = statistics::levy_stable_ppf(0.4, 1.5, 0.0);
    let r1, kr1 = statistics::levy_stable_cdf(q1, 1.5, 0.0);
    chk_val(bad, "levy cdf of the quantile", r1, kq1 + kr1, 0.4, 1e-9);
    let q2, kq2 = statistics::levy_stable_ppf(0.6, 0.8, 0.5);
    let r2, kr2 = statistics::levy_stable_cdf(q2, 0.8, 0.5);
    chk_val(bad, "levy cdf of the quantile", r2, kq2 + kr2, 0.6, 1e-9);
    // The quantile of the end points gives the bounds of the law.
    let e1, ce1 = statistics::levy_stable_ppf(0.0, 1.5, 0.0);
    chk_val(bad, "levy ppf 0", e1, ce1, (-1.0 / 0.0), 1e-12);
    let e2, ce2 = statistics::levy_stable_ppf(1.0, 1.5, 0.0);
    chk_val(bad, "levy ppf 1", e2, ce2, (1.0 / 0.0), 1e-12);
    // The known cases of the law have a closed form.
    let v1, c1 = statistics::levy_stable_pdf(1.0, 2.0, 0.0);
    chk_val(bad, "levy pdf of the normal law", v1, c1, 0.21969564473386122, 1e-12);
    let v2, c2 = statistics::levy_stable_cdf(1.0, 2.0, 0.0);
    chk_val(bad, "levy cdf of the normal law", v2, c2, 0.7602499389065233, 1e-12);
    let v3, c3 = statistics::levy_stable_pdf(1.0, 1.0, 0.0);
    chk_val(bad, "levy pdf of the Cauchy law", v3, c3, 0.15915494309189535, 1e-12);
    let v4, c4 = statistics::levy_stable_cdf(1.0, 1.0, 0.0);
    chk_val(bad, "levy cdf of the Cauchy law", v4, c4, 0.75, 1e-12);
    let v5, c5 = statistics::levy_stable_pdf(1.0, 0.5, 1.0);
    chk_val(bad, "levy pdf of the Levy law", v5, c5, 0.24197072451914337, 1e-12);
    let v6, c6 = statistics::levy_stable_cdf(1.0, 0.5, 1.0);
    chk_val(bad, "levy cdf of the Levy law", v6, c6, 0.31731050786291415, 1e-12);
    // A point that lies just outside the bound of the point near zeta
    // keeps its own value. The bound uses the shape of the law.
    let s1, cs1 = statistics::levy_stable_pdf(0.0018, 0.5, 1.0);
    chk_val(bad, "levy pdf below the bound near zeta", s1, cs1, 1.2040519459848489e-117, 1e-9);
    // A shape outside its range gives the domain fault.
    let _d1, cd1 = statistics::levy_stable_pdf(1.0, 0.0, 0.5);
    chk_code(bad, "levy shape 0", cd1, 5);
    let _d2, cd2 = statistics::levy_stable_pdf(1.0, 2.5, 0.5);
    chk_code(bad, "levy shape above 2", cd2, 5);
    let _d3, cd3 = statistics::levy_stable_cdf(1.0, 1.5, 1.5);
    chk_code(bad, "levy skewness above 1", cd3, 5);
    return bad;
}
