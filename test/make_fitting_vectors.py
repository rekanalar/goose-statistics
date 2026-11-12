#!/usr/bin/env python3
"""Writes test/suite/fitting_vectors.goose: the numbers of SciPy 1.18.1 for
logrank and goodness_of_fit, and the checks that use them.

The logrank rows are deterministic and match at 1e-12. The goodness_of_fit
statistic is deterministic and matches at 1e-12. Its p-value is a Monte
Carlo value: the generator measures it with a large SciPy run, and the
test accepts the value of the package within 5 combined standard errors.
SciPy fits the location of the exponential and the lognormal family, so
the generator fixes it at 0 with known_params, as the package does.
Run it from any directory: python3 test/make_fitting_vectors.py.
"""
import os
import warnings
import numpy as np
from scipy import stats

warnings.simplefilter("ignore")

MC_SCIPY = 20000
MC_GOOSE = 3000


def arr(v):
    return "[" + ", ".join(repr(float(x)) for x in v) + "]"


def lr_row(name, ux, rx, uy, ry, alt):
    x = stats.CensoredData(uncensored=ux, right=rx)
    y = stats.CensoredData(uncensored=uy, right=ry)
    r = stats.logrank(x, y, alternative=alt)
    code = {"two-sided": "statistics::ALT_TWO_SIDED", "less": "statistics::ALT_LESS",
            "greater": "statistics::ALT_GREATER"}[alt]
    return (f'    lrow(bad, {arr(ux)}, {arr(rx)}, {arr(uy)}, {arr(ry)}, {code}, '
            f'{float(r.statistic)!r}, {float(r.pvalue)!r}, "{name}");')


def gof_row(name, data, fam, famc, stat, statc):
    known = {"loc": 0} if fam in (stats.expon, stats.lognorm) else {}
    r = stats.goodness_of_fit(fam, data, known_params=known, statistic=stat,
                              n_mc_samples=MC_SCIPY, rng=np.random.default_rng(11))
    return (f'    grow(bad, {arr(data)}, {famc}, {statc}, '
            f'{float(r.statistic)!r}, {float(r.pvalue)!r}, "{name}");')


def main():
    o = []
    o.append("// The answers of SciPy 1.18.1 for the logrank test and the")
    o.append("// goodness-of-fit test. The generator test/make_fitting_vectors.py")
    o.append("// writes this file. It holds the numbers and the checks that use them.")
    o.append("namespace suite;")
    o.append("")
    o.append("import std;")
    o.append("import math;")
    o.append("import statistics from \"../..\";")
    o.append("")
    o.append("""fn lrow(bad: i64&, ux: const f64[:], rx: const f64[:], uy: const f64[:], ry: const f64[:], alt: i64, stat: f64, p: f64, name: u8[:]) {
    var none: f64[>..] = [];
    let x = statistics::CensoredData { uncensored: ux, left: none[..], right: rx, interval: none[..] };
    let y = statistics::CensoredData { uncensored: uy, left: none[..], right: ry, interval: none[..] };
    let r, c = statistics::logrank(x, y, alt);
    if c == 0 && near_tol(r.statistic, stat, 1e-12) && near_tol(r.pvalue, p, 1e-12) { return; }
    bad += 1;
    print("FAIL ", name, " got ", r.statistic, " ", r.pvalue, " code ", c, " want ", stat, " ", p, "\\n");
}

// The statistic matches at 1e-12. The p-value is a Monte Carlo value: it
// must lie within 5 combined standard errors of the SciPy value.
fn grow(bad: i64&, data: const f64[:], family: i64, stat: i64, want_stat: f64, want_p: f64, name: u8[:]) {
    var out: f64[>..] = [];
    var g = rng { s: 12345 };
    let r, c = statistics::goodness_of_fit(out, data, family, stat, %d, g);
    let se = sqrt(want_p * (1.0 - want_p) * (1.0 / %d.0 + 1.0 / %d.0)) + 1e-4;
    if c == 0 && near_tol(r.statistic, want_stat, 1e-12) && abs(r.pvalue - want_p) <= 5.0 * se { return; }
    bad += 1;
    print("FAIL ", name, " got ", r.statistic, " ", r.pvalue, " code ", c, " want ", want_stat, " ", want_p, "\\n");
}
""" % (MC_GOOSE, MC_GOOSE, MC_SCIPY))
    o.append("fn run_fitting_vectors() -> i64 {")
    o.append("    var bad: i64 = 0;")
    # logrank rows
    o.append(lr_row("logrank, the measured pair", [2, 3, 5, 8], [4, 6, 9], [1, 2, 2, 4], [3, 5, 7], "two-sided"))
    o.append(lr_row("logrank, less", [2, 3, 5, 8], [4, 6, 9], [1, 2, 2, 4], [3, 5, 7], "less"))
    o.append(lr_row("logrank, greater", [2, 3, 5, 8], [4, 6, 9], [1, 2, 2, 4], [3, 5, 7], "greater"))
    ux = [6, 13, 21, 30, 37, 38, 49, 50, 63, 79, 86, 98, 202, 219]
    rx = [31, 47, 80, 82, 82, 149]
    uy = [10, 10, 12, 13, 14, 15, 16, 17, 18, 20, 24, 24, 25, 28, 30, 33, 35, 37, 40, 40, 46, 48, 76, 81, 82, 91, 112, 181]
    ry = [34, 40, 70]
    o.append(lr_row("logrank, the glioma data", ux, rx, uy, ry, "two-sided"))
    o.append(lr_row("logrank, swapped groups", uy, ry, ux, rx, "two-sided"))
    o.append(lr_row("logrank, no censoring", [1.5, 2.5, 3.5, 6.0], [], [2.0, 4.5, 7.0, 8.0, 9.5], [], "two-sided"))
    o.append(lr_row("logrank, a censored time equals an event", [3, 5, 5], [5, 8], [2, 5, 9], [5], "greater"))
    # goodness of fit rows
    rng = np.random.default_rng(2026)
    nrm = list(np.round(rng.normal(10, 2, 30), 3))
    skew = list(np.round(rng.exponential(2.0, 30), 3))
    lgn = list(np.round(rng.lognormal(0.5, 0.6, 30), 3))
    uni = list(np.round(rng.uniform(1, 4, 30), 3))
    o.append(gof_row("gof, normal, anderson", nrm, stats.norm, "statistics::FAMILY_NORMAL", "ad", "statistics::GOF_ANDERSON"))
    o.append(gof_row("gof, normal, kolmogorov", nrm, stats.norm, "statistics::FAMILY_NORMAL", "ks", "statistics::GOF_KOLMOGOROV"))
    o.append(gof_row("gof, normal, cramer", nrm, stats.norm, "statistics::FAMILY_NORMAL", "cvm", "statistics::GOF_CRAMER"))
    o.append(gof_row("gof, skewed data against normal", skew, stats.norm, "statistics::FAMILY_NORMAL", "ad", "statistics::GOF_ANDERSON"))
    o.append(gof_row("gof, exponential, anderson", skew, stats.expon, "statistics::FAMILY_EXPONENTIAL", "ad", "statistics::GOF_ANDERSON"))
    o.append(gof_row("gof, exponential, cramer", skew, stats.expon, "statistics::FAMILY_EXPONENTIAL", "cvm", "statistics::GOF_CRAMER"))
    o.append(gof_row("gof, lognormal, anderson", lgn, stats.lognorm, "statistics::FAMILY_LOGNORMAL", "ad", "statistics::GOF_ANDERSON"))
    o.append(gof_row("gof, lognormal, kolmogorov", lgn, stats.lognorm, "statistics::FAMILY_LOGNORMAL", "ks", "statistics::GOF_KOLMOGOROV"))
    o.append(gof_row("gof, uniform, kolmogorov", uni, stats.uniform, "statistics::FAMILY_UNIFORM", "ks", "statistics::GOF_KOLMOGOROV"))
    o.append(gof_row("gof, uniform, cramer", uni, stats.uniform, "statistics::FAMILY_UNIFORM", "cvm", "statistics::GOF_CRAMER"))
    o.append("    return bad;")
    o.append("}")
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "suite", "fitting_vectors.goose")
    with open(path, "w") as f:
        f.write("\n".join(o) + "\n")


main()
