#!/usr/bin/env python3
"""Writes test/suite/multicomp_vectors.goose: the numbers of SciPy 1.18.1
for tukey_hsd, dunnett, bayes_mvs and mvsdist, and the checks that use them.

SciPy computes the Dunnett p-value with a randomized quasi Monte Carlo
that has an error of about 1.5e-4, so the checks use a bound of 1e-3 for
the p-value and 5e-3 for the interval. SciPy computes the Tukey p-value as
one minus the cdf, so a p-value below 1e-12 holds an absolute error of
about 1e-15; the checks add that absolute bound.
Run it from any directory: python3 test/make_multicomp_vectors.py
"""
import os
import random
import warnings

import numpy as np
from scipy import stats

warnings.simplefilter("ignore")
random.seed(20261112)

A = [2.1, 3.4, 1.9, 5.6, 4.4, 3.8, 2.7, 4.9]
B = [3.9, 5.2, 4.8, 6.1, 5.5, 4.1, 6.6, 5.0]
C = [1.2, 2.5, 2.2, 3.1, 1.8, 2.9, 2.4, 3.3]


def num(x):
    x = float(x)
    if x == float("inf"):
        return "(1.0 / 0.0)"
    if x == float("-inf"):
        return "(-1.0 / 0.0)"
    return repr(x)


def arr(v):
    return "[" + ", ".join(num(x) for x in v) + "]"


def flat(groups):
    return [x for g in groups for x in g]


def make_groups(sizes, shift):
    return [[round(random.gauss(shift * i, 1.0 + 0.3 * i), 3) for _ in range(n)]
            for i, n in enumerate(sizes)]


def tukey_rows(groups, level, equal_var):
    r = stats.tukey_hsd(*groups, equal_var=equal_var)
    ci = r.confidence_interval(level)
    k = len(groups)
    out = []
    for i in range(k - 1):
        for j in range(i + 1, k):
            out += [r.statistic[i][j], r.pvalue[i][j], ci.low[i][j], ci.high[i][j]]
    return out


def dunnett_rows(groups, control, alt, level, seed):
    others = [g for i, g in enumerate(groups) if i != control]
    r = stats.dunnett(*others, control=groups[control], alternative=alt,
                      rng=np.random.default_rng(seed))
    ci = r.confidence_interval(level)
    out = []
    for i in range(len(others)):
        out += [r.statistic[i], r.pvalue[i], ci.low[i], ci.high[i]]
    return out


def mvs_numbers(xs, level):
    m, v, s = stats.bayes_mvs(xs, level)
    d = stats.mvsdist(xs)
    return [m.statistic, m.minmax[0], m.minmax[1], v.statistic, v.minmax[0],
            v.minmax[1], s.statistic, s.minmax[0], s.minmax[1]], d


def dist_fields(xs):
    """gaussian, mean, mean_scale, mean_df, variance, variance_scale,
    variance_shape, deviation, deviation_scale, deviation_shape."""
    md, vd, sd = stats.mvsdist(xs)
    if len(xs) > 1000:
        return [True, md.mean(), md.std(), 0.0, vd.mean(), vd.std(), 0.0,
                sd.mean(), sd.std(), 0.0]
    return [False, md.kwds["loc"], md.kwds["scale"], md.args[0],
            vd.mean(), vd.kwds["scale"], vd.args[0],
            sd.mean(), sd.kwds["scale"], sd.args[0]]


TUKEY = [
    ("oracle equal n", [A, B, C], 0.95, True),
    ("oracle welch", [A, B, C], 0.95, False),
    ("oracle level 0.99", [A, B, C], 0.99, True),
    ("two groups", [A, B], 0.95, True),
    ("unequal n", make_groups([5, 9, 12, 7], 0.8), 0.95, True),
    ("unequal n welch", make_groups([5, 9, 12, 7], 0.8), 0.9, False),
    ("six groups", make_groups([6] * 6, 0.5), 0.95, True),
    ("small df", make_groups([3, 3, 4], 1.5), 0.95, True),
    ("large df", make_groups([40, 45, 38], 0.4), 0.95, True),
    ("far groups", make_groups([8, 8, 8], 6.0), 0.95, True),
]

DUNNETT = [
    ("oracle two-sided", [A, B, C], 0, 0, 0.95, 1),
    ("oracle greater", [A, B, C], 0, 2, 0.95, 2),
    ("oracle less", [A, B, C], 0, 1, 0.95, 3),
    ("control last", [A, B, C], 2, 0, 0.95, 4),
    ("unequal n", make_groups([6, 9, 12, 7], 0.9), 1, 0, 0.95, 5),
    ("five groups", make_groups([8] * 5, 0.7), 0, 0, 0.9, 6),
    ("one group", [A, C], 1, 0, 0.95, 7),
]

BAYES = [
    ("oracle", A, 0.9),
    ("level 0.95", A, 0.95),
    ("n 4", [1.0, 2.5, 2.0, 4.1], 0.9),
    ("n 30", [round(random.gauss(5, 2), 3) for _ in range(30)], 0.9),
    ("n 200", [round(random.gauss(-3, 0.5), 3) for _ in range(200)], 0.99),
    ("n 1000", [round(random.gauss(0, 1), 3) for _ in range(1000)], 0.9),
    ("n 1001", [round(random.gauss(10, 3), 3) for _ in range(1001)], 0.95),
]


def main():
    o = []
    w = o.append
    w("// The answers of SciPy 1.18.1 for tukey_hsd, dunnett, bayes_mvs and")
    w("// mvsdist. The generator test/make_multicomp_vectors.py writes this file.")
    w("// SciPy draws the Dunnett p-value with a randomized quasi Monte Carlo of")
    w("// about 1.5e-4 error, so the Dunnett rows use looser bounds.")
    w("namespace suite;")
    w("")
    w("import std;")
    w("import statistics from \"../..\";")
    w("")
    w("// A p-value or a statistic within a relative bound and an absolute bound.")
    w("fn mrow_near(a: f64, e: f64, rel: f64, ab: f64) -> bool {")
    w("    return abs(a - e) <= rel * abs(e) + ab;")
    w("}")
    w("")
    w("fn trow(bad: i64&, groups: const f64[:], sizes: const i64[:], level: f64, equal_var: bool, want: const f64[:], name: u8[:]) {")
    w("    var out: f64[>..] = [];")
    w("    let n, c = statistics::tukey_hsd(out, groups, sizes, level, equal_var);")
    w("    var ok = c == 0 && out.len == want.len && n * 4 == want.len;")
    w("    if ok {")
    w("        for i in 0..want.len {")
    w("            var rel: f64 = 1e-10;")
    w("            var ab: f64 = 1e-12;")
    w("            if i % 4 == 1 { ab = 1e-14; }")
    w("            if i % 4 != 1 { ab = 1e-12 * (1.0 + abs(want[i])); }")
    w("            if !mrow_near(out[i], want[i], rel, ab) { ok = false; }")
    w("        }")
    w("    }")
    w("    if ok { return; }")
    w("    bad += 1;")
    w("    print(\"FAIL tukey \", name, \" code \", c, \"\\n\");")
    w("}")
    w("")
    w("fn drow(bad: i64&, groups: const f64[:], sizes: const i64[:], control: i64, alt: i64, level: f64, want: const f64[:], name: u8[:]) {")
    w("    var out: f64[>..] = [];")
    w("    let n, c = statistics::dunnett(out, groups, sizes, control, alt, level);")
    w("    var ok = c == 0 && out.len == want.len && n * 4 == want.len;")
    w("    if ok {")
    w("        for i in 0..want.len {")
    w("            let col = i % 4;")
    w("            if col == 0 && !mrow_near(out[i], want[i], 1e-12, 1e-12) { ok = false; }")
    w("            if col == 1 && !mrow_near(out[i], want[i], 0.0, 1e-3) { ok = false; }")
    w("            if col >= 2 {")
    w("                if abs(want[i]) > 1.0e300 {")
    w("                    if out[i] != want[i] { ok = false; }")
    w("                } else if !mrow_near(out[i], want[i], 0.0, 5e-3 * (1.0 + abs(want[i]))) { ok = false; }")
    w("            }")
    w("        }")
    w("    }")
    w("    if ok { return; }")
    w("    bad += 1;")
    w("    print(\"FAIL dunnett \", name, \" code \", c, \"\\n\");")
    w("}")
    w("")
    w("fn brow(bad: i64&, xs: const f64[:], level: f64, want: const f64[:], name: u8[:]) {")
    w("    let r, c = statistics::bayes_mvs(xs, level);")
    w("    var ok = c == 0;")
    w("    if ok {")
    w("        let got: f64[9] = [r.mean, r.mean_interval.low, r.mean_interval.high, r.variance, r.variance_interval.low, r.variance_interval.high, r.deviation, r.deviation_interval.low, r.deviation_interval.high];")
    w("        for i in 0..9 { if !near_tol(got[i], want[i], 1e-9) { ok = false; } }")
    w("    }")
    w("    if ok { return; }")
    w("    bad += 1;")
    w("    print(\"FAIL bayes_mvs \", name, \" code \", c, \"\\n\");")
    w("}")
    w("")
    w("fn mrow(bad: i64&, xs: const f64[:], gaussian: bool, want: const f64[:], name: u8[:]) {")
    w("    let d, c = statistics::mvsdist(xs);")
    w("    var ok = c == 0 && d.gaussian == gaussian;")
    w("    if ok {")
    w("        let got: f64[9] = [d.mean, d.mean_scale, d.mean_df, d.variance, d.variance_scale, d.variance_shape, d.deviation, d.deviation_scale, d.deviation_shape];")
    w("        for i in 0..9 { if !near_tol(got[i], want[i], 1e-9) { ok = false; } }")
    w("    }")
    w("    if ok { return; }")
    w("    bad += 1;")
    w("    print(\"FAIL mvsdist \", name, \" code \", c, \"\\n\");")
    w("}")
    w("")
    w("fn run_multicomp_vectors() -> i64 {")
    w("    var bad: i64 = 0;")
    for name, groups, level, ev in TUKEY:
        sizes = [len(g) for g in groups]
        want = tukey_rows(groups, level, ev)
        w("    {")
        w("        let g = %s;" % arr(flat(groups)))
        w("        let s: i64[%d] = %s;" % (len(sizes), sizes))
        w("        let want = %s;" % arr(want))
        w("        trow(bad, g[..], s[..], %r, %s, want[..], \"%s\");" % (level, "true" if ev else "false", name))
        w("    }")
    for name, groups, control, alt, level, seed in DUNNETT:
        sizes = [len(g) for g in groups]
        alts = {0: "two-sided", 1: "less", 2: "greater"}[alt]
        want = dunnett_rows(groups, control, alts, level, seed)
        w("    {")
        w("        let g = %s;" % arr(flat(groups)))
        w("        let s: i64[%d] = %s;" % (len(sizes), sizes))
        w("        let want = %s;" % arr(want))
        w("        drow(bad, g[..], s[..], %d, %d, %r, want[..], \"%s\");" % (control, alt, level, name))
        w("    }")
    for name, xs, level in BAYES:
        nums, _ = mvs_numbers(xs, level)
        w("    {")
        w("        let x = %s;" % arr(xs))
        w("        let want = %s;" % arr(nums))
        w("        brow(bad, x[..], %r, want[..], \"%s\");" % (level, name))
        w("    }")
    for name, xs, level in BAYES:
        f = dist_fields(xs)
        w("    {")
        w("        let x = %s;" % arr(xs))
        w("        let want = %s;" % arr(f[1:]))
        w("        mrow(bad, x[..], %s, want[..], \"%s\");" % ("true" if f[0] else "false", name))
        w("    }")
    w("    return bad;")
    w("}")
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "suite", "multicomp_vectors.goose"), "w") as f:
        f.write("\n".join(o) + "\n")


if __name__ == "__main__":
    main()
