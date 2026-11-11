#!/usr/bin/env python3
"""Writes test/suite/normality_vectors.goose: numbers only, no algorithm.

Each row holds one sample and the answers of SciPy 1.18.1. The
Anderson-Darling statistic is the formula of SciPy on an exact fit. Its p-value is the
asymptotic rule applied to that statistic; the row comment names the
rule. The script pins the anchors of the acceptance sample first.
"""
import math
import random

import numpy as np
from scipy import stats

random.seed(20261111)

ACC = [1.1, 2.3, 3.7, 4.2, 5.9, 6.1, 7.4, 8.8, 2.0, 2.9, 3.1, 5.5,
       6.6, 7.2, 9.1, 9.9]


def arr(v):
    return "[" + ", ".join(repr(float(x)) for x in v) + "]"


def ad_p(a2, n):
    adj = a2 * (1.0 + 0.75 / n + 2.25 / (n * n))
    if adj >= 0.6:
        return math.exp(1.2937 - 5.709 * adj + 0.0186 * adj * adj)
    if adj >= 0.34:
        return math.exp(0.9177 - 4.279 * adj - 1.38 * adj * adj)
    if adj >= 0.2:
        return 1.0 - math.exp(-8.318 + 42.796 * adj - 59.938 * adj * adj)
    return 1.0 - math.exp(-13.436 + 101.14 * adj - 223.73 * adj * adj)


def ad_stat(x):
    # The weighted logarithms cancel against n.
    # The statistic carries noise of about 1e-10.
    # The rows of this test allow that noise.
    # SciPy's formula with logcdf and logsf. The fit uses exact sums, so
    # the statistic differs from stats.anderson by its float mean (1e-11).
    xs = np.sort(np.asarray(x, dtype=float))
    n = len(xs)
    mu = math.fsum(x) / n
    sd = math.sqrt(math.fsum((float(v) - mu) ** 2 for v in x) / (n - 1))
    i = np.arange(1, n + 1)
    a2 = -n - math.fsum((2 * i - 1.0) / n
                        * (stats.norm.logcdf(xs, mu, sd)
                           + stats.norm.logsf(xs, mu, sd)[::-1]))
    return float(a2), ad_p(float(a2), n)


def call_row(name, sample, stat, p, alt=None, tol=1e-13):
    if alt is None:
        call = "statistics::%s(%s)" % (name, arr(sample))
    else:
        call = "statistics::%s(%s, %d)" % (name, arr(sample), alt)
    return ("    { let r, c = %s; "
            "nrow(bad, r, c, %r, %r, %r, \"%s\"); }"
            % (call, float(stat), float(p), tol, name))


def emit(out, sample, tag):
    x = np.asarray(sample, dtype=float)
    n = len(x)
    if n >= 8:
        st, ps = stats.skewtest(x)
        out.append(call_row("skewtest", sample, st, ps, 0) + " // two-sided")
        st, ps = stats.skewtest(x, alternative="less")
        out.append(call_row("skewtest", sample, st, ps, 1) + " // less")
        st, ps = stats.skewtest(x, alternative="greater")
        out.append(call_row("skewtest", sample, st, ps, 2) + " // greater")
    st, ps = stats.kurtosistest(x)
    out.append(call_row("kurtosistest", sample, st, ps, 0) + " // two-sided")
    st, ps = stats.kurtosistest(x, alternative="less")
    out.append(call_row("kurtosistest", sample, st, ps, 1) + " // less")
    st, ps = stats.kurtosistest(x, alternative="greater")
    out.append(call_row("kurtosistest", sample, st, ps, 2) + " // greater")
    if n >= 8:
        st, ps = stats.normaltest(x)
        out.append(call_row("normaltest", sample, st, ps))
        a2, p = ad_stat(x)
        out.append(call_row("anderson_darling", sample, a2, p, None, 1e-10)
                   + " // asymptotic rule")
    st, ps = stats.jarque_bera(x)
    out.append(call_row("jarque_bera", sample, st, ps))
    _ = tag


def main():
    out = []
    out.append("// The answers of SciPy 1.18.1 for the normality tests. The")
    out.append("// generator test/make_normality_vectors.py writes this file.")
    out.append("// It holds numbers only. The Anderson-Darling p-value is the")
    out.append("// asymptotic rule applied to the sequential sum.")
    out.append("namespace suite;")
    out.append("")
    out.append("import std;")
    out.append("import statistics from \"../..\";")
    out.append("")
    out.append("fn nrow(bad: i64&, r: statistics::TestResult, c: i64, "
               "stat: f64, p: f64, tol: f64, name: u8[:]) {")
    out.append("    var ok = c == 0 && near_tol(r.statistic, stat, tol) "
               "&& near_tol(r.pvalue, p, tol);")
    out.append("    if ok { return; }")
    out.append("    bad += 1;")
    out.append("    if bad <= 12 { print(\"FAIL \", name, \" got \", "
               "r.statistic, \" \", r.pvalue, \" code \", c, \" want \", "
               "stat, \" \", p, \"\\n\"); }")
    out.append("}")
    out.append("")
    out.append("fn run_normality_vectors() -> i64 {")
    out.append("    var bad: i64 = 0;")
    emit(out, ACC, "acc")
    samples = {}
    for n in (5, 6, 7, 8, 20, 100, 500, 2000):
        samples[n] = [round(random.gauss(5.0, 2.0), 4) + random.random() * 1e-6
                      for _ in range(n)]
    for n in sorted(samples):
        emit(out, samples[n], "n%d" % n)
    # Skewed samples reach the upper branches of the Anderson-Darling rule.
    for n in (20, 100, 500):
        emit(out, [round(random.expovariate(0.5), 4) for _ in range(n)],
             "e%d" % n)
    # A symmetric sample makes the skewness zero; a two-point sample makes
    # the kurtosis denominator negative. Both are valid inputs.
    emit(out, [1.0, 2, 3, 4, 4, 3, 2, 1], "sym")
    emit(out, [1.0] * 50 + [2.0] * 50, "two")
    # The smallest sample of jarque_bera runs and the others refuse it.
    two = [1.5, 2.5]
    st, ps = stats.jarque_bera(np.asarray(two))
    out.append(call_row("jarque_bera", two, st, ps))
    out.append("    return bad;")
    out.append("}")
    text = "\n".join(out) + "\n"
    with open("test/suite/normality_vectors.goose", "w") as f:
        f.write(text)
    print("wrote test/suite/normality_vectors.goose: %d lines" % len(out))


if __name__ == "__main__":
    main()
