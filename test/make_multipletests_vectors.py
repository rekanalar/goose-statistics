#!/usr/bin/env python3
"""Writes test/suite/multipletests_vectors.goose: the numbers of the
reference libraries and the checks that use them.

Each list of p-values carries the adjusted values of statsmodels 0.15.0
(Hochberg, Hommel, Benjamini-Yekutieli, Benjamini-Hochberg, Holm,
Bonferroni) and the combined values of SciPy 1.18.1 (Fisher, Stouffer
with and without weights, Tippett). Run it from the package directory.
Fisher and Stouffer skip a list that holds 0 or 1, which they refuse.
"""
import random
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

random.seed(20261112)


def arr(v):
    return "[" + ", ".join(repr(float(x)) for x in v) + "]"


def make_list(n, kind):
    if kind == "ties":
        return [round(random.random(), 2) for _ in range(n)]
    if kind == "edge":
        v = [random.random() for _ in range(n)]
        v[0], v[-1] = 0.0, 1.0
        return v
    if kind == "small":
        return [random.random() ** 4 for _ in range(n)]
    return [random.random() for _ in range(n)]


METHODS = [("ADJUST_HOCHBERG", "simes-hochberg"), ("ADJUST_HOMMEL", "hommel"),
           ("ADJUST_BY", "fdr_by"), ("ADJUST_BH", "fdr_bh"),
           ("ADJUST_HOLM", "holm"), ("ADJUST_BONFERRONI", "bonferroni")]


def main():
    out = []
    out.append("// The answers of statsmodels 0.15.0 and SciPy 1.18.1 for the")
    out.append("// multiple-test corrections. The generator")
    out.append("// test/make_multipletests_vectors.py writes this file.")
    out.append("// It holds the numbers and the checks that use them.")
    out.append("namespace suite;")
    out.append("")
    out.append("import std;")
    out.append("import statistics from \"../..\";")
    out.append("")
    out.append("fn arow(bad: i64&, ps: const f64[:], method: i64, want: const f64[:], name: u8[:]) {")
    out.append("    var out: f64[>..] = [];")
    out.append("    let n, c = statistics::multipletests(out, ps, method);")
    out.append("    var ok = c == 0 && n == want.len && out.len == n;")
    out.append("    if ok {")
    out.append("        for i in 0..n { if !near_tol(out[i], want[i], 1e-13) { ok = false; } }")
    out.append("    }")
    out.append("    if ok { return; }")
    out.append("    bad += 1;")
    out.append("    if bad <= 12 { print(\"FAIL \", name, \" code \", c, \" n \", n, \"\\n\"); }")
    out.append("}")
    out.append("")
    out.append("fn frow(bad: i64&, ps: const f64[:], method: i64, want: const f64[:], name: u8[:]) {")
    out.append("    var out: f64[>..] = [];")
    out.append("    let n, c = statistics::false_discovery_control(out, ps, method);")
    out.append("    var ok = c == 0 && n == want.len && out.len == n;")
    out.append("    if ok {")
    out.append("        for i in 0..n { if !near_tol(out[i], want[i], 1e-13) { ok = false; } }")
    out.append("    }")
    out.append("    if ok { return; }")
    out.append("    bad += 1;")
    out.append("    if bad <= 12 { print(\"FAIL \", name, \" code \", c, \" n \", n, \"\\n\"); }")
    out.append("}")
    out.append("")
    out.append("fn crow(bad: i64&, ps: const f64[:], method: i64, w: const f64[:], stat: f64, p: f64, df: f64, tol: f64, name: u8[:]) {")
    out.append("    let r, c = statistics::combine_pvalues(ps, method, w);")
    out.append("    var ok = c == 0 && near_tol(r.statistic, stat, tol) && near_tol(r.pvalue, p, tol);")
    out.append("    if method == statistics::COMBINE_FISHER { ok = ok && r.df == df; }")
    out.append("    if ok { return; }")
    out.append("    bad += 1;")
    out.append("    if bad <= 12 { print(\"FAIL \", name, \" got \", r.statistic, \" \", r.pvalue, \" code \", c, \" want \", stat, \" \", p, \"\\n\"); }")
    out.append("}")
    out.append("")
    out.append("fn run_multipletests_vectors() -> i64 {")
    out.append("    var bad: i64 = 0;")
    cases = []
    for n in (3, 5, 10, 50, 200):
        for kind in ("plain", "ties", "small", "edge"):
            cases.append((n, kind))
    cases += [(2, "plain"), (2, "ties"), (1, "plain")]
    for idx, (n, kind) in enumerate(cases):
        ps = make_list(n, kind)
        tag = "%s%d" % (kind, n)
        out.append("    {")
        out.append("        let w0: f64[0] = [];")
        out.append("        let ps: f64[%d] = %s;" % (n, arr(ps)))
        for const, name in METHODS:
            want = multipletests(ps, method=name)[1]
            want = np.minimum(want, 1.0)
            out.append("        { let w: f64[%d] = %s; arow(bad, ps[..], statistics::%s, w[..], \"%s %s\"); }"
                       % (n, arr(want), const, name, tag))
        out.append("        arow(bad, ps[..], statistics::ADJUST_NONE, ps[..], \"none %s\");" % tag)
        # The two-step form of SciPy equals the BY engine.
        fdc = stats.false_discovery_control(ps, method="by")
        out.append("        { let w: f64[%d] = %s; frow(bad, ps[..], statistics::ADJUST_BY, w[..], \"fdc by %s\"); }"
                   % (n, arr(fdc), tag))
        fdb = stats.false_discovery_control(ps, method="bh")
        out.append("        { let w: f64[%d] = %s; frow(bad, ps[..], statistics::ADJUST_BH, w[..], \"fdc bh %s\"); }"
                   % (n, arr(fdb), tag))
        tp = stats.combine_pvalues(ps, method="tippett")
        out.append("        crow(bad, ps[..], statistics::COMBINE_TIPPETT, %s, %r, %r, 0.0, 1e-13, \"tippett %s\");"
                   % ("w0[..]", float(tp.statistic), float(tp.pvalue), tag))
        if 0.0 not in ps and 1.0 not in ps:
            f = stats.combine_pvalues(ps, method="fisher")
            # The tail of a large statistic keeps about 12 digits, so the
            # p-value row allows a wider bound than the other rows.
            out.append("        crow(bad, ps[..], statistics::COMBINE_FISHER, w0[..], %r, %r, %r, 1e-11, \"fisher %s\");"
                       % (float(f.statistic), float(f.pvalue), float(2 * n), tag))
            s = stats.combine_pvalues(ps, method="stouffer")
            out.append("        crow(bad, ps[..], statistics::COMBINE_STOUFFER, w0[..], %r, %r, 0.0, 1e-12, \"stouffer %s\");"
                       % (float(s.statistic), float(s.pvalue), tag))
            wt = [round(random.random() * 3 + 0.1, 3) for _ in range(n)]
            sw = stats.combine_pvalues(ps, method="stouffer", weights=wt)
            out.append("        { let wt: f64[%d] = %s; crow(bad, ps[..], statistics::COMBINE_STOUFFER, wt[..], %r, %r, 0.0, 1e-12, \"stouffer weighted %s\"); }"
                       % (n, arr(wt), float(sw.statistic), float(sw.pvalue), tag))
        out.append("    }")
    out.append("    return bad;")
    out.append("}")
    text = "\n".join(out) + "\n"
    with open("test/suite/multipletests_vectors.goose", "w") as f:
        f.write(text)
    print("wrote test/suite/multipletests_vectors.goose: %d lines" % text.count("\n"))


if __name__ == "__main__":
    main()
