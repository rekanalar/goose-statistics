#!/usr/bin/env python3
"""Writes test/suite/variance_vectors.goose: the numbers of SciPy 1.18.1 for
the tests of variance, scale and rank, and the checks that use them.

The first rows are the samples a, b and c of the measured oracle. The other
rows are seeded random samples of several sizes, with ties and without.
A row that SciPy refuses is skipped. For quantile_test, a bound that SciPy
gives as NaN becomes an infinite bound, as the package documents.
Run it from any directory with the system python3. It needs SciPy.
"""
import math
import os
import random
import warnings

import numpy as np
from scipy import stats

warnings.simplefilter("ignore")
random.seed(20261113)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "suite", "variance_vectors.goose")
ALT = {"two-sided": "statistics::ALT_TWO_SIDED", "less": "statistics::ALT_LESS",
       "greater": "statistics::ALT_GREATER"}
CEN = {"mean": "statistics::CENTER_MEAN", "median": "statistics::CENTER_MEDIAN",
       "trimmed": "statistics::CENTER_TRIMMED"}
TIES = {"below": "statistics::MEDIAN_BELOW", "above": "statistics::MEDIAN_ABOVE",
        "ignore": "statistics::MEDIAN_IGNORE"}
METH = {"auto": "statistics::MW_AUTO", "exact": "statistics::MW_EXACT",
        "asymptotic": "statistics::MW_ASYMPTOTIC"}

A = [2.1, 3.4, 1.9, 5.6, 4.4, 3.8, 2.7, 4.9]
B = [3.9, 5.2, 4.8, 6.1, 5.5, 4.1, 6.6, 5.0]
C = [1.2, 2.5, 2.2, 3.1, 1.8, 2.9, 2.4, 3.3]

rows = []
count = [0]


def num(v):
    v = float(v)
    if math.isnan(v):
        return "statistics::nan()"
    if math.isinf(v):
        return "vv_inf()" if v > 0 else "-vv_inf()"
    return repr(v)


def arr(v):
    return "[" + ", ".join(num(x) for x in v) + "]"


def arr_i(v):
    return "[" + ", ".join(str(int(x)) for x in v) + "]"


def flat(groups):
    return [x for g in groups for x in g]


def sizes(groups):
    return [len(g) for g in groups]


def tol_for(p):
    return "1e-12" if p > 1e-9 else "1e-9"


def block(lines):
    rows.append("    {\n" + "\n".join("        " + s for s in lines) + "\n    }")
    count[0] += 1


def groups_decl(groups):
    g, s = flat(groups), sizes(groups)
    return [f"let g: f64[{len(g)}] = {arr(g)};", f"let s: i64[{len(s)}] = {arr_i(s)};"]


def res_line(name, stat, p, df, df2, tol=None):
    tol = tol or tol_for(p)
    return (f'vv_res(bad, "{name}", r, c, {num(stat)}, {num(p)}, {num(df)}, {num(df2)}, {tol});')


def rnd_sample(n, mu, sd, ties):
    out = []
    for _ in range(n):
        v = random.gauss(mu, sd)
        out.append(round(v, 1) if ties else v)
    return out


def groups_case(k, lo, hi, ties, scales=None):
    gs = []
    for i in range(k):
        sd = (scales or [1.0] * k)[i]
        gs.append(rnd_sample(random.randint(lo, hi), random.uniform(-1, 1), sd, ties))
    return gs


# ---- levene, fligner, bartlett, alexandergovern ----------------------------

skipped = {"fligner": 0}


def fligner_exact_center(gs, center, trim):
    """Fligner with centers rounded once from the exact sum, as the package
    does. When the result differs from SciPy, SciPy sits on a tie of the
    deviations that its own rounded center breaks the other way."""
    from scipy.special import ndtri
    cs = []
    for g in gs:
        sg = sorted(g)
        if center == "mean":
            cs.append(math.fsum(g) / len(g))
        elif center == "median":
            cs.append(float(np.median(g)))
        else:
            cut = int(trim * len(g))
            kept = sg[cut:len(sg) - cut]
            cs.append(math.fsum(kept) / len(kept))
    z = [abs(x - c) for g, c in zip(gs, cs) for x in g]
    n = len(z)
    a = ndtri(stats.rankdata(z) / (2 * (n + 1.0)) + 0.5)
    i, tot = 0, 0.0
    for g in gs:
        ai = a[i:i + len(g)]
        i += len(g)
        tot += len(g) * (ai.mean() - a.mean()) ** 2
    return tot / np.var(a, ddof=1)

def add_groups_tests(tag, gs):
    for center in ["mean", "median", "trimmed"]:
        for trim in ([0.05, 0.25] if center == "trimmed" else [0.05]):
            try:
                r = stats.levene(*gs, center=center, proportiontocut=trim)
                lines = groups_decl(gs) + [
                    f"let r, c = statistics::levene(g[..], s[..], {CEN[center]}, {trim});",
                    res_line(f"levene {center} {trim} {tag}", r.statistic, r.pvalue,
                             len(gs) - 1, len(flat(gs)) - len(gs))]
                block(lines)
            except Exception:
                pass
            try:
                r = stats.fligner(*gs, center=center, proportiontocut=trim)
                alt_stat = fligner_exact_center(gs, center, trim)
                if abs(alt_stat - r.statistic) > 1e-12 * abs(r.statistic):
                    skipped["fligner"] += 1
                    continue
                lines = groups_decl(gs) + [
                    f"let r, c = statistics::fligner(g[..], s[..], {CEN[center]}, {trim});",
                    res_line(f"fligner {center} {trim} {tag}", r.statistic, r.pvalue,
                             len(gs) - 1, float("nan"))]
                block(lines)
            except Exception:
                pass
    for name, fn in [("bartlett", stats.bartlett), ("alexandergovern", stats.alexandergovern)]:
        try:
            r = fn(*gs)
            if not (math.isfinite(r.statistic) and math.isfinite(r.pvalue)):
                continue
            lines = groups_decl(gs) + [
                f"let r, c = statistics::{name}(g[..], s[..]);",
                res_line(f"{name} {tag}", r.statistic, r.pvalue, len(gs) - 1, float("nan"))]
            block(lines)
        except Exception:
            pass


def add_median_test(tag, gs):
    for ties in ["below", "above", "ignore"]:
        for corr in [True, False]:
            try:
                r = stats.median_test(*gs, ties=ties, correction=corr)
            except Exception:
                continue
            k = len(gs)
            table = [int(x) for x in r.table.flatten()]
            lines = groups_decl(gs) + [
                "var out: f64[>..] = [];",
                f"let r, c = statistics::median_test(out, g[..], s[..], {TIES[ties]}, {'true' if corr else 'false'});",
                res_line(f"median_test {ties} {corr} {tag}", r.statistic, r.pvalue, k - 1, r.median),
                f"let t: f64[{len(table)}] = {arr(table)};",
                f'vv_list(bad, "median_test table {ties} {corr} {tag}", out[..], t[..]);']
            block(lines)


# ---- two-sample tests -------------------------------------------------------

def exact_ansari(nx, ny, ab, alt):
    """The exact distribution by integer counts. SciPy keeps the counts in
    single precision, so its exact p-value drifts for samples near 54."""
    from fractions import Fraction
    n = nx + ny
    vals = [min(r, n + 1 - r) for r in range(1, n + 1)]
    dp = [dict() for _ in range(nx + 1)]
    dp[0][0] = 1
    for v in vals:
        for j in range(nx, 0, -1):
            for s_, c_ in list(dp[j - 1].items()):
                dp[j][s_ + v] = dp[j].get(s_ + v, 0) + c_
    tot = sum(dp[nx].values())
    le = Fraction(sum(c_ for s_, c_ in dp[nx].items() if s_ <= ab), tot)
    ge = Fraction(sum(c_ for s_, c_ in dp[nx].items() if s_ >= ab), tot)
    if alt == "greater":
        return float(le)
    if alt == "less":
        return float(ge)
    return float(min(Fraction(1), 2 * min(le, ge)))


def add_two_sample(tag, x, y):
    decl = [f"let x: f64[{len(x)}] = {arr(x)};", f"let y: f64[{len(y)}] = {arr(y)};"]
    for alt in ALT:
        try:
            r = stats.ranksums(x, y, alternative=alt)
            block(decl + [f"let r, c = statistics::ranksums(x[..], y[..], {ALT[alt]});",
                          res_line(f"ranksums {alt} {tag}", r.statistic, r.pvalue,
                                   float("nan"), float("nan"))])
        except Exception:
            pass
        try:
            r = stats.mood(x, y, alternative=alt)
            if math.isfinite(r.statistic):
                block(decl + [f"let r, c = statistics::mood(x[..], y[..], {ALT[alt]});",
                              res_line(f"mood {alt} {tag}", r.statistic, r.pvalue,
                                       float("nan"), float("nan"))])
        except Exception:
            pass
        for dist in ["t", "normal"]:
            try:
                r = stats.brunnermunzel(x, y, alternative=alt, distribution=dist)
                if not (math.isfinite(r.statistic) and math.isfinite(r.pvalue)):
                    continue
                df = float("nan")
                if dist == "t":
                    nx, ny = len(x), len(y)
                    rc = stats.rankdata(list(x) + list(y))
                    rx, ry = stats.rankdata(x), stats.rankdata(y)
                    tx = rc[:nx] - rx - rc[:nx].mean() + rx.mean()
                    ty = rc[nx:] - ry - rc[nx:].mean() + ry.mean()
                    sx, sy = (tx @ tx) / (nx - 1), (ty @ ty) / (ny - 1)
                    df = (nx * sx + ny * sy) ** 2 / ((nx * sx) ** 2 / (nx - 1) + (ny * sy) ** 2 / (ny - 1))
                block(decl + [f"let r, c = statistics::brunnermunzel(x[..], y[..], {ALT[alt]}, {'true' if dist == 't' else 'false'});",
                              res_line(f"brunnermunzel {alt} {dist} {tag}", r.statistic, r.pvalue,
                                       df, float("nan"))])
            except Exception:
                pass
        for meth in ["auto", "exact", "asymptotic"]:
            ties = len(set(list(x) + list(y))) < len(x) + len(y)
            if meth == "exact" and ties:
                continue
            try:
                r = stats.ansari(x, y, alternative=alt, method=meth)
                uses_exact = meth == "exact" or (meth == "auto" and not ties and len(x) < 55 and len(y) < 55)
                if uses_exact:
                    pe = exact_ansari(len(x), len(y), int(r.statistic), alt)
                    assert abs(pe - r.pvalue) <= 1e-5 * max(pe, 1e-300), (tag, pe, r.pvalue)
                    r = type(r)(r.statistic, pe)
                block(decl + [f"let r, c = statistics::ansari(x[..], y[..], {ALT[alt]}, {METH[meth]});",
                              res_line(f"ansari {alt} {meth} {tag}", r.statistic, r.pvalue,
                                       float("nan"), float("nan"), "1e-11")])
            except Exception:
                pass


def add_friedman(tag, cols):
    try:
        r = stats.friedmanchisquare(*cols)
    except Exception:
        return
    if not math.isfinite(r.statistic):
        return
    k, n = len(cols), len(cols[0])
    data = [cols[j][i] for i in range(n) for j in range(k)]
    block([f"let d: f64[{len(data)}] = {arr(data)};",
           f"let r, c = statistics::friedmanchisquare(d[..], {k});",
           res_line(f"friedmanchisquare {tag}", r.statistic, r.pvalue, k - 1, float("nan"))])


def add_quantile(tag, x):
    srt = sorted(x)
    for alt in ALT:
        for level in [0.9, 0.95, 0.99]:
            q = random.choice([srt[random.randrange(len(srt))], random.uniform(srt[0] - 0.5, srt[-1] + 0.5)])
            p = random.choice([0.1, 0.25, 0.5, 0.75, 0.9, random.uniform(0.05, 0.95)])
            r = stats.quantile_test(x, q=q, p=p, alternative=alt)
            ci = r.confidence_interval(level)
            lo = -math.inf if math.isnan(ci.low) else ci.low
            hi = math.inf if math.isnan(ci.high) else ci.high
            block([f"let x: f64[{len(x)}] = {arr(x)};", "var out: f64[>..] = [];",
                   f"let r, c = statistics::quantile_test(out, x[..], {num(q)}, {num(p)}, {ALT[alt]}, {level});",
                   res_line(f"quantile_test {alt} {level} {tag}", r.statistic, r.pvalue,
                            r.statistic_type, float("nan"), "1e-12"),
                   f"let w: f64[2] = {arr([lo, hi])};",
                   f'vv_list(bad, "quantile_test interval {alt} {level} {tag}", out[..], w[..]);'])


def add_from_stats(tag, m1, s1, n1, m2, s2, n2):
    for eq in [True, False]:
        for alt in ALT:
            r = stats.ttest_ind_from_stats(m1, s1, n1, m2, s2, n2, equal_var=eq, alternative=alt)
            if eq:
                df = n1 + n2 - 2
            else:
                v1, v2 = s1 ** 2 / n1, s2 ** 2 / n2
                df = (v1 + v2) ** 2 / (v1 ** 2 / (n1 - 1) + v2 ** 2 / (n2 - 1))
            block([f"let r, c = statistics::ttest_ind_from_stats({num(m1)}, {num(s1)}, {num(n1)}, {num(m2)}, {num(s2)}, {num(n2)}, {'true' if eq else 'false'}, {ALT[alt]});",
                   res_line(f"ttest_ind_from_stats {eq} {alt} {tag}", r.statistic, r.pvalue,
                            df, float("nan"))])


def main():
    # The measured oracle: the samples a, b and c.
    add_groups_tests("abc", [A, B, C])
    add_median_test("abc", [A, B, C])
    add_two_sample("ab", A, B)
    add_friedman("abc", [A, B, C])
    add_quantile("a", A)
    add_from_stats("oracle", 3.0, 1.2, 8, 4.5, 1.4, 9)
    # The tied example of the research notes.
    tx, ty = [1, 2, 2, 3, 4], [2, 3, 5, 6, 6, 7]
    add_two_sample("tied", tx, ty)
    add_friedman("tied", [[1, 2, 2, 4], [2, 2, 3, 5], [3, 1, 1, 2]])
    add_groups_tests("g3", [[1, 2, 3, 4, 5], [3, 3, 6, 7, 8], [3, 9, 10, 11, 12]])
    add_median_test("g3", [[1, 2, 3, 4, 5], [3, 3, 6, 7, 8], [3, 9, 10, 11, 12]])
    # Seeded random samples.
    for i, (lo, hi) in enumerate([(3, 6), (6, 12), (12, 30), (30, 80)]):
        for ties in [False, True]:
            for k in [2, 3, 5]:
                tag = f"r{i}{'t' if ties else 'c'}{k}"
                scales = [random.choice([0.5, 1.0, 2.0, 4.0]) for _ in range(k)]
                gs = groups_case(k, lo, hi, ties, scales)
                add_groups_tests(tag, gs)
                add_median_test(tag, gs)
            x = rnd_sample(random.randint(lo, hi), 0.0, 1.0, ties)
            y = rnd_sample(random.randint(lo, hi), random.uniform(-0.5, 0.5), random.choice([0.7, 1.0, 2.0]), ties)
            add_two_sample(f"s{i}{'t' if ties else 'c'}", x, y)
            add_quantile(f"s{i}{'t' if ties else 'c'}", x)
            for k in [3, 4, 6]:
                n = random.randint(max(lo, 3), hi)
                add_friedman(f"f{i}{'t' if ties else 'c'}{k}",
                             [rnd_sample(n, j * 0.3, 1.0, ties) for j in range(k)])
    # The sizes at the edge of the exact Ansari-Bradley form.
    for (m, n) in [(54, 54), (54, 20), (20, 54), (55, 30), (2, 2), (1, 5)]:
        add_two_sample(f"edge{m}x{n}", rnd_sample(m, 0.0, 1.0, False), rnd_sample(n, 0.0, 2.0, False))
    for i in range(6):
        m1, m2 = random.uniform(-3, 3), random.uniform(-3, 3)
        s1, s2 = random.uniform(0.1, 5), random.uniform(0.1, 5)
        add_from_stats(f"r{i}", m1, s1, random.randint(2, 200), m2, s2, random.randint(2, 200))

    head = '''// Generated by test/make_variance_vectors.py. Do not edit by hand.
// The measured rows of the variance calls: the answers of SciPy 1.18.1.
// The generator test/make_variance_vectors.py writes this file.
// It holds the numbers and the checks that use them.
namespace suite;

import std;
import math;
import statistics from "../..";

fn vv_inf() -> f64 { return 1.0 / 0.0; }

// One result against the row of SciPy. A want of NaN skips that field.
fn vv_res(bad: i64&, name: u8[:], r: statistics::TestResult, c: i64, stat: f64, p: f64, df: f64, df2: f64, tol: f64) {
    var ok = c == 0 && near_tol(r.statistic, stat, tol) && near_tol(r.pvalue, p, tol);
    if !is_nan(df) { ok = ok && near_tol(r.df, df, 1e-12); }
    if !is_nan(df2) { ok = ok && near_tol(r.df2, df2, 1e-12); }
    if ok { return; }
    bad += 1;
    if bad <= 12 { print("FAIL ", name, " got ", r.statistic, " ", r.pvalue, " ", r.df, " ", r.df2, " code ", c, " want ", stat, " ", p, " ", df, " ", df2, "\\n"); }
}

fn vv_list(bad: i64&, name: u8[:], got: const f64[:], want: const f64[:]) {
    if near_list(got, want, 1e-12) { return; }
    bad += 1;
    if bad <= 12 { print("FAIL ", name, " ", got, " want ", want, "\\n"); }
}

fn run_variance_vectors() -> i64 {
    var bad: i64 = 0;
'''
    with open(OUT, "w") as f:
        f.write(head)
        f.write("\n".join(rows))
        f.write("\n    return bad;\n}\n")
    print("rows:", count[0], "fligner rows skipped on a tie of the center:", skipped["fligner"])


main()
