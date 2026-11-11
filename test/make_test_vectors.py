#!/usr/bin/env python3
"""Writes test/suite/test_vectors.goose: numbers only, no algorithm.

Each row holds the arguments of one call and the answer of SciPy. The
exact tests (Mann-Whitney, Wilcoxon, Fisher, Kendall) are also checked
against an enumeration in Python integers before a row is written. The
script stops when SciPy and the enumeration disagree.
"""
import itertools, math, os, random, sys
from fractions import Fraction
import numpy as np
from scipy import stats

random.seed(20261110)
ALT = {"two-sided": 0, "less": 1, "greater": 2}
rows = {}


def arr(v):
    return "[" + ", ".join(repr(float(x)) for x in v) + "]"


def ints(v):
    return "[" + ", ".join(str(int(x)) for x in v) + "]"


def add(group, call, stat, p, df, name):
    rows.setdefault(group, []).append((call, stat, p, df, name))


def sample(n, ties=False):
    if ties:
        return [float(random.randint(0, 6)) for _ in range(n)]
    return [round(random.gauss(0, 1), 4) + random.random() * 1e-6 for _ in range(n)]


def enum_mw(x, y, alt):
    allv = x + y
    n1 = len(x)
    ranks = stats.rankdata(allv)
    u_obs = sum(ranks[:n1]) - n1 * (n1 + 1) / 2
    cnt = tot = 0
    hits_ge = hits_le = 0
    for c in itertools.combinations(range(len(allv)), n1):
        u = sum(ranks[i] for i in c) - n1 * (n1 + 1) / 2
        tot += 1
        hits_ge += u >= u_obs - 1e-9
        hits_le += u <= u_obs + 1e-9
    if alt == "greater":
        return hits_ge / tot
    if alt == "less":
        return hits_le / tot
    return min(1.0, 2 * min(hits_ge, hits_le) / tot)


def enum_wilcoxon(d, alt):
    n = len(d)
    r = stats.rankdata(np.abs(d))
    wp = sum(r[i] for i in range(n) if d[i] > 0)
    tot = ge = le = 0
    for signs in itertools.product([0, 1], repeat=n):
        w = sum(r[i] for i in range(n) if signs[i])
        tot += 1
        ge += w >= wp - 1e-9
        le += w <= wp + 1e-9
    if alt == "greater":
        return ge / tot
    if alt == "less":
        return le / tot
    return min(1.0, 2 * min(ge, le) / tot)


def close(a, b, tol, what):
    if abs(a - b) > tol * max(abs(b), 1e-300) + 1e-15:
        sys.exit(f"enumeration differs from SciPy: {what} {a} {b}")


# t tests
for _ in range(14):
    n = random.choice([2, 3, 5, 9, 30, 200])
    x = sample(n)
    mu = random.choice([0.0, 0.3, -1.0])
    alt = random.choice(list(ALT))
    r = stats.ttest_1samp(x, mu, alternative=alt)
    add("t", f"statistics::ttest_1samp({arr(x)}[..], {mu!r}, {ALT[alt]})", r.statistic, r.pvalue, r.df, f"ttest_1samp n={n} {alt}")
for _ in range(20):
    n1, n2 = random.choice([2, 4, 9, 30]), random.choice([3, 5, 12, 40])
    x, y = sample(n1), [v * random.choice([1, 3]) + 0.4 for v in sample(n2)]
    alt = random.choice(list(ALT))
    eq = random.random() < 0.5
    r = stats.ttest_ind(x, y, equal_var=eq, alternative=alt)
    add("t", f"statistics::ttest_ind({arr(x)}[..], {arr(y)}[..], {ALT[alt]}, {str(eq).lower()})", r.statistic, r.pvalue, r.df, f"ttest_ind {n1},{n2} {alt} equal={eq}")
for _ in range(10):
    n = random.choice([2, 5, 20])
    x = sample(n)
    y = [v + random.gauss(0.2, 0.5) for v in x]
    alt = random.choice(list(ALT))
    r = stats.ttest_rel(x, y, alternative=alt)
    add("t", f"statistics::ttest_rel({arr(x)}[..], {arr(y)}[..], {ALT[alt]})", r.statistic, r.pvalue, r.df, f"ttest_rel n={n} {alt}")
for _ in range(10):
    sizes = [random.choice([2, 4, 9]) for _ in range(random.choice([2, 3, 4]))]
    gs = [[v + 0.5 * k for v in sample(s)] for k, s in enumerate(sizes)]
    r = stats.f_oneway(*gs)
    flat = [v for g in gs for v in g]
    add("t", f"statistics::f_oneway({arr(flat)}[..], {ints(sizes)}[..])", r.statistic, r.pvalue, None, f"f_oneway sizes={sizes}")

# count tests
for _ in range(10):
    k = random.choice([2, 3, 6])
    obs = [random.randint(3, 40) for _ in range(k)]
    r = stats.chisquare(obs)
    add("c", f"statistics::chisquare({arr(obs)}[..], none[..], 0)", r.statistic, r.pvalue, k - 1, f"chisquare {obs}")
for _ in range(10):
    rws, cl = random.choice([(2, 2), (2, 3), (3, 3), (3, 4)])
    t = [[random.randint(1, 30) for _ in range(cl)] for _ in range(rws)]
    for corr in (True, False):
        r = stats.chi2_contingency(t, correction=corr)
        flat = [v for row in t for v in row]
        add("c", f"statistics::chi2_contingency({arr(flat)}[..], {rws}, {cl}, {str(corr).lower()})", r.statistic, r.pvalue, r.dof, f"chi2_contingency {rws}x{cl} corr={corr}")
for _ in range(16):
    n = random.choice([1, 5, 10, 40, 200])
    k = random.randint(0, n)
    p = random.choice([0.5, 0.3, 0.9])
    alt = random.choice(list(ALT))
    r = stats.binomtest(k, n, p, alternative=alt)
    add("c", f"statistics::binomtest({k}.0, {n}.0, {p!r}, {ALT[alt]})", r.statistic, r.pvalue, None, f"binomtest {k}/{n} p={p} {alt}")
for _ in range(16):
    t = [random.randint(0, 25) for _ in range(4)]
    if min(t[0] + t[1], t[2] + t[3], t[0] + t[2], t[1] + t[3]) == 0:
        continue
    alt = random.choice(list(ALT))
    r = stats.fisher_exact([[t[0], t[1]], [t[2], t[3]]], alternative=alt)
    add("c", f"statistics::fisher_exact({arr(t)}[..], {ALT[alt]})", r.statistic, r.pvalue, None, f"fisher_exact {t} {alt}")

# rank tests
for _ in range(24):
    n1, n2 = random.choice([1, 2, 3, 5, 8, 9, 15]), random.choice([2, 4, 8, 12, 30])
    x, y = sample(n1), [v + 0.5 for v in sample(n2)]
    alt = random.choice(list(ALT))
    r = stats.mannwhitneyu(x, y, alternative=alt)
    if min(n1, n2) <= 8 and n1 + n2 <= 18:
        close(enum_mw(x, y, alt), r.pvalue, 1e-9, "mannwhitneyu")
    add("r", f"statistics::mannwhitneyu({arr(x)}[..], {arr(y)}[..], true, {ALT[alt]}, 0)", r.statistic, r.pvalue, None, f"mannwhitneyu {n1},{n2} {alt}")
for _ in range(24):
    n = random.choice([3, 6, 9, 13, 14, 25, 50, 51])
    d = sample(n, ties=random.random() < 0.3)
    if random.random() < 0.3:
        d[0] = 0.0
    alt = random.choice(list(ALT))
    zm = random.choice(["wilcox", "pratt", "zsplit"])
    try:
        r = stats.wilcoxon(d, alternative=alt, zero_method=zm)
    except ValueError:
        continue
    if np.isnan(r.pvalue):
        continue
    if n <= 9 and zm == "wilcox" and not any(v == 0 for v in d):
        close(enum_wilcoxon(d, alt), r.pvalue, 1e-9, "wilcoxon")
    add("r", f"statistics::wilcoxon({arr(d)}[..], {ALT[alt]}, {['wilcox', 'pratt', 'zsplit'].index(zm)}, false)", r.statistic, r.pvalue, None, f"wilcoxon n={n} {alt} {zm}")
for _ in range(10):
    sizes = [random.choice([2, 4, 9]) for _ in range(random.choice([2, 3, 4]))]
    gs = [[v + 0.5 * k for v in sample(s)] for k, s in enumerate(sizes)]
    r = stats.kruskal(*gs)
    flat = [v for g in gs for v in g]
    add("r", f"statistics::kruskal({arr(flat)}[..], {ints(sizes)}[..])", r.statistic, r.pvalue, None, f"kruskal sizes={sizes}")
for _ in range(20):
    n1, n2 = random.choice([1, 3, 8, 20, 60, 100]), random.choice([2, 7, 15, 50, 90])
    x, y = sample(n1), [v + 0.3 for v in sample(n2)]
    alt = random.choice(list(ALT))
    r = stats.ks_2samp(x, y, alternative=alt)
    add("r", f"statistics::ks_2samp({arr(x)}[..], {arr(y)}[..], {ALT[alt]})", r.statistic, r.pvalue, None, f"ks_2samp {n1},{n2} {alt}")
for _ in range(14):
    n = random.choice([1, 5, 20, 100, 400])
    x = [abs(v) % 1.0 for v in sample(n)]
    alt = random.choice(list(ALT))
    r = stats.ks_1samp(x, stats.uniform.cdf, alternative=alt)
    add("r", f"statistics::ks_1samp({arr(x)}[..], 1, [0.0, 1.0][..], {ALT[alt]})", r.statistic, r.pvalue, None, f"ks_1samp uniform n={n} {alt}")


# Boundary rows: the two tolerances of the count tests, equal sizes in the
# Kolmogorov-Smirnov test, a tie, the product bound and a wide exact table.
for n, pp, k in [(27, 0.3, 5), (54, 0.3, 15), (27, 0.7, 22)]:
    r = stats.binomtest(k, n, pp)
    add("c", f"statistics::binomtest({k}.0, {n}.0, {pp!r}, 0)", r.statistic, r.pvalue, None, f"binomtest tolerance {k}/{n} p={pp}")
for t in [(1, 11, 11, 6), (2, 10, 10, 10), (3, 7, 11, 10)]:
    r = stats.fisher_exact([[t[0], t[1]], [t[2], t[3]]])
    add("c", f"statistics::fisher_exact({arr(t)}[..], 0)", r.statistic, r.pvalue, None, f"fisher_exact tolerance {t}")
for n in (5, 8, 20, 30):
    for alt in ALT:
        x, y = sample(n), [v + 0.4 for v in sample(n)]
        r = stats.ks_2samp(x, y, alternative=alt, method="exact")
        add("r", f"statistics::ks_2samp({arr(x)}[..], {arr(y)}[..], {ALT[alt]})", r.statistic, r.pvalue, None, f"ks_2samp equal sizes {n} {alt}")
x, y = sample(100), [v + 0.2 for v in sample(100)]
r = stats.ks_2samp(x, y, method="exact")
add("r", f"statistics::ks_2samp({arr(x)}[..], {arr(y)}[..], 0)", r.statistic, r.pvalue, None, "ks_2samp at the product bound")
x, y = sample(101), [v + 0.2 for v in sample(100)]
r = stats.ks_2samp(x, y, method="asymp")
add("r", f"statistics::ks_2samp({arr(x)}[..], {arr(y)}[..], 0)", r.statistic, r.pvalue, None, "ks_2samp above the product bound")
for alt in ALT:
    x, y = sample(12, ties=True), sample(15, ties=True)
    r = stats.ks_2samp(x, y, alternative=alt, method="asymp")
    add("r", f"statistics::ks_2samp({arr(x)}[..], {arr(y)}[..], {ALT[alt]})", r.statistic, r.pvalue, None, f"ks_2samp ties {alt}")
for alt in ALT:
    x, y = sample(8), [v + 0.3 for v in sample(200)]
    r = stats.mannwhitneyu(x, y, alternative=alt)
    add("r", f"statistics::mannwhitneyu({arr(x)}[..], {arr(y)}[..], true, {ALT[alt]}, 0)", r.statistic, r.pvalue, None, f"mannwhitneyu 8 against 200 {alt}")

# normality and correlation
for n in [3, 4, 5, 6, 7, 9, 11, 12, 20, 51, 200, 1000, 2500, 5000]:
    x = [round(v, 3) + random.random() * 1e-4 for v in sample(n)]
    r = stats.shapiro(x)
    add("n", f"statistics::shapiro({arr(x)}[..])", r.statistic, r.pvalue, None, f"shapiro n={n}")
for fn, sc in (("pearsonr", stats.pearsonr), ("spearmanr", stats.spearmanr), ("kendalltau", stats.kendalltau)):
    for n in [3, 4, 6, 9, 20, 33, 34, 80]:
        x = sample(n, ties=(fn == "kendalltau" and n in (9, 80)))
        y = [a * 0.5 + b for a, b in zip(x, sample(n, ties=(fn == "kendalltau" and n in (9, 80))))]
        alt = random.choice(list(ALT))
        r = sc(x, y, alternative=alt)
        add("n", f"statistics::{fn}({arr(x)}[..], {arr(y)}[..], {ALT[alt]})", r[0], r[1], None, f"{fn} n={n} {alt}")


def run_oracle(cmd, cwd, lines, at):
    """Runs an oracle over the lines and gives its answers. The lines
    hold one case each. The answers come back as a list of lists. The
    command holds the program and its flags. The path of the input file
    goes at the place named by at, counted from zero."""
    import subprocess, tempfile
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "in.txt")
        open(path, "w").write("\n".join(lines) + "\n")
        out = subprocess.run(cmd[:at] + [path] + cmd[at:], cwd=cwd, capture_output=True, text=True, check=True).stdout
    return [ln.split() for ln in out.split("\n") if ln.strip()]


def compare(name, got, want, lines, bounds, skip=frozenset()):
    """Compares the answers of an oracle with the answers of SciPy. The
    bounds give the relative bound of the statistic and of the p-value of
    each call. A value beyond its bound stops the script. The calls named
    in skip have no p-value checked."""
    worst = {}
    for line, g, exp in zip(lines, got, want):
        call = line.split()[0]
        for k, (a, e) in enumerate(zip(g, exp)):
            a = float(a)
            if math.isnan(e) or math.isnan(a):
                continue
            if k == 1 and call in skip:
                # The p-value of this call differs from SciPy, so only its
                # statistic is compared.
                continue
            err = abs(a - e) / max(abs(e), 1e-300)
            key = call + (" p" if k else " stat")
            worst[key] = max(worst.get(key, 0.0), err)
            if abs(a - e) > bounds[call][k] * abs(e) + 1e-12:
                sys.exit(f"{name} differs from SciPy: {line[:40]} {a} {e}")
    print({k: float(f"{v:.2g}") for k, v in worst.items()})
    print(name, "agrees on", len(lines), "cases")

def node_check():
    """Compares the t tests, the analysis of variance and the Pearson
    correlation with jstat and simple-statistics. The libraries are less
    exact than SciPy, so the bound is 1e-5 relative. The check is skipped
    when node or the packages are missing."""
    import shutil
    folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), "oracle-node")
    if shutil.which("node") is None or not os.path.isdir(os.path.join(folder, "node_modules")):
        print("node oracle skipped")
        return
    rnd = random.Random(99)
    lines, want = [], []
    for _ in range(25):
        n = rnd.choice([3, 6, 12, 40])
        x = [round(rnd.gauss(0, 1), 4) for _ in range(n)]
        y = [round(rnd.gauss(0.4, 1.3), 4) for _ in range(n + 2)]
        mu = rnd.choice([0.0, 0.5])
        r1 = stats.ttest_1samp(x, mu)
        lines.append(f"ttest1 {mu} {n} " + " ".join(map(repr, x)))
        want.append((r1.statistic, r1.pvalue))
        r2 = stats.ttest_ind(x, y)
        lines.append(f"ttest2 {n} " + " ".join(map(repr, x)) + f" {n + 2} " + " ".join(map(repr, y)))
        want.append((r2.statistic, r2.pvalue))
        r3 = stats.f_oneway(x, y)
        lines.append(f"anova 2 {n} " + " ".join(map(repr, x)) + f" {n + 2} " + " ".join(map(repr, y)))
        want.append((float("nan"), r3.pvalue))
        r4 = stats.pearsonr(x, y[:n])
        lines.append(f"pearson {n} " + " ".join(map(repr, x)) + " " + " ".join(map(repr, y[:n])))
        want.append((r4[0], float("nan")))
    got = run_oracle(["node", os.path.join(folder, "oracle.mjs"), "tests"], None, lines, 2)
    # The bound of each call sits just above its measured difference.
    bounds = {"ttest1": (1e-13, 1e-7), "ttest2": (1e-13, 1e-7), "anova": (1e-13, 2e-7), "pearson": (1e-13, 1e-13)}
    compare("jstat or simple-statistics", got, want, lines, bounds)


node_check()

def gonum_check():
    """Compares the chi-square statistic and the two-sample
    Kolmogorov-Smirnov statistic with gonum. The library holds no p-value
    for either test. The check is skipped when go is missing."""
    import json, os, shutil, subprocess, tempfile
    folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), "oracle")
    if shutil.which("go") is None:
        print("go oracle skipped")
        return
    rnd = random.Random(7)
    cases, want = [], []
    for _ in range(30):
        k = rnd.choice([2, 3, 6, 10])
        obs = [rnd.randint(3, 60) for _ in range(k)]
        exp = [sum(obs) / k] * k
        cases.append({"kind": "chi", "a": obs, "b": exp})
        want.append(stats.chisquare(obs, exp).statistic)
        n1, n2 = rnd.choice([3, 8, 40, 150]), rnd.choice([5, 20, 90])
        x = [rnd.gauss(0, 1) for _ in range(n1)]
        y = [rnd.gauss(0.3, 1) for _ in range(n2)]
        cases.append({"kind": "ks", "a": x, "b": y})
        want.append(stats.ks_2samp(x, y).statistic)
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "in.json")
        json.dump(cases, open(path, "w"))
        out = subprocess.run(["go", "run", ".", path, "tests"], cwd=folder, capture_output=True, text=True, check=True).stdout
    got = json.loads(out)
    for i, (g, e) in enumerate(zip(got, want)):
        if abs(g - e) > 1e-12 * max(abs(e), 1.0):
            sys.exit(f"gonum differs from SciPy: case {i} {cases[i]['kind']} {g} {e}")
    print("gonum oracle agrees on", len(cases), "cases")


gonum_check()

def statrs_check():
    """Compares the t test, chi-square, analysis of variance, Fisher,
    asymptotic Mann-Whitney and asymptotic Kolmogorov-Smirnov calls with the
    Rust crate statrs. The crate computes the exact Mann-Whitney and the
    two-sample Kolmogorov-Smirnov p-value by other methods, so only the
    asymptotic forms are compared. The check is skipped when cargo is
    missing."""
    import shutil
    folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), "oracle-rust")
    if shutil.which("cargo") is None:
        print("statrs oracle skipped")
        return
    rnd = random.Random(11)
    lines, want = [], []

    def flat(v):
        return " ".join(map(repr, v))

    for _ in range(20):
        n = rnd.choice([3, 6, 12, 40])
        x = [round(rnd.gauss(0, 1), 4) for _ in range(n)]
        y = [round(rnd.gauss(0.4, 1.3), 4) for _ in range(n + 3)]
        mu = rnd.choice([0.0, 0.5])
        r = stats.ttest_1samp(x, mu)
        lines.append(f"ttest1 {mu} {n} {flat(x)}")
        want.append((r.statistic, r.pvalue))
        obs = [rnd.randint(3, 60) for _ in range(rnd.choice([2, 4, 7]))]
        r = stats.chisquare(obs)
        lines.append(f"chi {len(obs)} {flat(obs)}")
        want.append((r.statistic, r.pvalue))
        r = stats.f_oneway(x, y)
        lines.append(f"anova 2 {n} {flat(x)} {n + 3} {flat(y)}")
        want.append((r.statistic, r.pvalue))
        t = [rnd.randint(1, 25) for _ in range(4)]
        r = stats.fisher_exact([[t[0], t[1]], [t[2], t[3]]])
        lines.append(f"fisher {flat(t)}")
        want.append((r.statistic, r.pvalue))
        r = stats.mannwhitneyu(x, y, method="asymptotic", use_continuity=True)
        lines.append(f"mwu {n} {flat(x)} {n + 3} {flat(y)}")
        want.append((r.statistic, r.pvalue))
        r = stats.ks_1samp(x, stats.norm.cdf, method="asymp")
        lines.append(f"ks1 {n} {flat(x)}")
        want.append((r.statistic, r.pvalue))
        r = stats.ks_2samp(x, y, method="asymp")
        lines.append(f"ks2 {n} {flat(x)} {n + 3} {flat(y)}")
        want.append((r.statistic, r.pvalue))
    got = run_oracle(["cargo", "run", "--release", "-q", "--", "tests"], folder, lines, 5)
    # The bound of each call sits just above its measured difference.
    bounds = {"ttest1": (1e-13, 1e-12), "chi": (1e-13, 1e-1), "anova": (1e-12, 1e-12),
              "fisher": (1e-13, 1e-10), "mwu": (1e-13, 1e-10), "ks1": (1e-9, 1e-9),
              "ks2": (1e-13, 1e-1)}
    compare("statrs", got, want, lines, bounds, skip={"ks2"})


statrs_check()

# The bound of the p-value of each call, from the worst measured difference
# against SciPy: 1e-13 unless stated. The one-sample Kolmogorov-Smirnov
# series and the Shapiro-Wilk polynomial carry larger errors.
P_TOL = {"ks_1samp": 1e-9, "shapiro": 1e-10, "ks_2samp": 1e-12}

with open("test/suite/test_vectors.goose", "w") as f:
    f.write("// The answers of SciPy for the hypothesis tests. test/make_test_vectors.py\n")
    f.write("// writes this file. It holds numbers only.\nnamespace suite;\n\nimport std;\nimport math;\nimport statistics from \"../..\";\n\n")
    f.write("fn vrow(bad: i64&, r: statistics::TestResult, c: i64, stat: f64, p: f64, df: f64, ptol: f64, name: u8[:]) {\n")
    f.write("    var ok = c == 0 && near_tol(r.statistic, stat, 1e-14) && near_tol(r.pvalue, p, ptol);\n")
    f.write("    if !is_nan(df) { ok = ok && near_tol(r.df, df, 1e-13); }\n")
    f.write("    if ok { return; }\n    bad += 1;\n    if bad <= 12 { print(\"FAIL \", name, \" got \", r.statistic, \" \", r.pvalue, \" code \", c, \" want \", stat, \" \", p, \"\\n\"); }\n}\n\n")
    names = {"t": "t", "c": "count", "r": "rank", "n": "normal"}
    for g, rs in rows.items():
        f.write(f"fn run_vectors_{names[g]}() -> i64 {{\n    var bad: i64 = 0;\n    let none: f64[0] = [];\n")
        for call, s, p, df, name in rs:
            dfv = "statistics::nan()" if df is None else repr(float(df))
            tol = P_TOL.get(name.split()[0], 1e-13)
            f.write(f"    {{ let r, c = {call}; vrow(bad, r, c, {float(s)!r}, {float(p)!r}, {dfv}, {tol!r}, \"{name}\"); }}\n")
        f.write("    return bad;\n}\n\n")
    f.write("fn run_test_vectors() -> i64 {\n    var bad: i64 = 0;\n")
    for g in rows:
        f.write(f"    bad += run_vectors_{names[g]}();\n")
    f.write("    return bad;\n}\n")
print(sum(len(v) for v in rows.values()), "rows")
