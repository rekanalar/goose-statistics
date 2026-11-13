#!/usr/bin/env python3
"""Make the reference vectors of the statistics package.

The script writes numbers only into test/suite/vectors.goose. The exact
value of each statistic comes from rational and high precision decimal
arithmetic. The script then compares NumPy, SciPy, Python statistics and
gonum with the exact value and stops on a difference beyond the stated
tolerance. Run it from the package directory:

    python3 test/make_vectors.py

It needs Python 3 with NumPy and SciPy, and Go with gonum in the module
cache for the gonum check. It uses the crate statrs of Rust and the
package simple-statistics of Node for two more checks. Use --no-go to
skip the gonum check, --no-rust to skip the statrs check and --no-node
to skip the simple-statistics check.
"""
import json, math, os, random, statistics as st, subprocess, sys
from decimal import Decimal, getcontext
from fractions import Fraction

import numpy as np
from scipy import stats

getcontext().prec = 80
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "suite", "vectors.goose")
METHODS = ["inverted_cdf", "averaged_inverted_cdf", "closest_observation",
           "interpolated_inverted_cdf", "hazen", "weibull", "linear",
           "median_unbiased", "normal_unbiased"]
QS = [0.0, 0.1, 0.25, 0.3, 0.5, 0.75, 0.9, 1.0]
EPS = 2.220446049250313e-16


def dsqrt(f):
    """Square root of a Fraction as a Decimal."""
    return (Decimal(f.numerator) / Decimal(f.denominator)).sqrt()


def dfloat(d):
    return float(d)


def frac_to_decimal(f):
    return Decimal(f.numerator) / Decimal(f.denominator)


def exact(xs):
    """The exact statistics of one sample, rounded once to a float."""
    n = len(xs)
    fx = [Fraction(x) for x in xs]
    s = sum(fx)
    m = s / n
    d = [x - m for x in fx]
    cm = {k: sum(v ** k for v in d) / n for k in (2, 3, 4)}
    r = {"sum": float(s), "mean": float(m)}
    r["pvar"] = float(cm[2])
    r["var"] = float(sum(v * v for v in d) / (n - 1)) if n > 1 else None
    r["pstdev"] = dfloat(dsqrt(cm[2]))
    r["stdev"] = dfloat(dsqrt(Fraction(sum(v * v for v in d)) / (n - 1))) if n > 1 else None
    r["sem"] = dfloat(dsqrt(Fraction(sum(v * v for v in d)) / (n - 1)) / Decimal(n).sqrt()) if n > 1 else None
    if m != 0 and n > 1:
        r["variation"] = dfloat(dsqrt(Fraction(sum(v * v for v in d)) / (n - 1)) / frac_to_decimal(m))
    else:
        r["variation"] = None
    if cm[2] != 0:
        m2 = frac_to_decimal(cm[2])
        g1 = frac_to_decimal(cm[3]) / (m2 * m2.sqrt())
        g2 = frac_to_decimal(cm[4] / (cm[2] ** 2)) - 3
        r["skew_b"] = dfloat(g1)
        r["kurt_f_b"] = dfloat(g2)
        r["kurt_p_b"] = dfloat(g2 + 3)
        if n > 2:
            r["skew_u"] = dfloat(g1 * Decimal(n * (n - 1)).sqrt() / (n - 2))
        if n > 3:
            fu = Decimal(n - 1) / Decimal((n - 2) * (n - 3)) * (Decimal(n + 1) * g2 + 6)
            r["kurt_f_u"] = dfloat(fu)
            r["kurt_p_u"] = dfloat(fu + 3)
        r["cm3"] = float(cm[3])
        r["cm4"] = float(cm[4])
    r["raw3"] = float(sum(x ** 3 for x in fx) / n)
    sf = sorted(fx)
    if n % 2:
        r["med"] = float(sf[n // 2])
    else:
        r["med"] = float((sf[n // 2 - 1] + sf[n // 2]) / 2)
    return r


def exact_geo(xs):
    """The geometric mean from high precision logarithms."""
    logs = sum(Decimal(x).ln() for x in xs) / len(xs)
    return dfloat(logs.exp())


def exact_harm(xs):
    """The harmonic mean from rational arithmetic. Python statistics gives 0
    for values near 5e-324, because its reciprocal overflows."""
    return float(Fraction(len(xs)) / sum(1 / Fraction(x) for x in xs))


def srt(xs):
    return sorted(xs)


def sample(rng, kind, n):
    if kind == 0:
        return [rng.uniform(-100, 100) for _ in range(n)]
    if kind == 1:
        return [float(rng.randint(-20, 60)) for _ in range(n)]
    if kind == 2:
        off = rng.choice([1e6, 1e9, 1e12])
        return [off + rng.randint(-50, 50) for _ in range(n)]
    if kind == 3:
        return [rng.choice([1e16, -1e16, 1.0, -1.0, 0.5, 3.0]) for _ in range(n)]
    if kind == 4:
        return [rng.uniform(0.125, 8.0) for _ in range(n)]
    if kind == 5:
        return [rng.choice([1e-310, 2e-310, 5e-324, 1e-300, 3.0]) for _ in range(n)]
    if kind == 6:
        v = rng.uniform(-5, 5)
        return [v] * n
    return [rng.gauss(0, 1) * 10 ** rng.randint(-8, 8) for _ in range(n)]


def make_samples(rng, count):
    out = []
    for k in range(count):
        kind = k % 8
        n = rng.randint(4, 40)
        out.append(sample(rng, kind, n))
    return out


def close(a, e, ulps, label, xs):
    if a == e:
        return
    tol = ulps * EPS * max(abs(e), 1e-320)
    if abs(a - e) > tol:
        sys.exit("DIFFERENCE %s: got %r expected %r sample %r" % (label, a, e, xs[:6]))


def check_second_sources(xs, ex):
    """Compare the second sources with the exact value."""
    n = len(xs)
    a = np.array(xs, dtype=float)
    close(float(st.fmean(xs)), ex["mean"], 2, "python fmean", xs) if max(map(abs, xs)) < 1e300 else None
    if n > 1 and max(map(abs, xs)) < 1e150 and min(map(abs, xs)) > 1e-150:
        v = st.variance(xs)
        close(float(v), ex["var"], 4, "python variance", xs)
    spread = ex["pstdev"]
    wellcond = spread > 0 and abs(ex["mean"]) < 1e4 * spread and 1e-100 < max(map(abs, xs)) < 1e100
    if all(x > 1e-150 for x in xs) and max(xs) < 1e150:
        g = st.geometric_mean(xs)
        if abs(g - exact_geo(xs)) > 64 * EPS * g:
            sys.exit("DIFFERENCE python geometric_mean %r" % g)
        h = st.harmonic_mean(xs)
        if abs(h - exact_harm(xs)) > 8 * EPS * h:
            sys.exit("DIFFERENCE python harmonic_mean %r" % h)
    if "skew_b" in ex and wellcond:
        sb = float(stats.skew(a))
        if abs(sb - ex["skew_b"]) > 1e-6 * max(1.0, abs(ex["skew_b"])):
            sys.exit("DIFFERENCE scipy skew %r %r" % (sb, ex["skew_b"]))
        kb = float(stats.kurtosis(a))
        if abs(kb - ex["kurt_f_b"]) > 1e-6 * max(1.0, abs(ex["kurt_f_b"])):
            sys.exit("DIFFERENCE scipy kurtosis %r %r" % (kb, ex["kurt_f_b"]))


def fl(x):
    r = repr(float(x))
    if r in ("inf", "-inf", "nan"):
        sys.exit("non finite vector value")
    return r


def arr(name, vals, typ, comment):
    lines = ["// " + comment, "const %s: %s[%d] = [" % (name, typ, len(vals))]
    row = []
    for v in vals:
        row.append(str(v))
        if len(row) == 8:
            lines.append("    " + ", ".join(row) + ",")
            row = []
    if row:
        lines.append("    " + ", ".join(row) + ",")
    lines.append("];")
    return "\n".join(lines) + "\n"


def numpy_quantile(xs, q, m):
    return float(np.quantile(np.array(xs, dtype=float), q, method=m))


def main():
    use_go = "--no-go" not in sys.argv
    use_rust = "--no-rust" not in sys.argv
    use_node = "--no-node" not in sys.argv
    rng = random.Random(20261109)
    # The generator compares the second sources on many samples. The
    # package keeps the first KEEP of them.
    COMPARE = 10000
    KEEP = 1200
    allsamples = make_samples(rng, COMPARE)
    kept = allsamples[:KEEP]
    cols = {k: [] for k in ["sum", "mean", "median", "mlow", "mhigh", "pvar", "var", "pstdev", "stdev",
                            "sem", "variation", "skew_b", "skew_u", "kurt_f_b", "kurt_p_b", "kurt_f_u",
                            "kurt_p_u", "cm3", "cm4", "raw3", "geo", "harm", "iqr", "tmean", "trim",
                            "grouped", "mode"]}
    flags = {"const": [], "pos": [], "zerosum": []}
    quant = []
    q99 = []
    for idx, xs in enumerate(allsamples):
        ex = exact(xs)
        check_second_sources(xs, ex)
        s = srt(xs)
        n = len(xs)
        # NumPy quantiles against the package formulas are bit equal.
        for m in METHODS:
            for q in QS:
                v = numpy_quantile(xs, q, m)
                if idx < KEEP:
                    quant.append(v)
        if idx < 200:
            for m in METHODS:
                for i in range(1, 100):
                    q99.append(numpy_quantile(xs, i / 100, m))
        if idx >= KEEP:
            continue
        const = len(set(xs)) == 1
        pos = all(x > 0 for x in xs)
        flags["const"].append(1 if const else 0)
        flags["pos"].append(1 if pos else 0)
        cols["sum"].append(ex["sum"])
        cols["mean"].append(ex["mean"])
        cols["median"].append(float(np.median(np.array(xs))) if max(map(abs, xs)) < 1e300 else st.median(xs))
        cols["mlow"].append(st.median_low(xs))
        cols["mhigh"].append(st.median_high(xs))
        cols["pvar"].append(ex["pvar"])
        cols["var"].append(ex["var"])
        cols["pstdev"].append(ex["pstdev"])
        cols["stdev"].append(ex["stdev"])
        cols["sem"].append(ex["sem"])
        cols["variation"].append(ex["variation"] if ex["variation"] is not None else 0.0)
        for k in ["skew_b", "skew_u", "kurt_f_b", "kurt_p_b", "kurt_f_u", "kurt_p_u", "cm3", "cm4"]:
            cols[k].append(ex.get(k, 0.0))
        cols["raw3"].append(ex["raw3"])
        cols["geo"].append(exact_geo(xs) if pos else 0.0)
        cols["harm"].append(exact_harm(xs) if pos else 0.0)
        qa = numpy_quantile(xs, 0.75, "linear")
        qb = numpy_quantile(xs, 0.25, "linear")
        cols["iqr"].append(qa - qb)
        lo, hi = s[n // 4], s[3 * n // 4]
        sel = [x for x in xs if lo <= x <= hi]
        cols["tmean"].append(float(sum(map(Fraction, sel)) / len(sel)))
        cut = int(n * 0.1)
        cols["trim"].append(float(sum(map(Fraction, s[cut:n - cut])) / (n - 2 * cut)))
        cols["grouped"].append(st.median_grouped(xs, 2.0) if max(map(abs, xs)) < 1e300 else 0.0)
        cols["mode"].append(st.mode(xs))
    # Weighted samples with whole number weights. The repeated sample is
    # the second source.
    wrng = random.Random(77)
    wxs, wws, wcols = [], [], {k: [] for k in ["wmean", "wpop", "wfreq", "wrel", "wmed"]}
    for k in range(10000):
        n = wrng.randint(2, 12)
        xs = [float(wrng.randint(-30, 30)) + wrng.choice([0.0, 0.5, 0.25]) for _ in range(n)]
        ws = [float(wrng.randint(0, 6)) for _ in range(n)]
        if sum(ws) < 2:
            ws[0] += 2.0
        rep = [x for x, w in zip(xs, ws) for _ in range(int(w))]
        W = sum(ws)
        fx = [Fraction(x) for x in xs]
        fw = [Fraction(w) for w in ws]
        mean = sum(a * b for a, b in zip(fx, fw)) / sum(fw)
        S = sum(w * (x - mean) ** 2 for x, w in zip(fx, fw))
        W2 = sum(w * w for w in fw)
        nonzero = sum(1 for w in ws if w > 0)
        if nonzero < 2 or W - W2 / W <= 0:
            continue
        wmean = float(mean)
        close(float(np.average(xs, weights=ws)), wmean, 4, "numpy average", xs)
        pop = float(S / sum(fw))
        freq = float(S / (sum(fw) - 1)) if W > 1 else None
        rel = float(S / (sum(fw) - W2 / sum(fw)))
        if freq is None:
            continue
        close(float(np.cov(xs, fweights=[int(w) for w in ws])), freq, 8, "numpy cov fweights", xs)
        med = float(np.quantile(np.array(rep), 0.5, method="inverted_cdf"))
        if k >= 600:
            continue
        wxs.append(xs)
        wws.append(ws)
        wcols["wmean"].append(wmean)
        wcols["wpop"].append(pop)
        wcols["wfreq"].append(freq)
        wcols["wrel"].append(rel)
        wcols["wmed"].append(med)
    print("weighted compared %d samples, kept %d" % (k + 1, len(wxs)))
    if use_go:
        run_gonum(allsamples[:2000])
    if use_rust:
        run_rust(allsamples[:2000])
    if use_node:
        run_node(allsamples[:2000])
    write(kept, cols, flags, quant, wxs, wws, wcols, q99)


def run_gonum(samples):
    """Compare gonum with the exact value. The oracle program is in test/oracle."""
    inp = os.path.join(HERE, "..", "..", "..", "tmp", "stat", "gonum_in.json")
    inp = os.path.abspath(inp)
    os.makedirs(os.path.dirname(inp), exist_ok=True)
    usable = [xs for xs in samples if max(map(abs, xs)) < 1e100 and min(map(abs, xs)) > 1e-100 or max(map(abs, xs)) < 1e100 and min(map(abs, xs)) == 0]
    json.dump(usable, open(inp, "w"))
    env = dict(os.environ, GOFLAGS="-mod=mod", GOPROXY="off")
    res = subprocess.run(["go", "run", ".", inp], cwd=os.path.join(HERE, "oracle"), env=env,
                         capture_output=True, text=True)
    if res.returncode != 0:
        sys.exit("gonum oracle failed:\n" + res.stderr)
    rows = json.loads(res.stdout)
    for xs, g in zip(usable, rows):
        ex = exact(xs)
        # gonum sums with a plain loop. Its error follows the size of the
        # values and not the size of the mean.
        top = max(map(abs, xs))
        if abs(g["mean"] - ex["mean"]) > 64 * EPS * top:
            sys.exit("DIFFERENCE gonum mean %r %r" % (g["mean"], ex["mean"]))
        if ex["var"] is not None and abs(g["variance"] - ex["var"]) > 1e-9 * max(ex["var"], top * top * 1e-12):
            sys.exit("DIFFERENCE gonum variance %r %r" % (g["variance"], ex["var"]))
    print("gonum checked %d samples" % len(rows))


def run_rust(samples):
    """Compare statrs of Rust with the exact value.

    The oracle program is in test/oracle-rust. Each line of the input
    file holds one sample of space separated values.
    """
    inp = os.path.abspath(os.path.join(HERE, "..", "..", "..", "tmp", "stat", "rust_in.txt"))
    os.makedirs(os.path.dirname(inp), exist_ok=True)
    usable = [xs for xs in samples if max(map(abs, xs)) < 1e100 and min(map(abs, xs)) > 1e-100]
    with open(inp, "w") as f:
        for xs in usable:
            f.write(" ".join(repr(float(x)) for x in xs) + "\n")
    res = subprocess.run(["cargo", "run", "--release", "--quiet", "--", inp],
                         cwd=os.path.join(HERE, "oracle-rust"), capture_output=True, text=True)
    if res.returncode != 0:
        sys.exit("statrs oracle failed:\n" + res.stderr)
    for xs, line in zip(usable, res.stdout.strip().split("\n")):
        if not line:
            continue
        g = [float(w) for w in line.split()]
        ex = exact(xs)
        top = max(map(abs, xs))
        if abs(g[0] - ex["mean"]) > 64 * EPS * top:
            sys.exit("DIFFERENCE statrs mean %r %r" % (g[0], ex["mean"]))
        if ex["var"] is not None and abs(g[1] - ex["var"]) > 1e-9 * max(ex["var"], top * top * 1e-12):
            sys.exit("DIFFERENCE statrs variance %r %r" % (g[1], ex["var"]))
        if abs(g[2] - ex["pvar"]) > 1e-9 * max(ex["pvar"], top * top * 1e-12):
            sys.exit("DIFFERENCE statrs population variance %r %r" % (g[2], ex["pvar"]))
        if all(x > 0 for x in xs):
            if abs(g[4] - exact_geo(xs)) > 1e-9 * max(exact_geo(xs), 1e-300):
                sys.exit("DIFFERENCE statrs geometric_mean %r %r" % (g[4], exact_geo(xs)))
            if abs(g[5] - exact_harm(xs)) > 1e-9 * max(exact_harm(xs), 1e-300):
                sys.exit("DIFFERENCE statrs harmonic_mean %r %r" % (g[5], exact_harm(xs)))
        if abs(g[7] - ex["med"]) > 64 * EPS * top:
            sys.exit("DIFFERENCE statrs median %r %r" % (g[7], ex["med"]))
    print("statrs checked %d samples" % len(usable))


def run_node(samples):
    """Compare simple-statistics of Node with the exact value.

    The oracle program is in test/oracle-node. Each line of the input
    file holds one sample of space separated values.
    """
    inp = os.path.abspath(os.path.join(HERE, "..", "..", "..", "tmp", "stat", "node_in.txt"))
    os.makedirs(os.path.dirname(inp), exist_ok=True)
    usable = [xs for xs in samples if max(map(abs, xs)) < 1e100 and min(map(abs, xs)) > 1e-100]
    with open(inp, "w") as f:
        for xs in usable:
            f.write(" ".join(repr(float(x)) for x in xs) + "\n")
    res = subprocess.run(["node", "oracle.mjs", inp],
                         cwd=os.path.join(HERE, "oracle-node"), capture_output=True, text=True)
    if res.returncode != 0:
        sys.exit("simple-statistics oracle failed:\n" + res.stderr)
    for xs, line in zip(usable, res.stdout.strip().split("\n")):
        if not line:
            continue
        g = [float(w) for w in line.split()]
        ex = exact(xs)
        top = max(map(abs, xs))
        if abs(g[0] - ex["mean"]) > 64 * EPS * top:
            sys.exit("DIFFERENCE simple-statistics mean %r %r" % (g[0], ex["mean"]))
        if ex["var"] is not None and abs(g[1] - ex["var"]) > 1e-9 * max(ex["var"], top * top * 1e-12):
            sys.exit("DIFFERENCE simple-statistics variance %r %r" % (g[1], ex["var"]))
        if abs(g[2] - ex["pvar"]) > 1e-9 * max(ex["pvar"], top * top * 1e-12):
            sys.exit("DIFFERENCE simple-statistics population variance %r %r" % (g[2], ex["pvar"]))
        if all(x > 0 for x in xs):
            if abs(g[4] - exact_geo(xs)) > 1e-9 * max(exact_geo(xs), 1e-300):
                sys.exit("DIFFERENCE simple-statistics geometricMean %r %r" % (g[4], exact_geo(xs)))
            if abs(g[5] - exact_harm(xs)) > 1e-9 * max(exact_harm(xs), 1e-300):
                sys.exit("DIFFERENCE simple-statistics harmonicMean %r %r" % (g[5], exact_harm(xs)))
        if abs(g[7] - ex["med"]) > 64 * EPS * top:
            sys.exit("DIFFERENCE simple-statistics median %r %r" % (g[7], ex["med"]))
    print("simple-statistics checked %d samples" % len(usable))


def write(kept, cols, flags, quant, wxs, wws, wcols, q99):
    flat, offs = [], [0]
    for xs in kept:
        flat.extend(xs)
        offs.append(len(flat))
    out = ["// Generated by test/make_vectors.py. Do not edit by hand.",
           "// Reference vectors of the statistics package.",
           "// The exact value of each statistic comes from rational arithmetic.",
           "// NumPy, SciPy, Python statistics, gonum, statrs and",
           "// simple-statistics agree with it within the",
           "// tolerance of the test. A value of 0 in a column marks a sample that",
           "// the statistic does not cover; the flag columns tell which.",
           "namespace suite;", "",
           "// The fractions of the quantile vectors.",
           "const QS: f64[8] = [0.0, 0.1, 0.25, 0.3, 0.5, 0.75, 0.9, 1.0];", ""]
    out.append(arr("V_FLAT", [fl(x) for x in flat], "f64", "All samples, one after the other."))
    out.append(arr("V_OFFS", offs, "i64", "The start of each sample in V_FLAT, and the end."))
    out.append(arr("V_CONST", flags["const"], "i64", "1 when every value of the sample is equal."))
    out.append(arr("V_POS", flags["pos"], "i64", "1 when every value of the sample is above zero."))
    for k, v in cols.items():
        out.append(arr("E_" + k.upper(), [fl(x) for x in v], "f64", "Expected " + k + " of each sample."))
    out.append(arr("E_QUANT", [fl(x) for x in quant], "f64",
                   "Expected quantile of each sample: 9 methods by 8 fractions (0, 0.1, 0.25, 0.3, 0.5, 0.75, 0.9, 1)."))
    out.append(arr("E_Q99", [fl(x) for x in q99], "f64",
                   "Expected quantile of each of the first 200 samples: 9 methods by 99 fractions (0.01 to 0.99)."))
    wflat, woffs = [], [0]
    wwflat = []
    for xs, ws in zip(wxs, wws):
        wflat.extend(xs)
        wwflat.extend(ws)
        woffs.append(len(wflat))
    out.append(arr("W_FLAT", [fl(x) for x in wflat], "f64", "Weighted samples, one after the other."))
    out.append(arr("W_WEIGHTS", [fl(x) for x in wwflat], "f64", "The whole number weights of W_FLAT."))
    out.append(arr("W_OFFS", woffs, "i64", "The start of each weighted sample, and the end."))
    for k, v in wcols.items():
        out.append(arr("E_" + k.upper(), [fl(x) for x in v], "f64", "Expected " + k + " of each weighted sample."))
    with open(OUT, "w") as f:
        f.write("\n".join(out))
    print("wrote", OUT, "samples", len(kept), "weighted", len(wxs))


main()
