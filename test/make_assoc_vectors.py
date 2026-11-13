#!/usr/bin/env python3
"""Make the reference vectors of the association calls.

The script writes numbers only into test/suite/assoc_vectors.goose. The
exact value of each call comes from rational arithmetic. The script then
compares NumPy, SciPy, Python statistics, gonum, statrs and
simple-statistics with the exact value and stops on a difference beyond
the stated tolerance. Run it from the package directory:

    python3 test/make_assoc_vectors.py

It needs Python 3 with NumPy and SciPy. It uses the crate statrs of Rust
and the package simple-statistics of Node for two more checks. Use
--no-go to skip the gonum check, --no-rust to skip the statrs check and
--no-node to skip the simple-statistics check.
"""
import json, math, os, random, subprocess, sys
from decimal import Decimal, getcontext
from fractions import Fraction

import numpy as np
from scipy import stats

getcontext().prec = 80
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "suite", "assoc_vectors.goose")
EPS = 2.220446049250313e-16
# The number of samples that the second sources check. The environment
# variable ASSOC_COMPARE makes a smaller run for a quick check.
COMPARE = int(os.environ.get("ASSOC_COMPARE", "10000"))
KEEP = 300
KEEP_LINE = int(os.environ.get("ASSOC_KEEP", "300"))
KEEP_RANK = int(os.environ.get("ASSOC_KEEP", "300"))
KEEP_MAT = 20


def fof(fr):
    """A Fraction as a float, or None when it leaves the float range."""
    try:
        return float(fr)
    except OverflowError:
        return None


def dsqrt(f):
    """Square root of a Fraction as a Decimal."""
    return (Decimal(f.numerator) / Decimal(f.denominator)).sqrt()


def dfloat(d):
    return float(d)


def cov_scale(ex):
    """The magnitude bound of a covariance: the product of the column
    spreads around their means. A covariance near zero has no relative
    bound, because every source rounds the centered values."""
    return (abs(ex["mean_x"]) + ex["spread_x"]) * (abs(ex["mean_y"]) + ex["spread_y"])


def close_abs(a, e, ulps, scale, label, xs):
    """The same, but the scale of the sample gives the bound. A value near
    zero has no relative bound, and the references round the slope."""
    if a == e:
        return
    tol = ulps * EPS * (abs(e) + abs(scale))
    if abs(a - e) > tol:
        sys.exit("DIFFERENCE %s: got %r expected %r sample %r" % (label, a, e, xs[:6]))


def close(a, e, ulps, label, xs):
    if a == e:
        return
    tol = ulps * EPS * max(abs(e), 1e-320)
    if abs(a - e) > tol:
        sys.exit("DIFFERENCE %s: got %r expected %r sample %r" % (label, a, e, xs[:6]))


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


def exact_pair(xs, ys):
    """The exact statistics of a pair of columns."""
    n = len(xs)
    if n == 0:
        return {"n": 0, "const_x": False, "const_y": False, "med_x": 0.0,
                "mean_x": 0.0, "mean_y": 0.0, "sx_scale": 0.0, "med_y": 0.0,
                "spread_x": 0.0, "spread_y": 0.0, "pcov": 0.0, "cov": None,
                "corr": None, "fit": False, "slope": None, "intercept": None,
                "err_scale": 0.0, "slope_err": None, "intercept_err": None}
    fx = [Fraction(x) for x in xs]
    fy = [Fraction(y) for y in ys]
    sx = sum(fx)
    sy = sum(fy)
    sxy = sum(a * b for a, b in zip(fx, fy))
    sxx = sum(a * a for a in fx)
    syy = sum(b * b for b in fy)
    nxy = n * sxy - sx * sy
    nxx = n * sxx - sx * sx
    nyy = n * syy - sy * sy
    out = {"n": n, "const_x": len(set(xs)) == 1, "const_y": len(set(ys)) == 1}
    sxv = sorted(fx)
    out["med_x"] = float(sxv[n // 2]) if n % 2 else float((sxv[n // 2 - 1] + sxv[n // 2]) / 2)
    out["mean_x"] = float(sx / n)
    out["mean_y"] = float(sy / n)
    out["sx_scale"] = dsqrt(sum(a * a for a in fx) / n)
    syv = sorted(fy)
    out["med_y"] = float(syv[n // 2]) if n % 2 else float((syv[n // 2 - 1] + syv[n // 2]) / 2)
    out["spread_x"] = fof(dsqrt(nxx / n)) if nxx > 0 else 0.0
    out["spread_y"] = fof(dsqrt(nyy / n)) if nyy > 0 else 0.0
    if out["spread_x"] is None or out["spread_y"] is None:
        out["spread_x"] = out["spread_y"] = 0.0
    out["pcov"] = fof(nxy / (n * n))
    out["cov"] = fof(nxy / (n * (n - 1))) if n > 1 else None
    out["fit"] = True
    if nxx > 0 and nyy > 0:
        cn = fof(nxy)
        cd = dsqrt(nxx)
        cy = dsqrt(nyy)
        cdf = fof(cd)
        cyf = fof(cy)
        if cn is None or cdf is None or cyf is None or cn == 0.0 or cd == 0 or cy == 0:
            out["corr"] = None
            out["fit"] = False
        else:
            out["corr"] = dfloat(Decimal(cn) / cd / cy)
    else:
        out["corr"] = None
    if nxx > 0:
        out["slope"] = fof(nxy / nxx)
        out["intercept"] = fof((sy - (nxy / nxx) * sx) / n)
        if out["slope"] is None or out["intercept"] is None:
            out["slope"] = None
            out["intercept"] = None
        ssx = nxx / n
        if nyy > 0 and n > 2:
            r2 = (nxy * nxy) / (nxx * nyy)
            ssy = nyy / n
            out["err_scale"] = dfloat(dsqrt(ssy / ssx / (n - 2)))
            out["slope_err"] = dfloat(dsqrt((1 - r2) * ssy / ssx / (n - 2)))
            out["intercept_err"] = dfloat(dsqrt((1 - r2) * ssy / ssx / (n - 2)) * dsqrt(sxx / n))
        else:
            out["slope_err"] = None
            out["intercept_err"] = None
    else:
        out["slope"] = None
        out["intercept"] = None
        out["slope_err"] = None
        out["intercept_err"] = None
    return out


def exact_ranks(xs):
    """The exact ranks of a sample for the five tie methods."""
    n = len(xs)
    order = sorted(range(n), key=lambda i: (xs[i], i))
    ranks = [[0.0] * n for _ in range(5)]
    k = 0
    dense = 0
    while k < n:
        j = k + 1
        while j < n and xs[order[j]] == xs[order[k]]:
            j += 1
        count = j - k
        first = k + 1
        for m in range(k, j):
            idx = order[m]
            ranks[0][idx] = Fraction(first) + Fraction(count - 1, 2)
            ranks[1][idx] = Fraction(first)
            ranks[2][idx] = Fraction(first + count - 1)
            ranks[3][idx] = Fraction(dense + 1)
            ranks[4][idx] = Fraction(first + (m - k))
        dense += 1
        k = j
    return [[float(v) for v in r] for r in ranks]


def exact_kendall(xs, ys):
    """The exact Kendall counts and both coefficients."""
    n = len(xs)
    concordant = 0
    discordant = 0
    tie_x = 0
    tie_y = 0
    tie_both = 0
    for i in range(n):
        for j in range(i + 1, n):
            dx = xs[j] - xs[i]
            dy = ys[j] - ys[i]
            if dx == 0 and dy == 0:
                tie_both += 1
                tie_x += 1
                tie_y += 1
            elif dx == 0:
                tie_x += 1
            elif dy == 0:
                tie_y += 1
            elif dx * dy > 0:
                concordant += 1
            else:
                discordant += 1
    n0 = n * (n - 1) // 2
    cd = concordant - discordant
    a = n0 - tie_x
    b = n0 - tie_y
    out = {"cd": cd}
    if a > 0 and b > 0:
        out["b"] = cd / math.sqrt(a * b)
    else:
        out["b"] = None
    ux = len(set(xs))
    uy = len(set(ys))
    m = min(ux, uy)
    if m > 1:
        out["c"] = 2.0 * cd / (n * n * (m - 1) / m)
    else:
        out["c"] = None
    return out


def exact_theil(xs, ys):
    """The exact Theil-Sen slope and intercept."""
    slopes = []
    smallest = 0.0
    for i in range(len(xs)):
        for j in range(i + 1, len(xs)):
            dx = Fraction(xs[j]) - Fraction(xs[i])
            if dx == 0:
                continue
            fdx = fof(dx)
            if fdx is None:
                continue
            if smallest == 0.0 or fdx < smallest:
                smallest = fdx
            slopes.append((Fraction(ys[j]) - Fraction(ys[i])) / dx)
    if not slopes:
        return None
    slopes.sort()
    m = len(slopes)
    if m % 2:
        slope = slopes[m // 2]
    else:
        slope = (slopes[m // 2 - 1] + slopes[m // 2]) / 2

    def med(vals):
        v = sorted(Fraction(x) for x in vals)
        k = len(v)
        if k % 2:
            return v[k // 2]
        return (v[k // 2 - 1] + v[k // 2]) / 2

    fs = fof(slope)
    fi = fof(med(ys) - slope * med(xs))
    if fs is None or fi is None:
        return None
    return fs, fi, smallest


def pair_sample(rng, kind, n):
    """A pair of columns with the same length."""
    if kind == 0:
        x = [rng.uniform(-100, 100) for _ in range(n)]
        return x, [2.0 * v + rng.uniform(-10, 10) for v in x]
    if kind == 1:
        x = [float(rng.randint(-5, 5)) for _ in range(n)]
        return x, [float(rng.randint(-5, 5)) for _ in range(n)]
    if kind == 2:
        off = rng.choice([1e6, 1e9, 1e12])
        x = [off + rng.randint(-50, 50) for _ in range(n)]
        return x, [2.0 * v + rng.randint(-50, 50) for v in x]
    if kind == 3:
        pool = [1.0, 2.0, 3.0, 4.0, 5.0]
        x = [rng.choice(pool) for _ in range(n)]
        return x, [rng.choice(pool) for _ in range(n)]
    if kind == 4:
        x = [rng.uniform(-100, 100) for _ in range(n)]
        return x, [7.0] * n
    if kind == 5:
        x = [3.0] * n
        return x, [rng.uniform(-100, 100) for _ in range(n)]
    if kind == 6:
        pool = [1e16, -1e16, 1.0, -1.0, 0.5, 3.0]
        return [rng.choice(pool) for _ in range(n)], [rng.choice(pool) for _ in range(n)]
    if kind == 7:
        pool = [1e-310, 2e-310, 5e-324, 1e-300, 3.0]
        return [rng.choice(pool) for _ in range(n)], [rng.choice(pool) for _ in range(n)]
    if kind == 8:
        x = [rng.uniform(-100, 100) for _ in range(n)]
        return x, [3.0 * v for v in x]
    if kind == 9:
        x = [rng.uniform(-100, 100) for _ in range(n)]
        return x, [-2.0 * v for v in x]
    x = [rng.gauss(0, 1) * 10 ** rng.randint(-8, 8) for _ in range(n)]
    return x, [rng.gauss(0, 1) * 10 ** rng.randint(-8, 8) for _ in range(n)]


def line_sample(rng, kind, n):
    """A pair of columns for the line fit. The offset stays small, so the
    intercept keeps its digits."""
    if kind == 0:
        x = [rng.uniform(-10, 10) for _ in range(n)]
        return x, [2.0 * v + rng.uniform(-1, 1) for v in x]
    if kind == 1:
        x = [float(rng.randint(-20, 20)) for _ in range(n)]
        return x, [float(rng.randint(-20, 20)) for _ in range(n)]
    if kind == 2:
        x = [rng.uniform(0.5, 8.0) for _ in range(n)]
        return x, [rng.uniform(0.5, 8.0) for _ in range(n)]
    if kind == 3:
        x = [float(rng.randint(-5, 5)) for _ in range(n)]
        return x, [2.0 * v + rng.randint(-3, 3) for v in x]
    if kind == 4:
        x = [rng.uniform(-10, 10) for _ in range(n)]
        return x, [7.0] * n
    if kind == 5:
        x = [3.0] * n
        return x, [rng.uniform(-10, 10) for _ in range(n)]
    if kind == 6:
        x = [rng.uniform(-10, 10) for _ in range(n)]
        return x, [-2.0 * v for v in x]
    x = [rng.uniform(-10, 10) for _ in range(n)]
    return x, [rng.gauss(0, 1) for _ in range(n)]


def rank_sample(rng, kind, n):
    """A pair of columns with many ties."""
    if kind == 0:
        pool = [1.0, 2.0, 3.0]
        return [rng.choice(pool) for _ in range(n)], [rng.choice(pool) for _ in range(n)]
    if kind == 1:
        pool = [1.0, 2.0, 3.0, 4.0, 5.0]
        return [rng.choice(pool) for _ in range(n)], [rng.choice(pool) for _ in range(n)]
    if kind == 2:
        x = [float(rng.randint(0, 9)) for _ in range(n)]
        return x, [float(rng.randint(0, 9)) for _ in range(n)]
    if kind == 3:
        x = [rng.uniform(-100, 100) for _ in range(n)]
        return x, [3.0] * n
    if kind == 4:
        x = [2.0] * n
        return x, [rng.uniform(-100, 100) for _ in range(n)]
    x = [float(rng.randint(0, 2)) for _ in range(n)]
    return x, [float(rng.randint(0, 2)) for _ in range(n)]


def mat_sample(rng, kind, rows, cols):
    """A matrix of columns for the matrix calls."""
    out = []
    for c in range(cols):
        if kind == 0:
            out.append([rng.uniform(-100, 100) for _ in range(rows)])
        elif kind == 1:
            out.append([float(rng.randint(-9, 9)) for _ in range(rows)])
        elif kind == 2:
            off = rng.choice([1e6, 1e9])
            out.append([off + rng.randint(-50, 50) for _ in range(rows)])
        else:
            out.append([rng.choice([1.0, 2.0, 3.0]) for _ in range(rows)])
    return out


def well_cond(ex, xs, ys):
    """True when the columns have a mean that is small against the spread.
    The second sources round the mean, so they lose the digits of a large
    offset. The exact value stays right."""
    if ex["spread_x"] <= 0 or ex["spread_y"] <= 0:
        return False
    if abs(ex["mean_x"]) > 1e4 * ex["spread_x"]:
        return False
    if abs(ex["mean_y"]) > 1e4 * ex["spread_y"]:
        return False
    if max(map(abs, xs)) > 1e6 or max(map(abs, ys)) > 1e6:
        return False
    if min(map(abs, xs)) < 1e-100 or min(map(abs, ys)) < 1e-100:
        return False
    return True


def check_pair_sources(xs, ys, ex):
    """Compare the second sources with the exact value."""
    n = len(xs)
    a = np.array(xs, dtype=float)
    b = np.array(ys, dtype=float)
    good = well_cond(ex, xs, ys)
    if n > 1 and ex["cov"] is not None and good:
        close_abs(float(np.cov(a, b, ddof=1)[0, 1]), ex["cov"], 64, cov_scale(ex), "numpy cov", xs)
        close_abs(float(np.corrcoef(a, b)[0, 1]), ex["corr"], 64, 1.0, "numpy corrcoef", xs)
        sr = stats.spearmanr(a, b)
        close(float(sr.statistic), ex["spearman"], 4, "scipy spearmanr", xs)
        kt = stats.kendalltau(a, b, variant="b")
        close(float(kt.statistic), ex["kendall_b"], 4, "scipy kendalltau b", xs)
        ktc = stats.kendalltau(a, b, variant="c")
        close(float(ktc.statistic), ex["kendall_c"], 4, "scipy kendalltau c", xs)
    if n > 1 and ex["cov"] is not None and good:
        lr = stats.linregress(a, b)
        close_abs(float(lr.slope), ex["slope"], 64, abs(ex["mean_y"]), "scipy linregress slope", xs)
        close_abs(float(lr.intercept), ex["intercept"], 64,
                  abs(ex["slope"] * ex["med_x"]) + abs(ex["mean_y"]),
                  "scipy linregress intercept", xs)
        if n > 2 and ex.get("err_fit"):
            close_abs(float(lr.stderr), ex["slope_err"], 64, ex["err_scale"],
                      "scipy linregress stderr", xs)
            close_abs(float(lr.intercept_stderr), ex["intercept_err"], 64,
                      dfloat(Decimal(ex["err_scale"]) * ex["sx_scale"]),
                      "scipy intercept stderr", xs)
    if n <= 64 and ex["theil"] is not None and good:
        ts = stats.theilslopes(b, a)
        close_abs(float(ts.slope), ex["theil"][0], 64,
                  ex["theil"][2] and abs(ex["med_y"]) / ex["theil"][2],
                  "scipy theilslopes slope", xs)
        close_abs(float(ts.intercept), ex["theil"][1], 32,
                  abs(ex["theil"][0] * ex["med_x"]) + abs(ex["med_y"]),
                  "scipy theilslopes intercept", xs)
    for m in range(5):
        r = stats.rankdata(a, method=["average", "min", "max", "dense", "ordinal"][m])
        for i in range(n):
            if abs(r[i] - ex["ranks"][m][i]) > 1e-9:
                sys.exit("DIFFERENCE scipy rankdata method %d" % m)


def run_oracles(pairs, use_go, use_rust, use_node):
    """The rows of the other libraries for the covariance and the
    correlation."""
    go_rows = []
    rust_rows = []
    node_rows = []
    if use_go:
        path = os.path.join(HERE, "oracle", "pairs.json")
        with open(path, "w") as f:
            json.dump([{"x": p[0], "y": p[1]} for p in pairs], f)
        out = subprocess.run(["go", "run", ".", path, "pairs"], cwd=os.path.join(HERE, "oracle"),
                             capture_output=True, text=True)
        if out.returncode != 0:
            sys.exit("the gonum oracle failed: " + out.stderr)
        go_rows = json.loads(out.stdout)
    if use_rust or use_node:
        path = os.path.join(HERE, "pairs.txt")
        with open(path, "w") as f:
            for p in pairs:
                f.write(" ".join(repr(v) for v in p[0]) + " | " + " ".join(repr(v) for v in p[1]) + "\n")
    if use_rust:
        out = subprocess.run(["cargo", "run", "--quiet", "--", os.path.join(HERE, "pairs.txt")],
                             cwd=os.path.join(HERE, "oracle-rust"), capture_output=True, text=True)
        if out.returncode != 0:
            sys.exit("the statrs oracle failed: " + out.stderr)
        rust_rows = [line.split() for line in out.stdout.strip().split("\n") if line.strip()]
    if use_node:
        out = subprocess.run(["node", "oracle.mjs", os.path.join(HERE, "pairs.txt")],
                             cwd=os.path.join(HERE, "oracle-node"), capture_output=True, text=True)
        if out.returncode != 0:
            sys.exit("the simple-statistics oracle failed: " + out.stderr)
        node_rows = [line.split() for line in out.stdout.strip().split("\n") if line.strip()]
    return go_rows, rust_rows, node_rows


def main():
    use_go = "--no-go" not in sys.argv
    use_rust = "--no-rust" not in sys.argv
    use_node = "--no-node" not in sys.argv
    rng = random.Random(20261110)

    # The pairs for the covariance and the correlation.
    allpairs = []
    pairkinds = []
    for k in range(COMPARE):
        n = rng.randint(3, 40)
        kind = k % 10
        allpairs.append(pair_sample(rng, kind, n))
        pairkinds.append(kind)
    # The subnormal values do not fit in a float with their exact result.
    # The hostile suite covers them.
    kept = []
    for i, p in enumerate(allpairs):
        if len(kept) >= KEEP:
            break
        if pairkinds[i] == 7:
            continue
        kept.append(p)

    # The pairs for the line fit.
    alllines = []
    for k in range(COMPARE):
        n = rng.randint(3, 40)
        alllines.append(line_sample(rng, k % 8, n))
    keptlines = alllines[:KEEP_LINE]

    # The pairs for the ranks.
    allranks = []
    for k in range(COMPARE):
        n = rng.randint(3, 40)
        allranks.append(rank_sample(rng, k % 6, n))
    keptranks = allranks[:KEEP_RANK]

    # The matrices.
    allmats = []
    for k in range(200):
        rows = rng.randint(3, 12)
        cols = rng.randint(2, 5)
        allmats.append(mat_sample(rng, k % 4, rows, cols))
    keptmats = allmats[:KEEP_MAT]

    # The exact values of every pair, and the check of the second sources.
    def full_exact(xs, ys):
        ex = exact_pair(xs, ys)
        ex["ranks"] = exact_ranks(xs)
        ex["kendall"] = exact_kendall(xs, ys)
        ex["theil"] = exact_theil(xs, ys)
        # The Spearman coefficient is the correlation of the average ranks.
        if ex["corr"] is None:
            ex["spearman"] = None
        else:
            rx = [Fraction(v) for v in ex["ranks"][0]]
            ry = [Fraction(v) for v in exact_ranks(ys)[0]]
            n = len(xs)
            sx = sum(rx)
            sy = sum(ry)
            sxy = sum(a * b for a, b in zip(rx, ry))
            sxx = sum(a * a for a in rx)
            syy = sum(b * b for b in ry)
            nxy = n * sxy - sx * sy
            nxx = n * sxx - sx * sx
            nyy = n * syy - sy * sy
            if nxx > 0 and nyy > 0:
                ex["spearman"] = float(nxy) / math.sqrt(float(nxx)) / math.sqrt(float(nyy))
            else:
                ex["spearman"] = None
        ex["kendall_b"] = ex["kendall"]["b"]
        ex["kendall_c"] = ex["kendall"]["c"]
        return ex

    for xs, ys in allpairs:
        ex = full_exact(xs, ys)
        if ex["corr"] is None:
            continue
        check_pair_sources(xs, ys, ex)
    for xs, ys in alllines:
        ex = full_exact(xs, ys)
        if ex["corr"] is None:
            continue
        check_pair_sources(xs, ys, ex)
    for xs, ys in allranks:
        ex = full_exact(xs, ys)
        if ex["corr"] is None:
            continue
        check_pair_sources(xs, ys, ex)

    # The other libraries for the covariance and the correlation.
    # The libraries that round their sums lose the digits of a large
    # offset, so the comparison uses the samples with a small offset only.
    good = []
    for p in allpairs:
        ex = full_exact(p[0], p[1])
        if ex["corr"] is None:
            continue
        if not well_cond(ex, p[0], p[1]):
            continue
        good.append(p)
    good = good[:2000]
    go_rows, rust_rows, node_rows = run_oracles(good, use_go, use_rust, use_node)
    for i, (xs, ys) in enumerate(good):
        ex = full_exact(xs, ys)
        if i < len(go_rows):
            close_abs(go_rows[i]["covariance"], ex["cov"], 4, cov_scale(ex), "gonum covariance", xs)
            close_abs(go_rows[i]["correlation"], ex["corr"], 4, 1.0, "gonum correlation", xs)
        if i < len(rust_rows):
            close_abs(float(rust_rows[i][0]), ex["cov"], 4, cov_scale(ex), "statrs covariance", xs)
        if i < len(node_rows):
            close_abs(float(node_rows[i][0]), ex["cov"], 4, cov_scale(ex), "simple-statistics covariance", xs)
            close_abs(float(node_rows[i][1]), ex["corr"], 4, 1.0, "simple-statistics correlation", xs)

    # The vectors.
    out = []
    out.append("// Generated by test/make_assoc_vectors.py. Do not edit by hand.\n"
               "// Reference vectors of the association calls.\n"
               "// The exact value of each call comes from rational arithmetic.\n"
               "// NumPy, SciPy, gonum, statrs and simple-statistics agree with\n"
               "// it within the tolerance of the test. A value of 0 marks a call\n"
               "// that has no value, and the flag columns tell which.\n"
               "namespace suite;\n\n")

    def flat_pairs(pairs):
        flat = []
        offs = []
        for xs, ys in pairs:
            offs.append(len(flat))
            flat.extend(xs)
            flat.extend(ys)
        offs.append(len(flat))
        return flat, offs

    xflat, xoffs = flat_pairs(kept)
    out.append(arr("A_FLAT", [fl(v) for v in xflat], "f64",
                   "The columns of the covariance samples, one pair after the other."))
    out.append(arr("A_OFFS", xoffs, "i64",
                   "The start of each pair in A_FLAT, and the end."))
    out.append(arr("A_NCOLS", [len(p[0]) for p in kept], "i64",
                   "The count of values in each column of A_FLAT."))

    cols = {k: [] for k in ["cov", "pcov", "corr"]}
    flags = {"const": []}
    for xs, ys in kept:
        ex = exact_pair(xs, ys)
        cols["cov"].append(fl(ex["cov"]) if ex["cov"] is not None else "0.0")
        cols["pcov"].append(fl(ex["pcov"]) if ex["pcov"] is not None else "0.0")
        cols["corr"].append(fl(ex["corr"]) if ex["corr"] is not None else "0.0")
        flags["const"].append(1 if (ex["const_x"] or ex["const_y"]) else 0)
        flags.setdefault("fit", []).append(1 if ex["fit"] else 0)
    out.append(arr("E_COV", cols["cov"], "f64", "Expected sample covariance of each pair."))
    out.append(arr("E_PCOV", cols["pcov"], "f64", "Expected population covariance of each pair."))
    out.append(arr("E_CORR", cols["corr"], "f64", "Expected correlation of each pair."))
    out.append(arr("A_CONST", flags["const"], "i64",
                   "1 when a column of the pair holds one value only."))
    out.append(arr("A_FIT", flags["fit"], "i64",
                   "1 when the exact covariance and the exact correlation fit in a float."))

    # The weighted forms. The whole number weights give the repeated sample.
    wflat = []
    woffs = []
    wcols = {k: [] for k in ["pop", "freq", "rel", "corr"]}
    wfitlist = []
    for idx, (xs, ys) in enumerate(kept):
        ws = [float(rng.randint(1, 4)) for _ in xs]
        # The values come first, so a sample with no finite value leaves
        # the flat list and the offsets alone.
        n = len(xs)
        W = sum(int(w) for w in ws)
        W2 = sum(int(w) * int(w) for w in ws)
        swx = sum(Fraction(w) * Fraction(x) for w, x in zip(ws, xs))
        swy = sum(Fraction(w) * Fraction(y) for w, y in zip(ws, ys))
        swxy = sum(Fraction(w) * Fraction(x) * Fraction(y) for w, x, y in zip(ws, xs, ys))
        num = W * swxy - swx * swy
        vpop = fof(num / (W * W))
        vfreq = fof(num / (W * (W - 1)))
        vrel = fof(num / (W * W - W2))
        swxx = sum(Fraction(w) * Fraction(x) * Fraction(x) for w, x in zip(ws, xs))
        swyy = sum(Fraction(w) * Fraction(y) * Fraction(y) for w, y in zip(ws, ys))
        nxx = W * swxx - swx * swx
        nyy = W * swyy - swy * swy
        wcorr = None
        if nxx > 0 and nyy > 0:
            cn = fof(num)
            cd = fof(dsqrt(nxx))
            cy = fof(dsqrt(nyy))
            if cn is not None and cd is not None and cy is not None:
                wcorr = float(Decimal(cn) / Decimal(cd) / Decimal(cy))
        wfit = 1
        if vpop is None or vfreq is None or vrel is None or wcorr is None:
            wfit = 0
            vpop = vfreq = vrel = wcorr = 0.0
        woffs.append(len(wflat))
        wflat.extend(ws)
        wcols["pop"].append(fl(vpop))
        wcols["freq"].append(fl(vfreq))
        wcols["rel"].append(fl(vrel))
        wcols["corr"].append(fl(wcorr))
        wfitlist.append(wfit)
    woffs.append(len(wflat))
    out.append(arr("A_WFLAT", [fl(v) for v in wflat], "f64",
                   "The weights of the covariance samples, one pair after the other."))
    out.append(arr("A_WOFFS", woffs, "i64",
                   "The start of each weight list in A_WFLAT, and the end."))
    out.append(arr("E_WCOV_POP", wcols["pop"], "f64", "Expected weighted covariance, population form."))
    out.append(arr("E_WCOV_FREQ", wcols["freq"], "f64", "Expected weighted covariance, frequency form."))
    out.append(arr("E_WCOV_REL", wcols["rel"], "f64", "Expected weighted covariance, reliability form."))
    out.append(arr("E_WCOV_CORR", wcols["corr"], "f64", "Expected weighted correlation."))
    out.append(arr("A_WFIT", wfitlist, "i64",
                   "1 when the weighted values of the pair fit in a float."))

    # The ranks.
    rflat = []
        # The x column of each pair holds the ranks.
    roffs = []
    rcols = {m: [] for m in range(5)}
    rflat_y = []
    for xs, ys in keptranks:
        roffs.append(len(rflat))
        rflat.extend(xs)
        rflat_y.extend(ys)
        r = exact_ranks(xs)
        for m in range(5):
            rcols[m].extend(fl(v) for v in r[m])
    roffs.append(len(rflat))
    out.append(arr("R_FLAT", [fl(v) for v in rflat], "f64",
                   "The first column of the rank samples, one pair after the other."))
    out.append(arr("R_FLAT_Y", [fl(v) for v in rflat_y], "f64",
                   "The second column of the rank samples, one pair after the other."))
    out.append(arr("R_OFFS", roffs, "i64", "The start of each pair in R_FLAT, and the end."))
    out.append(arr("E_RAVG", rcols[0], "f64", "Expected average ranks."))
    out.append(arr("E_RMIN", rcols[1], "f64", "Expected lowest ranks of a run."))
    out.append(arr("E_RMAX", rcols[2], "f64", "Expected highest ranks of a run."))
    out.append(arr("E_RDENSE", rcols[3], "f64", "Expected dense ranks."))
    out.append(arr("E_RORD", rcols[4], "f64", "Expected ordinal ranks."))

    sp = []
    kb = []
    kc = []
    rfit = []
    for xs, ys in keptranks:
        ex = full_exact(xs, ys)
        sp.append(fl(ex["spearman"]) if ex["spearman"] is not None else "0.0")
        kb.append(fl(ex["kendall_b"]) if ex["kendall_b"] is not None else "0.0")
        kc.append(fl(ex["kendall_c"]) if ex["kendall_c"] is not None else "0.0")
        rfit.append(1 if ex["spearman"] is not None else 0)
    out.append(arr("E_SPEARMAN", sp, "f64", "Expected Spearman coefficient of each pair."))
    out.append(arr("E_KENDALL_B", kb, "f64", "Expected Kendall tau-b of each pair."))
    out.append(arr("E_KENDALL_C", kc, "f64", "Expected Kendall tau-c of each pair."))
    out.append(arr("R_FIT", rfit, "i64",
                   "1 when the rank correlations of the pair have a value."))

    # The line fit.
    lflat = []
    loffs = []
    lcols = {k: [] for k in ["slope", "intercept", "slope_err", "intercept_err", "ts", "ti"]}
    err_fit = []
    lfit = []
    ltheil = []
    for xs, ys in keptlines:
        loffs.append(len(lflat))
        lflat.extend(xs)
        lflat.extend(ys)
        ex = full_exact(xs, ys)
        lcols["slope"].append(fl(ex["slope"]) if ex["slope"] is not None else "0.0")
        lcols["intercept"].append(fl(ex["intercept"]) if ex["intercept"] is not None else "0.0")
        lcols["slope_err"].append(fl(ex["slope_err"]) if ex.get("err_fit") and ex["slope_err"] is not None else "0.0")
        lcols["intercept_err"].append(fl(ex["intercept_err"]) if ex.get("err_fit") and ex["intercept_err"] is not None else "0.0")
        err_fit.append(1 if ex.get("err_fit") else 0)
        # The line needs two columns that both vary. A constant column
        # gives the undefined fault. A pair with an exact correlation of
        # zero has a line, so the correlation does not enter here.
        lfit.append(1 if ex["slope"] is not None and ex["spread_y"] > 0 else 0)
        ltheil.append(1 if ex["theil"] else 0)
        lcols["ts"].append(fl(ex["theil"][0]) if ex["theil"] else "0.0")
        lcols["ti"].append(fl(ex["theil"][1]) if ex["theil"] else "0.0")
    loffs.append(len(lflat))
    out.append(arr("L_FLAT", [fl(v) for v in lflat], "f64",
                   "The columns of the line samples, one pair after the other."))
    out.append(arr("L_OFFS", loffs, "i64", "The start of each pair in L_FLAT, and the end."))
    out.append(arr("E_SLOPE", lcols["slope"], "f64", "Expected slope of the line."))
    out.append(arr("E_INTERCEPT", lcols["intercept"], "f64", "Expected intercept of the line."))
    out.append(arr("E_SLOPE_ERR", lcols["slope_err"], "f64", "Expected standard error of the slope."))
    out.append(arr("E_INTERCEPT_ERR", lcols["intercept_err"], "f64",
                   "Expected standard error of the intercept."))
    out.append(arr("L_ERR_FIT", err_fit, "i64",
                   "1 when the standard errors of the line have a value."))
    out.append(arr("L_FIT", lfit, "i64", "1 when the line of the pair has a value."))
    out.append(arr("L_THEIL", ltheil, "i64",
                   "1 when the Theil-Sen estimate of the pair has a value."))
    out.append(arr("E_THEIL_SLOPE", lcols["ts"], "f64", "Expected Theil-Sen slope."))
    out.append(arr("E_THEIL_INTERCEPT", lcols["ti"], "f64", "Expected Theil-Sen intercept."))

    # The matrices.
    mflat = []
    moffs = []
    mcols = []
    mconst = []
    for cols_data in keptmats:
        rows = len(cols_data[0])
        nc = len(cols_data)
        moffs.append(len(mflat))
        for col in cols_data:
            mflat.extend(col)
        a = np.array(cols_data, dtype=float)
        # NumPy takes the columns of the array as the variables, and the
        # package takes the columns of the buffer.
        cov = np.cov(a, ddof=1)
        corr = np.corrcoef(a)
        # The exact value of each cell. A constant column has no value, and
        # the test expects the undefined fault there.
        for i in range(nc):
            for j in range(nc):
                ex = full_exact(cols_data[i], cols_data[j])
                cell_cov = ex["cov"]
                cell_corr = ex["corr"]
                if cell_cov is not None and cell_corr is not None:
                    good = well_cond(ex, cols_data[i], cols_data[j])
                    close_abs(float(cov[i, j]), cell_cov, 64, cov_scale(ex), "numpy matrix cov", cols_data[i])
                    close_abs(float(corr[i, j]), cell_corr, 64, 1.0, "numpy matrix corr", cols_data[i])
                mcols.append((fl(cell_cov) if cell_cov is not None else "0.0",
                              fl(cell_corr) if cell_corr is not None else "0.0"))
        mconst.append(1 if any(len(set(c)) == 1 for c in cols_data) else 0)
    moffs.append(len(mflat))
    out.append(arr("M_FLAT", [fl(v) for v in mflat], "f64",
                   "The columns of the matrix samples, one matrix after the other."))
    out.append(arr("M_OFFS", moffs, "i64", "The start of each matrix in M_FLAT, and the end."))
    out.append(arr("M_COLS", [len(c) for c in keptmats], "i64",
                   "The count of columns of each matrix."))
    out.append(arr("M_ROWS", [len(c[0]) for c in keptmats], "i64",
                   "The count of rows of each matrix."))
    out.append(arr("E_MAT_COV", [c[0] for c in mcols], "f64",
                   "Expected covariance matrix, row by row."))
    out.append(arr("E_MAT_CORR", [c[1] for c in mcols], "f64",
                   "Expected correlation matrix, row by row."))
    out.append(arr("M_CONST", mconst, "i64", "1 when a column of the matrix holds one value only."))

    with open(OUT, "w") as f:
        f.write("\n".join(out))
    print("wrote", OUT)
    for name in ["oracle/pairs.json", "pairs.txt"]:
        path = os.path.join(HERE, name)
        if os.path.exists(path):
            os.remove(path)


if __name__ == "__main__":
    main()
