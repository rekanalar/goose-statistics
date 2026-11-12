#!/usr/bin/env python3
"""Writes test/suite/exact_vectors.goose: the numbers of SciPy 1.18.1 for the
exact and distance tests and the checks that use them.

The file holds the documented rows, seeded random tables and samples (with
ties), and the one-sample and two-sample Kolmogorov-Smirnov dispatchers.
The p-values of barnard_exact come from the optimizer of SciPy. The
generator checks each of them against a dense grid of the nuisance
parameter, and it drops and prints the rows where SciPy ends below the
true maximum. Run it from a scratch directory with the system python3.
The output path is the first argument, test/suite/exact_vectors.goose by
default.
"""
import math
import random
import sys
import itertools
import numpy as np
from scipy import stats
from scipy.special import gammaln
from scipy.stats import PermutationMethod

random.seed(20261112)
ALT = {"two-sided": "statistics::ALT_TWO_SIDED", "less": "statistics::ALT_LESS",
       "greater": "statistics::ALT_GREATER"}
FAM = {
    "normal": ("statistics::FAMILY_NORMAL", lambda a: stats.norm(a[0], a[1])),
    "uniform": ("statistics::FAMILY_UNIFORM", lambda a: stats.uniform(a[0], a[1] - a[0])),
    "expon": ("statistics::FAMILY_EXPONENTIAL", lambda a: stats.expon(scale=1.0 / a[0])),
    "lognorm": ("statistics::FAMILY_LOGNORMAL", lambda a: stats.lognorm(a[1], scale=math.exp(a[0]))),
    "weibull": ("statistics::FAMILY_WEIBULL", lambda a: stats.weibull_min(a[0], scale=a[1])),
    "gamma": ("statistics::FAMILY_GAMMA", lambda a: stats.gamma(a[0], scale=a[1])),
    "beta": ("statistics::FAMILY_BETA", lambda a: stats.beta(a[0], a[1])),
    "chi2": ("statistics::FAMILY_CHI2", lambda a: stats.chi2(a[0])),
    "t": ("statistics::FAMILY_T", lambda a: stats.t(a[0])),
    "f": ("statistics::FAMILY_F", lambda a: stats.f(a[0], a[1])),
}
body = []
skipped = []


def arr(v):
    return "[" + ", ".join(repr(float(x)) for x in v) + "]"


def ilist(v):
    return "[" + ", ".join(str(int(x)) for x in v) + "]"


def draw(kind, n):
    if kind == "ties":
        return [round(random.gauss(5, 2), 0) for _ in range(n)]
    if kind == "shift":
        return [random.gauss(6.5, 2) for _ in range(n)]
    if kind == "heavy":
        return [random.gauss(5, 1) * random.choice([1, 1, 1, 3]) for _ in range(n)]
    return [random.gauss(5, 2) for _ in range(n)]


def emit(text):
    body.append("        " + text)


def trow(call, st, p, ts, tp, name, pre=""):
    emit("{ %slet r, c = %s; xtrow(bad, r, c, %r, %r, %r, %r, \"%s\"); }"
         % (pre, call, float(st), float(p), ts, tp, name))


def dense_barnard(table, alt, pooled):
    t = np.asarray(table)
    c1, c2 = t.sum(axis=0)
    x1 = np.arange(c1 + 1).reshape(-1, 1)
    x2 = np.arange(c2 + 1).reshape(1, -1)
    p1, p2 = x1 / c1, x2 / c2
    if pooled:
        p = (x1 + x2) / (c1 + c2)
        v = p * (1 - p) * (1 / c1 + 1 / c2)
    else:
        v = p1 * (1 - p1) / c1 + p2 * (1 - p2) / c2
    with np.errstate(all="ignore"):
        T = np.divide(p1 - p2, np.sqrt(v))
    T[p1 == p2] = 0
    T0 = T[t[0, 0], t[0, 1]]
    idx = {"two-sided": np.abs(T) >= abs(T0), "less": T <= T0, "greater": T >= T0}[alt]

    def lc(n):
        return gammaln(n + 1) - gammaln(np.arange(n + 1) + 1) - gammaln(n - np.arange(n + 1) + 1)
    w1 = np.exp(lc(c1)[x1] + lc(c2)[x2])
    n = c1 + c2
    A = np.array([w1[idx & (x1 + x2 == s)].sum() for s in range(n + 1)])
    s = np.arange(n + 1)
    best = 0.0
    for pi in np.linspace(1e-6, 1 - 1e-6, 4001):
        best = max(best, float((A * pi ** s * (1 - pi) ** (n - s)).sum()))
    return best


def table_rows():
    rows = [(7, 3, 2, 8), (7, 12, 8, 3), (74, 31, 43, 32), (0, 5, 3, 4), (5, 0, 0, 6),
            (1, 1, 1, 1), (10, 0, 0, 10), (3, 9, 12, 4), (20, 15, 11, 25)]
    for _ in range(14):
        tot = random.choice([8, 12, 20, 30, 45, 60])
        a, b, c, d = [random.randint(0, max(1, tot // 3)) for _ in range(4)]
        rows.append((a, b, c, d))
    for (a, b, c, d) in rows:
        if a + c == 0 or b + d == 0:
            continue
        for alt in ("two-sided", "less", "greater"):
            for pooled in (True, False):
                r = stats.barnard_exact([[a, b], [c, d]], alternative=alt, pooled=pooled)
                st = float(r.statistic)
                if math.isinf(st) or math.isnan(st):
                    continue
                dm = dense_barnard([[a, b], [c, d]], alt, pooled)
                if dm > float(r.pvalue) * (1 + 2e-9) + 1e-300:
                    skipped.append(("barnard", a, b, c, d, alt, pooled, float(r.pvalue), dm))
                    continue
                trow("statistics::barnard_exact(%r, %r, %r, %r, %s, %s)"
                     % (float(a), float(b), float(c), float(d), ALT[alt], "true" if pooled else "false"),
                     st, r.pvalue, 1e-12, 1e-9, "barnard %d %d %d %d %s %s" % (a, b, c, d, alt, pooled))
            r = stats.boschloo_exact([[a, b], [c, d]], alternative=alt)
            trow("statistics::boschloo_exact(%r, %r, %r, %r, %s)" % (float(a), float(b), float(c), float(d), ALT[alt]),
                 r.statistic, r.pvalue, 1e-12, 1e-9, "boschloo %d %d %d %d %s" % (a, b, c, d, alt))


def sample_for(kind, n):
    if kind == "unif":
        return [random.uniform(0.5, 4) for _ in range(n)]
    if kind == "unit":
        return [random.betavariate(2, 3) for _ in range(n)]
    if kind == "pos":
        return [random.expovariate(1.0) * 2 + 0.05 for _ in range(n)]
    return draw(kind, n)


PARAMS = {"normal": [5.0, 2.0], "uniform": [0.0, 4.5], "expon": [0.5], "lognorm": [0.5, 0.8],
          "weibull": [1.5, 2.0], "gamma": [2.0, 1.5], "beta": [2.0, 3.0], "chi2": [3.0],
          "t": [4.0], "f": [3.0, 5.0]}


def cvm_rows():
    for fam, n, kind in [("normal", 8, "plain"), ("normal", 20, "plain"), ("normal", 40, "shift"),
                         ("expon", 15, "pos"), ("uniform", 12, "unif"), ("lognorm", 25, "pos"),
                         ("weibull", 18, "pos"), ("gamma", 30, "pos"), ("beta", 14, "unit"),
                         ("chi2", 16, "pos"), ("t", 22, "plain"), ("f", 12, "pos"), ("normal", 3, "plain"),
                         ("normal", 60, "plain"), ("normal", 25, "heavy")]:
        xs = sample_for(kind, n)
        params = PARAMS[fam]
        const, mk = FAM[fam]
        r = stats.cramervonmises(xs, mk(params).cdf)
        trow("statistics::cramervonmises(xs[..], %s, pp[..])" % const, r.statistic, r.pvalue, 1e-12, 1e-5,
             "cramervonmises %s %d" % (fam, n),
             "let xs: f64[%d] = %s; let pp: f64[%d] = %s; " % (n, arr(xs), len(params), arr(params)))


def cvm2_rows():
    for (nx, ny, kind) in [(2, 2, "plain"), (3, 5, "plain"), (7, 6, "plain"), (8, 8, "ties"), (10, 12, "shift"),
                           (12, 12, "plain"), (15, 9, "ties"), (20, 17, "plain"), (20, 19, "shift"),
                           (4, 20, "plain"), (25, 30, "plain"), (40, 35, "shift"), (60, 60, "ties")]:
        xs = draw("ties" if kind == "ties" else "plain", nx)
        ys = draw(kind, ny)
        for method in ("auto", "exact", "asymptotic"):
            if method == "exact" and max(nx, ny) > 20:
                continue
            r = stats.cramervonmises_2samp(xs, ys, method=method)
            const = {"auto": "statistics::MW_AUTO", "exact": "statistics::MW_EXACT",
                     "asymptotic": "statistics::MW_ASYMPTOTIC"}[method]
            tp = 1e-12 if (method == "exact" or (method == "auto" and max(nx, ny) <= 20)) else 1e-9
            trow("statistics::cramervonmises_2samp(xs[..], ys[..], %s)" % const, r.statistic, r.pvalue, 1e-12, tp,
                 "cvm2 %d %d %s %s" % (nx, ny, kind, method),
                 "let xs: f64[%d] = %s; let ys: f64[%d] = %s; " % (nx, arr(xs), ny, arr(ys)))


def epps_rows():
    for (nx, ny, kind, t) in [(8, 8, "plain", None), (5, 5, "plain", None), (6, 9, "ties", None),
                              (12, 14, "shift", None), (24, 24, "plain", None), (25, 30, "plain", None),
                              (40, 35, "shift", None), (60, 80, "heavy", None), (10, 10, "plain", [0.3, 0.6, 1.0]),
                              (15, 11, "ties", [0.5]), (30, 20, "shift", [0.2, 0.4, 0.8, 1.2])]:
        xs = draw("ties" if kind == "ties" else "plain", nx)
        ys = draw(kind, ny)
        r = stats.epps_singleton_2samp(xs, ys) if t is None else stats.epps_singleton_2samp(xs, ys, t=t)
        tl = "tt[..]" if t else "tt[..]"
        trow("statistics::epps_singleton_2samp(xs[..], ys[..], %s)" % tl, r.statistic, r.pvalue, 1e-7, 1e-6,
             "epps %d %d %s" % (nx, ny, kind),
             "let xs: f64[%d] = %s; let ys: f64[%d] = %s; let tt: f64[%d] = %s; "
             % (nx, arr(xs), ny, arr(ys), len(t) if t else 0, arr(t) if t else "[]"))


def bws_stat_ranks(rx, ry, alt):
    """The statistic of two sorted lists of ranks, as SciPy computes it."""
    n = len(rx)
    m = len(ry)
    bx = 0.0
    by = 0.0
    for i in range(n):
        k = i + 1
        num = rx[i] - (m + n) / n * k
        if alt == "two-sided":
            num = num * num
        else:
            num = num * abs(num)
        den = k / (n + 1) * (1 - k / (n + 1)) * m * (m + n) / n
        bx += num / den
    for j in range(m):
        k = j + 1
        num = ry[j] - (m + n) / m * k
        if alt == "two-sided":
            num = num * num
        else:
            num = num * abs(num)
        den = k / (m + 1) * (1 - k / (m + 1)) * n * (m + n) / m
        by += num / den
    bx = bx / n
    by = by / m
    if alt == "two-sided":
        return (bx + by) / 2
    return (bx - by) / 2


def bws_ranks(x, y):
    z = np.concatenate([x, y])
    r = stats.rankdata(z, method="average")
    return np.sort(r[:len(x)]), np.sort(r[len(x):]), np.sort(r)


def bws_exact(x, y, alt):
    """The exact count of every split of the rank set, as SciPy does it when
    the count of splits fits the number of resamples."""
    rx, ry, zs = bws_ranks(x, y)
    obs = bws_stat_ranks(rx, ry, alt)
    gamma = abs(100.0 * np.finfo(float).eps * obs)
    nx = len(rx)
    total = len(zs)
    count = 0
    seen = 0
    for idx in itertools.combinations(range(total), nx):
        mask = set(idx)
        ax = []
        ay = []
        for i in range(total):
            if i in mask:
                ax.append(zs[i])
            else:
                ay.append(zs[i])
        v = bws_stat_ranks(ax, ay, alt)
        seen += 1
        if alt == "less":
            if v <= obs + gamma:
                count += 1
        else:
            if v >= obs - gamma:
                count += 1
    return obs, count / seen


def bws_rows():
    for (nx, ny, kind) in [(2, 2, "plain"), (3, 4, "plain"), (5, 6, "ties"), (7, 7, "plain"), (8, 8, "shift"),
                           (6, 11, "plain"), (9, 9, "ties"), (4, 17, "shift")]:
        xs = draw("ties" if kind == "ties" else "plain", nx)
        ys = draw(kind, ny)
        for alt in ("two-sided", "less", "greater"):
            st, p = bws_exact(xs, ys, alt)
            trow("statistics::bws_test(xs[..], ys[..], rg, %s)" % ALT[alt], st, p, 1e-12, 1e-12,
                 "bws %d %d %s %s" % (nx, ny, kind, alt),
                 "let xs: f64[%d] = %s; let ys: f64[%d] = %s; var rg = rng { s: 7 }; " % (nx, arr(xs), ny, arr(ys)))
    xs = draw("plain", 10)
    ys = draw("shift", 10)
    st, p = bws_exact(xs, ys, "two-sided")
    emit("{ let xs: f64[10] = %s; let ys: f64[10] = %s; xbwsmc(bad, xs[..], ys[..], %r, %r); }"
         % (arr(xs), arr(ys), float(st), float(p)))


def anderson_rows():
    for sizes, kind in [([5, 6], "plain"), ([8, 8, 8], "plain"), ([10, 7, 12, 5], "shift"), ([6, 6], "ties"),
                        ([9, 9, 9], "ties"), ([20, 15], "shift"), ([4, 4, 4, 4, 4], "plain"), ([30, 40, 25], "shift"),
                        ([3, 3], "plain"), ([15, 12], "plain")]:
        samples = []
        flat = []
        for i, n in enumerate(sizes):
            s = draw("ties" if kind == "ties" else ("shift" if (kind == "shift" and i % 2) else "plain"), n)
            samples.append(s)
            flat += s
        k = len(sizes)
        m = k - 1
        crit = (np.array([0.675, 1.281, 1.645, 1.96, 2.326, 2.573, 3.085])
                + np.array([-0.245, 0.25, 0.678, 1.149, 1.822, 2.364, 3.615]) / math.sqrt(m)
                + np.array([-0.105, -0.305, -0.362, -0.391, -0.396, -0.345, -0.154]) / m)
        for variant, const in [("midrank", "statistics::KSAMP_MIDRANK"), ("right", "statistics::KSAMP_RIGHT"),
                               ("continuous", "statistics::KSAMP_CONTINUOUS")]:
            r = stats.anderson_ksamp(samples, variant=variant)
            emit("{ let g: f64[%d] = %s; let sz: i64[%d] = %s; let cr: f64[7] = %s; xarow(bad, g[..], sz[..], %s, %r, %r, cr[..], \"anderson %s %s\"); }"
                 % (len(flat), arr(flat), k, ilist(sizes), arr(crit), const, float(r.statistic), float(r.pvalue),
                    "x".join(str(x) for x in sizes), variant))


def ks_rows():
    for (nx, ny, kind) in [(8, 8, "plain"), (5, 9, "ties"), (12, 7, "shift"), (20, 20, "plain"), (30, 41, "shift"),
                           (3, 3, "ties")]:
        xs = draw("ties" if kind == "ties" else "plain", nx)
        ys = draw(kind, ny)
        # A tie inside a sample leaves no exact path. The package then uses
        # the asymptotic form, and so does SciPy under the name asymp.
        tied = len(set(xs)) != nx or len(set(ys)) != ny
        for alt in ("two-sided", "less", "greater"):
            r = stats.ks_2samp(xs, ys, alternative=alt, method="asymp" if tied else "auto")
            emit("{ let xs: f64[%d] = %s; let ys: f64[%d] = %s; xksrow(bad, xs[..], ys[..], %s, %r, %r, %r, %r, \"kstest_2samp %d %d %s %s\"); }"
                 % (nx, arr(xs), ny, arr(ys), ALT[alt], float(r.statistic), float(r.pvalue),
                    float(r.statistic_location), float(r.statistic_sign), nx, ny, kind, alt))
    cases = [("normal", 15, "plain"), ("normal", 40, "plain"), ("uniform", 10, "unif"), ("expon", 20, "pos"),
             ("lognorm", 12, "pos"), ("weibull", 17, "pos"), ("gamma", 25, "pos"), ("beta", 9, "unit"),
             ("chi2", 14, "pos"), ("t", 19, "plain"), ("f", 11, "pos")]
    for fam, n, kind in cases:
        params = PARAMS[fam]
        xs = sample_for(kind, n)
        const, mk = FAM[fam]
        for alt in ("two-sided", "less", "greater"):
            r = stats.ks_1samp(xs, mk(params).cdf, alternative=alt)
            emit("{ let xs: f64[%d] = %s; let pp: f64[%d] = %s; xks1row(bad, xs[..], %s, pp[..], %s, %r, %r, %r, %r, \"kstest %s %d %s\"); }"
                 % (n, arr(xs), len(params), arr(params), const, ALT[alt], float(r.statistic), float(r.pvalue),
                    float(r.statistic_location), float(r.statistic_sign), fam, n, alt))


def poisson_rows():
    rows = [(10, 2.0, 18, 2.5, 0.0), (0, 100.0, 3, 100.0, 0.0), (5, 1.0, 5, 1.0, 0.0), (0, 1.0, 0, 1.0, 0.0),
            (3, 1.0, 0, 2.0, 0.5), (12, 3.0, 30, 3.0, 0.0), (1, 0.5, 9, 1.5, 0.0)]
    for _ in range(14):
        rows.append((random.randint(0, 40), round(random.uniform(0.5, 5), 2), random.randint(0, 40),
                     round(random.uniform(0.5, 5), 2), random.choice([0.0, 0.0, 0.2, 1.0])))
    for (k1, n1, k2, n2, diff) in rows:
        for alt in ("two-sided", "less", "greater"):
            with np.errstate(all="ignore"):
                try:
                    r = stats.poisson_means_test(k1, n1, k2, n2, diff=diff, alternative=alt)
                except Exception:
                    continue
            st, p = float(r.statistic), float(r.pvalue)
            if math.isnan(st) or math.isinf(st) or math.isnan(p):
                continue
            trow("statistics::poisson_means_test(%r, %r, %r, %r, %r, %s)"
                 % (float(k1), float(n1), float(k2), float(n2), float(diff), ALT[alt]), st, p, 1e-12, 1e-8,
                 "poisson %d %r %d %r %r %s" % (k1, n1, k2, n2, diff, alt))


def page_rows():
    for (m, n, mode) in [(8, 3, "plain"), (5, 4, "plain"), (6, 5, "ties"), (12, 3, "plain"), (4, 8, "plain"),
                         (3, 6, "ties"), (20, 3, "plain"), (2, 7, "plain"), (7, 4, "ranked"), (6, 5, "pred"),
                         (15, 4, "plain"), (30, 5, "ties"), (3, 10, "plain"), (9, 12, "plain"), (2, 12, "ranked"),
                         (5, 9, "pred")]:
        rows = []
        for i in range(m):
            if mode in ("ranked", "pred"):
                perm = list(range(1, n + 1))
                random.shuffle(perm)
                rows.append([float(v) for v in perm])
            elif mode == "ties":
                rows.append([float(round(random.gauss(5 + 0.3 * j, 1.5), 0)) for j in range(n)])
            else:
                rows.append([random.gauss(5 + 0.4 * j, 1.5) for j in range(n)])
        data = np.array(rows)
        flat = [v for r in rows for v in r]
        pred = None
        if mode == "pred":
            pred = list(range(1, n + 1))
            random.shuffle(pred)
        for method, const in [("auto", "statistics::PAGE_AUTO"), ("asymptotic", "statistics::PAGE_ASYMPTOTIC"),
                              ("exact", "statistics::PAGE_EXACT")]:
            if method == "exact" and (n > 8 or m > 12):
                continue
            r = stats.page_trend_test(data, ranked=(mode in ("ranked", "pred")), predicted_ranks=pred, method=method)
            used = 1.0 if r.method == "exact" else 2.0
            emit("{ let d: f64[%d] = %s; let pr: f64[%d] = %s; xprow(bad, d[..], %d, %s, pr[..], %s, %r, %r, %r, \"page %d %d %s %s\"); }"
                 % (m * n, arr(flat), len(pred) if pred else 0, arr(pred) if pred else "[]", n,
                    "true" if mode in ("ranked", "pred") else "false", const, float(r.statistic), float(r.pvalue), used,
                    m, n, mode, method))


ORACLE_A = [2.1, 3.4, 1.9, 5.6, 4.4, 3.8, 2.7, 4.9]
ORACLE_B = [3.9, 5.2, 4.8, 6.1, 5.5, 4.1, 6.6, 5.0]
ORACLE_C = [1.2, 2.5, 2.2, 3.1, 1.8, 2.9, 2.4, 3.3]


def ksamp_crit(m):
    """The seven critical values of the k-sample Anderson-Darling test."""
    return (np.array([0.675, 1.281, 1.645, 1.96, 2.326, 2.573, 3.085])
            + np.array([-0.245, 0.25, 0.678, 1.149, 1.822, 2.364, 3.615]) / math.sqrt(m)
            + np.array([-0.105, -0.305, -0.362, -0.391, -0.396, -0.345, -0.154]) / m)


def oracle_rows():
    """The measured rows of the documentation, one per call."""
    a = np.array(ORACLE_A)
    b = np.array(ORACLE_B)
    c = np.array(ORACLE_C)
    np2 = [3.5, 1.5]
    # The exact tests of the 2 by 2 table.
    for alt in ("two-sided", "less", "greater"):
        r = stats.barnard_exact([[7, 3], [2, 8]], alternative=alt)
        trow("statistics::barnard_exact(7.0, 3.0, 2.0, 8.0, %s)" % ALT[alt],
             float(r.statistic), float(r.pvalue), 1e-12, 1e-9, "barnard oracle %s" % alt)
        r = stats.boschloo_exact([[7, 3], [2, 8]], alternative=alt)
        trow("statistics::boschloo_exact(7.0, 3.0, 2.0, 8.0, %s)" % ALT[alt],
             float(r.statistic), float(r.pvalue), 1e-12, 1e-9, "boschloo oracle %s" % alt)
    # The one-sample Cramer-von Mises test.
    r = stats.cramervonmises(a, stats.norm(np2[0], np2[1]).cdf)
    trow("statistics::cramervonmises(xs[..], statistics::FAMILY_NORMAL, pp[..])",
         float(r.statistic), float(r.pvalue), 1e-12, 1e-9, "cramervonmises oracle",
         "let xs: f64[%d] = %s; let pp: f64[2] = %s; " % (len(a), arr(a), arr(np2)))
    # The two-sample Cramer-von Mises test.
    for method, const in [("exact", "statistics::MW_EXACT"), ("asymptotic", "statistics::MW_ASYMPTOTIC")]:
        r = stats.cramervonmises_2samp(a, b, method=method)
        trow("statistics::cramervonmises_2samp(xs[..], ys[..], %s)" % const,
             float(r.statistic), float(r.pvalue), 1e-12, 1e-12, "cvm2 oracle %s" % method,
             "let xs: f64[%d] = %s; let ys: f64[%d] = %s; " % (len(a), arr(a), len(b), arr(b)))
    # The Epps-Singleton test.
    r = stats.epps_singleton_2samp(a, b)
    trow("statistics::epps_singleton_2samp(xs[..], ys[..], tt[..])",
         float(r.statistic), float(r.pvalue), 1e-12, 1e-9, "epps oracle",
         "let xs: f64[%d] = %s; let ys: f64[%d] = %s; let tt: f64[0] = []; "
         % (len(a), arr(a), len(b), arr(b)))
    # The Baumgartner-Weiss-Schindler test. The exact count needs more
    # resamples than the default, so the enumeration runs.
    for alt in ("two-sided", "less", "greater"):
        st, p = bws_exact(a, b, alt)
        trow("statistics::bws_test(xs[..], ys[..], rg, %s)" % ALT[alt], st, p, 1e-12, 1e-12,
             "bws oracle %s" % alt,
             "let xs: f64[%d] = %s; let ys: f64[%d] = %s; var rg = rng { s: 7 }; "
             % (len(a), arr(a), len(b), arr(b)))
    # The k-sample Anderson-Darling test.
    samples = [a, b, c]
    flat = list(a) + list(b) + list(c)
    for variant, const in [("midrank", "statistics::KSAMP_MIDRANK"), ("right", "statistics::KSAMP_RIGHT"),
                           ("continuous", "statistics::KSAMP_CONTINUOUS")]:
        r = stats.anderson_ksamp(samples, variant=variant)
        emit("{ let g: f64[%d] = %s; let sz: i64[3] = [8, 8, 8]; let cr: f64[7] = %s; "
             "xarow(bad, g[..], sz[..], %s, %r, %r, cr[..], \"anderson oracle %s\"); }"
             % (len(flat), arr(flat), arr(ksamp_crit(2)), const, float(r.statistic), float(r.pvalue), variant))
    # The Kolmogorov-Smirnov dispatchers.
    r = stats.ks_2samp(a, b)
    emit("{ let xs: f64[%d] = %s; let ys: f64[%d] = %s; xksrow(bad, xs[..], ys[..], statistics::ALT_TWO_SIDED, "
         "%r, %r, %r, %r, \"kstest_2samp oracle\"); }"
         % (len(a), arr(a), len(b), arr(b), float(r.statistic), float(r.pvalue),
            float(r.statistic_location), float(r.statistic_sign)))
    r = stats.ks_1samp(a, stats.norm(np2[0], np2[1]).cdf)
    emit("{ let xs: f64[%d] = %s; let pp: f64[2] = %s; xks1row(bad, xs[..], statistics::FAMILY_NORMAL, pp[..], "
         "statistics::ALT_TWO_SIDED, %r, %r, %r, %r, \"kstest oracle\"); }"
         % (len(a), arr(a), arr(np2), float(r.statistic), float(r.pvalue),
            float(r.statistic_location), float(r.statistic_sign)))
    # The Poisson means test.
    r = stats.poisson_means_test(10, 2.0, 18, 2.5)
    trow("statistics::poisson_means_test(10.0, 2.0, 18.0, 2.5)",
         float(r.statistic), float(r.pvalue), 1e-12, 1e-9, "poisson_means oracle")
    # The Page trend test. The rows hold the three samples.
    data = np.array([a, b, c]).T
    flat = [v for row in data for v in row]
    for method, const in [("exact", "statistics::PAGE_EXACT"), ("asymptotic", "statistics::PAGE_ASYMPTOTIC")]:
        r = stats.page_trend_test(data, method=method)
        used = 1.0 if r.method == "exact" else 2.0
        emit("{ let d: f64[%d] = %s; let pr: f64[0] = []; xprow(bad, d[..], 3, false, pr[..], %s, %r, %r, %r, "
             "\"page oracle %s\"); }"
             % (len(flat), arr(flat), const, float(r.statistic), float(r.pvalue), used, method))


def main():
    import time; t0 = time.time()
    oracle_rows()
    table_rows(); print("table", time.time() - t0, flush=True)
    cvm_rows(); print("cvm_rows", time.time() - t0, flush=True)
    cvm2_rows(); print("cvm2_rows", time.time() - t0, flush=True)
    epps_rows(); print("epps_rows", time.time() - t0, flush=True)
    bws_rows(); print("bws_rows", time.time() - t0, flush=True)
    anderson_rows(); print("anderson_rows", time.time() - t0, flush=True)
    ks_rows(); print("ks_rows", time.time() - t0, flush=True)
    poisson_rows(); print("poisson_rows", time.time() - t0, flush=True)
    page_rows(); print("page_rows", time.time() - t0, flush=True)
    head = '''// The answers of SciPy 1.18.1 for the exact and distance tests. The
// generator test/make_exact_vectors.py writes this file. It holds the
// documented rows, seeded random tables and samples with ties, and the
// answers of the Kolmogorov-Smirnov dispatchers.
namespace suite;

import std;
import statistics from "../..";

fn xtrow(bad: i64&, rc: statistics::TestResult, c: i64, stat: f64, p: f64, ts: f64, tp: f64, name: u8[:]) {
    if c == 0 && near_tol(rc.statistic, stat, ts) && near_tol(rc.pvalue, p, tp) { return; }
    bad += 1;
    if bad <= 12 { print("FAIL ", name, " got ", rc.statistic, " ", rc.pvalue, " code ", c, " want ", stat, " ", p, "\\n"); }
}

// The seeded draw of bws_test above the exact bound: the same seed gives the same
// p-value, and the p-value lies within 0.03 of the exact one.
fn xbwsmc(bad: i64&, xs: const f64[:], ys: const f64[:], stat: f64, p: f64) {
    var r1 = rng { s: 11 };
    var r2 = rng { s: 11 };
    let a, ca = statistics::bws_test(xs, ys, r1);
    let b, cb = statistics::bws_test(xs, ys, r2);
    if ca == 0 && cb == 0 && a.pvalue == b.pvalue && near_tol(a.statistic, stat, 1e-12) && abs(a.pvalue - p) < 0.03 { return; }
    bad += 1;
    if bad <= 12 { print("FAIL bws seeded draw got ", a.pvalue, " ", b.pvalue, " want ", p, "\\n"); }
}

fn xarow(bad: i64&, g: const f64[:], sz: const i64[:], variant: i64, stat: f64, p: f64, crit: const f64[:], name: u8[:]) {
    var out: f64[>..] = [];
    let r, c = statistics::anderson_ksamp(out, g, sz, variant);
    if c == 0 && near_tol(r.statistic, stat, 1e-11) && near_tol(r.pvalue, p, 1e-9) && near_list(out[..], crit, 1e-13) { return; }
    bad += 1;
    if bad <= 12 { print("FAIL ", name, " got ", r.statistic, " ", r.pvalue, " code ", c, " want ", stat, " ", p, "\\n"); }
}

fn xksrow(bad: i64&, xs: const f64[:], ys: const f64[:], alt: i64, stat: f64, p: f64, loc: f64, sign: f64, name: u8[:]) {
    var out: f64[>..] = [];
    let r, c = statistics::kstest_2samp(out, xs, ys, alt);
    if c == 0 && near_tol(r.statistic, stat, 1e-12) && near_tol(r.pvalue, p, 1e-12) && out.len == 2
        && out[0] == loc && out[1] == sign && is_nan(r.df) && is_nan(r.df2) { return; }
    bad += 1;
    if bad <= 12 { print("FAIL ", name, " got ", r.statistic, " ", r.pvalue, " code ", c, " want ", stat, " ", p, "\\n"); }
}

fn xks1row(bad: i64&, xs: const f64[:], family: i64, pp: const f64[:], alt: i64, stat: f64, p: f64, loc: f64, sign: f64, name: u8[:]) {
    var out: f64[>..] = [];
    let r, c = statistics::kstest(out, xs, family, pp, alt);
    if c == 0 && near_tol(r.statistic, stat, 1e-12) && near_tol(r.pvalue, p, 1e-12) && out.len == 2
        && out[0] == loc && out[1] == sign && is_nan(r.df) && is_nan(r.df2) { return; }
    bad += 1;
    if bad <= 12 { print("FAIL ", name, " got ", r.statistic, " ", r.pvalue, " code ", c, " want ", stat, " ", p, "\\n"); }
}

fn xprow(bad: i64&, d: const f64[:], columns: i64, ranked: bool, pr: const f64[:], method: i64, stat: f64, p: f64, used: f64, name: u8[:]) {
    let r, c = statistics::page_trend_test(d, columns, ranked, pr, method);
    if c == 0 && near_tol(r.statistic, stat, 1e-12) && near_tol(r.pvalue, p, 1e-11) && r.df == used { return; }
    bad += 1;
    if bad <= 12 { print("FAIL ", name, " got ", r.statistic, " ", r.pvalue, " df ", r.df, " code ", c, " want ", stat, " ", p, "\\n"); }
}

fn run_exact_vectors() -> i64 {
    var bad: i64 = 0;
    {
'''
    # the helper block groups, in blocks of rows, so that no function is too large
    text = head + "\n".join(body) + "\n    }\n    return bad;\n}\n"
    path = sys.argv[1] if len(sys.argv) > 1 else "test/suite/exact_vectors.goose"
    with open(path, "w") as f:
        f.write(text)
    print("wrote %s: %d rows, %d lines" % (path, len(body), text.count("\n")))
    for s in skipped:
        print("SciPy below the dense grid, row dropped:", s)


if __name__ == "__main__":
    main()
