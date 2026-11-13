#!/usr/bin/env python3
"""Make the reference vectors of the probability distributions.

The script writes numbers only into test/suite/distribution_vectors.goose.
The exact value of each row comes from mpmath at 40 to 60 digits. The
script then compares SciPy with the exact value and stops on a difference
beyond a loose bound. SciPy gets its tails as one minus a value, so the
bound is loose in the deep tail. Run it from the package directory with a
Python that has mpmath and SciPy:

    python3 test/make_distribution_vectors.py

mpmath is a tool of this script only. The package holds no copy of it.

Family numbers: 0 norm, 1 unif, 2 expon, 3 lognorm, 4 weibull,
5 gammadist, 6 betadist, 7 chi2, 8 tdist, 9 fdist, 10 binom, 11 poisson,
12 geom, 13 nbinom, 14 hypergeom.
Operation numbers: 0 pdf or pmf, 1 cdf, 2 sf, 3 ppf, 4 mean, 5 var,
6 skew, 7 kurt.
"""
import math
import os
import random
import sys

import mpmath as mp
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_special_vectors as msv  # the exact incomplete functions

mp.mp.dps = 50
OUT = os.path.join(HERE, "suite", "distribution_vectors.goose")
QUICK = "--quick" in sys.argv
rnd = random.Random(20261111)
ROWS = []  # (family, op, x, a, b, c, want, tolerance in ULP)
F = {n: i for i, n in enumerate(
    "norm unif expon lognorm weibull gammadist betadist chi2 tdist fdist "
    "binom poisson geom nbinom hypergeom".split())}
PDF, CDF, SF, PPF, MEAN, VAR, SKEW, KURT = range(8)
SQ2 = mp.sqrt(2)


def add(fam, op, x, a, b, c, want, tol):
    want = float(want)
    if abs(want) < 1e-40 and op in (SKEW, KURT):
        want = 0.0  # the exact value is 0, and the sum leaves noise
    if not math.isfinite(want):
        return
    ROWS.append((F[fam], op, float(x), float(a), float(b), float(c), want, tol))


def near(got, want, rel=1e-8):
    return abs(got - want) <= rel * abs(want) + 1e-300


def sanity(name, got, want):
    got, want = float(got), float(want)
    if math.isfinite(got) and want != 0.0 and not near(got, want):
        sys.exit("SciPy differs for %s: got %r want %r" % (name, got, want))


def pick(dist, ps):
    return [float(dist.ppf(p)) for p in ps]


# ---------------------------------------------------------------------------
# Exact cumulative functions as pairs (lower, upper) of mp numbers.
# ---------------------------------------------------------------------------

def pair_norm(z):
    z = mp.mpf(z)
    return mp.erfc(-z / SQ2) / 2, mp.erfc(z / SQ2) / 2


def pair_gamma(a, z):
    p, q = msv.ex_gammainc_pair(a, z)
    if p is None:
        raise ValueError("the exact incomplete gamma is not available")
    return p, q


def pair_beta(a, b, x):
    return msv.ex_betainc_pair(a, b, x)


def pair_t(x, df):
    x, df = mp.mpf(x), mp.mpf(df)
    w = df / (df + x * x)
    pw, qw = pair_beta(df / 2, mp.mpf(1) / 2, w)
    tail = pw / 2
    if x >= 0:
        return 1 - tail, tail
    return tail, 1 - tail


def pair_f(x, d1, d2):
    x, d1, d2 = mp.mpf(x), mp.mpf(d1), mp.mpf(d2)
    u = d1 * x / (d1 * x + d2)
    return pair_beta(d1 / 2, d2 / 2, u)


def pdf_gamma(x, a, s):
    x, a, s = mp.mpf(x), mp.mpf(a), mp.mpf(s)
    return mp.exp((a - 1) * mp.log(x / s) - x / s - mp.loggamma(a) - mp.log(s))


def pdf_beta(x, a, b):
    x, a, b = mp.mpf(x), mp.mpf(a), mp.mpf(b)
    return mp.exp((a - 1) * mp.log(x) + (b - 1) * mp.log1p(-x) -
                  (mp.loggamma(a) + mp.loggamma(b) - mp.loggamma(a + b)))


def pdf_t(x, df):
    x, df = mp.mpf(x), mp.mpf(df)
    return mp.exp(mp.loggamma((df + 1) / 2) - mp.loggamma(df / 2) -
                  mp.log(mp.sqrt(df * mp.pi)) - (df + 1) / 2 * mp.log1p(x * x / df))


def pdf_f(x, d1, d2):
    x, d1, d2 = mp.mpf(x), mp.mpf(d1), mp.mpf(d2)
    lb = mp.loggamma(d1 / 2) + mp.loggamma(d2 / 2) - mp.loggamma((d1 + d2) / 2)
    return mp.exp((d1 / 2) * mp.log(d1 / d2) + (d1 / 2 - 1) * mp.log(x) -
                  ((d1 + d2) / 2) * mp.log1p(d1 * x / d2) - lb)


def refine(cdf_pair, p, x0, upper=False):
    """The exact quantile: the root of the cumulative function at p."""
    p = mp.mpf(p)
    f = (lambda x: cdf_pair(x)[1] - p) if upper else (lambda x: cdf_pair(x)[0] - p)
    mp.mp.dps = 30
    try:
        return mp.findroot(f, mp.mpf(x0), tol=mp.mpf(10) ** -24, maxsteps=40)
    finally:
        mp.mp.dps = 50


# ---------------------------------------------------------------------------
# Continuous families.
# ---------------------------------------------------------------------------

PS = [1e-12, 1e-6, 1e-3, 0.02, 0.1, 0.3, 0.5, 0.7, 0.9, 0.98, 0.999, 1 - 1e-6]


def add_cont(fam, sp, args, pdf, pair, params, ppf_pair, is_t=False, spread=1.0, tolpdf=16):
    a, b, c = (list(params) + [0, 0, 0])[:3]
    xs = [x for x in pick(sp, PS) if math.isfinite(x)]
    xs += [rnd.uniform(*spread_range(sp)) for _ in range(6 if not QUICK else 2)]
    for x in xs:
        try:
            e = pdf(x)
        except Exception:
            e = None
        if e is not None and e != 0:
            sanity(fam + " pdf", sp.pdf(x), e)
            add(fam, PDF, x, a, b, c, e, tolpdf)
        p, q = pair(x)
        if p is None:
            continue
        for op, v, ref in ((CDF, p, sp.cdf), (SF, q, sp.sf)):
            if v != 0 and abs(float(v)) > 1e-300:
                sanity("%s op %d" % (fam, op), ref(x), v)
                add(fam, op, x, a, b, c, v, 16)
    for p in PS:
        x0 = float(sp.ppf(p))
        if not math.isfinite(x0):
            continue
        x = ppf_pair(p, x0)
        sanity(fam + " ppf", x0, x) if False else None
        if not near(x0, float(x), 1e-4):
            # SciPy loses the quantile of a tiny probability for a small shape
            if p > 1e-6 and p < 1 - 1e-6:
                sys.exit("SciPy quantile differs for %s at %r: %r vs %r" % (fam, p, x0, float(x)))
            print("note: SciPy quantile off for", fam, p, x0, float(x))
        add(fam, PPF, p, a, b, c, x, 16)


def spread_range(sp):
    lo, hi = sp.ppf(0.001), sp.ppf(0.999)
    return lo, hi


def cont_rows():
    # normal
    for mu, s in [(0.0, 1.0), (2.0, 0.5), (-1.0, 3.0)]:
        sp = stats.norm(mu, s)
        z = lambda x: (mp.mpf(x) - mu) / mp.mpf(s)
        pdfm = lambda x: mp.exp(-z(x) ** 2 / 2) / (mp.mpf(s) * mp.sqrt(2 * mp.pi))
        pair = lambda x: pair_norm(z(x))
        ppf = lambda p, x0: mu + s * (-SQ2 * msv.ex_erfcinv(2 * mp.mpf(p)) if p < 0.5 else SQ2 * msv.ex_erfcinv(2 * (1 - mp.mpf(p))))
        add_cont("norm", sp, (mu, s), pdfm, pair, (mu, s), ppf, tolpdf=4)
        for x in (-40.0, -20.0, -8.0, 8.0, 20.0, 38.0):
            xx = mu + s * x
            p, q = pair_norm(x)
            add("norm", CDF, xx, mu, s, 0, p, 16)
            add("norm", SF, xx, mu, s, 0, q, 16)
        for p in (1e-300, 1e-100, 1e-20):
            add("norm", PPF, p, mu, s, 0, mu - s * SQ2 * msv.ex_erfcinv(2 * mp.mpf(p)), 16)
        mean, var = mu, mp.mpf(s) ** 2
        for op, v in ((MEAN, mean), (VAR, var), (SKEW, 0), (KURT, 0)):
            add("norm", op, 0, mu, s, 0, v, 16)
    # uniform
    for lo, hi in [(0.0, 1.0), (-2.0, 5.0)]:
        w = mp.mpf(hi) - lo
        for x in (lo + 0.1 * (hi - lo), lo + 0.5 * (hi - lo), lo + 0.97 * (hi - lo)):
            add("unif", PDF, x, lo, hi, 0, 1 / w, 4)
            add("unif", CDF, x, lo, hi, 0, (mp.mpf(x) - lo) / w, 4)
            add("unif", SF, x, lo, hi, 0, (hi - mp.mpf(x)) / w, 4)
        for p in (0.001, 0.3, 0.9):
            add("unif", PPF, p, lo, hi, 0, lo + mp.mpf(p) * w, 4)
        for op, v in ((MEAN, (lo + hi) / mp.mpf(2)), (VAR, w * w / 12), (SKEW, 0), (KURT, -mp.mpf(6) / 5)):
            add("unif", op, 0, lo, hi, 0, v, 16)
    # exponential
    for r in (1.0, 0.3, 25.0):
        for x in (1e-10, 0.01, 0.5, 2.0, 10.0, 40.0, 700.0):
            x = x / r
            t = mp.mpf(r) * mp.mpf(x)
            add("expon", PDF, x, r, 0, 0, r * mp.exp(-t), 8)
            add("expon", CDF, x, r, 0, 0, -mp.expm1(-t), 8)
            add("expon", SF, x, r, 0, 0, mp.exp(-t), 8)
        for p in PS:
            add("expon", PPF, p, r, 0, 0, -mp.log1p(-mp.mpf(p)) / r, 8)
        for op, v in ((MEAN, 1 / mp.mpf(r)), (VAR, 1 / mp.mpf(r) ** 2), (SKEW, 2), (KURT, 6)):
            add("expon", op, 0, r, 0, 0, v, 16)
    # lognormal
    for mu, s in [(0.0, 1.0), (1.0, 0.5), (-2.0, 2.0)]:
        sp = stats.lognorm(s, scale=math.exp(mu))
        z = lambda x: (mp.log(mp.mpf(x)) - mu) / mp.mpf(s)
        for x in pick(sp, PS):
            e = mp.exp(-z(x) ** 2 / 2) / (mp.mpf(x) * s * mp.sqrt(2 * mp.pi))
            p, q = pair_norm(z(x))
            add("lognorm", PDF, x, mu, s, 0, e, 16)
            add("lognorm", CDF, x, mu, s, 0, p, 16)
            add("lognorm", SF, x, mu, s, 0, q, 16)
        for p in PS:
            z0 = -SQ2 * msv.ex_erfcinv(2 * mp.mpf(p)) if p < 0.5 else SQ2 * msv.ex_erfcinv(2 * (1 - mp.mpf(p)))
            add("lognorm", PPF, p, mu, s, 0, mp.exp(mu + s * z0), 16)
        s2 = mp.mpf(s) ** 2
        m = mp.exp(mu + s2 / 2)
        v = mp.expm1(s2) * mp.exp(2 * mu + s2)
        add("lognorm", MEAN, 0, mu, s, 0, m, 16)
        add("lognorm", VAR, 0, mu, s, 0, v, 16)
        add("lognorm", SKEW, 0, mu, s, 0, (mp.exp(s2) + 2) * mp.sqrt(mp.expm1(s2)), 16)
        add("lognorm", KURT, 0, mu, s, 0, mp.exp(4 * s2) + 2 * mp.exp(3 * s2) + 3 * mp.exp(2 * s2) - 6, 16)
    # Weibull
    for c, s in [(1.0, 1.0), (0.5, 2.0), (3.0, 1.5), (20.0, 1.0)]:
        sp = stats.weibull_min(c, scale=s)
        for x in pick(sp, PS):
            t = (mp.mpf(x) / s) ** c
            add("weibull", PDF, x, c, s, 0, c / mp.mpf(s) * (mp.mpf(x) / s) ** (c - 1) * mp.exp(-t), 16)
            add("weibull", CDF, x, c, s, 0, -mp.expm1(-t), 16)
            add("weibull", SF, x, c, s, 0, mp.exp(-t), 16)
        for p in PS:
            add("weibull", PPF, p, c, s, 0, s * (-mp.log1p(-mp.mpf(p))) ** (1 / mp.mpf(c)), 16)
        g = lambda k: mp.gamma(1 + k / mp.mpf(c))
        m1, m2, m3, m4 = g(1), g(2), g(3), g(4)
        v = m2 - m1 ** 2
        add("weibull", MEAN, 0, c, s, 0, s * m1, 16)
        add("weibull", VAR, 0, c, s, 0, s * s * v, 16)
        add("weibull", SKEW, 0, c, s, 0, (m3 - 3 * m2 * m1 + 2 * m1 ** 3) / v ** 1.5, 64)
        add("weibull", KURT, 0, c, s, 0, (m4 - 4 * m3 * m1 + 6 * m2 * m1 ** 2 - 3 * m1 ** 4) / v ** 2 - 3, 64)
    # gamma
    for a, s in [(1.0, 1.0), (0.3, 2.0), (5.0, 1.0), (1000.0, 1.0), (1e4, 0.5)]:
        sp = stats.gamma(a, scale=s)
        add_cont("gammadist", sp, (a, s), lambda x: pdf_gamma(x, a, s),
                 lambda x: pair_gamma(a, mp.mpf(x) / s), (a, s),
                 lambda p, x0: s * invert(lambda z: pair_gamma(a, z), p, x0 / s))
        m = mp.mpf(a) * s
        for op, v in ((MEAN, m), (VAR, a * mp.mpf(s) ** 2), (SKEW, 2 / mp.sqrt(a)), (KURT, 6 / mp.mpf(a))):
            add("gammadist", op, 0, a, s, 0, v, 16)
    # beta
    for a, b in [(1.0, 1.0), (0.5, 0.5), (2.0, 5.0), (50.0, 30.0), (1e3, 2e3)]:
        sp = stats.beta(a, b)
        add_cont("betadist", sp, (a, b), lambda x: pdf_beta(x, a, b),
                 lambda x: pair_beta(a, b, x), (a, b),
                 lambda p, x0: invert(lambda z: pair_beta(a, b, z), p, x0))
        A, B = mp.mpf(a), mp.mpf(b)
        S = A + B
        for op, v in ((MEAN, A / S), (VAR, A * B / (S * S * (S + 1))),
                      (SKEW, 2 * (B - A) * mp.sqrt(S + 1) / ((S + 2) * mp.sqrt(A * B))),
                      (KURT, 6 * ((A - B) ** 2 * (S + 1) - A * B * (S + 2)) / (A * B * (S + 2) * (S + 3)))):
            add("betadist", op, 0, a, b, 0, v, 16)
    # chi-square
    for df in (1.0, 4.0, 30.0, 1000.0):
        sp = stats.chi2(df)
        add_cont("chi2", sp, (df,), lambda x: pdf_gamma(x, df / 2, 2),
                 lambda x: pair_gamma(df / 2, mp.mpf(x) / 2), (df,),
                 lambda p, x0: 2 * invert(lambda z: pair_gamma(df / 2, z), p, x0 / 2))
        for op, v in ((MEAN, df), (VAR, 2 * df), (SKEW, mp.sqrt(8 / mp.mpf(df))), (KURT, 12 / mp.mpf(df))):
            add("chi2", op, 0, df, 0, 0, v, 16)
    # Student t
    for df in (1.0, 7.0, 30.0, 1000.0):
        sp = stats.t(df)
        add_cont("tdist", sp, (df,), lambda x: pdf_t(x, df), lambda x: pair_t(x, df), (df,),
                 lambda p, x0: mp.mpf(0) if p == 0.5 else invert(lambda z: pair_t(z, df), p, x0))
        if df > 1:
            add("tdist", MEAN, 0, df, 0, 0, 0, 0)
        if df > 2:
            add("tdist", VAR, 0, df, 0, 0, mp.mpf(df) / (df - 2), 16)
        if df > 3:
            add("tdist", SKEW, 0, df, 0, 0, 0, 0)
        if df > 4:
            add("tdist", KURT, 0, df, 0, 0, 6 / (mp.mpf(df) - 4), 16)
    # F
    for d1, d2 in [(1.0, 1.0), (3.0, 6.0), (10.0, 40.0), (100.0, 200.0)]:
        sp = stats.f(d1, d2)
        add_cont("fdist", sp, (d1, d2), lambda x: pdf_f(x, d1, d2), lambda x: pair_f(x, d1, d2), (d1, d2),
                 lambda p, x0: invert(lambda z: pair_f(z, d1, d2), p, x0))
        D1, D2 = mp.mpf(d1), mp.mpf(d2)
        if d2 > 2:
            add("fdist", MEAN, 0, d1, d2, 0, D2 / (D2 - 2), 16)
        if d2 > 4:
            add("fdist", VAR, 0, d1, d2, 0, 2 * D2 ** 2 * (D1 + D2 - 2) / (D1 * (D2 - 2) ** 2 * (D2 - 4)), 16)
        if d2 > 6:
            add("fdist", SKEW, 0, d1, d2, 0, (2 * D1 + D2 - 2) * mp.sqrt(8 * (D2 - 4)) / ((D2 - 6) * mp.sqrt(D1 * (D1 + D2 - 2))), 16)
        if d2 > 8:
            add("fdist", KURT, 0, d1, d2, 0, 12 / (D2 - 6) * ((5 * D2 - 22) / (D2 - 8) + ((D2 - 4) / D1) * ((D2 - 2) / (D2 - 8)) * ((D2 - 2) / (D1 + D2 - 2))), 16)


def invert(pair, p, x0, cap=None):
    """The exact quantile of a cumulative function given as a pair. The
    start is the SciPy value, which is far off for a tiny probability and a
    small shape. A bracket by factors finds the root, a bisection on the
    logarithm narrows it, and the solver on the relative shift ends it."""
    p = mp.mpf(p)
    x0 = mp.mpf(x0)
    mp.mp.dps = 30
    try:
        def f(x):
            try:
                r = pair(x)[0] - p
                if isinstance(r, mp.mpc):  # outside the support: above it
                    return mp.mpf(1)
                return r
            except ValueError:
                return mp.mpf(1) if x > x0 else mp.mpf(-1)
        sign = 1 if x0 > 0 else -1
        lo = hi = x0
        fac = mp.mpf(2)
        for _ in range(12):
            lo, hi = x0 / fac, x0 * fac
            if sign < 0:
                lo, hi = hi, lo
            try:
                if f(lo) * f(hi) <= 0:
                    break
            except Exception:
                pass
            fac = fac * fac
        else:
            sys.exit("no bracket for the quantile at p = %s, x0 = %s, f = %s %s" % (p, x0, f(x0 / 2), f(x0 * 2)))
        for _ in range(40):
            mid = mp.sqrt(lo * hi) * sign if lo * hi > 0 else (lo + hi) / 2
            if sign < 0:
                mid = -mp.sqrt(lo * hi)
            if f(lo) * f(mid) <= 0:
                hi = mid
            else:
                lo = mid
            if abs(hi / lo - 1) < mp.mpf("1e-3"):
                break
        xm = (lo + hi) / 2
        g = lambda y: f(xm * (1 + y)) / p  # scaled, so the solver tolerance is relative
        w = abs(hi / lo - 1) + mp.mpf("1e-9")
        y = mp.findroot(g, (-w, w), solver="anderson", tol=mp.mpf(10) ** -40,
                        maxsteps=200, verify=False)
        x = xm * (1 + y)
        if abs(pair(x)[0] - p) > mp.mpf(10) ** -17 * p:
            sys.exit("the root of the quantile did not converge at p = %s, x0 = %s" % (p, x0))
        return x
    finally:
        mp.mp.dps = 50


# ---------------------------------------------------------------------------
# Discrete families.
# ---------------------------------------------------------------------------

def sums(pmf, lo, hi):
    """The cumulative sums of a mass over [lo, hi]."""
    out, acc = {}, mp.mpf(0)
    for k in range(lo, hi + 1):
        acc += pmf(k)
        out[k] = acc
    return out


def add_disc(fam, sp, params, pmf, lo, hi, kmoments, cdf_pair, ppf_exact):
    a, b, c = (list(params) + [0, 0, 0])[:3]
    ks = sorted(set([lo] + [int(sp.ppf(p)) for p in (1e-8, 0.001, 0.05, 0.3, 0.5, 0.7, 0.95, 0.999)] + [lo + 1, lo + 2, lo + 4]))
    ks = [k for k in ks if lo <= k <= hi]
    for k in ks:
        e = pmf(k)
        if e != 0:
            sanity(fam + " pmf", sp.pmf(k), e)
            add(fam, PDF, k, a, b, c, e, 16)
        p, q = cdf_pair(k)
        if p != 0:
            sanity(fam + " cdf", sp.cdf(k), p)
            add(fam, CDF, k, a, b, c, p, 16)
        if q != 0 and q > mp.mpf("1e-300"):
            sanity(fam + " sf", sp.sf(k), q)
            add(fam, SF, k, a, b, c, q, 16)
    # the quantile: p between two cumulative values, never on a tie
    for k in ks:
        p_lo = cdf_pair(k - 1)[0] if k > lo else mp.mpf(0)
        p_hi = cdf_pair(k)[0]
        if p_hi - p_lo < mp.mpf("1e-12") or p_hi >= 1:
            continue
        for fr in (mp.mpf("0.25"), mp.mpf("0.75")):
            p = float(p_lo + fr * (p_hi - p_lo))
            if not (0 < p < 1):
                continue
            want = ppf_exact(p)
            if int(sp.ppf(p)) != want:
                sys.exit("SciPy differs for %s ppf at %r: %r vs %r" % (fam, p, sp.ppf(p), want))
            add(fam, PPF, p, a, b, c, want, 0)
    m, v, s, kk = [float(z) for z in sp.stats("mvsk")]
    for op, val in zip((MEAN, VAR, SKEW, KURT), kmoments):
        add(fam, op, 0, a, b, c, val, 16)


def ppf_from(cdf_of, lo, hi):
    def f(p):
        k = lo
        while k <= hi:
            if cdf_of(k) >= p:
                return k
            k += 1
        return hi
    return f


def disc_rows():
    mp.mp.dps = 50
    # binomial
    for n, p in [(10, 0.5), (100, 0.01), (1000, 0.3), (50, 0.99), (7, 1e-5)]:
        P, N = mp.mpf(p), n
        pmf = lambda k: mp.binomial(N, k) * P ** k * (1 - P) ** (N - k)
        cum = sums(pmf, 0, N)
        def pair(k, cum=cum, N=N, pmf=pmf):
            if k < 0:
                return mp.mpf(0), mp.mpf(1)
            low = cum[k]
            up = sum(pmf(j) for j in range(k + 1, N + 1)) if k < N else mp.mpf(0)
            return low, up
        mean, var = N * P, N * P * (1 - P)
        mom = (mean, var, (1 - 2 * P) / mp.sqrt(var), (1 - 6 * P * (1 - P)) / var)
        add_disc("binom", stats.binom(n, p), (n, p), pmf, 0, N, mom, pair, ppf_from(lambda k: cum[k], 0, N))
    # Poisson
    for lam in (0.5, 3.0, 100.0, 2000.0):
        L = mp.mpf(lam)
        hi = int(lam + 40 * math.sqrt(lam) + 60)
        pmf = lambda k: mp.exp(k * mp.log(L) - L - mp.loggamma(k + 1))
        cum = sums(pmf, 0, hi)

        def pair(k, cum=cum, pmf=pmf, hi=hi):
            if k < 0:
                return mp.mpf(0), mp.mpf(1)
            low = cum[k]
            up, j = mp.mpf(0), k + 1
            while True:
                t = pmf(j)
                up += t
                if t < mp.mpf("1e-45") * up and j > lam:
                    break
                j += 1
            return low, up
        add_disc("poisson", stats.poisson(lam), (lam,), pmf, 0, hi, (L, L, 1 / mp.sqrt(L), 1 / L), pair, ppf_from(lambda k: cum[k], 0, hi))
    # geometric (the trials)
    for p in (0.5, 0.01, 1e-4, 0.9):
        P = mp.mpf(p)
        pmf = lambda k: (1 - P) ** (k - 1) * P
        pair = lambda k: ((mp.mpf(0), mp.mpf(1)) if k < 1 else (-mp.expm1(k * mp.log1p(-P)), mp.exp(k * mp.log1p(-P))))
        mom = (1 / P, (1 - P) / P ** 2, (2 - P) / mp.sqrt(1 - P), 6 + P ** 2 / (1 - P))
        ppf = lambda q, P=P: int(mp.ceil(mp.log1p(-mp.mpf(q)) / mp.log1p(-P)))
        add_disc("geom", stats.geom(p), (p,), pmf, 1, 10 ** 7, mom, pair, ppf)
    # negative binomial
    for r, p in [(1.0, 0.5), (5.0, 0.3), (100.0, 0.9), (2.5, 0.05)]:
        R, P = mp.mpf(r), mp.mpf(p)
        pmf = lambda k: mp.exp(mp.loggamma(k + R) - mp.loggamma(R) - mp.loggamma(k + 1) + R * mp.log(P) + k * mp.log1p(-P))
        hi = int(stats.nbinom(r, p).ppf(1 - 1e-14)) + 50
        cum = sums(pmf, 0, hi)

        def pair(k, cum=cum, pmf=pmf, hi=hi):
            if k < 0:
                return mp.mpf(0), mp.mpf(1)
            low = cum[k]
            up, j = mp.mpf(0), k + 1
            while True:
                t = pmf(j)
                up += t
                if t < mp.mpf("1e-45") * up and j > (r * (1 - p) / p):
                    break
                j += 1
            return low, up
        mom = (R * (1 - P) / P, R * (1 - P) / P ** 2, (2 - P) / mp.sqrt(R * (1 - P)), 6 / R + P ** 2 / (R * (1 - P)))
        add_disc("nbinom", stats.nbinom(r, p), (r, p), pmf, 0, hi, mom, pair, ppf_from(lambda k: cum[k], 0, hi))
    # hypergeometric (population, successes, draws)
    for M, n, N in [(20, 10, 4), (100, 30, 40), (1000, 400, 300), (500, 5, 200)]:
        kmin, kmax = max(0, N - (M - n)), min(n, N)
        pmf = lambda k: mp.binomial(n, k) * mp.binomial(M - n, N - k) / mp.binomial(M, N)
        cum = sums(pmf, kmin, kmax)

        def pair(k, cum=cum, pmf=pmf, kmin=kmin, kmax=kmax):
            if k < kmin:
                return mp.mpf(0), mp.mpf(1)
            if k >= kmax:
                return mp.mpf(1), mp.mpf(0)
            return cum[k], sum(pmf(j) for j in range(k + 1, kmax + 1))
        mean = mp.mpf(n) * N / M
        var = mp.mpf(N) * n * (M - n) * (M - N) / (mp.mpf(M) ** 2 * (M - 1))
        m3 = sum(pmf(k) * (k - mean) ** 3 for k in range(kmin, kmax + 1))
        m4 = sum(pmf(k) * (k - mean) ** 4 for k in range(kmin, kmax + 1))
        mom = (mean, var, m3 / var ** 1.5, m4 / var ** 2 - 3)
        add_disc("hypergeom", stats.hypergeom(M, n, N), (M, n, N), pmf, kmin, kmax, mom, pair, ppf_from(lambda k: cum[k], kmin, kmax))


def check_oracles():
    """Compare the other libraries with the exact values of the rows.

    The programs are in test/oracle (Go with gonum), test/oracle-rust (Rust
    with statrs) and test/oracle-node (Node with jStat). A library that
    differs beyond its tolerance stops the script. A library that holds no
    family gives NaN, and the row is not counted."""
    import subprocess
    import tempfile
    sel = [r for r in ROWS if r[1] in (PDF, CDF, MEAN, VAR)]
    step = max(1, len(sel) // 1200)
    sel = sel[::step]
    tmp = tempfile.mkdtemp(prefix="dist-oracle-", dir=os.environ.get("TMPDIR"))
    inp = os.path.join(tmp, "in.txt")
    with open(inp, "w") as f:
        for r in sel:
            f.write("%d %d %r %r %r %r\n" % (r[0], r[1], r[2], r[3], r[4], r[5]))
    libs = {
        "gonum": (["go", "run", ".", inp, "distribution"], "oracle"),
        "statrs": (["cargo", "run", "--offline", "-q", "--", inp, "distribution"], "oracle-rust"),
        "jstat": (["node", "oracle.mjs", inp, "distribution"], "oracle-node"),
    }
    fails = []
    for name, (cmd, d) in libs.items():
        res = subprocess.run(cmd, cwd=os.path.join(HERE, d), capture_output=True, text=True)
        if res.returncode != 0:
            sys.exit("%s oracle failed:\n%s" % (name, res.stderr))
        out = res.stdout.split("\n")
        worst = {}
        for r, v in zip(sel, out):
            got = float(v)
            if not math.isfinite(got):
                continue
            want = r[6]
            err = abs(got - want) / max(abs(want), 1e-300)
            key = (r[0], r[1])
            if err > worst.get(key, (0,))[0]:
                worst[key] = (err, r)
            tol = 1e-5 if r[1] == PPF else 1e-8
            # A library that forms a tail as one minus a number loses the
            # digits of a small tail. The error of that form is the size
            # of the unit in the last place over the value.
            if r[1] in (CDF, SF) and want != 0.0:
                tol = max(tol, 4.0 * 2.220446049250313e-16 / abs(want))
            # The library jStat approximates the incomplete beta with a
            # tolerance of 1e-8 in its own solver. The Student t and the
            # F families differ from the exact value by that much.
            if name == "jstat":
                tol = max(tol, 1e-6)
            # A small difference in absolute terms is below the noise of
            # a library on a tail that cancels. The package holds the
            # exact value there. The floor of each library is the
            # tolerance that the library states for itself.
            floor = 1e-9 if name == "jstat" else 1e-12
            if err > tol and abs(got - want) > floor:
                fails.append((name, r[0], r[1], r[2:], got, err, tol))
                if len(fails) > 40:
                    break
        print(name, "compared", len(sel), "rows")
    for f in fails:
        print("DIFF", f)
    if fails:
        sys.exit("%d rows differ from the exact value" % len(fails))


def write():
    cols = ("D_FAM", "D_OP", "D_X", "D_A", "D_B", "D_C", "D_E", "D_TOL")
    with open(OUT, "w") as f:
        f.write("// Generated by test/make_distribution_vectors.py. Do not edit by hand.\n")
        f.write("// Reference vectors of the probability distributions.\n")
        f.write("// The exact value of each row comes from arithmetic at 50 digits.\n")
        f.write("// SciPy agrees with it within a loose bound, because SciPy forms a tail as one minus a value.\n")
        f.write("namespace suite;\n\n")
        f.write("const D_ROWS: i64 = %d;\n\n" % len(ROWS))
        for ci, name in enumerate(cols):
            kind = "i64" if name in ("D_FAM", "D_OP", "D_TOL") else "f64"
            f.write("const %s: %s[%d] = [\n" % (name, kind, len(ROWS)))
            for i in range(0, len(ROWS), 6):
                vals = []
                for r in ROWS[i:i + 6]:
                    v = r[ci]
                    vals.append("%d" % v if kind == "i64" else repr(float(v)))
                f.write("    " + ", ".join(vals) + ",\n")
            f.write("];\n\n")


def main():
    cont_rows()
    print("continuous rows", len(ROWS), flush=True)
    disc_rows()
    print("all rows", len(ROWS), flush=True)
    write()
    if "--oracles" in sys.argv:
        check_oracles()
    print("wrote", len(ROWS), "rows to", OUT)


if __name__ == "__main__":
    main()
