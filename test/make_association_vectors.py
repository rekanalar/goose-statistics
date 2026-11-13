#!/usr/bin/env python3
"""Writes the vector files of the association and contingency groups:
test/suite/{contingency,divergence,clip}_vectors.goose and later groups.

Needs: Python 3, NumPy and SciPy 1.18.1. The files hold numbers only.
Run it from the package directory with the system python3.
"""
import warnings
import numpy as np
from scipy import stats
from scipy.stats import contingency as ct
from vecgen import G, f

warnings.simplefilter("ignore")
SCRIPT = "make_association_vectors.py"

OB = [[15, 20, 10], [30, 25, 12]]
T22 = [[10, 20], [30, 25]]
T43 = [[100, 150, 40], [203, 322, 60], [420, 700, 90], [320, 210, 70]]
T3D = (np.arange(24).reshape(2, 3, 4) + 1)


def flat(a):
    return [float(v) for v in np.asarray(a).ravel()]


def contingency_group():
    g = G("contingency", SCRIPT)
    g.decl("ob", flat(OB))
    g.decl("t22", flat(T22))
    g.decl("t43", flat(T43))
    g.decl("t3d", flat(T3D))
    g.decl_i("sh22", [2, 2])
    g.decl_i("sh23", [2, 3])
    g.decl_i("sh3d", [2, 3, 4])
    g.decl_i("sh43", [4, 3])
    g.fill("expected_freq 2x3", "statistics::expected_freq(out, ob[..], sh23[..])", flat(ct.expected_freq(OB)))
    g.fill("expected_freq 3d", "statistics::expected_freq(out, t3d[..], sh3d[..])", flat(ct.expected_freq(T3D)))
    for ax, m in enumerate(ct.margins(np.array(OB))):
        g.fill(f"margins 2x3 axis {ax}", f"statistics::margins(out, ob[..], sh23[..], {ax})", flat(m))
    for ax, m in enumerate(ct.margins(T3D)):
        g.fill(f"margins 3d axis {ax}", f"statistics::margins(out, t3d[..], sh3d[..], {ax})", flat(m))
    names = ("cramer", "tschuprow", "pearson")
    for tab, nm, rows, cols in ((OB, "ob", 2, 3), (T22, "t22", 2, 2), (T43, "t43", 4, 3)):
        a = np.array(tab)
        for mi, m in enumerate(names):
            for corr in (False, True):
                g.val(f"association {nm} {m} {corr}",
                      f"statistics::association({nm}[..], {rows}, {cols}, {mi}, {str(corr).lower()})",
                      ct.association(a, method=m, correction=corr))
    a = np.array(OB)
    for lam in (0, -0.5, 2.0 / 3.0):
        g.val(f"association lambda {lam}", f"statistics::association(ob[..], 2, 3, 0, false, {lam!r})",
              ct.association(a, method="cramer", lambda_=lam))
    g.raw("let neg: f64[4] = [-1.0, 2.0, 3.0, 4.0];")
    g.code("association negative", "statistics::association(neg[..], 2, 2)", 5)
    g.raw("let zrow: f64[4] = [0.0, 0.0, 3.0, 4.0];")
    g.code("association zero row", "statistics::association(zrow[..], 2, 2)", 5)
    g.raw("let frac: f64[4] = [1.5, 2.0, 3.0, 4.0];")
    g.code("association not whole", "statistics::association(frac[..], 2, 2)", 5)
    g.raw("let one: f64[3] = [1.0, 2.0, 3.0];")
    g.code("association one row", "statistics::association(one[..], 1, 3)", 2)
    # a shape with more axes than the bound
    g.raw("let many: i64[17] = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];")
    g.raw("let one_cell: f64[1] = [4.0];")
    g.codeonly("expected_freq too many axes", "statistics::expected_freq(out, one_cell[..], many[..])", 9)
    # a count above the exact limit
    g.raw("let huge: f64[4] = [1e17, 2.0, 3.0, 4.0];")
    g.codeonly("expected_freq count above the limit", "statistics::expected_freq(out, huge[..], sh22[..])", 5)
    # a table with two rows and one column
    g.raw("let twoone: f64[2] = [3.0, 5.0];")
    g.raw("let sh21: i64[2] = [2, 1];")
    g.code("association two rows one column", "statistics::association(twoone[..], 2, 1)", 2)
    g.codeonly("expected_freq bad shape", "statistics::expected_freq(out, ob[..], sh3d[..])", 6)
    g.codeonly("margins bad axis", "statistics::margins(out, ob[..], sh23[..], 2)", 5)
    g.codeonly("margins negative axis", "statistics::margins(out, ob[..], sh23[..], -1)", 5)
    g.raw("let zeros: f64[4] = [0.0, 0.0, 0.0, 0.0];")
    g.codeonly("expected_freq all zero", "statistics::expected_freq(out, zeros[..], sh22[..])", 7)
    g.write()


OBS = [16, 18, 16, 14, 12, 12]
EXP = [16, 16, 16, 16, 16, 8]
LAMS = [1.0, 0.0, -0.5, -1.0, -2.0, 2.0 / 3.0, 2.0, 3.0, -3.0, 0.25]


def divergence_group():
    g = G("divergence", SCRIPT)
    g.decl("obs", [float(v) for v in OBS])
    g.decl("exp", [float(v) for v in EXP])
    g.decl("small", [10.0, 10.0])
    g.decl("near", [12.0, 12.0])
    g.decl("zob", [0.0, 5.0, 7.0, 8.0])
    g.decl("zex", [4.0, 5.0, 5.0, 6.0])
    for lam in LAMS:
        r = stats.power_divergence(OBS, EXP, lambda_=lam)
        g.field(f"power_divergence lambda {lam}", f"statistics::power_divergence(obs[..], exp[..], {lam!r})",
                {"statistic": r.statistic, "pvalue": r.pvalue})
        r = stats.power_divergence(OBS, lambda_=lam)
        g.field(f"power_divergence uniform lambda {lam}", f"statistics::power_divergence_uniform(obs[..], {lam!r})",
                {"statistic": r.statistic, "pvalue": r.pvalue})
    for ddof in (1, 2, 3):
        r = stats.power_divergence(OBS, EXP, ddof=ddof)
        g.field(f"power_divergence ddof {ddof}", f"statistics::power_divergence(obs[..], exp[..], 1.0, {ddof})",
                {"statistic": r.statistic, "pvalue": r.pvalue, "df": len(OBS) - 1 - ddof})
    # a zero observed count
    for lam in (1.0, 0.0, -1.0, 0.5, -0.5):
        try:
            r = stats.power_divergence([0, 5, 7, 8], [4, 5, 5, 6], lambda_=lam)
        except Exception:
            continue
        if np.isfinite(r.statistic):
            g.field(f"power_divergence zero count {lam}", f"statistics::power_divergence(zob[..], zex[..], {lam!r})",
                    {"statistic": r.statistic, "pvalue": r.pvalue})
        else:
            g.code(f"power_divergence zero count {lam} undefined", f"statistics::power_divergence(zob[..], zex[..], {lam!r})", 7)
    g.code("power_divergence sums differ", "statistics::power_divergence(small[..], near[..])", 5)
    g.raw("let none: f64[0] = [];")
    g.raw("let one: f64[1] = [5.0];")
    g.raw("let neg: f64[3] = [-1.0, 2.0, 3.0];")
    g.raw("let zeroexp: f64[2] = [0.0, 3.0];")
    g.raw("let onetwo: f64[2] = [1.0, 2.0];")
    g.code("power_divergence one value", "statistics::power_divergence(one[..], none[..])", 2)
    g.code("power_divergence empty", "statistics::power_divergence(none[..], none[..])", 1)
    g.code("power_divergence negative count", "statistics::power_divergence(neg[..], none[..])", 5)
    g.code("power_divergence zero expected", "statistics::power_divergence(onetwo[..], zeroexp[..])", 7)
    g.code("power_divergence length", "statistics::power_divergence(obs[..], small[..])", 6)
    g.code("power_divergence ddof too big", "statistics::power_divergence(obs[..], exp[..], 1.0, 5)", 5)
    g.code("power_divergence negative ddof", "statistics::power_divergence(obs[..], exp[..], 1.0, -1)", 5)
    g.code("power_divergence lambda NaN", "statistics::power_divergence(obs[..], exp[..], 0.0 / 0.0)", 3)
    g.code("power_divergence lambda infinite", "statistics::power_divergence(obs[..], exp[..], 1.0 / 0.0)", 4)
    g.write()


A = [1.0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 100, 1000]


def clip_group():
    g = G("clip", SCRIPT)
    g.decl("a", A)
    g.decl("flat3", [3.0, 3.0, 3.0])
    rng = np.random.default_rng(1)
    b = np.r_[rng.normal(0, 1, 50), [20, -30]].tolist()
    g.decl("b", b)
    for nm, x, lo, hi in (("default", A, 4.0, 4.0), ("1.5", A, 1.5, 1.5), ("0.5 3", A, 0.5, 3.0),
                          ("zero low", A, 0.0, 4.0), ("50", b, 3.0, 3.0)):
        c, l, u = stats.sigmaclip(x, lo, hi)
        var = "a" if x is A else "b"
        g.field(f"sigmaclip {nm}", f"statistics::sigmaclip(out, {var}[..], {lo!r}, {hi!r})",
                {"low": l, "high": u, "kept": len(c)})
        g.n += 1
        k = g.n
        g.lines.append(f"    let _r{k}, c{k} = statistics::sigmaclip(out, {var}[..], {lo!r}, {hi!r});")
        g.lines.append(f"    let w{k}: f64[{len(c)}] = {'[' + ', '.join(f(v) for v in c) + ']'};")
        g.lines.append(f'    chk_list(bad, "sigmaclip {nm} kept", out[..], c{k}, w{k}[..], {g.tol});')
    c, l, u = stats.sigmaclip([3.0, 3, 3])
    g.field("sigmaclip constant", "statistics::sigmaclip(out, flat3[..])", {"low": l, "high": u, "kept": 3})
    g.raw("let none: f64[0] = [];")
    g.raw("let nan3: f64[3] = [1.0, 0.0 / 0.0, 2.0];")
    g.raw("let inf3: f64[3] = [1.0, 1.0 / 0.0, 2.0];")
    g.code("sigmaclip empty", "statistics::sigmaclip(out, none[..])", 1)
    g.code("sigmaclip NaN", "statistics::sigmaclip(out, nan3[..])", 3)
    g.code("sigmaclip infinity", "statistics::sigmaclip(out, inf3[..])", 4)
    g.code("sigmaclip negative low", "statistics::sigmaclip(out, a[..], -1.0, 4.0)", 5)
    g.code("sigmaclip negative high", "statistics::sigmaclip(out, a[..], 4.0, -0.5)", 5)
    g.code("sigmaclip NaN factor", "statistics::sigmaclip(out, a[..], 0.0 / 0.0, 4.0)", 3)
    g.code("sigmaclip infinite factor", "statistics::sigmaclip(out, a[..], 4.0, 1.0 / 0.0)", 4)
    g.write()


if __name__ == "__main__":
    contingency_group()
    divergence_group()
    clip_group()
