#!/usr/bin/env python3
"""Writes the vector file of the rank association group:
test/suite/assoc_rank_vectors.goose.

Needs: Python 3, NumPy and SciPy 1.18.1. The file holds numbers only.
Run it from the package directory with the system python3.
"""
import warnings
import numpy as np
from scipy import stats
from vecgen import G

warnings.simplefilter("ignore")
SCRIPT = "make_rank_vectors.py"

X = [2.1, 3.4, 1.9, 5.6, 4.4, 3.3, 6.1, 2.8, 4.9, 3.7]
Y = [1.2, 3.9, 2.2, 5.1, 4.8, 2.9, 6.6, 3.1, 4.1, 3.5]
XT = [1.0, 2.0, 2.0, 3.0, 4.0, 4.0, 4.0, 5.0, 6.0, 7.0]
YT = [2.0, 1.0, 3.0, 3.0, 5.0, 4.0, 6.0, 5.0, 7.0, 7.0]
B = [0, 0, 1, 1, 0, 1, 1, 0, 1, 1]
RANK = [3, 1, 2, 5, 4, 6, 9, 7, 8, 10]


def nm_x(x):
    return "xl[..]" if x is XL_REF[0] else ("xt[..]" if x is XT else "x[..]")


def nm_y(y):
    return "yl[..]" if y is YL_REF[0] else ("yt[..]" if y is YT else "y[..]")


XL_REF = [None]
YL_REF = [None]


def main():
    g = G("assoc_rank", SCRIPT)
    g.decl("x", X)
    g.decl("y", Y)
    g.decl("xt", XT)
    g.decl("yt", YT)
    g.decl("b", [float(v) for v in B])
    g.decl("rank", [float(v) for v in RANK])
    # spearmanrho
    for alt, k in (("two-sided", 0), ("less", 1), ("greater", 2)):
        r = stats.spearmanrho(X, Y, alternative=alt)
        g.field(f"spearmanrho {alt}", f"statistics::spearmanrho(x[..], y[..], {k})",
                {"statistic": r.statistic, "pvalue": r.pvalue})
    r = stats.spearmanrho(XT, YT)
    g.field("spearmanrho ties", "statistics::spearmanrho(xt[..], yt[..])",
            {"statistic": r.statistic, "pvalue": r.pvalue})
    r = stats.spearmanrho(X, Y)
    g.field("spearmanrho default", "statistics::spearmanrho(x[..], y[..])",
            {"statistic": r.statistic, "pvalue": r.pvalue})
    # pointbiserialr
    r = stats.pointbiserialr(B, Y)
    g.field("pointbiserialr", "statistics::pointbiserialr(b[..], y[..])",
            {"statistic": r.statistic, "pvalue": r.pvalue})
    r = stats.pointbiserialr(B, Y)
    g.field("pointbiserialr default", "statistics::pointbiserialr(b[..], y[..])",
            {"statistic": r.statistic, "pvalue": r.pvalue})
    # weightedtau
    for nm, x, y, call in (
        ("default", X, Y, "statistics::weightedtau(x[..], y[..])"),
        ("not additive", X, Y, "statistics::weightedtau(x[..], y[..], false)"),
        ("ties", XT, YT, "statistics::weightedtau(xt[..], yt[..])"),
        ("ties not additive", XT, YT, "statistics::weightedtau(xt[..], yt[..], false)"),
    ):
        add = nm != "not additive" and nm != "ties not additive"
        r = stats.weightedtau(x, y, additive=add)
        g.val(f"weightedtau {nm}", call, r.statistic)
    # the rank array of the caller
    r = stats.weightedtau(X, Y, rank=RANK)
    g.val("weightedtau rank array", "statistics::weightedtau_ranked(x[..], y[..], rank[..])", r.statistic)
    r = stats.weightedtau(XT, YT, rank=RANK)
    g.val("weightedtau rank array ties", "statistics::weightedtau_ranked(xt[..], yt[..], rank[..])", r.statistic)
    # chatterjeexi
    rng = np.random.default_rng(5)
    XL = rng.normal(size=40).tolist()
    YL = (np.array(XL) ** 2 + rng.normal(size=40) * 0.5).tolist()
    XL_REF[0] = XL
    YL_REF[0] = YL
    g.decl("xl", XL)
    g.decl("yl", YL)
    for nm, x, y, cont in (("default", X, Y, False), ("continuous", X, Y, True),
                           ("ties", XT, YT, False), ("ties continuous", XT, YT, True),
                           ("long", XL, YL, False), ("long continuous", XL, YL, True)):
        r = stats.chatterjeexi(x, y, y_continuous=cont)
        g.field(f"chatterjeexi {nm}",
                f"statistics::chatterjeexi({nm_x(x)}, {nm_y(y)}, {str(cont).lower()})",
                {"statistic": r.statistic, "pvalue": r.pvalue})
    # the fault edges
    g.raw("let none: f64[0] = [];")
    g.raw("let one: f64[1] = [2.0];")
    g.raw("let flat: f64[4] = [1.0, 1.0, 1.0, 1.0];")
    g.raw("let nan4: f64[4] = [1.0, 0.0 / 0.0, 3.0, 4.0];")
    g.raw("let inf4: f64[4] = [1.0, 1.0 / 0.0, 3.0, 4.0];")
    g.raw("let three: f64[3] = [1.0, 2.0, 3.0];")
    g.raw("let y4: f64[4] = [4.0, 3.0, 2.0, 1.0];")
    g.code("spearmanrho empty", "statistics::spearmanrho(none[..], none[..])", 1)
    g.code("spearmanrho one", "statistics::spearmanrho(one[..], one[..])", 2)
    g.code("spearmanrho length", "statistics::spearmanrho(x[..], three[..])", 6)
    g.code("spearmanrho NaN", "statistics::spearmanrho(nan4[..], y4[..])", 3)
    g.code("spearmanrho infinity", "statistics::spearmanrho(inf4[..], y4[..])", 4)
    g.code("spearmanrho unknown alternative", "statistics::spearmanrho(x[..], y[..], 9)", 5)
    g.code("pointbiserialr empty", "statistics::pointbiserialr(none[..], none[..])", 1)
    g.code("pointbiserialr one", "statistics::pointbiserialr(one[..], one[..])", 2)
    g.code("pointbiserialr constant", "statistics::pointbiserialr(flat[..], y4[..])", 7)
    g.code("pointbiserialr length", "statistics::pointbiserialr(x[..], three[..])", 6)
    g.code("pointbiserialr NaN", "statistics::pointbiserialr(nan4[..], y4[..])", 3)
    g.code("pointbiserialr unknown alternative", "statistics::pointbiserialr(b[..], y[..], 9)", 5)
    g.code("weightedtau empty", "statistics::weightedtau(none[..], none[..])", 1)
    g.code("weightedtau one", "statistics::weightedtau(one[..], one[..])", 2)
    g.code("weightedtau length", "statistics::weightedtau(x[..], three[..])", 6)
    g.code("weightedtau NaN", "statistics::weightedtau(nan4[..], y4[..])", 3)
    g.code("weightedtau rank length", "statistics::weightedtau_ranked(x[..], y[..], three[..])", 6)
    g.code("chatterjeexi empty", "statistics::chatterjeexi(none[..], none[..])", 1)
    g.code("chatterjeexi one", "statistics::chatterjeexi(one[..], one[..])", 2)
    g.code("chatterjeexi length", "statistics::chatterjeexi(x[..], three[..])", 6)
    g.code("chatterjeexi NaN", "statistics::chatterjeexi(nan4[..], y4[..])", 3)
    g.code("chatterjeexi infinity", "statistics::chatterjeexi(inf4[..], y4[..])", 4)
    g.write()


if __name__ == "__main__":
    main()
