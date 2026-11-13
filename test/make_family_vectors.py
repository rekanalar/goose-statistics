#!/usr/bin/env python3
"""Writes the vector files of the families and the laws:
test/suite/distributions_location_vectors.goose and later groups.

Needs: Python 3, NumPy and SciPy 1.18.1. The files hold numbers only.
Run it from the package directory with the system python3.
"""
import warnings
import numpy as np
from scipy import stats
from vecgen import G

warnings.simplefilter("ignore")
SCRIPT = "make_family_vectors.py"

XS = [-3.0, -1.0, -0.2, 0.0, 0.5, 1.0, 2.5, 8.0]
PS = [0.001, 0.1, 0.25, 0.5, 0.75, 0.9, 0.999]
SETS = [(0.0, 1.0, ""), (0.7, 2.3, "0.7, 2.3")]


def location_group():
    g = G("distributions_location", SCRIPT, tol="1e-11")
    for nm in ("cauchy", "laplace", "rayleigh", "gumbel_r", "gumbel_l"):
        d = getattr(stats, nm)
        for loc, scale, arg in SETS:
            tail = ", " + arg if arg else ""
            for x in XS:
                g.val(f"{nm} pdf {x} {loc}", f"statistics::{nm}_pdf({x!r}{tail})", d.pdf(x, loc, scale))
                g.val(f"{nm} cdf {x} {loc}", f"statistics::{nm}_cdf({x!r}{tail})", d.cdf(x, loc, scale))
                g.val(f"{nm} sf {x} {loc}", f"statistics::{nm}_sf({x!r}{tail})", d.sf(x, loc, scale))
            for p in PS:
                g.val(f"{nm} ppf {p} {loc}", f"statistics::{nm}_ppf({p!r}{tail})", d.ppf(p, loc, scale))
            if nm == "cauchy":
                for fn in ("mean", "var", "skew", "kurt"):
                    g.code(f"cauchy {fn} {loc}", f"statistics::cauchy_{fn}({arg})", 7)
            else:
                m, v, s, k = (float(q) for q in d.stats(loc, scale, moments="mvsk"))
                for fn, want in (("mean", m), ("var", v), ("skew", s), ("kurt", k)):
                    g.val(f"{nm} {fn} {loc}", f"statistics::{nm}_{fn}({arg})", want)
        # the far tails
        for x in (40.0, -40.0, 300.0, -300.0, 1e10, -1e10):
            g.val(f"{nm} sf {x}", f"statistics::{nm}_sf({x!r})", d.sf(x))
            g.val(f"{nm} cdf {x}", f"statistics::{nm}_cdf({x!r})", d.cdf(x))
            g.val(f"{nm} pdf {x}", f"statistics::{nm}_pdf({x!r})", d.pdf(x))
        for p in (1e-300, 1e-17, 1.0 - 1e-15):
            g.val(f"{nm} ppf {p}", f"statistics::{nm}_ppf({p!r})", d.ppf(p))
        # the score overflows when the location and the point have opposite signs
        for px, pl in ((1e308, 1e308), (-1e308, 1e308), (1e308, -1e308)):
            w = d.pdf(px, pl)
            label = f"{nm} pdf {px:g} {pl:g}"
            call = f"statistics::{nm}_pdf({px!r}, {pl!r})"
            if w != w:
                # SciPy has no density there. The package gives the undefined fault.
                g.code(label + " undefined", call, 7)
            else:
                g.val(label, call, w)
        g.code(f"{nm} scale 0", f"statistics::{nm}_pdf(1.0, 0.0, 0.0)", 5)
        g.code(f"{nm} scale negative", f"statistics::{nm}_cdf(1.0, 0.0, -1.0)", 5)
        g.code(f"{nm} ppf above 1", f"statistics::{nm}_ppf(1.5)", 5)
        g.code(f"{nm} ppf nan", f"statistics::{nm}_ppf(0.0 / 0.0)", 3)
        g.code(f"{nm} x infinity", f"statistics::{nm}_sf(1.0 / 0.0)", 4)
    # a result that does not fit in a double is a fault; SciPy gives infinity
    g.code("rayleigh ppf overflow", "statistics::rayleigh_ppf(0.5, 0.0, 1.7e308)", 8)
    g.code("gumbel_r ppf overflow", "statistics::gumbel_r_ppf(0.5, 0.0, 1e308)" if False else "statistics::gumbel_r_ppf(1e-300, 0.0, 1e308)", 8)
    g.code("gumbel_l ppf overflow", "statistics::gumbel_l_ppf(1.0 - 1e-16, 0.0, 1e308)", 8)
    g.code("rayleigh mean overflow", "statistics::rayleigh_mean(1e308, 1e308)", 8)
    g.val("rayleigh ppf 0", "statistics::rayleigh_ppf(0.0)", 0.0)
    g.code("rayleigh ppf 1", "statistics::rayleigh_ppf(1.0)", 5)
    g.code("gumbel_r ppf 0", "statistics::gumbel_r_ppf(0.0)", 5)
    g.write()


if __name__ == "__main__":
    location_group()
