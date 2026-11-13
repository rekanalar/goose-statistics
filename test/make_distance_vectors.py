#!/usr/bin/env python3
"""Writes the vector file of the distance group:
test/suite/distance_vectors.goose.

Needs: Python 3, NumPy and SciPy 1.18.1. The file holds numbers only.
Run it from the package directory with the system python3.
"""
import math
import warnings
import numpy as np
from scipy import stats
from vecgen import G, f, arr

warnings.simplefilter("ignore")
SCRIPT = "make_distance_vectors.py"

U = [0.0, 1.0, 3.0, 5.0, 8.0]
V = [0.5, 2.0, 2.5, 7.0]
UW = [1.0, 2.0, 1.0, 3.0, 1.0]
VW = [2.0, 1.0, 1.0, 1.0]
TIES_U = [1.0, 2.0, 2.0, 3.0]
TIES_V = [1.0, 1.0, 2.5]
G1 = [1.0, 2.0, 3.0, 4.0, 5.0]
G2 = [2.0, 4.0, 6.0, 8.0, 11.0, 13.0]
G3 = [5.0, 5.0, 6.0, 7.0]
ONES = [1.0] * 4

X = [2.1, 3.4, 1.9, 5.6, 4.4, 3.3, 6.1, 2.8, 4.9, 3.7]
Y = [1.2, 3.9, 2.2, 5.1, 4.8, 2.9, 6.6, 3.1, 4.1, 3.5]
XT = [1.0, 2.0, 2.0, 3.0, 4.0, 4.0, 4.0, 5.0, 6.0, 7.0]
YT = [2.0, 1.0, 3.0, 3.0, 5.0, 4.0, 6.0, 5.0, 7.0, 7.0]
RNG = np.random.default_rng(11)
XL = (RNG.normal(size=30) * 3.0 + 10.0).tolist()
YL = (0.8 * np.asarray(XL) + 1.5 + RNG.normal(size=30) * 0.4).tolist()


def main():
    g = G("distance", SCRIPT)
    g.decl("u", U)
    g.decl("v", V)
    g.decl("uw", UW)
    g.decl("vw", VW)
    g.decl("tu", TIES_U)
    g.decl("tv", TIES_V)
    g.raw("let vones: f64[4] = [1.0, 1.0, 1.0, 1.0];")
    # the groups of the O'Brien transform
    g.raw("let g1: f64[5] = [1.0, 2.0, 3.0, 4.0, 5.0];")
    g.raw("let g2: f64[6] = [2.0, 4.0, 6.0, 8.0, 11.0, 13.0];")
    g.raw("let g3: f64[4] = [5.0, 5.0, 6.0, 7.0];")
    g.raw("let ga: f64[15] = [1.0, 2.0, 3.0, 4.0, 5.0, 2.0, 4.0, 6.0, 8.0, 11.0, 13.0, 5.0, 5.0, 6.0, 7.0];")
    g.decl_i("sz1", [len(G1)])
    g.decl_i("sz3", [len(G1), len(G2), len(G3)])
    g.decl_i("szshort", [2])
    g.decl_i("szempty", [0])
    # the distances of two samples
    g.val("wasserstein", "statistics::wasserstein_distance(u[..], v[..])",
          stats.wasserstein_distance(U, V))
    g.val("wasserstein weights", "statistics::wasserstein_distance_weighted(u[..], v[..], uw[..], vw[..])",
          stats.wasserstein_distance(U, V, UW, VW))
    g.val("wasserstein one weight list", "statistics::wasserstein_distance_weighted(u[..], v[..], uw[..], vones[..])",
          stats.wasserstein_distance(U, V, UW))
    g.val("wasserstein ties", "statistics::wasserstein_distance(tu[..], tv[..])",
          stats.wasserstein_distance(TIES_U, TIES_V))
    g.val("wasserstein equal", "statistics::wasserstein_distance(u[..], u[..])",
          stats.wasserstein_distance(U, U))
    g.val("energy", "statistics::energy_distance(u[..], v[..])",
          stats.energy_distance(U, V))
    g.val("energy weights", "statistics::energy_distance_weighted(u[..], v[..], uw[..], vw[..])",
          stats.energy_distance(U, V, UW, VW))
    g.val("energy ties", "statistics::energy_distance(tu[..], tv[..])",
          stats.energy_distance(TIES_U, TIES_V))
    g.val("energy equal", "statistics::energy_distance(u[..], u[..])",
          stats.energy_distance(U, U))
    # the O'Brien transform of the groups
    for name, groups, sizes, call in (
            ("one", [G1], [len(G1)], "g1[..], sz1[..]"),
            ("three", [G1, G2, G3], [len(G1), len(G2), len(G3)], "ga[..], sz3[..]")):
        want = [v for a in stats.obrientransform(*groups) for v in np.asarray(a).tolist()]
        g.fill(f"obrientransform {name}",
               f"statistics::obrientransform(out, {call})", want)
    # the fault edges
    g.raw("let none: f64[0] = [];")
    g.raw("let nan4: f64[4] = [1.0, 0.0 / 0.0, 3.0, 4.0];")
    g.raw("let inf4: f64[4] = [1.0, 1.0 / 0.0, 3.0, 4.0];")
    g.raw("let negw: f64[5] = [-1.0, 2.0, 1.0, 3.0, 1.0];")
    g.raw("let zerow: f64[5] = [0.0, 0.0, 0.0, 0.0, 0.0];")
    g.raw("let shortw: f64[4] = [1.0, 2.0, 1.0, 3.0];")
    g.raw("let shortw3: f64[3] = [1.0, 2.0, 1.0];")
    g.raw("let negw4: f64[4] = [-1.0, 2.0, 1.0, 3.0];")
    g.code("wasserstein empty", "statistics::wasserstein_distance(none[..], v[..])", 1)
    g.code("wasserstein NaN", "statistics::wasserstein_distance(nan4[..], v[..])", 3)
    g.code("wasserstein infinity", "statistics::wasserstein_distance(inf4[..], v[..])", 4)
    g.code("wasserstein weight length", "statistics::wasserstein_distance_weighted(u[..], v[..], shortw[..], vw[..])", 6)
    g.code("wasserstein negative weight", "statistics::wasserstein_distance_weighted(u[..], v[..], negw[..], vw[..])", 5)
    g.code("wasserstein zero weight", "statistics::wasserstein_distance_weighted(u[..], v[..], zerow[..], vw[..])", 7)
    g.code("energy empty", "statistics::energy_distance(none[..], v[..])", 1)
    g.code("energy NaN", "statistics::energy_distance(nan4[..], v[..])", 3)
    g.code("energy infinity", "statistics::energy_distance(inf4[..], v[..])", 4)
    g.code("energy weight length", "statistics::energy_distance_weighted(u[..], v[..], uw[..], shortw3[..])", 6)
    g.code("energy negative weight", "statistics::energy_distance_weighted(u[..], v[..], uw[..], negw4[..])", 5)
    g.codeonly("obrientransform short group", "statistics::obrientransform(out, g1[..], szshort[..])", 2)
    g.codeonly("obrientransform empty", "statistics::obrientransform(out, none[..], szempty[..])", 1)
    g.write()


def line_group():

    g = G("line_siegel", SCRIPT)
    g.decl("x", X)
    g.decl("y", Y)
    g.decl("xt", XT)
    g.decl("yt", YT)
    g.decl("xl", XL)
    g.decl("yl", YL)
    for name, xs, ys in (("default", X, Y), ("ties", XT, YT), ("long", XL, YL)):
        for method, m in (("hierarchical", 0), ("separate", 1)):
            r = stats.siegelslopes(ys, xs, method=method)
            var = "x" if xs is X else ("xt" if xs is XT else "xl")
            vary = "y" if ys is Y else ("yt" if ys is YT else "yl")
            g.val2(f"siegelslopes {name} {method}",
                   f"statistics::siegelslopes({vary}[..], {var}[..], {m})", (r.slope, r.intercept))
    # the fault edges
    g.raw("let none: f64[0] = [];")
    g.raw("let one: f64[1] = [2.0];")
    g.raw("let two: f64[2] = [1.0, 2.0];")
    g.raw("let three: f64[3] = [1.0, 2.0, 3.0];")
    g.raw("let nan3: f64[3] = [1.0, 0.0 / 0.0, 3.0];")
    g.raw("let inf3: f64[3] = [1.0, 1.0 / 0.0, 3.0];")
    g.raw("let flat3: f64[3] = [1.0, 1.0, 1.0];")
    g.code3("siegelslopes empty", "statistics::siegelslopes(none[..], none[..])", 1)
    g.code3("siegelslopes one", "statistics::siegelslopes(one[..], one[..])", 2)
    g.code3("siegelslopes length", "statistics::siegelslopes(y[..], three[..])", 6)
    g.code3("siegelslopes NaN", "statistics::siegelslopes(nan3[..], three[..])", 3)
    g.code3("siegelslopes infinity", "statistics::siegelslopes(inf3[..], three[..])", 4)
    g.code3("siegelslopes constant x", "statistics::siegelslopes(three[..], flat3[..])", 7)
    g.code3("siegelslopes unknown method", "statistics::siegelslopes(y[..], x[..], 9)", 5)
    g.write()




def laws_group():
    """The rows of truncate, order_statistic and dirichlet_multinomial."""
    from scipy import stats as st
    g = G("laws", SCRIPT)
    val = g.val

    def near_zero_val(label, call, want, tol=None):
        """A reference value of rounding noise is checked as one plus the value."""
        if abs(want) >= 1e-12:
            return val(label, call, want, tol)
        g.n += 1
        k = g.n
        g.lines.append(f"    let v{k}, c{k} = {call};")
        g.lines.append(f'    chk_val(bad, "{label}", v{k} + 1.0, c{k}, 1.0, 1e-12);')
    g.val = near_zero_val
    inf = "inf()"
    fams = {"norm": (0, "0.0, 1.0", st.Normal()),
            "uniform": (1, "0.0, 1.0", st.Uniform(a=0.0, b=1.0)),
            "logistic": (2, "0.0, 1.0", st.Logistic())}
    def lit(v):
        return "(1.0 / 0.0)" if v == float("inf") else ("(-1.0 / 0.0)" if v == float("-inf") else repr(float(v)))
    for lb, ub in ((-1.0, 2.0), (0.5, 1.5), (float("-inf"), 0.0), (0.0, float("inf"))):
        T = st.truncate(st.Normal(), lb, ub)
        b = f"{lit(lb)}, {lit(ub)}"
        for x in (-0.5, 0.0, 0.7, 1.2):
            g.val(f"truncate norm {lb} {ub} pdf {x}", f"statistics::truncate_pdf({x}, 0, 0.0, 1.0, {b})", float(T.pdf(x)))
            g.val(f"truncate norm {lb} {ub} cdf {x}", f"statistics::truncate_cdf({x}, 0, 0.0, 1.0, {b})", float(T.cdf(x)))
        for p in (0.1, 0.5, 0.9):
            g.val(f"truncate norm {lb} {ub} ppf {p}", f"statistics::truncate_ppf({p}, 0, 0.0, 1.0, {b})", float(T.icdf(p)))
        g.val(f"truncate norm {lb} {ub} mean", f"statistics::truncate_mean(0, 0.0, 1.0, {b})", float(T.mean()), "1e-10")
        g.val(f"truncate norm {lb} {ub} var", f"statistics::truncate_var(0, 0.0, 1.0, {b})", float(T.variance()), "1e-10")
    T = st.truncate(st.Logistic(), 0.5, 2.0)
    g.val("truncate logistic pdf", "statistics::truncate_pdf(1.0, 2, 0.0, 1.0, 0.5, 2.0)", float(T.pdf(1.0)))
    g.val("truncate logistic cdf", "statistics::truncate_cdf(1.0, 2, 0.0, 1.0, 0.5, 2.0)", float(T.cdf(1.0)))
    g.val("truncate logistic ppf", "statistics::truncate_ppf(0.3, 2, 0.0, 1.0, 0.5, 2.0)", float(T.icdf(0.3)))
    g.val("truncate logistic mean", "statistics::truncate_mean(2, 0.0, 1.0, 0.5, 2.0)", float(T.mean()), "1e-10")
    g.val("truncate logistic var", "statistics::truncate_var(2, 0.0, 1.0, 0.5, 2.0)", float(T.variance()), "1e-10")
    # a Laplace law truncated over a range that holds its location: the
    # density has a kink there, so the integral splits at the location.
    # SciPy cannot truncate a Laplace law, so the reference comes from a
    # quadrature of the density.
    from scipy.integrate import quad
    def lap_d(x, loc, scale):
        return math.exp(-abs((x - loc) / scale)) / (2.0 * scale)
    def trunc_lap(loc, scale, lo, hi):
        mass = quad(lambda x: lap_d(x, loc, scale), lo, hi)[0]
        mean = quad(lambda x: x * lap_d(x, loc, scale), lo, hi)[0] / mass
        var = quad(lambda x: (x - mean) ** 2 * lap_d(x, loc, scale), lo, hi)[0] / mass
        return mass, mean, var
    mass, mref, vref = trunc_lap(1.0, 1.0, 0.0, 3.0)
    g.val("truncate laplace mean", "statistics::truncate_mean(4, 1.0, 1.0, 0.0, 3.0)", mref, "1e-12")
    g.val("truncate laplace var", "statistics::truncate_var(4, 1.0, 1.0, 0.0, 3.0)", vref, "1e-12")
    g.val("truncate laplace pdf", "statistics::truncate_pdf(0.4, 4, 1.0, 1.0, 0.0, 3.0)", lap_d(0.4, 1.0, 1.0) / mass)
    T = st.truncate(st.Uniform(a=0.0, b=1.0), 0.2, 0.6)
    g.val("truncate uniform pdf", "statistics::truncate_pdf(0.3, 1, 0.0, 1.0, 0.2, 0.6)", float(T.pdf(0.3)))
    g.val("truncate uniform cdf", "statistics::truncate_cdf(0.3, 1, 0.0, 1.0, 0.2, 0.6)", float(T.cdf(0.3)))
    g.val("truncate uniform ppf", "statistics::truncate_ppf(0.3, 1, 0.0, 1.0, 0.2, 0.6)", float(T.icdf(0.3)))
    g.val("truncate uniform mean", "statistics::truncate_mean(1, 0.0, 1.0, 0.2, 0.6)", float(T.mean()), "1e-10")
    g.val("truncate uniform var", "statistics::truncate_var(1, 0.0, 1.0, 0.2, 0.6)", float(T.variance()), "1e-10")
    g.code("truncate bounds out of order", "statistics::truncate_pdf(0.0, 0, 0.0, 1.0, 2.0, 1.0)", 5)
    g.code("truncate NaN bound", "statistics::truncate_pdf(0.0, 0, 0.0, 1.0, 0.0 / 0.0, 1.0)", 3)
    g.code("truncate NaN point", "statistics::truncate_pdf(0.0 / 0.0, 0, 0.0, 1.0, 0.0, 1.0)", 3)
    g.code("truncate infinite point", "statistics::truncate_cdf(1.0 / 0.0, 0, 0.0, 1.0, 0.0, 1.0)", 4)
    g.code("truncate unknown family", "statistics::truncate_pdf(0.0, 99, 0.0, 1.0, 0.0, 1.0)", 5)
    g.code("truncate bad scale", "statistics::truncate_pdf(0.0, 0, 0.0, -1.0, 0.0, 1.0)", 5)
    g.code("truncate no mass", "statistics::truncate_pdf(50.5, 0, 0.0, 1.0, 50.0, 51.0)", 7)
    g.code("truncate no mass quantile", "statistics::truncate_ppf(0.5, 0, 0.0, 1.0, 50.0, 51.0)", 7)
    g.code("truncate Cauchy infinite low bound", "statistics::truncate_mean(3, 0.0, 1.0, -1.0 / 0.0, 1.0)", 7)
    g.code("truncate Cauchy infinite mean", "statistics::truncate_mean(3, 0.0, 1.0, 0.0, 1.0 / 0.0)", 7)
    g.code("truncate quantile of an open end", "statistics::truncate_ppf(0.0, 0, 0.0, 1.0, -1.0 / 0.0, 0.0)", 5)
    g.code("truncate probability out of range", "statistics::truncate_ppf(1.5, 0, 0.0, 1.0, 0.0, 1.0)", 5)
    # the order statistic
    for r, n in ((1, 5), (3, 5), (5, 5), (2, 7)):
        for nm, (fam, par, base) in fams.items():
            O = st.order_statistic(base, r=r, n=n)
            call = f"{fam}, {par}, {r}, {n}"
            g.val(f"order {nm} {r} {n} pdf", f"statistics::order_statistic_pdf(0.4, {call})", float(O.pdf(0.4)))
            g.val(f"order {nm} {r} {n} pdf below the location", f"statistics::order_statistic_pdf(-0.4, {call})", float(O.pdf(-0.4)))
            g.val(f"order {nm} {r} {n} cdf", f"statistics::order_statistic_cdf(0.4, {call})", float(O.cdf(0.4)))
            g.val(f"order {nm} {r} {n} ppf", f"statistics::order_statistic_ppf(0.3, {call})", float(O.icdf(0.3)))
            g.val(f"order {nm} {r} {n} mean", f"statistics::order_statistic_mean({call})", float(O.mean()), "1e-9")
            g.val(f"order {nm} {r} {n} var", f"statistics::order_statistic_var({call})", float(O.variance()), "1e-9")
    g.code("order rank above n", "statistics::order_statistic_pdf(0.0, 0, 0.0, 1.0, 6, 5)", 5)
    g.code("order rank zero", "statistics::order_statistic_pdf(0.0, 0, 0.0, 1.0, 0, 5)", 5)
    g.code("order n too long", "statistics::order_statistic_cdf(0.0, 0, 0.0, 1.0, 1, 2000000)", 9)
    g.code("order NaN point", "statistics::order_statistic_pdf(0.0 / 0.0, 0, 0.0, 1.0, 1, 5)", 3)
    g.code("order Cauchy mean", "statistics::order_statistic_mean(3, 0.0, 1.0, 2, 5)", 7)
    g.code("order unknown family", "statistics::order_statistic_cdf(0.0, 99, 0.0, 1.0, 1, 5)", 5)
    # the Dirichlet multinomial
    al = [1.5, 2.0, 0.5]
    g.decl("al", al)
    for x in ([2, 1, 1], [0, 0, 4], [4, 0, 0]):
        g.raw(f"let cx{''.join(map(str, x))}: f64[3] = {arr_f(x)};")
        nm = "cx" + "".join(map(str, x))
        g.val(f"dirmult pmf {x}", f"statistics::dirichlet_multinomial_pmf({nm}[..], al[..], 4)", float(st.dirichlet_multinomial.pmf(x, al, 4)))
        g.val(f"dirmult logpmf {x}", f"statistics::dirichlet_multinomial_logpmf({nm}[..], al[..], 4)", float(st.dirichlet_multinomial.logpmf(x, al, 4)))
    g.fill("dirmult mean", "statistics::dirichlet_multinomial_mean(out, al[..], 4)", st.dirichlet_multinomial.mean(al, 4).tolist())
    g.fill("dirmult var", "statistics::dirichlet_multinomial_var(out, al[..], 4)", st.dirichlet_multinomial.var(al, 4).tolist())
    g.fill("dirmult cov", "statistics::dirichlet_multinomial_cov(out, al[..], 4)", st.dirichlet_multinomial.cov(al, 4).ravel().tolist())
    g.raw("let big: f64[3] = [500.0, 300.0, 200.0];")
    g.raw("let bigal: f64[3] = [0.5, 1.5, 3.0];")
    g.val("dirmult big", "statistics::dirichlet_multinomial_logpmf(big[..], bigal[..], 1000)",
          float(st.dirichlet_multinomial.logpmf([500, 300, 200], [0.5, 1.5, 3.0], 1000)), "1e-10")
    g.raw("let cbad: f64[3] = [2.0, 1.0, 0.0];")
    g.raw("let alzero: f64[3] = [1.5, 0.0, 0.5];")
    g.raw("let cz: f64[3] = [0.0, 0.0, 0.0];")
    g.raw("let cfrac: f64[3] = [1.5, 1.5, 1.0];")
    g.raw("let cshort: f64[2] = [2.0, 2.0];")
    g.raw("let cnan: f64[3] = [2.0, 0.0 / 0.0, 2.0];")
    g.raw("let alnan: f64[3] = [1.5, 0.0 / 0.0, 0.5];")
    g.raw("let alinf: f64[3] = [1.5, 1.0 / 0.0, 0.5];")
    g.raw("let none: f64[0] = [];")
    g.code("dirmult counts do not sum to n", "statistics::dirichlet_multinomial_pmf(cbad[..], al[..], 4)", 5)
    g.code("dirmult alpha zero", "statistics::dirichlet_multinomial_pmf(cx211[..], alzero[..], 4)", 5)
    g.val("dirmult no draws logpmf", "statistics::dirichlet_multinomial_logpmf(cz[..], al[..], 0)", 0.0)
    g.code("dirmult no draws", "statistics::dirichlet_multinomial_pmf(cz[..], al[..], 0)", 0)
    g.code("dirmult negative n", "statistics::dirichlet_multinomial_pmf(cz[..], al[..], -1)", 5)
    g.code("dirmult fractional count", "statistics::dirichlet_multinomial_pmf(cfrac[..], al[..], 4)", 5)
    g.code("dirmult length", "statistics::dirichlet_multinomial_pmf(cshort[..], al[..], 4)", 6)
    g.code("dirmult NaN count", "statistics::dirichlet_multinomial_pmf(cnan[..], al[..], 4)", 3)
    g.code("dirmult NaN alpha", "statistics::dirichlet_multinomial_pmf(cx211[..], alnan[..], 4)", 3)
    g.code("dirmult infinite alpha", "statistics::dirichlet_multinomial_pmf(cx211[..], alinf[..], 4)", 4)
    g.code("dirmult empty", "statistics::dirichlet_multinomial_pmf(none[..], none[..], 4)", 1)
    g.codeonly("dirmult mean fault", "statistics::dirichlet_multinomial_mean(out, alzero[..], 4)", 5)
    g.write()


PZ = [0.4, 0.9, 1.3, 1.8, 2.4, 2.9, 3.6, 4.1, 5.5, 7.2]


def probplot_group():
    """The rows of probplot, ppcc_plot and ppcc_max."""
    from scipy import stats as st
    g = G("probplot", SCRIPT)
    g.decl("z", PZ)
    g.decl("xt", XT)
    fams = (("norm", 0, 0.0, ()), ("expon", 1, 0.0, ()), ("uniform", 2, 0.0, ()),
            ("tukeylambda", 3, 0.3, (0.3,)), ("logistic", 4, 0.0, ()), ("laplace", 5, 0.0, ()),
            ("cauchy", 6, 0.0, ()), ("gumbel_r", 7, 0.0, ()), ("gumbel_l", 8, 0.0, ()),
            ("rayleigh", 9, 0.0, ()), ("weibull_min", 10, 1.5, (1.5,)))

    def plot(label, data, name, fam, shape, want):
        g.n += 1
        k = g.n
        (osm, osr), (slope, icpt, r) = want
        g.lines.append(f"    let p{k}, c{k} = statistics::probplot(out, out2, {name}[..], {fam}, {shape});")
        g.lines.append(f"    let a{k}: f64[{len(osm)}] = {arr(osm)};")
        g.lines.append(f"    let b{k}: f64[{len(osr)}] = {arr(osr)};")
        g.lines.append(f'    chk_list(bad, "{label} quantiles", out[..], c{k}, a{k}[..], 1e-12);')
        g.lines.append(f'    chk_list(bad, "{label} sorted", out2[..], c{k}, b{k}[..], 1e-12);')
        g.lines.append(f'    chk_val(bad, "{label} slope", p{k}.slope, c{k}, {f(slope)}, 1e-12);')
        g.lines.append(f'    chk_val(bad, "{label} intercept", p{k}.intercept, c{k}, {f(icpt)}, 1e-12);')
        g.lines.append(f'    chk_val(bad, "{label} correlation", p{k}.r, c{k}, {f(r)}, 1e-12);')
    for nm, fam, shape, sp in fams:
        plot(f"probplot {nm}", PZ, "z", fam, shape, st.probplot(PZ, sparams=sp, dist=nm, fit=True))
    plot("probplot ties", XT, "xt", 0, 0.0, st.probplot(XT))
    g.raw("let two: f64[2] = [1.0, 2.0];")
    plot("probplot two values", [1.0, 2.0], "two", 0, 0.0, st.probplot([1.0, 2.0]))
    # the plot of the correlations over a grid
    for label, a, b, n in (("a", -2, 2, 9), ("b", 0, 1, 5), ("one", -2, 2, 1)):
        svals, ppcc = st.ppcc_plot(PZ, a, b, N=n)
        g.fill(f"ppcc_plot {label}", f"statistics::ppcc_plot(out, z[..], {float(a)}, {float(b)}, {n})", ppcc.tolist(), "1e-12")
    # the search of the best shape; the tolerance of the search is 1e-8
    g.val("ppcc_max", "statistics::ppcc_max(z[..])", float(st.ppcc_max(PZ)), "1e-7")
    g.val("ppcc_max bracket", "statistics::ppcc_max(z[..], 3, -1.0, 1.0)", float(st.ppcc_max(PZ, brack=(-1.0, 1.0))), "1e-7")
    g.val("ppcc_max weibull", "statistics::ppcc_max(z[..], 10, 0.5, 3.0)", float(st.ppcc_max(PZ, dist="weibull_min", brack=(0.5, 3.0))), "1e-7")
    nz = np.random.default_rng(2).normal(size=40).tolist()
    g.decl("nz", nz)
    g.val("ppcc_max normal sample", "statistics::ppcc_max(nz[..])", float(st.ppcc_max(nz)), "1e-7")
    # the fault edges
    g.raw("let none: f64[0] = [];")
    g.raw("let one: f64[1] = [1.0];")
    g.raw("let flat: f64[3] = [2.0, 2.0, 2.0];")
    g.raw("let nan3: f64[3] = [1.0, 0.0 / 0.0, 3.0];")
    g.raw("let inf3: f64[3] = [1.0, 1.0 / 0.0, 3.0];")
    g.code("probplot empty", "statistics::probplot(out, out2, none[..])", 1)
    g.code("probplot one value", "statistics::probplot(out, out2, one[..])", 2)
    g.code("probplot NaN", "statistics::probplot(out, out2, nan3[..])", 3)
    g.code("probplot infinity", "statistics::probplot(out, out2, inf3[..])", 4)
    g.code("probplot constant sample", "statistics::probplot(out, out2, flat[..])", 7)
    g.code("probplot unknown family", "statistics::probplot(out, out2, z[..], 99)", 5)
    g.code("probplot bad weibull shape", "statistics::probplot(out, out2, z[..], 10, -1.0)", 5)
    g.code("probplot NaN shape", "statistics::probplot(out, out2, z[..], 3, 0.0 / 0.0)", 3)
    g.codeonly("ppcc_plot bounds out of order", "statistics::ppcc_plot(out, z[..], 2.0, -2.0, 9)", 5)
    g.codeonly("ppcc_plot grid too large", "statistics::ppcc_plot(out, z[..], 0.0, 1.0, 5000)", 9)
    g.codeonly("ppcc_plot empty grid", "statistics::ppcc_plot(out, z[..], 0.0, 1.0, 0)", 5)
    g.codeonly("ppcc_plot normal family", "statistics::ppcc_plot(out, z[..], 0.0, 1.0, 5, 0)", 5)
    g.codeonly("ppcc_plot constant sample", "statistics::ppcc_plot(out, flat[..], 0.0, 1.0, 5)", 7)
    g.code("ppcc_max one value", "statistics::ppcc_max(one[..])", 2)
    g.code("ppcc_max normal family", "statistics::ppcc_max(z[..], 0)", 5)
    g.code("ppcc_max constant sample", "statistics::ppcc_max(flat[..])", 7)
    g.write()


AS = [0.9, 1.4, 2.1, 0.3, 3.8, 1.1, 0.7, 2.6, 1.9, 5.2, 0.5, 1.6]


def anderson_group():
    """The rows of the Anderson-Darling test against four more families."""
    from scipy import stats as st
    g = G("anderson_families", SCRIPT)
    s2 = (np.random.default_rng(4).weibull(1.8, size=40) * 3.1).tolist()
    g.decl("s", AS)
    g.decl("s2", s2)
    g.raw("var crit: f64[>..] = [];")
    g.raw("var lev: f64[>..] = [];")
    g.raw("var par: f64[>..] = [];")

    def row(label, name, data, fam, dist, tol):
        r = st.anderson(data, dist=dist)
        pv = st.anderson(data, dist=dist, method="interpolate").pvalue
        fit = list(r.fit_result.params)
        g.n += 1
        k = g.n
        g.lines.append(f"    let t{k}, c{k} = statistics::anderson_family(crit, lev, par, {name}[..], {fam});")
        g.lines.append(f'    chk_val(bad, "{label} statistic", t{k}.statistic, c{k}, {f(r.statistic)}, {tol});')
        g.lines.append(f'    chk_val(bad, "{label} p-value", t{k}.pvalue, c{k}, {f(pv)}, {tol});')
        g.lines.append(f"    let wc{k}: f64[{len(r.critical_values)}] = {arr(r.critical_values.tolist())};")
        g.lines.append(f"    let wl{k}: f64[{len(r.critical_values)}] = {arr(r.significance_level.tolist())};")
        g.lines.append(f"    let wp{k}: f64[{len(fit)}] = {arr(fit)};")
        g.lines.append(f'    chk_list(bad, "{label} critical values", crit[..], c{k}, wc{k}[..], 1e-12);')
        g.lines.append(f'    chk_list(bad, "{label} levels", lev[..], c{k}, wl{k}[..], 1e-12);')
        g.lines.append(f'    chk_list(bad, "{label} parameters", par[..], c{k}, wp{k}[..], {tol});')
    out_hi = [0.1, 0.2, 0.2, 0.3, 0.3, 0.3, 0.4, 0.4, 0.5, 0.5, 90.0]
    g.decl("sh", out_hi)
    row("anderson exponential", "s", AS, 1, "expon", "1e-12")
    row("anderson exponential beyond the table", "sh", out_hi, 1, "expon", "1e-9")
    row("anderson logistic", "s", AS, 2, "logistic", "2e-5")
    row("anderson gumbel left", "s", AS, 4, "gumbel_l", "1e-6")
    row("anderson gumbel right", "s", AS, 5, "gumbel_r", "1e-6")
    row("anderson weibull", "s2", s2, 7, "weibull_min", "1e-5")
    row("anderson gumbel alias", "s", AS, 3, "gumbel", "1e-6")
    row("anderson extreme1 alias", "s", AS, 6, "extreme1", "1e-6")
    g.raw("let one: f64[1] = [1.0];")
    g.raw("let neg: f64[3] = [-1.0, 1.0, 2.0];")
    g.raw("let flat: f64[4] = [2.0, 2.0, 2.0, 2.0];")
    g.raw("let nan3: f64[3] = [1.0, 0.0 / 0.0, 3.0];")
    g.raw("let none: f64[0] = [];")
    g.code("anderson empty", "statistics::anderson_family(crit, lev, par, none[..], 1)", 1)
    g.code("anderson one value", "statistics::anderson_family(crit, lev, par, one[..], 2)", 2)
    g.code("anderson NaN", "statistics::anderson_family(crit, lev, par, nan3[..], 1)", 3)
    g.code("anderson negative value", "statistics::anderson_family(crit, lev, par, neg[..], 1)", 5)
    g.code("anderson normal family", "statistics::anderson_family(crit, lev, par, s[..], 0)", 5)
    g.code("anderson unknown family", "statistics::anderson_family(crit, lev, par, s[..], 99)", 5)
    g.code("anderson constant logistic", "statistics::anderson_family(crit, lev, par, flat[..], 2)", 7)
    g.code("anderson constant gumbel", "statistics::anderson_family(crit, lev, par, flat[..], 5)", 7)
    g.code("anderson weibull with no interior fit", "statistics::anderson_family(crit, lev, par, s[..], 7)", 10)
    g.write()


def fitting_mle_group():
    """The rows of the goodness-of-fit test against gamma, beta and Weibull."""
    from scipy import stats as st
    g = G("fitting_mle", SCRIPT)
    rng = np.random.default_rng(11)
    gm = rng.gamma(2.5, 1.7, size=60).tolist()
    bt = rng.beta(2.0, 3.5, size=60).tolist()
    wb = (rng.weibull(1.8, size=60) * 3.1).tolist()
    g.decl("gm", gm)
    g.decl("bt", bt)
    g.decl("wb", wb)
    g.raw("var rg = rng { s: 7 };")
    rows = (("gamma", 5, "gm", gm, st.gamma, {"loc": 0}, 1e-6),
            ("beta", 6, "bt", bt, st.beta, {"loc": 0, "scale": 1}, 2e-3),
            ("weibull", 4, "wb", wb, st.weibull_min, {"loc": 0}, 1e-4))
    for nm, fam, var, data, dist, fixed, tol in rows:
        for sname, sc, stat in (("ad", 0, "ad"), ("ks", 1, "ks"), ("cvm", 2, "cvm")):
            r = st.goodness_of_fit(dist, data, known_params=fixed, statistic=stat,
                                   n_mc_samples=99, rng=np.random.default_rng(3))
            params = [float(v) for v in r.fit_result.params]
            g.n += 1
            k = g.n
            g.lines.append(f"    let t{k}, c{k} = statistics::goodness_of_fit(out, {var}[..], {fam}, {sc}, 19, rg);")
            g.lines.append(f'    chk_val(bad, "goodness_of_fit {nm} {sname} statistic", t{k}.statistic, c{k}, {f(r.statistic)}, {tol});')
            g.lines.append(f'    check(bad, c{k} == 0 && t{k}.pvalue > 0.0 && t{k}.pvalue <= 1.0, "goodness_of_fit {nm} {sname} p-value");')
            if sname == "ad":
                g.lines.append(f"    let w{k}: f64[{len(params)}] = {arr(params)};")
                g.lines.append(f'    chk_list(bad, "goodness_of_fit {nm} parameters", out[..], c{k}, w{k}[..], {tol});')
    g.raw("let neg: f64[4] = [-1.0, 1.0, 2.0, 3.0];")
    g.raw("let big: f64[4] = [0.2, 0.4, 1.5, 0.6];")
    g.raw("let flat: f64[4] = [2.0, 2.0, 2.0, 2.0];")
    g.code("goodness_of_fit gamma negative value", "statistics::goodness_of_fit(out, neg[..], 5, 0, 19, rg)", 5)
    g.code("goodness_of_fit weibull negative value", "statistics::goodness_of_fit(out, neg[..], 4, 0, 19, rg)", 5)
    g.code("goodness_of_fit beta value above one", "statistics::goodness_of_fit(out, big[..], 6, 0, 19, rg)", 5)
    g.code("goodness_of_fit gamma constant sample", "statistics::goodness_of_fit(out, flat[..], 5, 0, 19, rg)", 7)
    g.code("goodness_of_fit chi-square family", "statistics::goodness_of_fit(out, gm[..], 7, 0, 19, rg)", 5)
    g.write()


def somers_group():
    """The rows of somersd on two samples and on a table."""
    from scipy import stats as st
    g = G("assoc_somers", SCRIPT)
    g.decl("x", X)
    g.decl("y", Y)
    g.decl("xt", XT)
    g.decl("yt", YT)
    tab = [[27, 25, 14, 7], [7, 9, 22, 12], [4, 11, 17, 28]]
    flat = [float(v) for row in tab for v in row]
    g.decl("tab", flat)
    tz = [[5, 0, 2], [0, 0, 0], [1, 4, 3]]
    g.decl("tz", [float(v) for row in tz for v in row])
    alts = (("two-sided", "two-sided", 0), ("less", "less", 1), ("greater", "greater", 2))
    for nm, xs, ys in (("ties", "xt", "yt"), ("distinct", "x", "y")):
        data = (XT, YT) if nm == "ties" else (X, Y)
        for label, alt, code in alts:
            r = st.somersd(data[0], data[1], alternative=alt)
            g.field(f"somersd {nm} {label}", f"statistics::somersd({xs}[..], {ys}[..], {code})",
                    {"statistic": float(r.statistic), "pvalue": float(r.pvalue)}, "1e-12")
    for label, alt, code in alts:
        r = st.somersd(tab, alternative=alt)
        g.field(f"somersd_table {label}", f"statistics::somersd_table(tab[..], 3, 4, {code})",
                {"statistic": float(r.statistic), "pvalue": float(r.pvalue)}, "1e-12")
    r = st.somersd(tz)
    g.field("somersd_table with an empty row", "statistics::somersd_table(tz[..], 3, 3)",
            {"statistic": float(r.statistic), "pvalue": float(r.pvalue)}, "1e-12")
    g.raw("let none: f64[0] = [];")
    g.raw("let one: f64[1] = [1.0];")
    g.raw("let ones: f64[4] = [1.0, 1.0, 1.0, 1.0];")
    g.raw("let four: f64[4] = [1.0, 2.0, 3.0, 4.0];")
    g.raw("let nan4: f64[4] = [1.0, 0.0 / 0.0, 3.0, 4.0];")
    g.raw("let inf4: f64[4] = [1.0, 1.0 / 0.0, 3.0, 4.0];")
    g.raw("let negt: f64[4] = [1.0, -1.0, 2.0, 3.0];")
    g.raw("let fract: f64[4] = [1.0, 1.5, 2.0, 3.0];")
    g.raw("let single: f64[2] = [3.0, 4.0];")
    g.raw("let lone: f64[4] = [4.0, 0.0, 0.0, 0.0];")
    g.code("somersd empty", "statistics::somersd(none[..], none[..])", 1)
    g.code("somersd one value", "statistics::somersd(one[..], one[..])", 2)
    g.code("somersd length", "statistics::somersd(four[..], single[..])", 6)
    g.code("somersd NaN", "statistics::somersd(nan4[..], four[..])", 3)
    g.code("somersd infinity", "statistics::somersd(four[..], inf4[..])", 4)
    g.code("somersd constant first sample", "statistics::somersd(ones[..], four[..])", 7)
    g.code("somersd bad alternative", "statistics::somersd(four[..], four[..], 9)", 5)
    g.code("somersd_table empty", "statistics::somersd_table(none[..], 0, 0)", 1)
    g.code("somersd_table negative count", "statistics::somersd_table(negt[..], 2, 2)", 5)
    g.code("somersd_table fractional count", "statistics::somersd_table(fract[..], 2, 2)", 5)
    g.code("somersd_table shape", "statistics::somersd_table(four[..], 3, 2)", 6)
    g.code("somersd_table NaN", "statistics::somersd_table(nan4[..], 2, 2)", 3)
    g.code("somersd_table single row", "statistics::somersd_table(four[..], 1, 4)", 7)
    g.code("somersd_table one nonzero cell", "statistics::somersd_table(lone[..], 2, 2)", 5)
    g.write()


def risk_group():
    """The rows of the odds ratio and the relative risk."""
    from scipy.stats.contingency import odds_ratio, relative_risk
    g = G("contingency_risk", SCRIPT)
    tables = {
        "drug": [7, 15, 58, 472],
        "small": [1, 2, 3, 4],
        "balanced": [5, 5, 5, 5],
        "big": [120, 80, 200, 300],
    }
    for nm, tab in tables.items():
        g.decl(nm, [float(v) for v in tab])
    alts = (("two-sided", 0, "two-sided"), ("less", 1, "less"), ("greater", 2, "greater"))
    for nm, tab in tables.items():
        for kind, kcode in (("sample", "statistics::ODDS_SAMPLE"), ("conditional", "statistics::ODDS_CONDITIONAL")):
            for label, acode, aname in alts:
                r = odds_ratio([tab[:2], tab[2:]], kind=kind)
                ci = r.confidence_interval(0.95, alternative=aname)
                fields = {"statistic": float(r.statistic), "low": float(ci.low)}
                if acode != 2:
                    fields["high"] = float(ci.high)
                tol = "1e-12" if kind == "sample" else "1e-9"
                g.field(f"odds {nm} {kind} {label}",
                        f"statistics::odds_ratio({nm}[..], {kcode}, {acode}, 0.95)", fields, tol)
    r = odds_ratio([[7, 15], [58, 472]], kind="sample")
    ci = r.confidence_interval(0.99)
    g.field("odds drug sample level 0.99", "statistics::odds_ratio(drug[..], statistics::ODDS_SAMPLE, 0, 0.99)",
            {"statistic": float(r.statistic), "low": float(ci.low), "high": float(ci.high)}, "1e-12")
    g.decl("zero_cell", [0.0, 5.0, 5.0, 0.0])
    g.decl("diagonal", [10.0, 0.0, 0.0, 10.0])
    r = odds_ratio([[0, 5], [5, 0]], kind="conditional")
    ci = r.confidence_interval(0.95)
    g.field("odds zero cell conditional",
            "statistics::odds_ratio(zero_cell[..], statistics::ODDS_CONDITIONAL, 0, 0.95)",
            {"statistic": float(r.statistic), "low": float(ci.low), "high": float(ci.high)}, "1e-9")
    g.code("odds zero cell sample", "statistics::odds_ratio(zero_cell[..], statistics::ODDS_SAMPLE, 0, 0.95)", 7)
    g.code("odds diagonal sample", "statistics::odds_ratio(diagonal[..], statistics::ODDS_SAMPLE, 0, 0.95)", 8)
    g.code("odds diagonal conditional", "statistics::odds_ratio(diagonal[..], statistics::ODDS_CONDITIONAL, 0, 0.95)", 8)
    for (ec, et, cc, ct) in ((27, 122, 44, 487), (10, 75, 12, 225), (2, 10, 4, 20)):
        r = relative_risk(ec, et, cc, ct)
        ci = r.confidence_interval(0.95)
        g.field(f"relative risk {ec} {et} {cc} {ct}",
                f"statistics::relative_risk({float(ec)}, {float(et)}, {float(cc)}, {float(ct)}, 0.95)",
                {"statistic": float(r.relative_risk), "low": float(ci.low), "high": float(ci.high)}, "1e-12")
    g.code("relative risk of two empty groups", "statistics::relative_risk(0.0, 10.0, 0.0, 20.0, 0.95)", 7)
    g.code("relative risk of an infinite ratio", "statistics::relative_risk(5.0, 10.0, 0.0, 20.0, 0.95)", 8)
    g.code("relative risk of no cases", "statistics::relative_risk(0.0, 10.0, 3.0, 20.0, 0.95)", 7)
    g.write()


def nd_group():
    """The rows of the distance of two sets of points."""
    g = G("distance_nd", SCRIPT)
    rng = np.random.default_rng(23)
    cases = {
        "p1": (np.array([[0, 0], [1, 0], [0, 2], [3, 3]], float), np.array([[1, 1], [2, 0], [4, 4]], float)),
        "line": (np.array([[0.0], [1.0], [3.0]]), np.array([[0.5], [2.0]])),
        "rand2": (rng.normal(size=(12, 2)) * 2.0, rng.normal(size=(9, 2)) + 1.0),
        "rand3": (rng.normal(size=(8, 3)), rng.normal(size=(8, 3)) * 1.5),
        "dup": (np.array([[1, 1], [1, 1], [2, 2], [5, 0]], float), np.array([[1, 1], [2, 2], [2, 2]], float)),
        "big": (rng.normal(size=(15, 4)) * 1000.0, rng.normal(size=(20, 4)) * 1000.0 + 300.0),
    }
    for nm, (a, b) in cases.items():
        d = a.shape[1]
        g.decl(nm + "u", a.ravel().tolist())
        g.decl(nm + "v", b.ravel().tolist())
        g.val(f"wasserstein_distance_nd {nm}", f"statistics::wasserstein_distance_nd({nm}u[..], {nm}v[..], {d})",
              float(stats.wasserstein_distance_nd(a, b)))
    a, b = cases["p1"]
    wu, wv = [1.0, 2.0, 1.0, 1.0], [2.0, 1.0, 1.0]
    g.decl("wu", wu)
    g.decl("wv", wv)
    g.val("wasserstein_distance_nd weights", "statistics::wasserstein_distance_nd_weighted(p1u[..], p1v[..], 2, wu[..], wv[..])",
          float(stats.wasserstein_distance_nd(a, b, wu, wv)))
    zu, zv = [1.0, 0.0, 3.0, 2.0], [0.0, 1.0, 5.0]
    g.decl("zu", zu)
    g.decl("zv", zv)
    g.val("wasserstein_distance_nd zero weights", "statistics::wasserstein_distance_nd_weighted(p1u[..], p1v[..], 2, zu[..], zv[..])",
          float(stats.wasserstein_distance_nd(a, b, zu, zv)))
    a, b = cases["rand2"]
    rw, rv = rng.random(12) + 0.1, rng.random(9) + 0.1
    g.decl("rwu", rw.tolist())
    g.decl("rwv", rv.tolist())
    g.val("wasserstein_distance_nd random weights", "statistics::wasserstein_distance_nd_weighted(rand2u[..], rand2v[..], 2, rwu[..], rwv[..])",
          float(stats.wasserstein_distance_nd(a, b, rw, rv)))
    g.raw("let none: f64[0] = [];")
    g.raw("let nan2: f64[4] = [0.0, 0.0 / 0.0, 1.0, 1.0];")
    g.raw("let inf2: f64[4] = [0.0, 1.0 / 0.0, 1.0, 1.0];")
    g.raw("let odd3: f64[3] = [0.0, 1.0, 2.0];")
    g.raw("let negw: f64[2] = [1.0, -1.0];")
    g.raw("let okw: f64[2] = [1.0, 1.0];")
    g.raw("let zerow: f64[2] = [0.0, 0.0];")
    g.raw("let nanw: f64[2] = [1.0, 0.0 / 0.0];")
    g.raw("let infw: f64[2] = [1.0, 1.0 / 0.0];")
    g.raw("let huge: f64[4] = [1e200, 0.0, -1e200, 0.0];")
    g.raw("let orig: f64[4] = [0.0, 0.0, 0.0, 0.0];")
    g.code("wasserstein_distance_nd empty", "statistics::wasserstein_distance_nd(none[..], nan2[..], 2)", 1)
    g.code("wasserstein_distance_nd NaN", "statistics::wasserstein_distance_nd(nan2[..], inf2[..], 2)", 3)
    g.code("wasserstein_distance_nd infinity", "statistics::wasserstein_distance_nd(orig[..], inf2[..], 2)", 4)
    g.code("wasserstein_distance_nd a length that dims does not divide", "statistics::wasserstein_distance_nd(odd3[..], orig[..], 2)", 6)
    g.code("wasserstein_distance_nd a second set that dims does not divide", "statistics::wasserstein_distance_nd(orig[..], odd3[..], 2)", 6)
    g.code("wasserstein_distance_nd zero dims", "statistics::wasserstein_distance_nd(orig[..], orig[..], 0)", 5)
    g.code("wasserstein_distance_nd too many dims", "statistics::wasserstein_distance_nd(orig[..], orig[..], 17)", 9)
    g.code("wasserstein_distance_nd overflow", "statistics::wasserstein_distance_nd(huge[..], orig[..], 2)", 8)
    g.code("wasserstein_distance_nd_weighted negative weight", "statistics::wasserstein_distance_nd_weighted(orig[..], orig[..], 2, negw[..], okw[..])", 5)
    g.code("wasserstein_distance_nd_weighted weight length", "statistics::wasserstein_distance_nd_weighted(orig[..], orig[..], 2, okw[..], wu[..])", 6)
    g.code("wasserstein_distance_nd_weighted all zero", "statistics::wasserstein_distance_nd_weighted(orig[..], orig[..], 2, zerow[..], okw[..])", 7)
    g.code("wasserstein_distance_nd_weighted NaN weight", "statistics::wasserstein_distance_nd_weighted(orig[..], orig[..], 2, nanw[..], okw[..])", 3)
    g.code("wasserstein_distance_nd_weighted infinite weight", "statistics::wasserstein_distance_nd_weighted(orig[..], orig[..], 2, infw[..], okw[..])", 4)
    g.write()


def mgc_group():
    """The rows of the multiscale graph correlation."""
    from scipy.spatial.distance import cdist
    g = G("assoc_mgc", SCRIPT)
    rng = np.random.default_rng(5)
    xl = rng.normal(size=40)
    yl = xl ** 2 + rng.normal(size=40) * 0.5
    g.raw("var map: f64[>..] = [];")
    g.raw("var rg = rng { s: 3 };")

    def rows(label, xn, yn, xd, yd, x, y, reps=0, **kw):
        r = stats.multiscale_graphcorr(np.asarray(x, float).reshape(len(x) // xd, xd) if np.ndim(x) == 1 else x,
                                       np.asarray(y, float).reshape(len(y) // yd, yd) if np.ndim(y) == 1 else y,
                                       reps=reps, **kw)
        k = g.n = g.n + 1
        g.lines.append(f"    let r{k}, c{k} = statistics::multiscale_graphcorr(map, {xn}[..], {yn}[..], {xd}, {yd}, {reps}, rg);")
        emit(label, k, r)

    def emit(label, k, r):
        m = np.asarray(r.mgc_dict["mgc_map"])
        sc = r.mgc_dict["opt_scale"]
        kx, ky = (int(sc[0]), int(sc[1])) if np.ndim(sc) else (m.shape[0], m.shape[1])
        g.lines.append(f'    chk_val(bad, "{label} statistic", r{k}.statistic, c{k}, {f(float(r.statistic))}, 1e-12);')
        g.lines.append(f'    chk_val(bad, "{label} p-value", r{k}.pvalue, c{k}, {f(float(r.pvalue))}, 1e-12);')
        g.lines.append(f'    chk_code(bad, "{label} scale x", r{k}.kx, {kx});')
        g.lines.append(f'    chk_code(bad, "{label} scale y", r{k}.ky, {ky});')
        g.lines.append(f'    chk_code(bad, "{label} rows", r{k}.rows, {m.shape[0]});')
        g.lines.append(f'    chk_code(bad, "{label} columns", r{k}.cols, {m.shape[1]});')
        g.lines.append(f"    let w{k}: f64[{m.size}] = {arr(m.ravel().tolist())};")
        g.lines.append(f'    chk_list(bad, "{label} map", map[..], c{k}, w{k}[..], 1e-12);')

    cases = {
        "ten": (X, Y, 1, 1),
        "twenty": (xl[:20].tolist(), yl[:20].tolist(), 1, 1),
        "forty": (xl.tolist(), yl.tolist(), 1, 1),
        "ties": ([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0], [1.0, 4.0, 9.0, 16.0, 25.0, 36.0, 49.0, 64.0], 1, 1),
        "repeats": ([1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0, 5.0, 5.0], [2.0, 1.0, 4.0, 3.0, 6.0, 5.0, 8.0, 7.0, 9.0, 10.0], 1, 1),
        # a sample whose statistic depends on the level of the threshold
        "level": ([-4.3304350561105, -6.889742312428999, -4.69365864339732, -5.197933378593724, -3.0397183648127886, -0.3806393574940867, -3.728583555963779, -5.362923650942799, -6.0047557787585415, -4.176482257730299, -4.374489517570227, -5.973533065149551, -3.8935778351601718, -5.97253951217548, -2.260040898901467, -2.3411748435334996, -2.304918584629819, -4.8621417885500415, -3.240385327512395, -4.3011200499615025, -4.7223076885062625, -4.640829666956407, -6.216649640824947, -6.453125956739864, -2.7813812896327685, -4.398182046719554], [-0.1300184447764965, 0.2850012919542796, 0.70459650949732, 0.7993904241925855, 0.46301571414167675, 1.078001535356247, 0.7199290642125717, 0.09727987519845249, 0.22087513261797215, 0.9801690065257207, 0.8307624556493003, -0.20258275242065027, 1.217303088631079, 0.8894491468657377, 1.215392290086881, 0.4601178295290064, 0.4985998393788531, 0.135294122846557, 1.7373740871209176, 0.5511886953978938, 1.3221820643512245, 0.34484340907276656, 0.6995700817872275, -0.1030006814202602, 0.4604824390972973, 1.058136474224566], 1, 1),
        # a sample whose region of large correlations is too small: the
        # statistic stays the global one
        "region": ([1.8738096985117645, 1.7067342334547426, 2.078696565792778, 2.0942927903661843, 1.5778924771693161, 1.774200054473722, 2.728654452182945, 1.9050610279702205, 1.653910543085179, 2.5781938699501796, 1.7346512683912105], [-3.27491314676858, -2.557543397043428, -1.0901240796991987, -1.8974558806870485, -3.0319797632943684, -2.6532881822108356, -2.7288722121312654, 0.40636246837172685, -2.9562763841105433, -2.691174216282487, 1.9604142071414552], 1, 1),
        # a sample whose statistic depends on the shape of the threshold
        "shape": ([1.6996752418028578, 1.7261057556395314, 1.6132857200779325, 1.978562500855361, 2.485647305996454, 3.327258751386191, 1.19946831592581, 1.7521444779818878], [4.159483638064041, 3.805858839494032, 3.910727159555692, 5.677836623010066, 6.543851248512019, 6.33057517025161, 5.528267523067718, 5.9463793763532795], 1, 1),
    }
    for nm, (x, y, xd, yd) in cases.items():
        g.decl(nm + "x", x)
        g.decl(nm + "y", y)
        rows(f"multiscale_graphcorr {nm}", nm + "x", nm + "y", xd, yd, x, y)
    # a point cloud of two dimensions against one dimension
    p2 = np.round(rng.normal(size=(14, 2)) * 2, 3)
    q1 = np.round(p2[:, 0] ** 2 + rng.normal(size=14) * 0.3, 3)
    g.decl("p2", p2.ravel().tolist())
    g.decl("q1", q1.tolist())
    rows("multiscale_graphcorr two dimensions", "p2", "q1", 2, 1, p2, q1.reshape(-1, 1))
    # the same sample reps 1 gives a p-value in the set {1/2, 1}
    # the pooled test of two samples of different sizes
    u = np.round(rng.normal(size=(9, 2)), 3)
    v = np.round(rng.normal(size=(12, 2)) + 0.8, 3)
    g.decl("tu", u.ravel().tolist())
    g.decl("tv", v.ravel().tolist())
    r = stats.multiscale_graphcorr(u, v, reps=0)
    k = g.n = g.n + 1
    g.lines.append(f"    let r{k}, c{k} = statistics::multiscale_graphcorr_twosample(map, tu[..], tv[..], 2, 0, rg);")
    emit("multiscale_graphcorr_twosample", k, r)
    # the form with distance tables
    a = cdist(np.asarray(X).reshape(-1, 1), np.asarray(X).reshape(-1, 1))
    b = cdist(np.asarray(Y).reshape(-1, 1), np.asarray(Y).reshape(-1, 1))
    g.decl("da", a.ravel().tolist())
    g.decl("db", b.ravel().tolist())
    r = stats.multiscale_graphcorr(a, b, compute_distance=None, reps=0)
    k = g.n = g.n + 1
    g.lines.append(f"    let r{k}, c{k} = statistics::multiscale_graphcorr_dist(map, da[..], db[..], 0, rg);")
    emit("multiscale_graphcorr_dist", k, r)
    # a distance table that is not symmetric
    asym = a + np.arange(10).reshape(-1, 1) * 0.0 + np.triu(np.ones((10, 10))) * 0.01 * a
    g.decl("dasym", asym.ravel().tolist())
    r = stats.multiscale_graphcorr(asym, b, compute_distance=None, reps=0)
    k = g.n = g.n + 1
    g.lines.append(f"    let r{k}, c{k} = statistics::multiscale_graphcorr_dist(map, dasym[..], db[..], 0, rg);")
    emit("multiscale_graphcorr_dist asymmetric", k, r)
    # a sample of one value: every distance is 0
    g.raw("let flat: f64[6] = [2.0, 2.0, 2.0, 2.0, 2.0, 2.0];")
    g.raw("let ramp: f64[6] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];")
    rows("multiscale_graphcorr a constant sample", "flat", "ramp", 1, 1, [2.0] * 6, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    # the faults
    g.raw("let none: f64[0] = [];")
    g.raw("let four: f64[4] = [1.0, 2.0, 3.0, 4.0];")
    g.raw("let five: f64[5] = [1.0, 2.0, 3.0, 4.0, 5.0];")
    g.raw("let six: f64[6] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];")
    g.raw("let nan6: f64[6] = [1.0, 0.0 / 0.0, 3.0, 4.0, 5.0, 6.0];")
    g.raw("let inf6: f64[6] = [1.0, 1.0 / 0.0, 3.0, 4.0, 5.0, 6.0];")
    g.raw("let huge: f64[6] = [1e200, -1e200, 0.0, 1.0, 2.0, 3.0];")
    calls = [
        ("empty", "none[..], none[..], 1, 1, 0", 1),
        ("a first sample of four points", "four[..], four[..], 1, 1, 0", 2),
        ("a length mismatch", "five[..], six[..], 1, 1, 0", 6),
        ("a width that does not divide", "six[..], six[..], 4, 1, 0", 6),
        ("a zero width", "six[..], six[..], 0, 1, 0", 5),
        ("a negative count of repetitions", "six[..], six[..], 1, 1, -1", 5),
        ("too many repetitions", "six[..], six[..], 1, 1, 10001", 9),
        ("NaN", "nan6[..], six[..], 1, 1, 0", 3),
        ("infinity", "six[..], inf6[..], 1, 1, 0", 4),
        ("overflow", "huge[..], six[..], 1, 1, 0", 8),
    ]
    for label, args, code in calls:
        g.code(f"multiscale_graphcorr {label}", f"statistics::multiscale_graphcorr(map, {args}, rg)", code)
    g.code("multiscale_graphcorr_dist a table that is not a square", "statistics::multiscale_graphcorr_dist(map, six[..], six[..], 0, rg)", 6)
    g.code("multiscale_graphcorr_dist empty", "statistics::multiscale_graphcorr_dist(map, none[..], none[..], 0, rg)", 1)
    g.code("multiscale_graphcorr_twosample a short sample", "statistics::multiscale_graphcorr_twosample(map, four[..], six[..], 1, 0, rg)", 2)
    g.write()


def arr_f(x):
    return "[" + ", ".join(repr(float(v)) for v in x) + "]"


if __name__ == "__main__":
    main()
    line_group()
    laws_group()
    probplot_group()
    anderson_group()
    fitting_mle_group()
    somers_group()
    risk_group()
    nd_group()
    mgc_group()
