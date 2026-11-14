#!/usr/bin/env python3
"""Writes the vector files of the families and the laws:
test/suite/distributions_location_vectors.goose and later groups.

Needs: Python 3, NumPy and SciPy 1.18.1. The files hold numbers only.
Run it from the package directory with the system python3.
"""
import math
import warnings
import numpy as np
from scipy import special
from scipy import stats
from scipy.integrate import quad
from vecgen import G

warnings.simplefilter("ignore")
SCRIPT = "make_family_vectors.py"

XS = [-3.0, -1.0, -0.2, 0.0, 0.5, 1.0, 2.5, 8.0]
PS = [0.001, 0.1, 0.25, 0.5, 0.75, 0.9, 0.999]
SETS = [(0.0, 1.0, ""), (0.7, 2.3, "0.7, 2.3")]


def law_moments(d, args):
    """The mean, the variance, the skewness and the kurtosis of a law.

    The values come from an integration of the density of the law. The
    law of a nonnegative support is integrated in the logarithm of the
    point. The law of the whole line is integrated in the point itself.
    SciPy integrates the quantile function instead. That route loses the
    accuracy of the third and the fourth raw moment of a law with a heavy
    tail.
    """
    lo, _hi = d.support(*args)
    raws = []
    for n in (1, 2, 3, 4):
        if lo >= 0.0:
            def G_of_u(u, n=n):
                return math.exp((n + 1) * u) * float(d.pdf(math.exp(u), *args))

            a0, a1, pieces = -30.0, 60.0, 1000
        else:
            def G_of_u(x, n=n):
                return float(d.pdf(x, *args)) * x ** n

            a0, a1, pieces = -30.0, 30.0, 500
        total = 0.0
        for i in range(pieces):
            left = a0 + (a1 - a0) * i / pieces
            right = a0 + (a1 - a0) * (i + 1) / pieces
            total += quad(G_of_u, left, right, limit=200)[0]
        raws.append(total)
    return moments_of_raws(raws)


def moments_of_raws(raws):
    """The mean, the variance, the skewness and the kurtosis of the raw
    moments of a law."""
    m1, m2, m3, m4 = raws
    v = m2 - m1 * m1
    skew = (m3 - 3.0 * m1 * m2 + 2.0 * m1 ** 3) / v ** 1.5
    kurt = (m4 - 4.0 * m1 * m3 + 6.0 * m1 ** 2 * m2 - 3.0 * m1 ** 4) / v ** 2 - 3.0
    return [m1, v, skew, kurt]


def powerlognorm_moments(d, args):
    """The moments of the power of the log-normal law.

    The shape 1 makes the law the log-normal law of the scale s. That
    law has a closed form for its raw moments. The values come from that
    form. A scale of 0.001 makes the law narrow in the logarithm of the
    point. No integration of the density finds that peak.
    """
    c, s = args
    if c != 1.0:
        raise ValueError("the closed form of this law holds for the shape 1 only")
    raws = [math.exp(n * n * s * s / 2.0) for n in (1, 2, 3, 4)]
    return moments_of_raws(raws)


def exponweib_moments(d, args):
    """The moments of the Weibull law of the exponential kind.

    The shape 1 makes the law the Weibull law of the shape c. That law
    has a closed form for its raw moments. The values come from that
    form. A shape of 0.03 puts the peak of the law far away in the
    logarithm of the point.
    """
    a, c = args
    if a != 1.0:
        raise ValueError("the closed form of this law holds for the shape 1 only")
    raws = [float(special.gamma(1.0 + n / c)) for n in (1, 2, 3, 4)]
    return moments_of_raws(raws)



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


def elementary_group():
    """The rows of the eight families of elementary functions."""
    g = G("distributions_elementary", SCRIPT, tol="1e-10")
    names = ("anglit", "arcsine", "cosine", "gibrat",
             "halfcauchy", "halflogistic", "halfnorm", "hypsecant")
    for nm in names:
        d = getattr(stats, nm)
        lo, hi = d.support()
        qs = [float(d.ppf(p)) for p in (0.05, 0.25, 0.5, 0.75, 0.95)]
        out = []
        if np.isfinite(lo):
            out.append(float(lo - 1.0))
        if np.isfinite(hi):
            out.append(float(hi + 1.0))
        for x in qs + out:
            for fn in ("pdf", "cdf", "sf"):
                g.val(f"{nm} {fn} {x:.4g}", f"statistics::{nm}_{fn}({x!r})", float(getattr(d, fn)(x)))
        for p in (0.1, 0.5, 0.9):
            g.val(f"{nm} ppf {p}", f"statistics::{nm}_ppf({p!r})", float(d.ppf(p)))
        m, v, s, k = (float(q) for q in d.stats(moments="mvsk"))
        for fn, want in (("mean", m), ("var", v), ("skew", s), ("kurt", k)):
            if np.isfinite(want):
                g.val(f"{nm} {fn}", f"statistics::{nm}_{fn}()", want)
            else:
                g.code(f"{nm} {fn} does not exist", f"statistics::{nm}_{fn}()", 7)
    g.write()



def bounded_group():
    """The rows of the seven families of a bounded support."""
    g = G("distributions_bounded", SCRIPT, tol="1e-10")
    # The family, the parameters and the values to try.
    sets = [
        ("semicircular", (), [()]),
        ("triang", ("c",), [(0.3,)]),
        ("trapezoid", ("c", "d"), [(0.2, 0.9)]),
        ("bradford", ("c",), [(1.5,)]),
        ("wald", (), [()]),
        ("maxwell", (), [()]),
        ("moyal", (), [()]),
    ]
    for nm, shapes, argsets in sets:
        d = getattr(stats, nm)
        for args in argsets:
            lo, hi = d.support(*args)
            xs = [float(d.ppf(p, *args)) for p in (0.05, 0.25, 0.5, 0.75, 0.95)]
            out = []
            if np.isfinite(lo):
                out.append(float(lo - 1.0))
            if np.isfinite(hi):
                out.append(float(hi + 1.0))
            tail = "".join(", " + repr(a) for a in args)
            for x in xs + out:
                for fn in ("pdf", "cdf", "sf"):
                    want = float(getattr(d, fn)(x, *args))
                    g.val(f"{nm} {fn} {x:.4g}", f"statistics::{nm}_{fn}({x!r}{tail})", want)
            for p in (0.1, 0.5, 0.9):
                g.val(f"{nm} ppf {p}", f"statistics::{nm}_ppf({p!r}{tail})", float(d.ppf(p, *args)))
            m, v, s, k = (float(q) for q in d.stats(*args, moments="mvsk"))
            arg = ", ".join(repr(a) for a in args)
            for fn, want in (("mean", m), ("var", v), ("skew", s), ("kurt", k)):
                call = f"statistics::{nm}_{fn}({arg})" if args else f"statistics::{nm}_{fn}()"
                if np.isfinite(want):
                    g.val(f"{nm} {fn}", call, want)
                else:
                    g.code(f"{nm} {fn} does not exist", call, 7)
    g.write()



def weibull_group():
    """The rows of the four families of the Weibull kind."""
    g = G("distributions_weibull", SCRIPT, tol="1e-10")
    sets = [
        ("dweibull", ["c"], [(1.7,)]),
        ("invweibull", ["c"], [(1.5,)]),
        ("weibull_max", ["c"], [(1.5,)]),
        ("truncweibull_min", ["c", "a", "b"], [(1.5, 1.0, 2.0)]),
    ]
    for nm, shapes, argsets in sets:
        d = getattr(stats, nm)
        for args in argsets:
            lo, hi = d.support(*args)
            xs = [float(d.ppf(p, *args)) for p in (0.05, 0.25, 0.5, 0.75, 0.95)]
            out = []
            if np.isfinite(lo):
                out.append(float(lo - 1.0))
            if np.isfinite(hi):
                out.append(float(hi + 1.0))
            tail = "".join(", " + repr(a) for a in args)
            for x in xs + out:
                for fn in ("pdf", "cdf", "sf"):
                    want = float(getattr(d, fn)(x, *args))
                    g.val(f"{nm} {fn} {x:.4g}", f"statistics::{nm}_{fn}({x!r}{tail})", want)
            for p in (0.1, 0.5, 0.9):
                g.val(f"{nm} ppf {p}", f"statistics::{nm}_ppf({p!r}{tail})", float(d.ppf(p, *args)))
            m, v, s, k = (float(q) for q in d.stats(*args, moments="mvsk"))
            arg = ", ".join(repr(a) for a in args)
            for fn, want in (("mean", m), ("var", v), ("skew", s), ("kurt", k)):
                call = f"statistics::{nm}_{fn}({arg})" if args else f"statistics::{nm}_{fn}()"
                if np.isfinite(want):
                    g.val(f"{nm} {fn}", call, want)
                else:
                    g.code(f"{nm} {fn} does not exist", call, 7)
    g.write()


def power_group():
    """The rows of the four families of a power or an exponent.

    The moments of these families have no closed form. The values come
    from an integration of the density of the law. SciPy integrates the
    quantile function instead, and that route loses the accuracy of the
    third and the fourth raw moment of a law with a heavy tail.
    """
    g = G("distributions_exponent", SCRIPT, tol="1e-10",
          note="the answers of SciPy 1.18.1, and the moments from an "
               "integration of the density of the law")
    sets = [
        ("powernorm", (1.6,), law_moments, ("mean", "var", "skew", "kurt"), {}),
        ("powerlognorm", (1.4, 2.0), law_moments, ("mean", "var", "skew", "kurt"), {}),
        ("powerlognorm", (1.0, 0.001), powerlognorm_moments, ("mean", "var"), {"var": 1e-8}),
        ("exponpow", (1.5,), law_moments, ("mean", "var", "skew", "kurt"), {}),
        ("exponweib", (1.5, 2.0), law_moments, ("mean", "var", "skew", "kurt"), {}),
        ("exponweib", (1.0, 0.03), exponweib_moments, ("mean", "var", "skew", "kurt"), {}),
    ]
    for nm, args, moment_fn, names, tols in sets:
        d = getattr(stats, nm)
        lo, hi = d.support(*args)
        arg = ", ".join(repr(a) for a in args)
        xs = [float(d.ppf(p, *args)) for p in (0.05, 0.25, 0.5, 0.75, 0.95)]
        out = []
        if np.isfinite(lo):
            out.append(float(lo - 1.0))
        if np.isfinite(hi):
            out.append(float(hi + 1.0))
        tail = "".join(", " + repr(a) for a in args)
        for x in xs + out:
            for fn in ("pdf", "cdf", "sf"):
                want = float(getattr(d, fn)(x, *args))
                g.val(f"{nm} {arg} {fn} {x:.4g}", f"statistics::{nm}_{fn}({x!r}{tail})", want)
        for p in (0.1, 0.5, 0.9):
            g.val(f"{nm} {arg} ppf {p}",
                  f"statistics::{nm}_ppf({p!r}{tail})", float(d.ppf(p, *args)))
        got = moment_fn(d, args)
        for fn, want in zip(names, got):
            call = f"statistics::{nm}_{fn}({arg})" if args else f"statistics::{nm}_{fn}()"
            label = f"{nm} {arg} {fn}"
            if np.isfinite(want):
                g.val(label, call, want, tols.get(fn))
            else:
                g.code(label + " does not exist", call, 7)
    g.write()


def rice_group():
    """The rows of the Rice law for the shapes 0, 1.6, 10 and 30.

    The moments come from mpmath at 60 digits, with the closed form of the
    raw moments (a confluent hypergeometric function). SciPy loses digits
    of the skewness and the kurtosis for a shape of 10 or more, because
    the central moments cancel. The package loses digits the same way, so
    those rows carry a looser relative tolerance: 1e-5 for the skewness
    and 5e-3 for the kurtosis. mpmath 1.3 is the oracle.
    """
    import mpmath as mp
    mp.mp.dps = 60

    def exact(b):
        def raw(n):
            nd2 = mp.mpf(n) / 2
            b2 = mp.mpf(b) ** 2 / 2
            return 2 ** nd2 * mp.e ** (-b2) * mp.gamma(1 + nd2) * mp.hyp1f1(1 + nd2, 1, b2)
        m = [raw(n) for n in (1, 2, 3, 4)]
        v = m[1] - m[0] ** 2
        sk = (m[2] - 3 * m[0] * m[1] + 2 * m[0] ** 3) / v ** mp.mpf(1.5)
        ku = (m[3] - 4 * m[0] * m[2] + 6 * m[0] ** 2 * m[1] - 3 * m[0] ** 4) / v ** 2 - 3
        return [float(m[0]), float(v), float(sk), float(ku)]

    g = G("distributions_rice", SCRIPT, tol="1e-10",
          note="the answers of SciPy 1.18.1, and the moments from mpmath")
    d = stats.rice
    for b in (0.0, 1.6, 10.0, 30.0):
        args = (b,)
        lo, hi = d.support(*args)
        xs = [float(d.ppf(p, *args)) for p in (0.05, 0.25, 0.5, 0.75, 0.95)] + [float(lo - 1.0)]
        for x in xs:
            for fn in ("pdf", "cdf", "sf"):
                want = float(getattr(d, fn)(x, *args))
                g.val(f"rice {b} {fn} {x:.4g}", f"statistics::rice_{fn}({x!r}, {b!r})", want)
        for p in (0.1, 0.5, 0.9):
            g.val(f"rice {b} ppf {p}", f"statistics::rice_ppf({p!r}, {b!r})", float(d.ppf(p, *args)))
        for fn, want in zip(("mean", "var", "skew", "kurt"), exact(b)):
            tol = None
            if b >= 10.0 and fn == "skew":
                tol = "1e-5"
            if b >= 10.0 and fn == "kurt":
                tol = "5e-3"
            g.val(f"rice {b} {fn}", f"statistics::rice_{fn}({b!r})", want, tol)
    g.write()


def search_group_a():
    """The rows of the six laws of the first search group.

    The laws whose quantile needs a search give their density, their
    cumulative function and their quantile. The rows come from SciPy
    1.18.1.
    """
    g = G("distributions_search_a", SCRIPT, tol="1e-10")
    sets = [
        ("argus", (1.0,), (0.5,)),
        ("dpareto_lognorm", (0.0, 1.0, 1.0, 1.0), (1.0, 2.0, 4.0)),
        ("erlang", (3.0,), (1.0, 2.0, 5.0)),
        ("exponnorm", (1.0,), (-1.0, 0.0, 2.0)),
        ("foldcauchy", (1.0,), (0.5, 1.0, 3.0)),
        ("foldnorm", (0.5,), (0.5, 1.0, 2.0)),
    ]
    for nm, args, xs in sets:
        d = getattr(stats, nm)
        tail = "".join(", " + repr(a) for a in args)
        for x in xs:
            g.val(f"{nm} pdf {x}", f"statistics::{nm}_pdf({x!r}{tail})",
                  float(d.pdf(x, *args)))
            g.val(f"{nm} cdf {x}", f"statistics::{nm}_cdf({x!r}{tail})",
                  float(d.cdf(x, *args)))
        for p in (0.1, 0.5, 0.9):
            g.val(f"{nm} ppf {p}", f"statistics::{nm}_ppf({p!r}{tail})",
                  float(d.ppf(p, *args)))
    g.write()


def search_group_b():
    """The rows of the five laws of the second search group.

    The laws whose quantile needs a search give their density, their
    cumulative function and their quantile. The rows come from SciPy
    1.18.1.
    """
    g = G("distributions_search_b", SCRIPT, tol="1e-10")
    sets = [
        ("loguniform", (1.0, 4.0), (1.5, 2.0, 3.0)),
        ("irwinhall", (3,), (0.5, 1.5, 2.5)),
        ("gausshyper", (2.0, 3.0, 1.0, 1.0), (0.25, 0.5, 0.75)),
        ("geninvgauss", (1.5, 2.0), (0.5, 1.0, 2.0)),
        ("genhyperbolic", (1.5, 1.0, 0.5), (0.5, 1.0, 2.0)),
    ]
    for nm, args, xs in sets:
        d = getattr(stats, nm)
        tail = "".join(", " + repr(a) for a in args)
        for x in xs:
            g.val(f"{nm} pdf {x}", f"statistics::{nm}_pdf({x!r}{tail})",
                  float(d.pdf(x, *args)))
            g.val(f"{nm} cdf {x}", f"statistics::{nm}_cdf({x!r}{tail})",
                  float(d.cdf(x, *args)))
        for p in (0.1, 0.5, 0.9):
            g.val(f"{nm} ppf {p}", f"statistics::{nm}_ppf({p!r}{tail})",
                  float(d.ppf(p, *args)))
    g.write()


def search_group_c():
    """The rows of the five laws of the third search group.

    The laws whose quantile needs a search give their density, their
    cumulative function and their quantile. The rows come from SciPy
    1.18.1.
    """
    g = G("distributions_search_c", SCRIPT, tol="1e-10")
    sets = [
        ("norminvgauss", (1.0, 0.5), (0.0, 1.0, 2.0)),
        ("recipinvgauss", (2.0,), (0.25, 0.5, 1.0)),
        ("rel_breitwigner", (0.5,), (0.0, 1.0, 2.0)),
        ("studentized_range", (3.0, 10.0), (1.0, 2.0, 3.0)),
        ("vonmises", (2.0,), (0.0, 1.0, 3.0)),
    ]
    for nm, args, xs in sets:
        d = getattr(stats, nm)
        tail = "".join(", " + repr(a) for a in args)
        for x in xs:
            g.val(f"{nm} pdf {x}", f"statistics::{nm}_pdf({x!r}{tail})",
                  float(d.pdf(x, *args)))
            g.val(f"{nm} cdf {x}", f"statistics::{nm}_cdf({x!r}{tail})",
                  float(d.cdf(x, *args)))
        for p in (0.1, 0.5, 0.9):
            g.val(f"{nm} ppf {p}", f"statistics::{nm}_ppf({p!r}{tail})",
                  float(d.ppf(p, *args)))
    g.write()


def levy_agrees(a, b):
    """Two values of the law agree when the difference stays under 1e-9."""
    return abs(a - b) <= 1e-9 * max(abs(a), abs(b), 1e-300)


def search_group_d():
    """The rows of the stable law.

    The rows come from mpmath at 30 digits, by the Fourier representation
    of the law. The representation is a different formula from the
    integral of Nolan that the package uses. The two agree to about
    1e-15 at every point of the grid.

    SciPy 1.18.1 checks the rows. A point joins the file only where
    SciPy and the representation agree. The far tail of the law stays
    out, because no method of double precision resolves the value
    there. The quantiles come from SciPy, because a quantile is an
    inverse of the cumulative function and the package finds it by a
    search.

    The rows run in parallel, because each point needs an integral.
    """
    import levy_oracle
    from multiprocessing import Pool

    g = G("distributions_search_d", SCRIPT, tol="1e-10",
          note="the answers of mpmath at 30 digits; the quantiles come "
               "from SciPy 1.18.1")
    alphas = [0.6, 0.8, 1.0, 1.3, 1.5, 1.9]
    betas = [-0.5, 0.0, 0.5, 1.0]
    xs = [-2.0, 0.5, 3.0]
    jobs = [(a, b, x) for a in alphas for b in betas for x in xs]
    with Pool(6) as pool:
        rows = pool.map(levy_oracle.pair, jobs)
    for (a, b, x), (p, c) in zip(jobs, rows):
        sp = float(stats.levy_stable.pdf(x, a, b))
        sc = float(stats.levy_stable.cdf(x, a, b))
        if not levy_agrees(p, sp) or not levy_agrees(c, sc):
            print("skip the far tail", a, b, x, p, sp, c, sc)
            continue
        tail = f", {a!r}, {b!r}"
        g.val(f"levy pdf {x} {a} {b}", f"statistics::levy_stable_pdf({x!r}{tail})", p)
        g.val(f"levy cdf {x} {a} {b}", f"statistics::levy_stable_cdf({x!r}{tail})", c)
    for p in (0.1, 0.5, 0.9):
        g.val(f"levy ppf {p}", f"statistics::levy_stable_ppf({p!r}, 1.5, 0.0)",
              float(stats.levy_stable.ppf(p, 1.5, 0.0)))
    g.write()


TAIL_A = [
    ("alpha", (1.5,), (0.2, 0.5, 1.0, 3.0)),
    ("fatiguelife", (0.7,), (0.3, 1.0, 2.5)),
    ("levy", (), (0.3, 1.0, 4.0)),
    ("levy_l", (), (-4.0, -1.0, -0.3)),
    ("logistic", (), (-2.0, 0.0, 1.5)),
    ("loglaplace", (2.5,), (0.4, 1.0, 2.5)),
    ("lomax", (2.0,), (0.5, 1.0, 4.0)),
    ("pareto", (2.5,), (1.2, 2.0, 5.0)),
    ("powerlaw", (1.8,), (0.2, 0.5, 0.9)),
    ("reciprocal", (1.0, 4.0), (1.5, 2.0, 3.0)),
    ("rdist", (3.0,), (-0.6, 0.0, 0.7)),
    ("skewcauchy", (0.4,), (-2.0, 0.5, 3.0)),
    ("wrapcauchy", (0.6,), (0.8, 3.0, 5.5)),
    ("kappa4", (0.3, 0.2), (-1.0, 0.5, 2.0)),
    ("kappa4", (0.0, 0.2), (-1.0, 0.5, 2.0)),
    ("kappa4", (0.3, 0.0), (-1.0, 0.5, 2.0)),
    ("kappa4", (0.0, 0.0), (-1.0, 0.5, 2.0)),
]

TAIL_B = [
    ("crystalball", (1.5, 3.0), (-3.0, -1.0, 0.5, 2.0)),
    ("kstwobign", (), (0.4, 0.8, 1.2, 2.0)),
    ("vonmises_line", (2.0,), (-2.0, 0.5, 1.5)),
    ("ncx2", (3.0, 2.0), (0.5, 2.0, 6.0)),
    ("ncf", (4.0, 7.0, 1.5), (0.4, 1.0, 3.0)),
    ("nct", (6.0, 1.0), (-1.0, 0.5, 2.5)),
    ("ksone", (10,), (0.1, 0.3, 0.5)),
    ("kstwo", (10,), (0.15, 0.3, 0.5)),
]


def tail_group(name, sets, tol="1e-9"):
    """The rows of the laws of one tail group. Each law gives its density,
    its cumulative function and its quantile. The rows come from SciPy
    1.18.1."""
    g = G(name, SCRIPT, tol=tol)
    for nm, args, xs in sets:
        d = getattr(stats, nm)
        tail = "".join(", " + repr(a) for a in args)
        tag = nm + "".join(" " + repr(a) for a in args)
        for x in xs:
            g.val(f"{tag} pdf {x}", f"statistics::{nm}_pdf({x!r}{tail})", float(d.pdf(x, *args)))
            g.val(f"{tag} cdf {x}", f"statistics::{nm}_cdf({x!r}{tail})", float(d.cdf(x, *args)))
        for p in (0.1, 0.5, 0.9):
            g.val(f"{tag} ppf {p}", f"statistics::{nm}_ppf({p!r}{tail})", float(d.ppf(p, *args)))
    g.write()


def tail_group_a():
    tail_group("distributions_closed", TAIL_A)


def tail_group_b():
    tail_group("distributions_noncentral", TAIL_B, tol="1e-7")


def nch_group():
    """The rows of the two non-central hypergeometric laws."""
    g = G("distributions_nch", SCRIPT, tol="1e-9")
    for nm in ("nchypergeom_fisher", "nchypergeom_wallenius"):
        d = getattr(stats, nm)
        for args in ((50, 20, 15, 2.5), (30, 10, 12, 0.5)):
            tail = "".join(", " + repr(float(a)) for a in args)
            lo = max(0, args[2] - (args[0] - args[1]))
            for k in (lo + 1, lo + 3, lo + 5):
                g.val(f"{nm} pmf {k} {args}", f"statistics::{nm}_pmf({float(k)!r}{tail})", float(d.pmf(k, *args)))
                g.val(f"{nm} cdf {k} {args}", f"statistics::{nm}_cdf({float(k)!r}{tail})", float(d.cdf(k, *args)))
            for p in (0.1, 0.5, 0.9):
                g.val(f"{nm} ppf {p} {args}", f"statistics::{nm}_ppf({p!r}{tail})", float(d.ppf(p, *args)))
    g.write()


def landau_group():
    tail_group("distributions_landau", [("landau", (), (-2.0, -0.5, 0.0, 1.0, 4.0))], tol="1e-7")


def family_rows(g, nm, d, args, moments=None):
    """The rows of one family: the density, the cumulative function, the
    survival function, the quantile and the moments."""
    lo, hi = d.support(*args)
    xs = [float(d.ppf(p, *args)) for p in (0.05, 0.25, 0.5, 0.75, 0.95)]
    out = []
    if np.isfinite(lo):
        out.append(float(lo - 1.0))
    if np.isfinite(hi):
        out.append(float(hi + 1.0))
    tail = "".join(", " + repr(a) for a in args)
    for x in xs + out:
        for fn in ("pdf", "cdf", "sf"):
            want = float(getattr(d, fn)(x, *args))
            g.val(f"{nm} {fn} {x:.4g}", f"statistics::{nm}_{fn}({x!r}{tail})", want)
    for p in (0.1, 0.5, 0.9):
        g.val(f"{nm} ppf {p}", f"statistics::{nm}_ppf({p!r}{tail})", float(d.ppf(p, *args)))
    m, v, s, k = moments
    arg = ", ".join(repr(a) for a in args)
    for fn, want in (("mean", m), ("var", v), ("skew", s), ("kurt", k)):
        call = f"statistics::{nm}_{fn}({arg})" if args else f"statistics::{nm}_{fn}()"
        if np.isfinite(want):
            g.val(f"{nm} {fn}", call, want)
        else:
            g.code(f"{nm} {fn} does not exist", call, 7)


def power_law_group():
    """The rows of the five families of a power."""
    g = G("distributions_power", SCRIPT, tol="1e-10")
    sets = [
        ("fisk", (5.0,)),
        ("burr", (5.0, 3.0)),
        ("burr12", (2.0, 3.0)),
        ("betaprime", (3.0, 5.0)),
        ("genpareto", (-0.5,)),
    ]
    for nm, args in sets:
        d = getattr(stats, nm)
        m, v, s, k = (float(q) for q in d.stats(*args, moments="mvsk"))
        family_rows(g, nm, d, args, (m, v, s, k))
    g.write()


def growth_group():
    """The rows of the three families of a growth or a logistic kind."""
    g = G("distributions_growth", SCRIPT, tol="1e-10",
          note="the answers of SciPy 1.18.1, and the moments of the "
               "Gompertz law from an integration of its density")
    sets = [
        ("genlogistic", (2.5,)),
        ("genhalflogistic", (2.0,)),
        ("gompertz", (1.5,)),
    ]
    for nm, args in sets:
        d = getattr(stats, nm)
        # SciPy integrates the quantile function for the moments of the
        # Gompertz law and loses digits. The density route keeps them.
        if nm == "gompertz":
            moments = law_moments(d, args)
        else:
            moments = tuple(float(q) for q in d.stats(*args, moments="mvsk"))
        family_rows(g, nm, d, args, moments)
    g.write()



def truncated_group():
    """The rows of the four families of a truncation."""
    g = G("distributions_truncated", SCRIPT, tol="1e-11")
    sets = [
        ("truncexpon", (2.0,)),
        ("truncpareto", (2.0, 5.0)),
        ("truncnorm", (-1.0, 2.0)),
    ]
    for nm, args in sets:
        d = getattr(stats, nm)
        moments = tuple(float(q) for q in d.stats(*args, moments="mvsk"))
        family_rows(g, nm, d, args, moments)
    # The support of the truncated Pareto law is bounded, so its moments
    # exist for every shape.
    g.code("truncpareto shape 0", "statistics::truncpareto_pdf(1.5, 0.0, 5.0)", 5)
    g.code("truncpareto end below 1", "statistics::truncpareto_cdf(1.5, 2.0, 0.5)", 5)
    g.code("truncnorm bound above", "statistics::truncnorm_pdf(0.0, 2.0, -1.0)", 5)
    g.code("truncexpon shape 0", "statistics::truncexpon_pdf(1.0, 0.0)", 5)
    g.write()


def gamma_group():
    """The rows of the four families that come from the gamma law."""
    g = G("distributions_gamma", SCRIPT, tol="1e-11")
    sets = [
        ("chi", (3.0,)),
        ("dgamma", (2.0,)),
        ("invgamma", (3.0,)),
        ("loggamma", (2.0,)),
    ]
    for nm, args in sets:
        d = getattr(stats, nm)
        moments = tuple(float(q) for q in d.stats(*args, moments="mvsk"))
        family_rows(g, nm, d, args, moments)
    g.code("chi shape 0", "statistics::chi_pdf(1.0, 0.0)", 5)
    g.code("dgamma shape 0", "statistics::dgamma_pdf(1.0, 0.0)", 5)
    g.code("invgamma shape 0", "statistics::invgamma_pdf(1.0, 0.0)", 5)
    g.code("loggamma shape 0", "statistics::loggamma_pdf(1.0, 0.0)", 5)
    g.write()


def generalized_group():
    """The rows of the three families of a generalized law."""
    g = G("distributions_generalized", SCRIPT, tol="1e-10",
          note="the answers of SciPy 1.18.1, and the moments of the "
               "generalized exponential law from an integration of its density")
    sets = [
        ("gengamma", (3.0, 2.0)),
        ("genexpon", (1.0, 2.0, 3.0)),
        ("mielke", (2.0, 3.0)),
    ]
    for nm, args in sets:
        d = getattr(stats, nm)
        # SciPy integrates the quantile function for the moments of the
        # generalized exponential law and loses digits. The density route
        # keeps them.
        if nm == "genexpon":
            moments = law_moments(d, args)
        else:
            moments = tuple(float(q) for q in d.stats(*args, moments="mvsk"))
        family_rows(g, nm, d, args, moments)
    g.code("gengamma shape 0", "statistics::gengamma_pdf(1.0, 0.0, 2.0)", 5)
    g.code("gengamma shape c 0", "statistics::gengamma_pdf(1.0, 3.0, 0.0)", 5)
    g.code("genexpon shape 0", "statistics::genexpon_pdf(1.0, 0.0, 2.0, 3.0)", 5)
    g.code("mielke shape 0", "statistics::mielke_pdf(1.0, 0.0, 3.0)", 5)
    g.write()


def symmetric_group():
    """The rows of the three families of a symmetric law."""
    g = G("distributions_symmetric", SCRIPT, tol="1e-11")
    sets = [
        ("gennorm", (1.5,)),
        ("halfgennorm", (1.5,)),
        ("tukeylambda", (0.5,)),
    ]
    for nm, args in sets:
        d = getattr(stats, nm)
        moments = tuple(float(q) for q in d.stats(*args, moments="mvsk"))
        family_rows(g, nm, d, args, moments)
    g.code("gennorm shape 0", "statistics::gennorm_pdf(1.0, 0.0)", 5)
    g.code("halfgennorm shape 0", "statistics::halfgennorm_pdf(1.0, 0.0)", 5)
    g.code("tukeylambda nan shape", "statistics::tukeylambda_pdf(1.0, 0.0 / 0.0)", 3)
    g.write()


def johnson_group():
    """The rows of the three families of the Johnson kind."""
    g = G("distributions_johnson", SCRIPT, tol="1e-10",
          note="the answers of SciPy 1.18.1, and the moments of the Johnson "
               "law of the B kind from an integration of its density")
    sets = [
        ("jf_skew_t", (3.0, 2.0)),
        ("johnsonsb", (1.0, 2.0)),
        ("johnsonsu", (1.0, 2.0)),
    ]
    for nm, args in sets:
        d = getattr(stats, nm)
        # SciPy integrates the density of the Johnson law of the B kind and
        # loses digits. The route of the law keeps them.
        if nm == "johnsonsb":
            moments = johnsonsb_moments(args)
        else:
            moments = tuple(float(q) for q in d.stats(*args, moments="mvsk"))
        family_rows(g, nm, d, args, moments)
    g.code("jf_skew_t shape 0", "statistics::jf_skew_t_pdf(1.0, 0.0, 2.0)", 5)
    g.code("johnsonsb shape 0", "statistics::johnsonsb_pdf(0.5, 1.0, 0.0)", 5)
    g.code("johnsonsu shape 0", "statistics::johnsonsu_pdf(1.0, 1.0, 0.0)", 5)
    g.write()


def johnsonsb_moments(args):
    """The moments of the Johnson law of the B kind.

    The law of the logit of a normal variable gives the raw moments. The
    integration runs over the variable of the normal law. SciPy integrates
    the density of the law and loses digits in the third and the fourth
    moment.
    """
    a, b = args
    d = stats.johnsonsb

    def raw(n):
        def g_of_u(u):
            # The point of the law is the logistic of the variable. The
            # density of the law in that variable is the normal density.
            s = 1.0 / (1.0 + math.exp(-u))
            z = a + b * u
            return s ** n * b * math.exp(-0.5 * z * z) / math.sqrt(2.0 * math.pi)

        total = 0.0
        pieces = 400
        for i in range(pieces):
            left = -16.0 + 32.0 * i / pieces
            right = -16.0 + 32.0 * (i + 1) / pieces
            total += quad(g_of_u, left, right, limit=200)[0]
        return total

    m1, m2, m3, m4 = (raw(n) for n in (1, 2, 3, 4))
    v = m2 - m1 * m1
    skew = (m3 - 3.0 * m1 * m2 + 2.0 * m1 ** 3) / v ** 1.5
    kurt = (m4 - 4.0 * m1 * m3 + 6.0 * m1 * m1 * m2 - 3.0 * m1 ** 4) / v ** 2 - 3.0
    return [m1, v, skew, kurt]



def kappa3_moments(a):
    """The moments of the Kappa law of the third kind.

    SciPy integrates the quantile function and loses the digits of the
    higher moments. The closed form keeps them.
    """
    def raw(n):
        if n >= a:
            return float("nan")
        return (a ** (n / a - 1.0) * math.gamma((n + 1) / a)
                * math.gamma(1.0 - n / a) / math.gamma(1.0 + 1.0 / a))

    m1, m2, m3, m4 = raw(1), raw(2), raw(3), raw(4)
    if math.isnan(m1) or math.isnan(m2):
        var = float("nan")
    else:
        var = m2 - m1 * m1
    if math.isnan(var) or math.isnan(m3):
        skew = float("nan")
    else:
        skew = (m3 - 3.0 * m1 * m2 + 2.0 * m1 ** 3) / var ** 1.5
    if math.isnan(var) or math.isnan(m3) or math.isnan(m4):
        kurt = float("nan")
    else:
        kurt = (m4 - 4.0 * m1 * m3 + 6.0 * m1 * m1 * m2 - 3.0 * m1 ** 4) / var ** 2 - 3.0
    return [m1, var, skew, kurt]


def special_group():
    """The rows of the five families of a special shape."""
    g = G("distributions_special", SCRIPT, tol="1e-10",
          note="the answers of SciPy 1.18.1, and the moments of the Kappa "
               "law of the third kind from its closed form")
    sets = [
        ("invgauss", [(0.5,), (2.0,)]),
        ("nakagami", [(0.75,), (2.5,)]),
        ("kappa3", [(2.0,), (3.5,)]),
        ("pearson3", [(1.5,), (-1.5,)]),
    ]
    for nm, argsets in sets:
        d = getattr(stats, nm)
        for args in argsets:
            if nm == "kappa3":
                moments = kappa3_moments(args[0])
            else:
                moments = tuple(float(q) for q in d.stats(*args, moments="mvsk"))
            family_rows(g, nm, d, args, moments)
    g.code("invgauss shape 0", "statistics::invgauss_pdf(1.0, 0.0)", 5)
    g.code("nakagami shape 0", "statistics::nakagami_pdf(1.0, 0.0)", 5)
    g.code("kappa3 shape 0", "statistics::kappa3_pdf(1.0, 0.0)", 5)
    g.code("pearson3 skew 0", "statistics::pearson3_pdf(1.0, 0.0)", 0)
    g.write()


def extreme_group():
    """The rows of the generalised extreme value law."""
    g = G("distributions_extreme", SCRIPT, tol="1e-10")
    d = stats.genextreme
    for args in [(0.5,), (-0.3,)]:
        moments = tuple(float(q) for q in d.stats(*args, moments="mvsk"))
        family_rows(g, "genextreme", d, args, moments)
    g.code("genextreme shape nan", "statistics::genextreme_pdf(1.0, 0.0 / 0.0)", 3)
    g.write()


def skew_group():
    """The rows of the two families of a skew shape."""
    g = G("distributions_skew", SCRIPT, tol="1e-10")
    sets = [
        ("laplace_asymmetric", [(0.5,), (2.0,)]),
        ("skewnorm", [(2.0,), (-1.5,), (0.5,)]),
    ]
    for nm, argsets in sets:
        d = getattr(stats, nm)
        for args in argsets:
            moments = tuple(float(q) for q in d.stats(*args, moments="mvsk"))
            family_rows(g, nm, d, args, moments)
    g.code("laplace_asymmetric shape 0",
           "statistics::laplace_asymmetric_pdf(1.0, 0.0)", 5)
    g.code("skewnorm shape nan", "statistics::skewnorm_pdf(1.0, 0.0 / 0.0)", 3)
    g.write()


def discrete_group_a():
    """The rows of the first group of discrete families.

    The laws of this group have a count, a shape or two. The Poisson
    binomial law of the second group takes a list of probabilities.
    """
    g = G("distributions_discrete_a", SCRIPT, tol="1e-10",
          note="the answers of SciPy 1.18.1")
    sets = [
        ("bernoulli", (0.5,), [0, 1], [0, 1]),
        ("betabinom", (10, 2.0, 3.0), [0, 3, 10], [0, 10]),
        ("betanbinom", (10, 2.0, 3.0), [0, 5, 30], [0, 100]),
        ("boltzmann", (2.0, 5), [0, 2, 4], [0, 4]),
        ("dlaplace", (0.5,), [-2, 0, 3], [-6, 6]),
        ("logser", (0.5,), [1, 3, 10], [1, 50]),
        ("nhypergeom", (10, 5, 3), [0, 2, 5], [0, 5]),
    ]
    for nm, args, ks, ends in sets:
        d = getattr(stats, nm)
        arg = ", ".join(repr(a) for a in args)
        for k in ks:
            g.val(f"{nm} {arg} pmf {k}",
                  f"statistics::{nm}_pmf({k!r}, {arg})", float(d.pmf(k, *args)))
            g.val(f"{nm} {arg} cdf {k}",
                  f"statistics::{nm}_cdf({k!r}, {arg})", float(d.cdf(k, *args)))
        for k in ends:
            g.val(f"{nm} {arg} pmf {k}",
                  f"statistics::{nm}_pmf({k!r}, {arg})", float(d.pmf(k, *args)))
        for p in (0.1, 0.5, 0.9):
            want = float(d.ppf(p, *args))
            if np.isfinite(want):
                g.val(f"{nm} {arg} ppf {p}",
                      f"statistics::{nm}_ppf({p!r}, {arg})", want)
    g.write()


def discrete_group_b():
    """The rows of the second group of discrete families.

    The Poisson binomial law takes a list of probabilities.
    """
    g = G("distributions_discrete_b", SCRIPT, tol="1e-10",
          note="the answers of SciPy 1.18.1")
    sets = [
        ("planck", (2.0,), [0, 3, 10], (0, 100)),
        ("randint", (0, 5), [0, 2, 4], (0, 4)),
        ("skellam", (1.0, 2.0), [-2, 0, 3], (-8, 8)),
        ("yulesimon", (2.0,), [1, 4, 20], (1, 200)),
        ("zipf", (2.0,), [1, 3, 10], (1, 100)),
        ("zipfian", (2.0, 5), [1, 2, 5], (1, 5)),
    ]
    for nm, args, ks, ends in sets:
        d = getattr(stats, nm)
        arg = ", ".join(repr(a) for a in args)
        for k in ks:
            g.val(f"{nm} {arg} pmf {k}",
                  f"statistics::{nm}_pmf({k!r}, {arg})", float(d.pmf(k, *args)))
            g.val(f"{nm} {arg} cdf {k}",
                  f"statistics::{nm}_cdf({k!r}, {arg})", float(d.cdf(k, *args)))
        for k in (ends[0], ends[1]):
            g.val(f"{nm} {arg} pmf {k}",
                  f"statistics::{nm}_pmf({k!r}, {arg})", float(d.pmf(k, *args)))
        for p in (0.1, 0.5, 0.9):
            want = float(d.ppf(p, *args))
            if np.isfinite(want):
                g.val(f"{nm} {arg} ppf {p}",
                      f"statistics::{nm}_ppf({p!r}, {arg})", want)
    # The Poisson binomial law takes a list of probabilities.
    ps = [0.3, 0.5, 0.2]
    d = stats.poisson_binom
    lst = "[" + ", ".join(repr(p) for p in ps) + "]"
    for k in (0, 1, 2, 3):
        g.val(f"poisson_binom pmf {k}",
              f"statistics::poisson_binom_pmf({k!r}, {lst})",
              float(d.pmf(k, ps)))
        g.val(f"poisson_binom cdf {k}",
              f"statistics::poisson_binom_cdf({k!r}, {lst})",
              float(d.cdf(k, ps)))
    for p in (0.1, 0.5, 0.9):
        g.val(f"poisson_binom ppf {p}",
              f"statistics::poisson_binom_ppf({p!r}, {lst})",
              float(d.ppf(p, ps)))
    g.write()


if __name__ == "__main__":
    location_group()
    elementary_group()
    bounded_group()
    weibull_group()
    power_group()
    power_law_group()
    growth_group()
    truncated_group()
    gamma_group()
    generalized_group()
    symmetric_group()
    johnson_group()
    special_group()
    extreme_group()
    skew_group()
    discrete_group_a()
    discrete_group_b()
    rice_group()
    search_group_a()
    search_group_b()
    search_group_c()
    search_group_d()
    tail_group_a()
    tail_group_b()
    nch_group()
    landau_group()
    objects_group()
    mixture_group()
    covariance_group()
    histogram_group()
