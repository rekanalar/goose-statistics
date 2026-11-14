"""The rows of the distribution objects, the mixture, the covariance
factories and the histogram law.

The rows come from SciPy 1.18.1. Needs Python 3 with SciPy. The script
lives in the package, under test/, and the user of the package never
runs it.
"""
import numpy as np
import scipy.stats as st

from vecgen import G

SCRIPT = "make_object_vectors.py"


def objects_group():
    """The rows of the distribution objects: Normal, Uniform, Logistic and
    Binomial, with the methods of SciPy 1.18.1."""
    g = G("distributions_objects", SCRIPT, tol="1e-10")
    sets = [
        ("normal", "statistics::Dist { kind: 1, a: 1.0, b: 2.0 }", st.Normal(mu=1.0, sigma=2.0)),
        ("normal", "statistics::Dist { kind: 1, a: -0.5, b: 0.75 }", st.Normal(mu=-0.5, sigma=0.75)),
        ("uniform", "statistics::Dist { kind: 2, a: 0.0, b: 2.0 }", st.Uniform(a=0.0, b=2.0)),
        ("uniform", "statistics::Dist { kind: 2, a: -1.5, b: 0.5 }", st.Uniform(a=-1.5, b=0.5)),
        ("logistic", "statistics::Dist { kind: 3, a: 0.0, b: 1.0 }", st.Logistic()),
        ("binomial", "statistics::Dist { kind: 4, a: 10.0, b: 0.5 }", st.Binomial(n=10, p=0.5)),
        ("binomial", "statistics::Dist { kind: 4, a: 20.0, b: 0.25 }", st.Binomial(n=20, p=0.25)),
    ]
    xs = [-1.0, 0.0, 0.5, 1.0, 2.0]
    for nm, ctor, d in sets:
        for x in xs:
            g.val(f"{nm} pdf {x}", f"statistics::dist_pdf({ctor}, {x!r})",
                  float(d.pdf(x)))
            g.val(f"{nm} cdf {x}", f"statistics::dist_cdf({ctor}, {x!r})",
                  float(d.cdf(x)))
            g.val(f"{nm} ccdf {x}", f"statistics::dist_ccdf({ctor}, {x!r})",
                  float(d.ccdf(x)))
        for p in (0.1, 0.5, 0.9):
            g.val(f"{nm} icdf {p}", f"statistics::dist_icdf({ctor}, {p!r})",
                  float(d.icdf(p)))
            g.val(f"{nm} iccdf {p}", f"statistics::dist_iccdf({ctor}, {p!r})",
                  float(d.iccdf(p)))
        g.val(f"{nm} logpdf 1", f"statistics::dist_logpdf({ctor}, 1.0)",
              float(d.logpdf(1.0)))
        g.val(f"{nm} logcdf 1", f"statistics::dist_logcdf({ctor}, 1.0)",
              float(d.logcdf(1.0)))
        g.val(f"{nm} logccdf 1", f"statistics::dist_logccdf({ctor}, 1.0)",
              float(d.logccdf(1.0)))
        g.val(f"{nm} ilogcdf -1", f"statistics::dist_ilogcdf({ctor}, -1.0)",
              float(d.ilogcdf(-1.0)))
        g.val(f"{nm} ilogccdf -1", f"statistics::dist_ilogccdf({ctor}, -1.0)",
              float(d.ilogccdf(-1.0)))
        g.val(f"{nm} mean", f"statistics::dist_mean({ctor})", float(d.mean()))
        g.val(f"{nm} median", f"statistics::dist_median({ctor})", float(d.median()))
        g.val(f"{nm} mode", f"statistics::dist_mode({ctor})", float(d.mode()))
        g.val(f"{nm} variance", f"statistics::dist_variance({ctor})",
              float(d.variance()))
        g.val(f"{nm} standard deviation",
              f"statistics::dist_standard_deviation({ctor})",
              float(d.standard_deviation()))
        g.val(f"{nm} skewness", f"statistics::dist_skewness({ctor})",
              float(d.skewness()))
        g.val(f"{nm} kurtosis", f"statistics::dist_kurtosis({ctor})",
              float(d.kurtosis()))
        g.val(f"{nm} entropy", f"statistics::dist_entropy({ctor})",
              float(d.entropy()))
        for o in (1, 2, 3, 4):
            g.val(f"{nm} moment {o}", f"statistics::dist_moment({ctor}, {o})",
                  float(d.moment(o)))
    g.write()


def mixture_group():
    """The rows of the mixture of two normal laws."""
    g = G("distributions_mixture", SCRIPT, tol="1e-10")
    comps = [st.Normal(mu=0.0, sigma=1.0), st.Normal(mu=3.0, sigma=1.0)]
    w = [0.3, 0.7]
    m = st.Mixture(comps, weights=w)
    for x in [-2.0, 0.0, 1.0, 3.0, 5.0]:
        g.val(f"mixture pdf {x}", f"statistics::mixture_pdf(comps, w, {x!r})",
              float(m.pdf(x)))
        g.val(f"mixture cdf {x}", f"statistics::mixture_cdf(comps, w, {x!r})",
              float(m.cdf(x)))
        g.val(f"mixture ccdf {x}", f"statistics::mixture_ccdf(comps, w, {x!r})",
              float(m.ccdf(x)))
    for p in (0.1, 0.5, 0.9):
        g.val(f"mixture icdf {p}", f"statistics::mixture_icdf(comps, w, {p!r})",
              float(m.icdf(p)))
    g.val("mixture mean", "statistics::mixture_mean(comps, w)", float(m.mean()))
    g.val("mixture variance", "statistics::mixture_variance(comps, w)",
          float(m.variance()))
    g.val("mixture standard deviation",
          "statistics::mixture_standard_deviation(comps, w)",
          float(m.standard_deviation()))
    g.val("mixture skewness", "statistics::mixture_skewness(comps, w)",
          float(m.skewness()))
    g.val("mixture kurtosis", "statistics::mixture_kurtosis(comps, w)",
          float(m.kurtosis()))
    for o in (1, 2, 3, 4):
        g.val(f"mixture moment {o}", f"statistics::mixture_moment(comps, w, {o})",
              float(m.moment(o)))
    g.write("    var comps: statistics::Dist[2] = [statistics::Dist { kind: 1, a: 0.0, b: 1.0 }, statistics::Dist { kind: 1, a: 3.0, b: 1.0 }];\n    let w: f64[2] = [0.3, 0.7];\n")


def covariance_group():
    """The rows of the four factories of a covariance matrix."""
    g = G("distributions_covariance", SCRIPT, tol="1e-9")
    d = np.diag([1.0, 2.0, 3.0])
    prec = np.linalg.inv(d)
    evals = np.array([1.0, 2.0, 3.0])
    evecs = np.eye(3)
    chol = np.linalg.cholesky(d)
    x = np.array([1.0, -2.0, 0.5])
    cases = (
        ("diagonal", 1, [1.0, 2.0, 3.0], st.Covariance.from_diagonal(np.array([1.0, 2.0, 3.0]))),
        ("precision", 2, [float(v) for v in prec.reshape(-1)],
         st.Covariance.from_precision(prec)),
        ("cholesky", 3, [float(v) for v in chol.reshape(-1)],
         st.Covariance.from_cholesky(chol)),
        ("eigen", 4, [float(v) for v in evals] + [float(v) for v in evecs.reshape(-1)],
         st.Covariance.from_eigendecomposition((evals, evecs))),
    )
    decl = ""
    for nm, kind, data, c in cases:
        var = f"data_{nm}"
        n = len(data)
        decl += f"    let {var}: f64[{n}] = [" + ", ".join(repr(float(v)) for v in data) + "];\n"
        g.fill(f"covariance {nm} matrix",
               f"statistics::covariance_matrix(out, {kind}, {var}, 3)",
               [float(v) for v in np.asarray(c.covariance).reshape(-1)])
        g.val(f"covariance {nm} log pdet",
              f"statistics::covariance_log_pdet({kind}, {var}, 3)", float(c.log_pdet))
        g.val(f"covariance {nm} rank",
              f"statistics::covariance_rank({kind}, {var}, 3)", float(c.rank))
        g.fill(f"covariance {nm} whiten",
               f"statistics::covariance_whiten(out, {kind}, {var}, 3, x)",
               [float(v) for v in np.asarray(c.whiten(x)).reshape(-1)])
        g.fill(f"covariance {nm} colorize",
               f"statistics::covariance_colorize(out, {kind}, {var}, 3, x)",
               [float(v) for v in np.asarray(c.colorize(x)).reshape(-1)])
    g.write(decl + "    let x: f64[3] = [1.0, -2.0, 0.5];\n")


def histogram_group():
    """The rows of the histogram law."""
    g = G("distributions_histogram", SCRIPT, tol="1e-10")
    edges = [0.0, 1.0, 2.0, 3.0]
    content = [1.0, 2.0, 1.0]
    d = st.rv_histogram((np.asarray(content), np.asarray(edges)), density=True)
    for x in [-1.0, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0]:
        g.val(f"histogram pdf {x}",
              f"statistics::histogram_pdf(content, edges, true, {x!r})",
              float(d.pdf(x)))
        g.val(f"histogram cdf {x}",
              f"statistics::histogram_cdf(content, edges, true, {x!r})",
              float(d.cdf(x)))
    for p in (0.1, 0.25, 0.5, 0.75, 0.9):
        g.val(f"histogram ppf {p}",
              f"statistics::histogram_ppf(content, edges, true, {p!r})",
              float(d.ppf(p)))
    g.val("histogram mean", "statistics::histogram_mean(content, edges, true)",
          float(d.mean()))
    g.val("histogram variance", "statistics::histogram_var(content, edges, true)",
          float(d.var()))
    g.val("histogram entropy", "statistics::histogram_entropy(content, edges, true)",
          float(d.entropy()))
    g.val("histogram moment 2", "statistics::histogram_moment(content, edges, true, 2)",
          float(d.moment(2)))
    c2 = st.rv_histogram((np.asarray([2.0, 4.0, 2.0]), np.asarray(edges)), density=False)
    g.val("histogram count pdf 1.5",
          "statistics::histogram_pdf(ccontent, edges, false, 1.5)", float(c2.pdf(1.5)))
    g.val("histogram count cdf 1.5",
          "statistics::histogram_cdf(ccontent, edges, false, 1.5)", float(c2.cdf(1.5)))
    g.write("    let content: f64[3] = [1.0, 2.0, 1.0];\n    let ccontent: f64[3] = [2.0, 4.0, 2.0];\n    let edges: f64[4] = [0.0, 1.0, 2.0, 3.0];\n")


if __name__ == "__main__":
    objects_group()
    mixture_group()
    covariance_group()
    histogram_group()
