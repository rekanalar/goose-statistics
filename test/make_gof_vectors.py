"""The rows of the goodness-of-fit test over the chi-square, student
and F families, and of the Filliben statistic over every family.

The rows come from SciPy 1.18.1. Needs Python 3 with SciPy. The script
lives in the package, under test/, and the user of the package never
runs it.
"""
import numpy as np
import scipy.stats as st

from vecgen import G, arr

SCRIPT = "make_gof_vectors.py"


def gof_group():
    g = G("fitting_shapes", SCRIPT, tol="5e-4",
          note="the answers of SciPy 1.18.1; the fits of the chi-square and "
               "the F family fix the location at 0, so the values agree to about 1e-4")
    rng = np.random.default_rng(20261114)
    seeded = lambda: np.random.default_rng(20261114)
    new = [
        ("chi2", st.chi2, {'loc': 0}, 7, st.chi2.rvs(df=6.0, scale=0.7, size=60, random_state=rng)),
        ("t", st.t, {}, 8, st.t.rvs(df=6.0, loc=0.5, scale=1.2, size=60, random_state=rng)),
        ("f", st.f, {'loc': 0}, 9, st.f.rvs(dfn=5.0, dfd=9.0, scale=1.1, size=60, random_state=rng)),
    ]
    old = [
        ("norm", 0, st.norm, {}, st.norm.rvs(loc=1.0, scale=2.0, size=40, random_state=seeded())),
        ("uniform", 1, st.uniform, {}, st.uniform.rvs(loc=-1.0, scale=3.0, size=40, random_state=seeded())),
        ("expon", 2, st.expon, {'loc': 0}, st.expon.rvs(scale=1.5, size=40, random_state=seeded())),
        ("lognorm", 3, st.lognorm, {'loc': 0}, st.lognorm.rvs(s=0.6, size=40, random_state=seeded())),
        ("weibull", 4, st.weibull_min, {'loc': 0}, st.weibull_min.rvs(c=1.8, size=40, random_state=seeded())),
        ("gamma", 5, st.gamma, {'loc': 0}, st.gamma.rvs(a=2.0, size=40, random_state=seeded())),
        ("beta", 6, st.beta, {'loc': 0, 'scale': 1}, st.beta.rvs(a=2.0, b=3.0, size=40, random_state=seeded())),
    ]
    decl = ""
    for name, _d, _kp, _fam, x in new:
        decl += f"    let xs_{name}: f64[{len(x)}] = {arr(x)};\n"
    for name, _fam, _d, _kp, x in old:
        decl += f"    let xs_{name}: f64[{len(x)}] = {arr(x)};\n"
    decl += "    var rg = rng { s: 11 };\n"
    for name, d, kp, fam, x in new:
        for sn, si in (("ad", 0), ("ks", 1), ("cvm", 2), ("filliben", 3)):
            r = st.goodness_of_fit(d, x, known_params=kp or None, statistic=sn,
                                   n_mc_samples=19, rng=np.random.default_rng(99))
            g.field(f"{name} {sn}",
                    f"statistics::goodness_of_fit(out, xs_{name}[..], {fam}, {si}, 1, rg)",
                    {"statistic": float(r.statistic)})
    for name, fam, d, kp, x in old:
        r = st.goodness_of_fit(d, x, known_params=kp or None, statistic="filliben",
                               n_mc_samples=19, rng=np.random.default_rng(99))
        g.field(f"{name} filliben",
                f"statistics::goodness_of_fit(out, xs_{name}[..], {fam}, 3, 1, rg)",
                {"statistic": float(r.statistic)})
    g.write(decl)


if __name__ == "__main__":
    gof_group()
