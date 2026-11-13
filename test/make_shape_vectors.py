#!/usr/bin/env python3
"""Writes test/suite/{descriptive,frequency,entropy,trimmed,circular,binned}_vectors.goose:
the answers of SciPy 1.18.1 for the descriptive and shape calls.

Each file holds one function run_<group>_vectors() with one check per
measured row. Run it from the package directory with the system python3.
"""
import math
import warnings
import numpy as np
from scipy import stats

warnings.simplefilter("ignore")

S1 = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
S2 = [2.5, -1.0, 4.0, 3.25, 10.5, 0.5, 7.75, 3.0, 8.0, 6.0]
TIES = [1.0, 2.5, 2.5, 4.0]
POS = [0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0]
ANG = [0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0]
DE = [0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3]
XB = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
VB = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0]


def f(v):
    v = float(v)
    if math.isnan(v):
        return "(0.0 / 0.0)"
    if math.isinf(v):
        return "(1.0 / 0.0)" if v > 0 else "(-1.0 / 0.0)"
    return repr(v)


def arr(v):
    return "[" + ", ".join(f(x) for x in v) + "]"


class G:
    """One output file."""
    def __init__(self, name, tol="1e-12"):
        self.name = name
        self.tol = tol
        self.lines = []
        self.n = 0

    def decl(self, nm, v):
        self.lines.append(f"    let {nm}: f64[{len(v)}] = {arr(v)};")

    def val(self, label, call, want, tol=None):
        self.n += 1
        k = self.n
        self.lines.append(f"    let v{k}, c{k} = {call};")
        self.lines.append(f'    chk_val(bad, "{label}", v{k}, c{k}, {f(want)}, {tol or self.tol});')

    def lst(self, label, call, want, tol=None, name="out", dummy=True):
        """call fills the list `name` and returns a value then a code."""
        self.n += 1
        k = self.n
        self.lines.append(f"    let _r{k}, c{k} = {call};")
        self.lines.append(f"    let w{k}: f64[{len(want)}] = {arr(want)};")
        self.lines.append(f'    chk_list(bad, "{label}", {name}[..], c{k}, w{k}[..], {tol or self.tol});')

    def code(self, label, call, want):
        self.n += 1
        k = self.n
        self.lines.append(f"    let _v{k}, c{k} = {call};")
        self.lines.append(f'    chk_code(bad, "{label}", c{k}, {want});')

    def raw(self, text):
        self.lines.append("    " + text)

    def write(self):
        head = (f"// The measured rows of the {self.name} calls: the answers of SciPy 1.18.1.\n"
                f"// The generator test/make_shape_vectors.py writes this file.\n"
                "namespace suite;\n\nimport std;\nimport math;\nimport statistics from \"../..\";\n\n"
                f"fn run_{self.name}_vectors() -> i64 {{\n    var bad: i64 = 0;\n"
                "    var out: f64[>..] = [];\n    var out2: f64[>..] = [];\n"
                + ("    var bst: f64[>..] = [];\n    var bshape: f64[>..] = [];\n    var bnum: f64[>..] = [];\n    var bempty: f64[>..] = [];\n" if self.name == "binned" else ""))
        body = "\n".join(self.lines)
        open(f"test/suite/{self.name}_vectors.goose", "w").write(head + body + "\n    return bad;\n}\n")


def shared(g, *names):
    d = {"s1": S1, "s2": S2, "ties": TIES, "pos": POS, "ang": ANG, "de": DE}
    for n in names:
        g.decl(n, d[n])


# ---------------------------------------------------------------- descriptive
g = G("descriptive")
shared(g, "s1", "s2", "ties", "pos")
for p in (25, 50, 75):
    g.val(f"scoreatpercentile {p}", f"statistics::scoreatpercentile(s2[..], {float(p)})", stats.scoreatpercentile(S2, p))
    g.val(f"scoreatpercentile {p} limit", f"statistics::scoreatpercentile_limit(s2[..], {float(p)}, 2.0, 6.0)", stats.scoreatpercentile(S2, p, limit=(2.0, 6.0)))
for k, m in enumerate(("fraction", "lower", "higher")):
    g.val(f"scoreatpercentile 25 {m}", f"statistics::scoreatpercentile_method(s2[..], 25.0, {k})", stats.scoreatpercentile(S2, 25, interpolation_method=m))
    g.val(f"scoreatpercentile limit 30 {m}", f"statistics::scoreatpercentile_limit_method(s2[..], 30.0, 0.0, 10.5, {k})", stats.scoreatpercentile(S2, 30, limit=(0.0, 10.5), interpolation_method=m))
    g.val(f"scoreatpercentile limit 60 {m}", f"statistics::scoreatpercentile_limit_method(s2[..], 60.0, 0.0, 10.5, {k})", stats.scoreatpercentile(S2, 60, limit=(0.0, 10.5), interpolation_method=m))
g.val("scoreatpercentile limit 50", "statistics::scoreatpercentile_limit(s2[..], 50.0, 0.0, 10.5)", stats.scoreatpercentile(S2, 50, limit=(0.0, 10.5)))
g.val("scoreatpercentile limit 80", "statistics::scoreatpercentile_limit(s2[..], 80.0, 0.0, 10.5)", stats.scoreatpercentile(S2, 80, limit=(0.0, 10.5)))
for k, kind in enumerate(("rank", "weak", "strict", "mean")):
    g.val(f"percentileofscore {kind} 3", f"statistics::percentileofscore(s2[..], 3.0, {k})", stats.percentileofscore(S2, 3.0, kind=kind))
    g.val(f"percentileofscore {kind} ties", f"statistics::percentileofscore(ties[..], 2.5, {k})", stats.percentileofscore(TIES, 2.5, kind=kind))
    g.val(f"percentileofscore {kind} above", f"statistics::percentileofscore(s2[..], 99.0, {k})", stats.percentileofscore(S2, 99.0, kind=kind))
g.val("median_abs_deviation", "statistics::median_abs_deviation(s2[..], 1.0, false)", stats.median_abs_deviation(S2))
g.val("median_abs_deviation normal", "statistics::median_abs_deviation_normal(s2[..])", stats.median_abs_deviation(S2, scale="normal"))
g.val("median_abs_deviation mean", "statistics::median_abs_deviation(s2[..], 1.0, true)", stats.median_abs_deviation(S2, center=np.mean))
for p in (0.0, 1.0, 2.0, -1.0, 0.5, 1e-6, -1e-6):
    g.val(f"pmean {p}", f"statistics::pmean(pos[..], {p})", stats.pmean(POS, p))
g.raw("let w: f64[8] = [1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0];")
g.val("pmean weighted", "statistics::pmean_weighted(pos[..], w[..], 1.0)", stats.pmean(POS, 1, weights=[1, 1, 2, 2, 3, 3, 4, 4]))
for d in (0, 1):
    g.val(f"gstd {d}", f"statistics::gstd(pos[..], {d})", stats.gstd(POS, ddof=d))
    g.lst(f"gzscore {d}", f"statistics::gzscore(out, pos[..], {d})", stats.gzscore(POS, ddof=d))
    g.lst(f"zmap {d}", f"statistics::zmap(out, s2[..], s1[..], {d})", stats.zmap(S2, S1, ddof=d))
for bias in (True, False):
    for d in (0, 1):
        r = stats.describe(S2, ddof=d, bias=bias)
        g.raw(f"let ds{int(bias)}{d}, dc{int(bias)}{d} = statistics::describe(s2[..], {d}, {str(bias).lower()});")
        want = [int(r.nobs), r.minmax[0], r.minmax[1], r.mean, r.variance, r.skewness, r.kurtosis]
        g.raw(f"let dw{int(bias)}{d}: f64[7] = {arr(want)};")
        g.raw(f"let dg{int(bias)}{d}: f64[7] = [ds{int(bias)}{d}.n, ds{int(bias)}{d}.min, ds{int(bias)}{d}.max, ds{int(bias)}{d}.mean, ds{int(bias)}{d}.variance, ds{int(bias)}{d}.skewness, ds{int(bias)}{d}.kurtosis];")
        g.raw(f'chk_list(bad, "describe bias={bias} ddof={d}", dg{int(bias)}{d}[..], dc{int(bias)}{d}, dw{int(bias)}{d}[..], 1e-12);')
g.lst("lmoment default", "statistics::lmoment(out, s2[..], 4, true)", stats.lmoment(S2))
g.lst("lmoment raw", "statistics::lmoment(out, s2[..], 4, false)", stats.lmoment(S2, standardize=False))
g.lst("lmoment order 3", "statistics::lmoment(out, s2[..], 3, true)", [stats.lmoment(S2, order=k) for k in (1, 2, 3)])
for a in (0.5, 0.2, 0.8, 0.0, 1.0):
    g.raw(f"let ex{int(a*10)}: f64[4] = [1.0, 4.0, 2.0, -1.0];")
    g.val(f"expectile {a}", f"statistics::expectile(ex{int(a*10)}[..], {a})", stats.expectile([1.0, 4.0, 2.0, -1.0], a))
g.raw("let ew: f64[4] = [1.0, 3.0, 1.0, 1.0];")
g.val("expectile weighted", "statistics::expectile_weighted(ex8[..], ew[..], 0.8)", stats.expectile([1.0, 4.0, 2.0, -1.0], 0.8, weights=[1, 3, 1, 1]))
g.val("expectile mean", "statistics::expectile(s2[..], 0.5)", stats.expectile(S2, 0.5))
g.write()

# ---------------------------------------------------------------- frequency
g = G("frequency")
shared(g, "s2", "ties")
e = stats.ecdf(S2)
g.lst("ecdf quantiles", "statistics::ecdf(out, out2, s2[..])", e.cdf.quantiles)
g.lst("ecdf probabilities", "statistics::ecdf(out, out2, s2[..])", e.cdf.probabilities, name="out2")
et = stats.ecdf([2.5, -1.0, 4.0, 3.0, 2.5])
g.raw("let tt: f64[5] = [2.5, -1.0, 4.0, 3.0, 2.5];")
g.lst("ecdf ties quantiles", "statistics::ecdf(out, out2, tt[..])", et.cdf.quantiles)
g.lst("ecdf ties probabilities", "statistics::ecdf(out, out2, tt[..])", et.cdf.probabilities, name="out2")
g.val("ecdf_value 3", "statistics::ecdf_value(s2[..], 3.0)", float(e.cdf.evaluate(3.0)))
ci = e.cdf.confidence_interval(0.95)
g.lst("ecdf_interval low", "statistics::ecdf_interval(out, out2, s2[..], 0.95)", ci.low.probabilities)
g.lst("ecdf_interval high", "statistics::ecdf_interval(out, out2, s2[..], 0.95)", ci.high.probabilities, name="out2")
unc, rt = [1.0, 2.0, 3.0, 5.0, 8.0], [4.0, 6.0, 7.0]
km = stats.ecdf(stats.CensoredData(uncensored=unc, right=rt))
g.raw("let unc: f64[5] = [1.0, 2.0, 3.0, 5.0, 8.0];")
g.raw("let rt: f64[3] = [4.0, 6.0, 7.0];")
g.lst("ecdf_censored times", "statistics::ecdf_censored(out, out2, unc[..], rt[..])", km.sf.quantiles)
g.lst("ecdf_censored survival", "statistics::ecdf_censored(out, out2, unc[..], rt[..])", km.sf.probabilities, name="out2")
kci = km.sf.confidence_interval(0.95)
g.lst("ecdf_censored_interval low", "statistics::ecdf_censored_interval(out, out2, unc[..], rt[..], 0.95)", kci.low.probabilities)
g.lst("ecdf_censored_interval high", "statistics::ecdf_censored_interval(out, out2, unc[..], rt[..], 0.95)", kci.high.probabilities, name="out2")
for k, m in enumerate(("interpolated_inverted_cdf", "hazen", "weibull", "linear", "median_unbiased", "normal_unbiased", "inverted_cdf", "averaged_inverted_cdf", "closest_observation")):
    for y in (3.0, 5.5, 0.7):
        g.val(f"estimated_cdf {m} {y}", f"statistics::estimated_cdf(s2[..], {y}, {k})", stats.estimated_cdf(S2, y, method=m))
g.raw("let cf: f64[6] = [1.0, 4.0, 2.0, 1.0, 3.0, 1.0];")
c1 = stats.cumfreq([1, 4, 2, 1, 3, 1], numbins=4, defaultreallimits=(1.5, 5))
r1 = stats.relfreq([1, 4, 2, 1, 3, 1], numbins=4, defaultreallimits=(1.5, 5))
g.raw("let fq1, fc1 = statistics::cumfreq_limits(out, cf[..], 4, 1.5, 5.0);")
g.raw(f"let fw1: f64[4] = {arr(c1.cumcount)};")
g.raw('chk_list(bad, "cumfreq counts", out[..], fc1, fw1[..], 1e-12);')
g.raw(f'chk_val(bad, "cumfreq lowerlimit", fq1.lowerlimit, fc1, {f(c1.lowerlimit)}, 1e-12);')
g.raw(f'chk_val(bad, "cumfreq binsize", fq1.binsize, fc1, {f(c1.binsize)}, 1e-12);')
g.raw(f'chk_val(bad, "cumfreq extrapoints", fq1.extrapoints, fc1, {f(c1.extrapoints)}, 1e-12);')
g.raw("let fq2, fc2 = statistics::relfreq_limits(out, cf[..], 4, 1.5, 5.0);")
g.raw(f"let fw2: f64[4] = {arr(r1.frequency)};")
g.raw('chk_list(bad, "relfreq frequency", out[..], fc2, fw2[..], 1e-12);')
g.raw(f'chk_val(bad, "relfreq extrapoints", fq2.extrapoints, fc2, {f(r1.extrapoints)}, 1e-12);')
c2 = stats.cumfreq(S2, numbins=3)
g.raw("let fq3, fc3 = statistics::cumfreq(out, s2[..], 3);")
g.raw(f"let fw3: f64[3] = {arr(c2.cumcount)};")
g.raw('chk_list(bad, "cumfreq default counts", out[..], fc3, fw3[..], 1e-12);')
g.raw(f'chk_val(bad, "cumfreq default lowerlimit", fq3.lowerlimit, fc3, {f(c2.lowerlimit)}, 1e-12);')
g.raw(f'chk_val(bad, "cumfreq default binsize", fq3.binsize, fc3, {f(c2.binsize)}, 1e-12);')
g.raw("let cws: f64[4] = [1.0, 2.0, 3.0, 4.0];")
g.raw("let cw: f64[4] = [1.0, 2.0, 3.0, 4.0];")
cw1 = stats.cumfreq([1, 2, 3, 4], numbins=4, weights=[1, 2, 3, 4])
rw1 = stats.relfreq([1, 2, 3, 4], numbins=4, weights=[1, 2, 3, 4])
g.raw("let fqw, fqwc = statistics::cumfreq_weighted(out, cws[..], cw[..], 4);")
g.raw(f"let fww: f64[4] = {arr(cw1.cumcount)};")
g.raw('chk_list(bad, "cumfreq weighted counts", out[..], fqwc, fww[..], 1e-12);')
g.raw(f'chk_val(bad, "cumfreq weighted lowerlimit", fqw.lowerlimit, fqwc, {f(cw1.lowerlimit)}, 1e-12);')
g.raw(f'chk_val(bad, "cumfreq weighted binsize", fqw.binsize, fqwc, {f(cw1.binsize)}, 1e-12);')
g.raw(f'chk_val(bad, "cumfreq weighted extrapoints", fqw.extrapoints, fqwc, {f(cw1.extrapoints)}, 1e-12);')
g.raw("let fqw2, fqwc2 = statistics::relfreq_weighted(out, cws[..], cw[..], 4);")
g.raw(f"let fww2: f64[4] = {arr(rw1.frequency)};")
g.raw('chk_list(bad, "relfreq weighted frequency", out[..], fqwc2, fww2[..], 1e-12);')
g.raw(f'chk_val(bad, "relfreq weighted extrapoints", fqw2.extrapoints, fqwc2, {f(rw1.extrapoints)}, 1e-12);')
g.val("tiecorrect", "statistics::tiecorrect(rk[..])", stats.tiecorrect(stats.rankdata([1, 2.5, 2.5, 4])))
g.lines.insert(g.lines.index([l for l in g.lines if "tiecorrect" in l][0]) - 1, "    let rk: f64[4] = " + arr(stats.rankdata([1, 2.5, 2.5, 4])) + ";")
g.raw("let rk2: f64[9] = " + arr(stats.rankdata([1, 3, 2, 4, 5, 7, 2, 8, 4])) + ";")
g.val("tiecorrect nine", "statistics::tiecorrect(rk2[..])", stats.tiecorrect(stats.rankdata([1, 3, 2, 4, 5, 7, 2, 8, 4])))
g.raw("let rk3: f64[3] = [1.0, 2.0, 3.0];")
g.val("tiecorrect none", "statistics::tiecorrect(rk3[..])", 1.0)
g.write()

# ---------------------------------------------------------------- entropy
g = G("entropy")
shared(g, "s2", "de")
g.raw("let h2: f64[2] = [0.5, 0.5];")
g.raw("let h3: f64[2] = [0.9, 0.1];")
g.raw("let h4: f64[2] = [1.0, 3.0];")
g.val("entropy", "statistics::entropy(h2[..], 0.0)", stats.entropy([0.5, 0.5]))
g.val("entropy base 2", "statistics::entropy(h2[..], 2.0)", stats.entropy([0.5, 0.5], base=2))
g.val("entropy unnormalized", "statistics::entropy(h4[..], 0.0)", stats.entropy([1.0, 3.0]))
g.val("relative_entropy", "statistics::relative_entropy(h2[..], h3[..], 0.0)", stats.entropy([0.5, 0.5], [0.9, 0.1]))
for k, m in ((1, "vasicek"), (2, "van es"), (3, "ebrahimi"), (4, "correa")):
    g.val(f"differential_entropy {m}", f"statistics::differential_entropy(de[..], {k}, 0, 0.0)", stats.differential_entropy(DE, method=m))
    g.val(f"differential_entropy {m} window 3", f"statistics::differential_entropy(de[..], {k}, 3, 0.0)", stats.differential_entropy(DE, method=m, window_length=3))
g.val("differential_entropy auto", "statistics::differential_entropy(de[..], 0, 0, 0.0)", stats.differential_entropy(DE))
g.val("differential_entropy base 2", "statistics::differential_entropy(de[..], 1, 0, 2.0)", stats.differential_entropy(DE, method="vasicek", base=2))
for n in (1, 2, 3, 4):
    g.val(f"kstat {n}", f"statistics::kstat(s2[..], {n})", stats.kstat(S2, n))
for n in (1, 2):
    g.val(f"kstatvar {n}", f"statistics::kstatvar(s2[..], {n})", stats.kstatvar(S2, n))
g.write()

# ---------------------------------------------------------------- trimmed
g = G("trimmed")
shared(g, "s2")
g.val("tmin", "statistics::tmin(s2[..], -1.0e300, true)", stats.tmin(S2))
g.val("tmin lower 3", "statistics::tmin(s2[..], 3.0, true)", stats.tmin(S2, lowerlimit=3))
g.val("tmin lower 3 exclusive", "statistics::tmin(s2[..], 3.0, false)", stats.tmin(S2, lowerlimit=3, inclusive=False))
g.val("tmax", "statistics::tmax(s2[..], 1.0e300, true)", stats.tmax(S2))
g.val("tmax upper 4", "statistics::tmax(s2[..], 4.0, true)", stats.tmax(S2, upperlimit=4))
g.val("tmax upper 4 exclusive", "statistics::tmax(s2[..], 4.0, false)", stats.tmax(S2, upperlimit=4, inclusive=False))
g.val("tvar", "statistics::tvar(s2[..], 0.0, 0.0, false, false, true, true, 1)", stats.tvar(S2))
g.val("tstd", "statistics::tstd(s2[..], 0.0, 0.0, false, false, true, true, 1)", stats.tstd(S2))
g.val("tstd limits", "statistics::tstd(s2[..], 2.0, 6.0, true, true, true, true, 1)", stats.tstd(S2, limits=(2.0, 6.0)))
g.val("tsem", "statistics::tsem(s2[..], 0.0, 0.0, false, false, true, true, 1)", stats.tsem(S2))
g.val("tsem limits", "statistics::tsem(s2[..], 2.0, 6.0, true, true, true, true, 1)", stats.tsem(S2, limits=(2.0, 6.0)))
g.lst("trimboth 0.1", "statistics::trimboth(out, s2[..], 0.1)", sorted(stats.trimboth(S2, 0.1)))
g.lst("trimboth 0.25", "statistics::trimboth(out, s2[..], 0.25)", sorted(stats.trimboth(S2, 0.25)))
g.lst("trim1 right", "statistics::trim1(out, s2[..], 0.25, 1)", sorted(stats.trim1(S2, 0.25, tail="right")))
g.lst("trim1 right large", "statistics::trim1(out, s2[..], 0.6, 1)", sorted(stats.trim1(S2, 0.6, tail="right")))
g.lst("trim1 left", "statistics::trim1(out, s2[..], 0.25, 0)", sorted(stats.trim1(S2, 0.25, tail="left")))
g.write()

# ---------------------------------------------------------------- circular
g = G("circular")
shared(g, "ang")
g.val("circmean", "statistics::circmean(ang[..], 0.0, 6.283185307179586)", stats.circmean(ANG))
g.val("circmean low -pi", "statistics::circmean(ang[..], -3.141592653589793, 3.141592653589793)", stats.circmean(ANG, low=-math.pi, high=math.pi))
g.val("circvar", "statistics::circvar(ang[..], 0.0, 6.283185307179586)", stats.circvar(ANG))
g.val("circstd", "statistics::circstd(ang[..], 0.0, 6.283185307179586, false)", stats.circstd(ANG))
g.val("circstd normalized", "statistics::circstd(ang[..], 0.0, 6.283185307179586, true)", stats.circstd(ANG, normalize=True))
g.val("circmean range 10", "statistics::circmean(ang[..], 0.0, 10.0)", stats.circmean(ANG, high=10.0))
g.val("circstd range 10", "statistics::circstd(ang[..], 0.0, 10.0, false)", stats.circstd(ANG, high=10.0))
for nm, rows in (("a", [[3.0, 4.0], [6.0, -8.0]]), ("b", [[0.8660254, 0.5, 0.0], [0.8660254, -0.5, 0.0]]), ("c", [[1.0, 2.0, 2.0], [2.0, 0.5, 1.0], [-1.0, 3.0, 0.5]])):
    d = stats.directional_stats(np.array(rows))
    flat = [x for r in rows for x in r]
    dims = len(rows[0])
    g.raw(f"let dr{nm}: f64[{len(flat)}] = {arr(flat)};")
    g.lst(f"directional {nm} direction", f"statistics::directional_stats(out, dr{nm}[..], {dims}, true)", d.mean_direction)
    g.val(f"directional {nm} resultant", f"statistics::directional_stats(out, dr{nm}[..], {dims}, true)", d.mean_resultant_length)
g.write()

# ---------------------------------------------------------------- binned
g = G("binned")
g.raw("let xb: f64[8] = " + arr(XB) + ";")
g.raw("let vb: f64[8] = " + arr(VB) + ";")
names = ("mean", "std", "median", "count", "sum", "min", "max")
for k, s in enumerate(names):
    r = stats.binned_statistic(XB, VB, statistic=s, bins=4)
    g.raw(f"let bt{k}, bc{k} = statistics::binned_statistic(out, out2, bshape, bnum, bempty, xb[..], vb[..], {k}, 4);")
    g.raw(f"let bw{k}: f64[4] = {arr(r.statistic)};")
    g.raw(f'chk_list(bad, "binned_statistic {s}", out[..], bc{k}, bw{k}[..], 1e-12);')
    g.raw("out.resize(0); out2.resize(0); bshape.resize(0); bnum.resize(0); bempty.resize(0);")
edge_x = [1.0, 2.75, 4.5, 6.25, 8.0]
edge_v = [10.0, 20.0, 30.0, 40.0, 50.0]
re_ = stats.binned_statistic(edge_x, edge_v, statistic="count", bins=4)
g.raw("out.resize(0); out2.resize(0); bshape.resize(0); bnum.resize(0); bempty.resize(0);")
g.raw("let ex: f64[5] = " + arr(edge_x) + ";")
g.raw("let ev: f64[5] = " + arr(edge_v) + ";")
g.raw("let be2, bec2 = statistics::binned_statistic(out, out2, bshape, bnum, bempty, ex[..], ev[..], 3, 4);")
g.raw("let bew2: f64[4] = " + arr(re_.statistic) + ";")
g.raw('chk_list(bad, "binned_statistic on the edges", out[..], bec2, bew2[..], 1e-12);')
g.raw("let ben2: f64[5] = " + arr(re_.binnumber) + ";")
g.raw('chk_list(bad, "binned_statistic on the edges, bin numbers", bnum[..], bec2, ben2[..], 0.0);')
g.raw("out.resize(0); out2.resize(0); bshape.resize(0); bnum.resize(0); bempty.resize(0);")
g.raw("out.resize(0); out2.resize(0); bshape.resize(0); bnum.resize(0); bempty.resize(0);")
ge = [1.0, 2.75, 4.5, 6.25, 8.0]
for k, s in enumerate(names):
    r = stats.binned_statistic(XB, VB, statistic=s, bins=ge)
    g.raw("out.resize(0); out2.resize(0); bshape.resize(0); bnum.resize(0); bempty.resize(0);")
    g.raw("let ge: f64[5] = " + arr(ge) + ";")
    g.raw(f"let btge{k}, bcge{k} = statistics::binned_statistic_edges(out, out2, bshape, bnum, bempty, xb[..], vb[..], {k}, ge[..]);")
    g.raw(f"let bwge{k}: f64[4] = {arr(r.statistic)};")
    g.raw(f'chk_list(bad, "binned_statistic_edges {s}", out[..], bcge{k}, bwge{k}[..], 1e-12);')
    g.raw(f"let bwnge{k}: f64[8] = {arr(r.binnumber)};")
    g.raw(f'chk_list(bad, "binned_statistic_edges {s} bin numbers", bnum[..], bcge{k}, bwnge{k}[..], 0.0);')
r = stats.binned_statistic(XB, VB, statistic="count", bins=4)
g.raw("let be, bec = statistics::binned_statistic(bst, out2, bshape, bnum, bempty, xb[..], vb[..], 3, 4);")
g.raw(f"let bew: f64[5] = {arr(r.bin_edges)};")
g.raw('chk_list(bad, "binned_statistic edges", out2[..], bec, bew[..], 1e-12);')
g.raw(f"let bnw: f64[8] = {arr(r.binnumber)};")
g.raw('chk_list(bad, "binned_statistic binnumber", bnum[..], bec, bnw[..], 0.0);')
g.raw("out.resize(0); out2.resize(0); bshape.resize(0); bnum.resize(0); bempty.resize(0);")
r2 = stats.binned_statistic_2d(XB, XB, VB, statistic="mean", bins=[3, 3])
g.raw("let b2, b2c = statistics::binned_statistic_2d(out, out2, bshape, bnum, bempty, xb[..], xb[..], vb[..], 0, 3, 3);")
g.raw("let b2w: f64[9] = " + arr([0.0 if math.isnan(x) else x for x in r2.statistic.ravel()]) + ";")
g.raw('chk_list(bad, "binned_statistic_2d mean", out[..], b2c, b2w[..], 1e-12);')
g.raw("let b2e: f64[9] = " + arr([1.0 if math.isnan(x) else 0.0 for x in r2.statistic.ravel()]) + ";")
g.raw('chk_list(bad, "binned_statistic_2d empty", bempty[..], b2c, b2e[..], 0.0);')
g.raw("let b2n: f64[8] = " + arr(r2.binnumber) + ";")
g.raw('chk_list(bad, "binned_statistic_2d binnumber", bnum[..], b2c, b2n[..], 0.0);')
g.raw("let b2x: f64[8] = " + arr(list(r2.x_edge) + list(r2.y_edge)) + ";")
g.raw('chk_list(bad, "binned_statistic_2d edges", out2[..], b2c, b2x[..], 1e-12);')
samp = np.array([[1.0, 1.0], [2.0, 2.0], [3.0, 1.0], [4.0, 2.0], [5.0, 1.0]])
r3 = stats.binned_statistic_dd(samp, np.arange(5), statistic="count", bins=[2, 2])
g.raw("out.resize(0); out2.resize(0); bshape.resize(0); bnum.resize(0); bempty.resize(0);")
g.raw("let sd: f64[10] = " + arr(samp.ravel()) + ";")
g.raw("let vd: f64[5] = [0.0, 1.0, 2.0, 3.0, 4.0];")
g.raw("let b3, b3c = statistics::binned_statistic_dd(out, out2, bshape, bnum, bempty, sd[..], 2, vd[..], 3, 2);")
g.raw("let b3w: f64[4] = " + arr(r3.statistic.ravel()) + ";")
g.raw('chk_list(bad, "binned_statistic_dd count", out[..], b3c, b3w[..], 1e-12);')
g.raw("let b3n: f64[5] = " + arr(r3.binnumber) + ";")
g.raw('chk_list(bad, "binned_statistic_dd binnumber", bnum[..], b3c, b3n[..], 0.0);')
# The calls with a trailing block. The block gives the sum of the values of
# a bin, as the callable statistic of SciPy does.
g.raw("out.resize(0); out2.resize(0); bshape.resize(0); bnum.resize(0); bempty.resize(0);")
rw = stats.binned_statistic(XB, VB, statistic=np.sum, bins=4)
g.raw("let btw, btwc = statistics::binned_statistic_with(out, out2, bshape, bnum, bempty, xb[..], vb[..], 4) { v => sum_of(v) };")
g.raw(f"let btwr: f64[4] = {arr(rw.statistic)};")
g.raw('chk_list(bad, "binned_statistic_with sum", out[..], btwc, btwr[..], 1e-12);')
g.raw("out.resize(0); out2.resize(0); bshape.resize(0); bnum.resize(0); bempty.resize(0);")
rw2 = stats.binned_statistic_2d(XB, XB, VB, statistic=np.sum, bins=[3, 3])
g.raw("let btw2, btw2c = statistics::binned_statistic_2d_with(out, out2, bshape, bnum, bempty, xb[..], xb[..], vb[..], 3, 3) { v => sum_of(v) };")
g.raw(f"let btw2r: f64[9] = {arr(rw2.statistic.ravel())};")
g.raw('chk_list(bad, "binned_statistic_2d_with sum", out[..], btw2c, btw2r[..], 1e-12);')
g.raw("out.resize(0); out2.resize(0); bshape.resize(0); bnum.resize(0); bempty.resize(0);")
rw3 = stats.binned_statistic_dd(samp, np.arange(5), statistic=np.sum, bins=[2, 2])
g.raw("let btw3, btw3c = statistics::binned_statistic_dd_with(out, out2, bshape, bnum, bempty, sd[..], 2, vd[..], 2) { v => sum_of(v) };")
g.raw(f"let btw3r: f64[4] = {arr(rw3.statistic.ravel())};")
g.raw('chk_list(bad, "binned_statistic_dd_with sum", out[..], btw3c, btw3r[..], 1e-12);')
g.write()
