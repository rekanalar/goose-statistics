#!/usr/bin/env python3
"""Writes test/suite/resample_vectors.goose from SciPy.

Run it from a directory that holds no other Python file:
    python3 test/make_resample_vectors.py > test/suite/resample_vectors.goose

The file holds numbers and the checks that use them. Rows that depend on
a random stream use a large SciPy run as the oracle, and the Goose test
compares within a stated Monte Carlo bound.
"""
import math
import sys

import numpy as np
import scipy
import scipy.stats as st
from scipy.special import ndtri

A = [2.1, 3.4, 1.9, 5.6, 4.4, 3.8, 2.7, 4.9]
B = [3.9, 5.2, 4.8, 6.1, 5.5, 4.1, 6.6, 5.0]


def g(x):
    return repr(float(x))


def lst(xs):
    return "[" + ", ".join(repr(float(x)) for x in xs) + "]"


def fixed(name, xs):
    return "    let %s: f64[%d] = %s;\n" % (name, len(xs), lst(xs))


out = []
w = out.append

# ---- the bootstrap pieces from one fixed list of bootstrap values ----
a = np.array(A)
theta = a.mean()
res = {}
for m in ("percentile", "basic", "BCa"):
    rng = np.random.default_rng(3)
    r = st.bootstrap((a,), np.mean, method=m, n_resamples=999, rng=rng)
    res[m] = (r.confidence_interval.low, r.confidence_interval.high)
dist = r.bootstrap_distribution
p0 = ((dist < theta).sum() + (dist <= theta).sum()) / (2 * len(dist))
z0 = float(ndtri(p0))
n = len(a)
jk = np.array([np.delete(a, i).mean() for i in range(n)])
U = (n - 1) * (jk.mean() - jk)
ahat = float((U ** 3).sum() / n ** 3 / 6 / ((U ** 2).sum() / n ** 2) ** 1.5)
# a one-sided interval from the same list
rng = np.random.default_rng(3)
rl = st.bootstrap((a,), np.mean, method="BCa", n_resamples=999, rng=rng,
                  alternative="greater", confidence_level=0.9)
rng = np.random.default_rng(3)
rp = st.bootstrap((a,), np.mean, method="percentile", n_resamples=999, rng=rng,
                  alternative="less", confidence_level=0.9)

# ---- the exact permutation tests ----
pe = st.permutation_test((A, B), lambda x, y: np.mean(x) - np.mean(y),
                         n_resamples=np.inf, rng=1)
diff = np.array(A) - np.array(B)
ps = st.permutation_test((diff,), lambda x: np.mean(x), permutation_type="samples",
                         n_resamples=np.inf, rng=1)
pl = st.permutation_test((A, B), lambda x, y: np.mean(x) - np.mean(y),
                         n_resamples=np.inf, alternative="less", rng=1)
pp = st.permutation_test((A[:5], B[:5]), lambda x, y: st.pearsonr(x, y)[0],
                         permutation_type="pairings", n_resamples=np.inf, rng=1)
# paired samples of three values: 2^3 swaps of the two samples
pw = st.permutation_test((A[:4], B[:4]), lambda x, y: np.mean(x) - np.mean(y),
                         permutation_type="samples", n_resamples=np.inf, rng=1)

# ---- the large SciPy runs that give the oracles of the seeded rows ----
BIG = 400000
rng = np.random.default_rng(21)
big_b = st.bootstrap((a,), np.mean, method="percentile", n_resamples=BIG, rng=rng)
rng = np.random.default_rng(22)
big_bca = st.bootstrap((a,), np.mean, method="BCa", n_resamples=BIG, rng=rng)
rng = np.random.default_rng(23)
big_basic = st.bootstrap((a,), np.mean, method="basic", n_resamples=BIG, rng=rng)

# A normal with mean 3 and sd 1.5 as the null. The two samples of the
# draws are independent of the stream of Goose.
rg = np.random.default_rng(24)
mcs = st.monte_carlo_test(A, lambda size: rg.normal(3.0, 1.5, size),
                          lambda x, axis: np.mean(x, axis=axis),
                          n_resamples=BIG, alternative="two-sided")
rg = np.random.default_rng(25)
mcl = st.monte_carlo_test(A, lambda size: rg.normal(3.0, 1.5, size),
                          lambda x, axis: np.mean(x, axis=axis),
                          n_resamples=BIG, alternative="greater")
# Power of the one-sample t test, exact by the noncentral t distribution.
nn, delta, level = 20, 0.5, 0.05
tcrit = st.t.ppf(1 - level / 2, nn - 1)
nc = delta * math.sqrt(nn)
power_exact = float(st.nct.sf(tcrit, nn - 1, nc) + st.nct.cdf(-tcrit, nn - 1, nc))

# ---- the Sobol estimators on fixed evaluation lists ----
srng = np.random.default_rng(11)
sn, sd = 8, 3
fA = srng.normal(size=sn)
fB = srng.normal(size=sn)
fAB = srng.normal(size=(sd, sn))
si = st.sobol_indices(func={"f_A": fA, "f_B": fB, "f_AB": fAB.reshape(sd, 1, sn)}, n=sn)

w("// The answers of SciPy %s for the resampling calls. The generator\n" % scipy.__version__)
w("// test/make_resample_vectors.py writes this file. It holds the numbers\n")
w("// and the checks that use them.\n")
w("namespace suite;\n\nimport std;\nimport math;\nimport statistics from \"../..\";\n\n")
w("// The list of a literal.\n")
w("fn rv_list(vals: const f64[:]) -> f64[>..] {\n    var out: f64[>..] = [];\n    for v in vals { out.push(v); }\n    return out;\n}\n\n")
w("// The sum divided by the count of the list.\n")
w("fn rv_mean(s: const f64[:]) -> f64 {\n    var t: f64 = 0.0;\n    for v in s { t += v; }\n    return t / s.len;\n}\n\n")
w("// The mean of the first half minus the mean of the second half.\n")
w("fn rv_diff(s: const f64[:]) -> f64 {\n    let h = s.len / 2;\n    return rv_mean(s[0..h]) - rv_mean(s[h..s.len]);\n}\n\n")
w("// The Pearson correlation of the two halves.\n")
w("fn rv_corr(s: const f64[:]) -> f64 {\n    let h = s.len / 2;\n    let mx = rv_mean(s[0..h]);\n    let my = rv_mean(s[h..s.len]);\n")
w("    var sxy: f64 = 0.0;\n    var sxx: f64 = 0.0;\n    var syy: f64 = 0.0;\n")
w("    for i in 0..h {\n        sxy += (s[i] - mx) * (s[h + i] - my);\n        sxx += (s[i] - mx) * (s[i] - mx);\n        syy += (s[h + i] - my) * (s[h + i] - my);\n    }\n")
w("    return sxy / sqrt(sxx * syy);\n}\n\n")
w("fn rv_ttest_p(s: const f64[:]) -> f64 {\n    let r, c = statistics::ttest_1samp(s, 0.0);\n    return r.pvalue;\n}\n\n")
w("fn run_resample_vectors() -> i64 {\n    var bad: i64 = 0;\n")
w("    let a = rv_list([%s]);\n" % ", ".join(repr(x) for x in A))
w("    let b = rv_list([%s]);\n" % ", ".join(repr(x) for x in B))
w("    var ab: f64[>..] = [];\n    for x in a[..] { ab.push(x); }\n    for x in b[..] { ab.push(x); }\n")
w("    let two: i64[2] = [8, 8];\n    let one: i64[1] = [8];\n")
w(fixed("boot", dist).replace("let boot: f64[999]", "let boot: f64[999]"))
w("    let theta = %r;\n" % float(theta))
w("    let ahat = %r;\n    let z0 = %r;\n" % (ahat, z0))

w("""
    // The jackknife acceleration and the bias term of the BCa interval.
    let ga, ca = statistics::rs_accel(a[..]) { s => rv_mean(s) };
    check(bad, ca == 0 && near_tol(ga, ahat, 1e-11), "bca, the acceleration of the jackknife");
    let gz, cz = statistics::rs_bias(boot[..], theta);
    check(bad, cz == 0 && near_tol(gz, z0, 1e-13), "bca, the bias term");

    // The three intervals of one fixed list of bootstrap values.
    let lo1, hi1, c1 = statistics::rs_interval(boot[..], theta, 0.0, 0.95, statistics::ALT_TWO_SIDED, statistics::INTERVAL_PERCENTILE);
    check(bad, c1 == 0 && near_tol(lo1, %r, 1e-13) && near_tol(hi1, %r, 1e-13), "percentile interval, the fixed list");
    let lo2, hi2, c2 = statistics::rs_interval(boot[..], theta, 0.0, 0.95, statistics::ALT_TWO_SIDED, statistics::INTERVAL_BASIC);
    check(bad, c2 == 0 && near_tol(lo2, %r, 1e-12) && near_tol(hi2, %r, 1e-12), "basic interval, the fixed list");
    let lo3, hi3, c3 = statistics::rs_interval(boot[..], theta, ahat, 0.95, statistics::ALT_TWO_SIDED, statistics::INTERVAL_BCA);
    check(bad, c3 == 0 && near_tol(lo3, %r, 1e-12) && near_tol(hi3, %r, 1e-12), "bca interval, the fixed list");
    let lo4, hi4, c4 = statistics::rs_interval(boot[..], theta, ahat, 0.9, statistics::ALT_GREATER, statistics::INTERVAL_BCA);
    check(bad, c4 == 0 && near_tol(lo4, %r, 1e-12) && hi4 > 1e300, "bca interval, one-sided greater");
    let lo5, hi5, c5 = statistics::rs_interval(boot[..], theta, 0.0, 0.9, statistics::ALT_LESS, statistics::INTERVAL_PERCENTILE);
    check(bad, c5 == 0 && near_tol(hi5, %r, 1e-13) && lo5 < -1e300, "percentile interval, one-sided less");
""" % (res["percentile"][0], res["percentile"][1], res["basic"][0], res["basic"][1],
       res["BCa"][0], res["BCa"][1], rl.confidence_interval.low, rp.confidence_interval.high))

w("""
    // The exact permutation tests. All splits are listed, so no draw is made.
    var r0 = rng { s: 1 };
    var nul: f64[>..] = [];
    let exact = statistics::ResampleMethod { kind: statistics::RESAMPLE_PERMUTATION, resamples: 1000000, .. };
    let t1, k1 = statistics::permutation_test(nul, ab[..], two[..], exact, statistics::ALT_TWO_SIDED, r0) { s => rv_diff(s) };
    check(bad, k1 == 0 && nul.len == %d && near_tol(t1.pvalue, %r, 1e-13) && near_tol(t1.statistic, %r, 1e-12), "permutation, all 12870 splits");
    let t2, k2 = statistics::permutation_test(nul, ab[..], two[..], exact, statistics::ALT_LESS, r0) { s => rv_diff(s) };
    check(bad, k2 == 0 && near_tol(t2.pvalue, %r, 1e-13), "permutation, the less side");
    var d: f64[8];
    for i in 0..8 { d[i] = a[i] - b[i]; }
    let sgn = statistics::ResampleMethod { kind: statistics::RESAMPLE_PERMUTATION, resamples: 1000, permutation: statistics::PERMUTE_SAMPLES, .. };
    let t3, k3 = statistics::permutation_test(nul, d[..], one[..], sgn, statistics::ALT_TWO_SIDED, r0) { s => rv_mean(s) };
    check(bad, k3 == 0 && nul.len == 256 && near_tol(t3.pvalue, %r, 1e-13), "permutation, the sign flips");
    var p5: f64[10];
    for i in 0..5 { p5[i] = a[i]; p5[5 + i] = b[i]; }
    let five: i64[2] = [5, 5];
    let prs = statistics::ResampleMethod { kind: statistics::RESAMPLE_PERMUTATION, resamples: 20000, permutation: statistics::PERMUTE_PAIRINGS, .. };
    let t4, k4 = statistics::permutation_test(nul, p5[..], five[..], prs, statistics::ALT_TWO_SIDED, r0) { s => rv_corr(s) };
    check(bad, k4 == 0 && nul.len == %d && near_tol(t4.pvalue, %r, 1e-12) && near_tol(t4.statistic, %r, 1e-12), "permutation, the pairings");
    var p4: f64[8];
    for i in 0..4 { p4[i] = a[i]; p4[4 + i] = b[i]; }
    let four: i64[2] = [4, 4];
    let swp = statistics::ResampleMethod { kind: statistics::RESAMPLE_PERMUTATION, resamples: 1000, permutation: statistics::PERMUTE_SAMPLES, .. };
    let t5, k5 = statistics::permutation_test(nul, p4[..], four[..], swp, statistics::ALT_TWO_SIDED, r0) { s => rv_diff(s) };
    check(bad, k5 == 0 && nul.len == 16 && near_tol(t5.pvalue, %r, 1e-12), "permutation, the swaps of two samples");
""" % (len(pe.null_distribution), pe.pvalue, pe.statistic, pl.pvalue, ps.pvalue,
       len(pp.null_distribution), pp.pvalue, pp.statistic, pw.pvalue))

w("""
    // The seeded bootstrap. The oracle is a SciPy run of %d resamples.
    // The Monte Carlo error of the endpoint of 9999 resamples stays
    // below 0.05, about 0.1 standard errors of the mean.
    var r1 = rng { s: 3 };
    var bo: f64[>..] = [];
    let bpm = statistics::ResampleMethod { resamples: 9999, interval: statistics::INTERVAL_PERCENTILE, .. };
    let br1, bc1 = statistics::bootstrap(bo, a[..], one[..], bpm, 0.95, statistics::ALT_TWO_SIDED, r1) { s => rv_mean(s) };
    check(bad, bc1 == 0 && bo.len == 9999 && abs(br1.interval.low - %r) < 0.05 && abs(br1.interval.high - %r) < 0.05
        && abs(br1.standard_error - %r) < 0.02, "bootstrap percentile, seeded against the large run");
    var r2 = rng { s: 4 };
    let bbm = statistics::ResampleMethod { resamples: 9999, interval: statistics::INTERVAL_BASIC, .. };
    let br2, bc2 = statistics::bootstrap(bo, a[..], one[..], bbm, 0.95, statistics::ALT_TWO_SIDED, r2) { s => rv_mean(s) };
    check(bad, bc2 == 0 && abs(br2.interval.low - %r) < 0.05 && abs(br2.interval.high - %r) < 0.05, "bootstrap basic, seeded against the large run");
    var r3 = rng { s: 5 };
    let bam = statistics::ResampleMethod { resamples: 9999, interval: statistics::INTERVAL_BCA, .. };
    let br3, bc3 = statistics::bootstrap(bo, a[..], one[..], bam, 0.95, statistics::ALT_TWO_SIDED, r3) { s => rv_mean(s) };
    check(bad, bc3 == 0 && abs(br3.interval.low - %r) < 0.05 && abs(br3.interval.high - %r) < 0.05, "bootstrap bca, seeded against the large run");
""" % (BIG, big_b.confidence_interval.low, big_b.confidence_interval.high, big_b.standard_error,
       big_basic.confidence_interval.low, big_basic.confidence_interval.high,
       big_bca.confidence_interval.low, big_bca.confidence_interval.high))

w("""
    // The seeded Monte Carlo test against a normal null of mean 3 and
    // deviation 1.5. The error of a p-value of 9999 draws stays below 0.02.
    var r4 = rng { s: 6 };
    let mcm = statistics::ResampleMethod { kind: statistics::RESAMPLE_MONTE_CARLO, resamples: 9999, .. };
    let np: f64[2] = [3.0, 1.5];
    let m1, mk1 = statistics::monte_carlo_test(nul, a[..], statistics::FAMILY_NORMAL, np[..], mcm, statistics::ALT_TWO_SIDED, r4) { s => rv_mean(s) };
    check(bad, mk1 == 0 && nul.len == 9999 && abs(m1.pvalue - %r) < 0.02 && near_tol(m1.statistic, theta, 1e-12), "monte carlo, the two-sided p-value");
    let m2, mk2 = statistics::monte_carlo_test(nul, a[..], statistics::FAMILY_NORMAL, np[..], mcm, statistics::ALT_GREATER, r4) { s => rv_mean(s) };
    check(bad, mk2 == 0 && abs(m2.pvalue - %r) < 0.02, "monte carlo, the greater side");
    let ue: f64[2] = [0.0, 6.0];
    let m3, mk3 = statistics::monte_carlo_test(nul, a[..], statistics::FAMILY_UNIFORM, ue[..], mcm, statistics::ALT_TWO_SIDED, r4) { s => rv_mean(s) };
    check(bad, mk3 == 0 && m3.pvalue > 0.1 && m3.pvalue <= 1.0, "monte carlo, a uniform null");

    // The power of the one-sample t test. The exact value comes from the
    // noncentral t distribution. The error of 20000 draws stays below 0.015.
    var r5 = rng { s: 7 };
    let pn: f64[2] = [0.5, 1.0];
    let pw, pk = statistics::power(statistics::FAMILY_NORMAL, pn[..], 20, 0.05, 20000, r5) { s => rv_ttest_p(s) };
    check(bad, pk == 0 && abs(pw - %r) < 0.015, "power, the one-sample t test");
""" % (mcs.pvalue, mcl.pvalue, power_exact))

w("""
    // The estimators of the Sobol indices on fixed evaluation lists.
    let fa = %s;
    let fb = %s;
    let fab = %s;
    let fi: f64[3] = %s;
    let ft: f64[3] = %s;
    var first: f64[>..] = [];
    var total: f64[>..] = [];
    let sn, sk = statistics::sobol_from_evaluations(first, total, fa[..], fb[..], fab[..], 3);
    check(bad, sk == 0 && sn == 3 && near_list(first[..], fi[..], 1e-12) && near_list(total[..], ft[..], 1e-12), "sobol, the estimators on fixed lists");
""" % (lst(fA), lst(fB), lst(fAB.reshape(-1)),
       lst(si.first_order), lst(si.total_order)))
w("    return bad;\n}\n")

import re
text = "".join(out)
text = re.sub(r"np\.float64\(([^)]*)\)", r"\1", text)
sys.stdout.write(text)
