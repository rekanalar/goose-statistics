#!/usr/bin/env python3
"""Writes test/suite/transform_vectors.goose: the answers of SciPy 1.18.1
for the Box-Cox and Yeo-Johnson transforms and the checks that use them.

The rows hold the transformed values, the log-likelihood, the lambda of
both searches, the confidence interval of the Box-Cox lambda and the normal
plot tables, for fixed samples and for seeded random samples.
Run it from the package directory with the system python3.
"""
import random
import warnings
import numpy as np
from scipy import stats

warnings.simplefilter("ignore")
random.seed(20261112)

POS = [0.8, 1.5, 2.2, 3.1, 4.7, 6.3, 9.8, 12.5]
MIXED = [-2.1, -0.4, 0.3, 1.1, 1.9, 2.6, 4.8, 7.5]
A = [2.1, 3.4, 1.9, 5.6, 4.4, 3.8, 2.7, 4.9]
B = [3.9, 5.2, 4.8, 6.1, 5.5, 4.1, 6.6, 5.0]
C = [1.2, 2.5, 2.2, 3.1, 1.8, 2.9, 2.4, 3.3]
NEG = [-x for x in POS]


def arr(v):
    return "[" + ", ".join(repr(float(x)) for x in v) + "]"


def samples():
    out = [("pos", POS, "pos"), ("mixed", MIXED, "mixed"), ("a", A, "pos"),
           ("b", B, "pos"), ("c", C, "pos"), ("neg", NEG, "neg")]
    for n in (5, 9, 20, 50, 120):
        v = [random.lognormvariate(0.5, 0.9) for _ in range(n)]
        out.append((f"lognormal{n}", v, "pos"))
        w = [random.gauss(1.0, 2.5) + random.expovariate(0.8) * random.choice((0, 1)) for _ in range(n)]
        out.append((f"mixed{n}", w, "mixed"))
        u = [-random.lognormvariate(0.2, 0.7) for _ in range(n)]
        out.append((f"neg{n}", u, "neg"))
        g = [random.gammavariate(2.0, 3.0) for _ in range(n)]
        out.append((f"gamma{n}", g, "pos"))
    return out


def main():
    rows = []
    for name, x, kind in samples():
        xa = np.array(x)
        tag = name
        if kind == "pos":
            for lmb in (0.0, 0.3, -0.7, 2.0, 5.0, -3.0, 40.0, -25.0):
                rows.append(f'    tfbcllf(bad, {arr(x)}[..], {lmb!r}, {float(stats.boxcox_llf(lmb, xa))!r}, "boxcox_llf {tag} {lmb}");')
            for lmb in (0.0, 0.3, -0.7, 2.0):
                rows.append(f'    tfbc(bad, {arr(x)}[..], {lmb!r}, {arr(stats.boxcox(xa, lmb))}[..], "boxcox {tag} {lmb}");')
            mle = float(stats.boxcox_normmax(xa, method="mle"))
            pr = float(stats.boxcox_normmax(xa, method="pearsonr"))
            rows.append(f'    tfbcmax(bad, {arr(x)}[..], {mle!r}, {pr!r}, "boxcox_normmax {tag}");')
            y, lam = stats.boxcox(xa)
            rows.append(f'    tfbcfit(bad, {arr(x)}[..], {float(lam)!r}, {arr(y)}[..], "boxcox_fit {tag}");')
            for alpha in (0.05, 0.1):
                _, _, ci = stats.boxcox(xa, alpha=alpha)
                rows.append(f'    tfbcci(bad, {arr(x)}[..], {alpha!r}, {float(ci[0])!r}, {float(ci[1])!r}, "boxcox interval {tag} {alpha}");')
            for la, lb, n in ((-1.0, 1.0, 5), (-2.0, 3.0, 17)):
                lm, r = stats.boxcox_normplot(xa, la, lb, N=n)
                pairs = [v for t in zip(lm, r) for v in t]
                rows.append(f'    tfbcplot(bad, {arr(x)}[..], {la!r}, {lb!r}, {n}, {arr(pairs)}[..], "boxcox_normplot {tag} {n}");')
        for lmb in (0.0, 0.4, 2.0, 1.7, -1.0, 3.5, -6.0, 30.0):
            try:
                v = float(stats.yeojohnson_llf(lmb, xa))
            except Exception:
                continue
            if np.isfinite(v):
                rows.append(f'    tfyjllf(bad, {arr(x)}[..], {lmb!r}, {v!r}, "yeojohnson_llf {tag} {lmb}");')
        for lmb in (0.0, 0.5, 2.0, 1.7, -1.0):
            rows.append(f'    tfyj(bad, {arr(x)}[..], {lmb!r}, {arr(stats.yeojohnson(xa, lmb))}[..], "yeojohnson {tag} {lmb}");')
        y, lam = stats.yeojohnson(xa)
        rows.append(f'    tfyjfit(bad, {arr(x)}[..], {float(lam)!r}, {arr(y)}[..], "yeojohnson_fit {tag}");')
        for la, lb, n in ((-1.0, 1.0, 5), (-2.0, 3.0, 17)):
            lm, r = stats.yeojohnson_normplot(xa, la, lb, N=n)
            pairs = [v for t in zip(lm, r) for v in t]
            rows.append(f'    tfyjplot(bad, {arr(x)}[..], {la!r}, {lb!r}, {n}, {arr(pairs)}[..], "yeojohnson_normplot {tag} {n}");')

    head = '''// The answers of SciPy 1.18.1 for the Box-Cox and Yeo-Johnson transforms.
// The generator test/make_transform_vectors.py writes this file.
// It holds the numbers and the checks that use them.
namespace suite;

import std;
import math;
import statistics from "../..";

// True when a equals e within tol, relative to 1 plus the size of e.
fn tfnear(a: f64, e: f64, tol: f64) -> bool {
    return abs(a - e) <= tol * (1.0 + abs(e));
}

fn tfbad(bad: i64&, ok: bool, what: u8[:]) {
    if ok { return; }
    bad += 1;
    if bad <= 12 { print("FAIL ", what, "\\n"); }
}

fn tfbcllf(bad: i64&, xs: const f64[:], l: f64, want: f64, name: u8[:]) {
    let v, c = statistics::boxcox_llf(l, xs);
    if c == 0 && near_tol(v, want, 1e-12) { return; }
    bad += 1;
    if bad <= 12 { print("FAIL ", name, " got ", v, " code ", c, " want ", want, "\\n"); }
}

fn tfyjllf(bad: i64&, xs: const f64[:], l: f64, want: f64, name: u8[:]) {
    let v, c = statistics::yeojohnson_llf(l, xs);
    if c == 0 && near_tol(v, want, 1e-12) { return; }
    bad += 1;
    if bad <= 12 { print("FAIL ", name, " got ", v, " code ", c, " want ", want, "\\n"); }
}

fn tfvalues(got: const f64[:], want: const f64[:], tol: f64) -> bool {
    if got.len != want.len { return false; }
    for i in 0..got.len {
        if !tfnear(got[i], want[i], tol) { return false; }
    }
    return true;
}

fn tfbc(bad: i64&, xs: const f64[:], l: f64, want: const f64[:], name: u8[:]) {
    var out: f64[>..] = [];
    let r, c = statistics::boxcox(out, xs, l);
    tfbad(bad, c == 0 && r == l && tfvalues(out[..], want, 1e-12), name);
}

fn tfyj(bad: i64&, xs: const f64[:], l: f64, want: const f64[:], name: u8[:]) {
    var out: f64[>..] = [];
    let r, c = statistics::yeojohnson(out, xs, l);
    tfbad(bad, c == 0 && r == l && tfvalues(out[..], want, 1e-12), name);
}

// The lambda of the searches agrees with SciPy within the tolerance of the
// Brent search, which is 1.48e-8 relative to the point. The maximum
// likelihood is a sharp optimum. The Pearson fit is a flat optimum, so it
// gets a wider bound.
fn tfbcmax(bad: i64&, xs: const f64[:], mle: f64, pr: f64, name: u8[:]) {
    let a, ca = statistics::boxcox_normmax(xs, statistics::NORMMAX_MLE);
    let b, cb = statistics::boxcox_normmax(xs);
    let b2, cb2 = statistics::boxcox_normmax(xs, statistics::NORMMAX_PEARSONR);
    var ok = ca == 0 && cb == 0 && cb2 == 0 && b == b2;
    if !(tfnear(a, mle, 1e-6) && tfnear(b, pr, 1e-5)) { ok = false; }
    if ok { return; }
    bad += 1;
    if bad <= 12 { print("FAIL ", name, " got ", a, " ", b, " codes ", ca, " ", cb, " want ", mle, " ", pr, "\\n"); }
}

fn tfbcfit(bad: i64&, xs: const f64[:], lam: f64, want: const f64[:], name: u8[:]) {
    var out: f64[>..] = [];
    var again: f64[>..] = [];
    let l, c = statistics::boxcox_fit(out, xs);
    let _r, c2 = statistics::boxcox(again, xs, l);
    tfbad(bad, c == 0 && c2 == 0 && tfnear(l, lam, 1e-6) && tfvalues(out[..], again[..], 0.0)
        && tfvalues(out[..], want, 1e-5), name);
}

fn tfyjfit(bad: i64&, xs: const f64[:], lam: f64, want: const f64[:], name: u8[:]) {
    var out: f64[>..] = [];
    var again: f64[>..] = [];
    let l, c = statistics::yeojohnson_fit(out, xs);
    let m, cm = statistics::yeojohnson_normmax(xs);
    let _r, c2 = statistics::yeojohnson(again, xs, l);
    tfbad(bad, c == 0 && c2 == 0 && cm == 0 && m == l && tfnear(l, lam, 1e-6)
        && tfvalues(out[..], again[..], 0.0) && tfvalues(out[..], want, 1e-5), name);
}

fn tfbcci(bad: i64&, xs: const f64[:], alpha: f64, lo: f64, hi: f64, name: u8[:]) {
    let a, b, m, c = statistics::boxcox_interval(xs, alpha);
    tfbad(bad, c == 0 && tfnear(a, lo, 1e-6) && tfnear(b, hi, 1e-6) && a < m && m < b, name);
}

fn tfbcplot(bad: i64&, xs: const f64[:], la: f64, lb: f64, n: i64, want: const f64[:], name: u8[:]) {
    var out: f64[>..] = [];
    let k, c = statistics::boxcox_normplot(out, xs, la, lb, n);
    tfbad(bad, c == 0 && k == n && tfvalues(out[..], want, 1e-10), name);
}

fn tfyjplot(bad: i64&, xs: const f64[:], la: f64, lb: f64, n: i64, want: const f64[:], name: u8[:]) {
    var out: f64[>..] = [];
    let k, c = statistics::yeojohnson_normplot(out, xs, la, lb, n);
    tfbad(bad, c == 0 && k == n && tfvalues(out[..], want, 1e-10), name);
}

fn run_transform_vectors() -> i64 {
    var bad: i64 = 0;
'''
    tail = "    return bad;\n}\n"
    path = "test/suite/transform_vectors.goose"
    with open(path, "w") as f:
        f.write(head + "\n".join(rows) + "\n" + tail)
    print(len(rows), "rows")


if __name__ == "__main__":
    main()
