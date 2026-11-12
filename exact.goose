// The exact and distance tests: barnard_exact, boschloo_exact,
// cramervonmises, cramervonmises_2samp, epps_singleton_2samp, bws_test,
// anderson_ksamp, kstest, kstest_2samp, poisson_means_test and
// page_trend_test.
//
// Every call returns a TestResult and then a fault code, like the tests of
// tests.goose. Several samples sit in one buffer and a list of sizes gives
// the length of each sample. The helpers of this file start with ex_.
namespace statistics;

import std;
import math;

// The variants of the k-sample Anderson-Darling statistic.
const KSAMP_MIDRANK: i64 = 0;      // ties share the mean of their positions
const KSAMP_RIGHT: i64 = 1;        // the right empirical distribution
const KSAMP_CONTINUOUS: i64 = 2;   // no ties: the shortest sum

// The method of the Page trend test.
const PAGE_AUTO: i64 = 0;          // the exact form for a small table
const PAGE_EXACT: i64 = 1;         // the exact form, or a fault
const PAGE_ASYMPTOTIC: i64 = 2;    // the normal form

const EX_PI: f64 = 3.141592653589793;
// The cells of the double sum of the Poisson means test.
const EX_POISSON_CELLS: i64 = 25000000;
// The most points of the characteristic function of the Epps-Singleton test.
const EX_MAX_POINTS: i64 = 16;

// ---- shared helpers ---------------------------------------------------


// Row n of the triangle of Pascal as doubles. A row up to 1000 stays finite.
fn ex_binom_row(out: f64[>..]&, n: i64) {
    out.resize(0);
    out.push(1.0);
    for r in 1..n + 1 {
        out.push(1.0);
        var k = r - 1;
        while k >= 1 {
            out[k] = out[k] + out[k - 1];
            k -= 1;
        }
    }
}

// The binomial coefficient as a double, from the product form.
fn ex_choose(n: i64, k: i64) -> f64 {
    if k < 0 || k > n { return 0.0; }
    let kk = min(k, n - k);
    var r: f64 = 1.0;
    for i in 1..kk + 1 { r = r * (n - kk + i) / i; }
    return r;
}

// ---- the nuisance parameter of the 2 by 2 tests -----------------------

// The log of the sum of a[s] pi^s (1 - pi)^(n - s), where la[s] holds the
// log of a[s] and -infinity for a zero. The end points use one term.
fn ex_logp(la: const f64[:], n: i64, pi: f64) -> f64 {
    if pi <= 0.0 { return la[0]; }
    if pi >= 1.0 { return la[n]; }
    let lp = log(pi);
    let lq = log(1.0 - pi);
    var mx: f64 = -inf();
    for s in 0..n + 1 {
        if la[s] > -inf() {
            let t = la[s] + s * lp + (n - s) * lq;
            if t > mx { mx = t; }
        }
    }
    if mx == -inf() { return mx; }
    var tot = acc_new();
    for s in 0..n + 1 {
        if la[s] > -inf() {
            acc_add(tot, exp(la[s] + s * lp + (n - s) * lq - mx));
        }
    }
    return mx + log(acc_total(tot));
}

// The largest value over pi from 0 to 1 of the polynomial with the
// coefficients a[0..n]. A grid of grid + 1 points finds each local
// maximum, and Brent refines the maxima near the best one.
fn ex_max_nuisance(a: const f64[:], n: i64, grid: i64) -> f64, i64 {
    var la: f64[>..] = [];
    for s in 0..n + 1 {
        if a[s] > 0.0 { la.push(log(a[s])); } else { la.push(-inf()); }
    }
    var fv: f64[>..] = [];
    var best: f64 = -inf();
    for i in 0..grid + 1 {
        let v = ex_logp(la[..], n, i / (grid as f64));
        fv.push(v);
        if v > best { best = v; }
    }
    if best == -inf() { return 0.0, OK; }
    var refined: i64 = 0;
    for i in 0..grid + 1 {
        if refined >= 64 { break; }
        if fv[i] < best - 0.05 { continue; }
        if i > 0 && fv[i] < fv[i - 1] { continue; }
        if i < grid && fv[i] < fv[i + 1] { continue; }
        refined += 1;
        let lo = max(0.0, (i - 1) / (grid as f64));
        let hi = min(1.0, (i + 1) / (grid as f64));
        var x0 = i / (grid as f64);
        if i == 0 || i == grid { x0 = 0.5 * (lo + hi); }
        let xm, cm = brent_min(lo, x0, hi, 1e-13) { p => -ex_logp(la[..], n, p) };
        if cm == OK {
            let v = ex_logp(la[..], n, xm);
            if v > best { best = v; }
        }
    }
    return exp(best), OK;
}

// The checks that the two 2 by 2 tests share. The counts come first, then
// the alternative, the grid and the total of the table.
fn ex_check_table(a: f64, b: f64, c: f64, d: f64, alt: i64, grid: i64) -> i64 {
    var code = chk_whole(a);
    if code != OK { return code; }
    code = chk_whole(b);
    if code != OK { return code; }
    code = chk_whole(c);
    if code != OK { return code; }
    code = chk_whole(d);
    if code != OK { return code; }
    if check_alt(alt) != OK { return ERR_DOMAIN; }
    if grid < 1 { return ERR_DOMAIN; }
    if grid > MAX_NUISANCE_GRID { return ERR_TOO_LONG; }
    if a + b + c + d > MAX_EXACT_TABLE { return ERR_TOO_LONG; }
    if a + c == 0.0 || b + d == 0.0 { return ERR_UNDEFINED; }
    return OK;
}

// The Wald statistic of the table with x1 in the first column and x2 in the
// second. It is 0 where the two proportions agree.
fn ex_wald(x1: i64, x2: i64, c1: i64, c2: i64, pooled: bool) -> f64 {
    let f1 = c1 as f64;
    let f2 = c2 as f64;
    let p1 = x1 / f1;
    let p2 = x2 / f2;
    if p1 == p2 { return 0.0; }
    var v: f64 = 0.0;
    if pooled {
        let p = (x1 + x2) / (f1 + f2);
        v = p * (1.0 - p) * (1.0 / f1 + 1.0 / f2);
    } else {
        v = p1 * (1.0 - p1) / f1 + p2 * (1.0 - p2) / f2;
    }
    return (p1 - p2) / sqrt(v);
}

// Signature: barnard_exact(a, b, c, d, alt = ALT_TWO_SIDED, pooled = true,
// grid = 1024) -> TestResult, code. The table is [[a, b], [c, d]]. The
// p-value is the largest value over the nuisance parameter. A column with
// no count, or an infinite statistic, gives ERR_UNDEFINED.
fn barnard_exact(a: f64, b: f64, c: f64, d: f64, alt: i64 = ALT_TWO_SIDED, pooled: bool = true, grid: i64 = 1024) -> TestResult, i64 {
    let code = ex_check_table(a, b, c, d, alt, grid);
    if code != OK { return fault_result(code); }
    let c1 = (a + c) as i64;
    let c2 = (b + d) as i64;
    let n = c1 + c2;
    let t0 = ex_wald(a as i64, b as i64, c1, c2, pooled);
    if is_bad(t0) { return fault_result(ERR_UNDEFINED); }
    var r1: f64[>..] = [];
    var r2: f64[>..] = [];
    ex_binom_row(r1, c1);
    ex_binom_row(r2, c2);
    var sums: f64[>..] = [];
    for s in 0..n + 1 { sums.push(0.0); }
    for x1 in 0..c1 + 1 {
        for x2 in 0..c2 + 1 {
            let t = ex_wald(x1, x2, c1, c2, pooled);
            var inside = false;
            if alt == ALT_TWO_SIDED { inside = abs(t) >= abs(t0); }
            else if alt == ALT_LESS { inside = t <= t0; }
            else { inside = t >= t0; }
            if inside { sums[x1 + x2] += r1[x1] * r2[x2]; }
        }
    }
    let p, cp = ex_max_nuisance(sums[..], n, grid);
    if cp != OK { return fault_result(cp); }
    return good_result(t0, p, nan(), nan());
}

// The one-sided Fisher p-values of one column count s: pv[i - lo] for the
// first column count i from lo = max(0, s - c2) up. The lower tail of
// the first column (less), or of the second column (greater).
fn ex_fisher_row(pv: f64[>..]&, r1: const f64[:], r2: const f64[:], rt: const f64[:], c1: i64, c2: i64, s: i64, less: bool) {
    pv.resize(0);
    let lo = max(0, s - c2);
    let hi = min(c1, s);
    for i in lo..hi + 1 { pv.push(0.0); }
    if less {
        var run: f64 = 0.0;
        for i in lo..hi + 1 {
            run += r1[i] * r2[s - i] / rt[s];
            pv[i - lo] = run;
        }
    } else {
        var run: f64 = 0.0;
        var i = hi;
        while i >= lo {
            run += r1[i] * r2[s - i] / rt[s];
            pv[i - lo] = run;
            i -= 1;
        }
    }
}

// One side of Boschloo's test: the statistic and the p-value.
fn ex_boschloo_side(a: i64, b: i64, c: i64, d: i64, less: bool, grid: i64) -> f64, f64, i64 {
    let c1 = a + c;
    let c2 = b + d;
    let n = c1 + c2;
    var r1: f64[>..] = [];
    var r2: f64[>..] = [];
    var rt: f64[>..] = [];
    ex_binom_row(r1, c1);
    ex_binom_row(r2, c2);
    ex_binom_row(rt, n);
    var pv: f64[>..] = [];
    ex_fisher_row(pv, r1[..], r2[..], rt[..], c1, c2, a + b, less);
    let stat = pv[a - max(0, a + b - c2)];
    let limit = stat * (1.0 + 1e-13);
    var sums: f64[>..] = [];
    for s in 0..n + 1 { sums.push(0.0); }
    for s in 0..n + 1 {
        ex_fisher_row(pv, r1[..], r2[..], rt[..], c1, c2, s, less);
        let lo = max(0, s - c2);
        for k in 0..pv.len {
            if pv[k] <= limit { sums[s] += r1[lo + k] * r2[s - lo - k]; }
        }
    }
    let p, cp = ex_max_nuisance(sums[..], n, grid);
    return stat, p, cp;
}

// Signature: boschloo_exact(a, b, c, d, alt = ALT_TWO_SIDED, grid = 1024)
// -> TestResult, code. The statistic is the one-sided Fisher p-value of the
// table. The two-sided p-value is twice the smaller one-sided p-value.
fn boschloo_exact(a: f64, b: f64, c: f64, d: f64, alt: i64 = ALT_TWO_SIDED, grid: i64 = 1024) -> TestResult, i64 {
    let code = ex_check_table(a, b, c, d, alt, grid);
    if code != OK { return fault_result(code); }
    let ia = a as i64;
    let ib = b as i64;
    let ic = c as i64;
    let id = d as i64;
    if alt == ALT_LESS || alt == ALT_GREATER {
        let st, p, cp = ex_boschloo_side(ia, ib, ic, id, alt == ALT_LESS, grid);
        if cp != OK { return fault_result(cp); }
        return good_result(st, p, nan(), nan());
    }
    let sl, pl, cl = ex_boschloo_side(ia, ib, ic, id, true, grid);
    if cl != OK { return fault_result(cl); }
    let sg, pg, cg = ex_boschloo_side(ia, ib, ic, id, false, grid);
    if cg != OK { return fault_result(cg); }
    if pl < pg { return good_result(sl, 2.0 * pl, nan(), nan()); }
    return good_result(sg, 2.0 * pg, nan(), nan());
}

// ---- the Kolmogorov-Smirnov dispatchers -------------------------------

// Signature: kstest(out, xs, family, params, alt = ALT_TWO_SIDED)
// -> TestResult, code. The result is that of ks_1samp. The list out holds
// the location of the statistic and its sign, +1 or -1. A tie of the two
// signs gives -1, as SciPy does.
fn kstest(out: f64[>..]&, xs: const f64[:], family: i64, params: const f64[:], alt: i64 = ALT_TWO_SIDED) -> TestResult, i64 {
    out.resize(0);
    let dplus, dminus, lplus, lminus, c = ks_1samp_sweep(xs, family, params, alt);
    if c != OK { return fault_result(c); }
    let r, cr = ks_1samp(xs, family, params, alt);
    if cr != OK { return fault_result(cr); }
    var plus = dplus > dminus;
    if alt == ALT_GREATER { plus = true; }
    if alt == ALT_LESS { plus = false; }
    if plus { out.push(lplus); out.push(1.0); } else { out.push(lminus); out.push(-1.0); }
    return r, OK;
}

// The value of the empirical function at count j of n values, as numpy
// linspace makes it. SciPy picks the location from these values.
fn ex_lin(j: i64, n: i64) -> f64 {
    if j == n { return 1.0; }
    return j * (1.0 / n);
}

// Signature: kstest_2samp(out, xs, ys, alt = ALT_TWO_SIDED) -> TestResult,
// code. The result is that of ks_2samp. The list out holds the location of
// the statistic and its sign.
fn kstest_2samp(out: f64[>..]&, xs: const f64[:], ys: const f64[:], alt: i64 = ALT_TWO_SIDED) -> TestResult, i64 {
    out.resize(0);
    if check_alt(alt) != OK { return fault_result(ERR_DOMAIN); }
    let r, cr = ks_2samp(xs, ys, alt);
    if cr != OK { return fault_result(cr); }
    let sx = sort_copy(xs);
    let sy = sort_copy(ys);
    let n1 = sx.len;
    let n2 = sy.len;
    var dmax: f64 = -inf();
    var dmin: f64 = inf();
    var lmax: f64 = 0.0;
    var lmin: f64 = 0.0;
    var i: i64 = 0;
    var j: i64 = 0;
    while i < n1 || j < n2 {
        var v: f64 = 0.0;
        if i < n1 && j < n2 {
            v = if sx[i] <= sy[j] { sx[i] } else { sy[j] };
        } else if i < n1 {
            v = sx[i];
        } else {
            v = sy[j];
        }
        while i < n1 && sx[i] == v { i += 1; }
        while j < n2 && sy[j] == v { j += 1; }
        let diff = ex_lin(i, n1) - ex_lin(j, n2);
        if diff > dmax { dmax = diff; lmax = v; }
        if diff < dmin { dmin = diff; lmin = v; }
    }
    let minus_side = clip01(-dmin);
    var minus = minus_side > dmax;
    if alt == ALT_GREATER { minus = false; }
    if alt == ALT_LESS { minus = true; }
    if minus { out.push(lmin); out.push(-1.0); } else { out.push(lmax); out.push(1.0); }
    return r, OK;
}

// ---- the Poisson means test -------------------------------------------

// Signature: poisson_means_test(k1, n1, k2, n2, diff = 0.0,
// alt = ALT_TWO_SIDED) -> TestResult, code. The counts k1 and k2 are whole
// numbers. The sizes n1 and n2 are above 0. The sum of the Poisson mass runs
// over the window from the quantile 1e-10 to the quantile 1 - 1e-16.
fn poisson_means_test(k1: f64, n1: f64, k2: f64, n2: f64, diff: f64 = 0.0, alt: i64 = ALT_TWO_SIDED) -> TestResult, i64 {
    var code = chk_whole(k1);
    if code != OK { return fault_result(code); }
    code = check_shape(n1);
    if code != OK { return fault_result(code); }
    code = chk_whole(k2);
    if code != OK { return fault_result(code); }
    code = check_shape(n2);
    if code != OK { return fault_result(code); }
    code = check_scalar(diff);
    if code != OK { return fault_result(code); }
    if diff < 0.0 || check_alt(alt) != OK { return fault_result(ERR_DOMAIN); }
    let lam2 = (k1 + k2) / (n1 + n2) - diff * n1 / (n1 + n2);
    if lam2 <= 0.0 { return good_result(0.0, 1.0, nan(), nan()); }
    let v0 = k1 / (n1 * n1) + k2 / (n2 * n2);
    let t0 = (k1 / n1 - k2 / n2 - diff) / sqrt(v0);
    if is_bad(t0) { return fault_result(ERR_UNDEFINED); }
    let mu1 = n1 * (lam2 + diff);
    let mu2 = n2 * lam2;
    let lo1, c1 = poisson_ppf(1e-10, mu1);
    let hi1, c2 = poisson_ppf(1.0 - 1e-16, mu1);
    let lo2, c3 = poisson_ppf(1e-10, mu2);
    let hi2, c4 = poisson_ppf(1.0 - 1e-16, mu2);
    if c1 != OK { return fault_result(c1); }
    if c2 != OK { return fault_result(c2); }
    if c3 != OK { return fault_result(c3); }
    if c4 != OK { return fault_result(c4); }
    let w1 = hi1 - lo1 + 1.0;
    let w2 = hi2 - lo2 + 1.0;
    if w1 * w2 > EX_POISSON_CELLS { return fault_result(ERR_TOO_LONG); }
    var p1: f64[>..] = [];
    var p2: f64[>..] = [];
    for k in 0..w1 as i64 {
        let m, cm = poisson_pmf(lo1 + k, mu1);
        if cm != OK { return fault_result(cm); }
        p1.push(m);
    }
    for k in 0..w2 as i64 {
        let m, cm = poisson_pmf(lo2 + k, mu2);
        if cm != OK { return fault_result(cm); }
        p2.push(m);
    }
    var tot = acc_new();
    for i in 0..p1.len {
        let x1 = lo1 + i;
        let l1 = x1 / n1;
        for j in 0..p2.len {
            let x2 = lo2 + j;
            let l2 = x2 / n2;
            let t = (l1 - l2 - diff) / sqrt(l1 / n1 + l2 / n2);
            var inside = false;
            if alt == ALT_TWO_SIDED { inside = abs(t) >= abs(t0); }
            else if alt == ALT_LESS { inside = t <= t0; }
            else { inside = t >= t0; }
            if inside { acc_add(tot, p1[i] * p2[j]); }
        }
    }
    return good_result(t0, acc_total(tot), nan(), nan());
}

// ---- the modified Bessel function K -----------------------------------

const EX_BESSEL_NODES: i64 = 24;     // nodes of the Gauss-Legendre rule
const EX_BESSEL_PANEL: f64 = 0.5;    // width of one panel of the integral
const EX_CVM_TERMS: i64 = 4000;      // terms of the series of the Cramer-von Mises law

// e^x K_v(x) for x above 0, by the integral of exp(-x (cosh t - 1))
// cosh(v t) over t from 0 up. The integral runs over panels, and the
// nodes gx and gw are the Gauss-Legendre rule. The upper end T is the point
// where the integrand has fallen below e^-50. The first term of the
// tail beyond T is added.
fn ex_bessel_ke(v: f64, x: f64, gx: const f64[:], gw: const f64[:]) -> f64 {
    var top: f64 = 0.5;
    var steps: i64 = 0;
    while x * (ex_cosh(top) - 1.0) - v * top < 50.0 && steps < 1400 {
        top += 0.5;
        steps += 1;
    }
    let panels = ceil(top / EX_BESSEL_PANEL) as i64;
    let h = top / panels;
    var s: f64 = 0.0;
    for p in 0..panels {
        let a = p * h;
        var part: f64 = 0.0;
        for i in 0..gx.len {
            let t = a + (gx[i] + 1.0) * h / 2.0;
            part += gw[i] * exp(-x * (ex_cosh(t) - 1.0)) * ex_cosh(v * t);
        }
        s += part * h / 2.0;
    }
    let ft = exp(-x * (ex_cosh(top) - 1.0)) * ex_cosh(v * top);
    let d = x * ex_sinh(top) - v;
    if d > 0.0 { s += ft / d; }
    return s;
}

fn ex_cosh(t: f64) -> f64 { return 0.5 * (exp(t) + exp(-t)); }
fn ex_sinh(t: f64) -> f64 { return 0.5 * (exp(t) - exp(-t)); }

// The function K of order v at x above 0. The scaled form above holds the
// work. A value that underflows gives 0.
fn ex_bessel_k(v: f64, x: f64) -> f64 {
    var gx: f64[>..] = [];
    var gw: f64[>..] = [];
    gauss_legendre(gx, gw, EX_BESSEL_NODES);
    return exp(-x) * ex_bessel_ke(v, x, gx[..], gw[..]);
}

// ---- the Cramer-von Mises law -----------------------------------------

// The terms of eq. 1.2 of Csorgo and Faraway: the limit law of the statistic
// for x above 0. The sum stops at a term below 1e-7, as SciPy does.
fn ex_cvm_inf(x: f64, gx: const f64[:], gw: const f64[:]) -> f64, i64 {
    var tot: f64 = 0.0;
    var ratio: f64 = sqrt(EX_PI);          // Gamma(k + 1/2) / k!
    for k in 0..EX_CVM_TERMS {
        let y = 4.0 * k + 1.0;
        let q = y * y / (16.0 * x);
        let u = ratio / (EX_PI * sqrt(EX_PI) * sqrt(x));
        let term = u * sqrt(y) * exp(-2.0 * q) * ex_bessel_ke(0.25, q, gx, gw);
        tot += term;
        if abs(term) < 1e-7 { return tot, OK; }
        ratio = ratio * (k + 0.5) / (k + 1.0);
    }
    return nan(), ERR_FIT;
}

// The helper functions of the correction of order 1/n. Each one holds
// the factor exp(-2 z) that the scaled function Bessel gives.
fn ex_ed2(y: f64, gx: const f64[:], gw: const f64[:]) -> f64 {
    let z = y * y / 4.0;
    let b = ex_bessel_ke(0.25, z, gx, gw) + ex_bessel_ke(0.75, z, gx, gw);
    return exp(-2.0 * z) * pow(y / 2.0, 1.5) * b / sqrt(EX_PI);
}

fn ex_ed3(y: f64, gx: const f64[:], gw: const f64[:]) -> f64 {
    let z = y * y / 4.0;
    let b = 2.0 * ex_bessel_ke(0.25, z, gx, gw) + 3.0 * ex_bessel_ke(0.75, z, gx, gw) - ex_bessel_ke(1.25, z, gx, gw);
    return exp(-2.0 * z) * pow(y / 2.0, 2.5) * b / sqrt(EX_PI);
}

// The correction term of eq. 1.10 without the part V(x) / 12. The sum
// stops at a term below 1e-7.
fn ex_psi1_mod(x: f64, gx: const f64[:], gw: const f64[:]) -> f64, i64 {
    var tot: f64 = 0.0;
    var ratio: f64 = sqrt(EX_PI);          // Gamma(k + 1/2) / k!
    let sx = 2.0 * sqrt(x);
    let y1 = pow(x, 0.75);
    let y2 = pow(x, 1.25);
    for k in 0..EX_CVM_TERMS {
        let m = 2.0 * k + 1.0;
        let g1 = ratio;                     // Gamma(k + 1/2) / k!
        let g2 = ratio * (k + 0.5);         // Gamma(k + 3/2) / k!
        let e1 = m * g1 * ex_ed2((4.0 * k + 3.0) / sx, gx, gw) / (9.0 * y1);
        let e2 = g1 * ex_ed3((4.0 * k + 1.0) / sx, gx, gw) / (72.0 * y2);
        let e3 = 2.0 * (m + 2.0) * g2 * ex_ed3((4.0 * k + 5.0) / sx, gx, gw) / (12.0 * y2);
        let e4 = 7.0 * m * g1 * ex_ed2((4.0 * k + 1.0) / sx, gx, gw) / (144.0 * y1);
        let e5 = 7.0 * m * g1 * ex_ed2((4.0 * k + 5.0) / sx, gx, gw) / (144.0 * y1);
        let z = -(e1 + e2 + e3 + e4 + e5) / EX_PI;
        tot += z;
        if abs(z) < 1e-7 { return tot, OK; }
        ratio = ratio * (k + 0.5) / (k + 1.0);
    }
    return nan(), ERR_FIT;
}

// The cumulative function of the statistic for n values, by the finite
// form of Csorgo and Faraway. The form is poor for x above 2 or 4.
fn ex_cvm_cdf(x: f64, n: i64, gx: const f64[:], gw: const f64[:]) -> f64, i64 {
    let nf = n as f64;
    if x >= nf / 3.0 { return 1.0, OK; }
    if !(1.0 / (12.0 * nf) < x) { return 0.0, OK; }
    let v, c1 = ex_cvm_inf(x, gx, gw);
    if c1 != OK { return nan(), c1; }
    let ps, c2 = ex_psi1_mod(x, gx, gw);
    if c2 != OK { return nan(), c2; }
    return v * (1.0 + 1.0 / (12.0 * nf)) + ps / nf, OK;
}

// Signature: cramervonmises(xs, family, params) -> TestResult, code.
// The family is one of the FAMILY_ constants and params holds its
// parameters, as in ks_1samp. The statistic is n times the Cramer-von Mises
// criterion. The p-value is one minus the finite law, cut at 0.
fn cramervonmises(xs: const f64[:], family: i64, params: const f64[:]) -> TestResult, i64 {
    var c = check_sample(xs, 2);
    if c != OK { return fault_result(c); }
    if family < 0 || family > FAMILY_F { return fault_result(ERR_DOMAIN); }
    if params.len != family_arity(family) { return fault_result(ERR_LENGTH); }
    for p in params {
        c = check_scalar(p);
        if c != OK { return fault_result(c); }
    }
    let s = sort_copy(xs);
    let n = s.len;
    let nf = n as f64;
    var tot = acc_new();
    for k in 0..n {
        let f, cf = family_cdf(family, s[k], params);
        if cf != OK { return fault_result(cf); }
        let u = (2.0 * (k + 1) - 1.0) / (2.0 * nf);
        acc_add(tot, (u - f) * (u - f));
    }
    let w = 1.0 / (12.0 * nf) + acc_total(tot);
    var gx: f64[>..] = [];
    var gw: f64[>..] = [];
    gauss_legendre(gx, gw, EX_BESSEL_NODES);
    let cdf, cc = ex_cvm_cdf(w, n, gx[..], gw[..]);
    if cc != OK { return fault_result(cc); }
    return good_result(w, max(1.0 - cdf, 0.0), nan(), nan());
}

// ---- the two-sample Cramer-von Mises test -----------------------------

// Merge two sorted lists of (value, count) pairs. A value in both lists
// adds its counts.
fn ex_merge(ov: i64[>..]&, of: i64[>..]&, av: const i64[:], af: const i64[:], bv: const i64[:], bf: const i64[:]) {
    ov.resize(0);
    of.resize(0);
    var i: i64 = 0;
    var j: i64 = 0;
    while i < av.len || j < bv.len {
        if j >= bv.len || (i < av.len && av[i] < bv[j]) {
            ov.push(av[i]); of.push(af[i]); i += 1;
        } else if i >= av.len || bv[j] < av[i] {
            ov.push(bv[j]); of.push(bf[j]); j += 1;
        } else {
            ov.push(av[i]); of.push(af[i] + bf[j]); i += 1; j += 1;
        }
    }
}

// The count of the assignments of m and n values with the Anderson sum at
// or above zeta, over the count of all assignments. The table follows the
// lattice paths of Xiao, Gordon and Yakovlev. The code is ERR_TOO_LONG
// when a list passes MAX_EXACT_CELLS or the work passes 5e8 steps.
fn ex_cvm_exact(m: i64, n: i64, zeta: i64) -> f64, i64 {
    let g = gcd_counts(m, n);
    let lcm = m / g * n;
    let a = lcm / m;
    let b = lcm / n;
    // The layer before u = 0 holds one list with the value 0, counted once.
    // Each later list of that layer is empty. ps holds the offsets.
    var pv: i64[>..] = [0];
    var pf: i64[>..] = [1];
    var ps: i64[>..] = [0];
    for v in 0..m + 1 { ps.push(1); }
    var work: i64 = 0;
    var tv: i64[>..] = [];
    var tf: i64[>..] = [];
    var mv: i64[>..] = [];
    var mf: i64[>..] = [];
    for u in 0..n + 1 {
        var nv: i64[>..] = [];
        var nf: i64[>..] = [];
        var ns: i64[>..] = [0];
        tv.resize(0);
        tf.resize(0);
        for v in 0..m + 1 {
            let lo = ps[v];
            let hi = ps[v + 1];
            ex_merge(mv, mf, tv[..], tf[..], pv[lo..hi], pf[lo..hi]);
            if mv.len > MAX_EXACT_CELLS { return nan(), ERR_TOO_LONG; }
            work += mv.len;
            if work > 500000000 { return nan(), ERR_TOO_LONG; }
            tv.resize(0);
            tf.resize(0);
            let res = (a * v - b * u) * (a * v - b * u);
            for k in 0..mv.len { tv.push(mv[k] + res); tf.push(mf[k]); }
            for k in 0..tv.len { nv.push(tv[k]); nf.push(tf[k]); }
            ns.push(nv.len);
        }
        pv.resize(0);
        pf.resize(0);
        for k in 0..nv.len { pv.push(nv[k]); pf.push(nf[k]); }
        ps.resize(0);
        for k in 0..ns.len { ps.push(ns[k]); }
    }
    let lo = ps[m];
    let hi = ps[m + 1];
    var above: i64 = 0;
    var total: i64 = 0;
    for k in lo..hi {
        total += pf[k];
        if pv[k] >= zeta { above += pf[k]; }
    }
    if total == 0 { return nan(), ERR_UNDEFINED; }
    return above / (total as f64), OK;
}

// The floor of num over den for a positive den.
fn ex_floor_div(num: i64, den: i64) -> i64 {
    if num >= 0 { return num / den; }
    return -((-num + den - 1) / den);
}

// Signature: cramervonmises_2samp(xs, ys, method = MW_AUTO) -> TestResult,
// code. MW_AUTO uses the exact form when both samples hold 20 values or
// fewer. MW_EXACT forces it up to 100 values. MW_ASYMPTOTIC uses the limit
// law. The exact form uses midranks, so a tie gives the value of SciPy.
fn cramervonmises_2samp(xs: const f64[:], ys: const f64[:], method: i64 = MW_AUTO) -> TestResult, i64 {
    var c = check_two(xs, ys, 2);
    if c != OK { return fault_result(c); }
    if method < MW_AUTO || method > MW_ASYMPTOTIC { return fault_result(ERR_DOMAIN); }
    let nx = xs.len;
    let ny = ys.len;
    let sx = sort_copy(xs);
    let sy = sort_copy(ys);
    var z: f64[>..] = [];
    for x in sx { z.push(x); }
    for y in sy { z.push(y); }
    var r: f64[>..] = [];
    let cnt, cr = rankdata(r, z[..], TIE_AVERAGE);
    if cr != OK { return fault_result(cr); }
    var sumx = acc_new();
    var sumy = acc_new();
    for i in 0..nx {
        let d = 2.0 * r[i] - 2.0 * (i + 1);
        acc_add(sumx, d * d);
    }
    for j in 0..ny {
        let d = 2.0 * r[nx + j] - 2.0 * (j + 1);
        acc_add(sumy, d * d);
    }
    // four times U, as a double
    let u4 = nx * acc_total(sumx) + ny * acc_total(sumy);
    let k = nx as f64 * ny;
    let nn = (nx + ny) as f64;
    let t = u4 / 4.0 / (k * nn) - (4.0 * k - 1.0) / (6.0 * nn);
    var use_exact = method == MW_EXACT;
    if method == MW_AUTO { use_exact = max(nx, ny) <= 20; }
    var gx: f64[>..] = [];
    var gw: f64[>..] = [];
    if use_exact {
        if max(nx, ny) > 100 { return fault_result(ERR_TOO_LONG); }
        let g = gcd_counts(nx, ny);
        let lcm = (nx / g * ny) as f64;
        let mn = k;
        let d = 3.0 * u4 - 2.0 * mn * (4.0 * mn - 1.0);
        let num = lcm * lcm * nn * d;
        if abs(num) > 9.0e18 { return fault_result(ERR_TOO_LONG); }
        let zeta = ex_floor_div(num as i64, 12 * nx * ny * nx * ny);
        let p, cp = ex_cvm_exact(nx, ny, zeta);
        if cp != OK { return fault_result(cp); }
        return good_result(t, p, nan(), nan());
    }
    let et = (1.0 + 1.0 / nn) / 6.0;
    var vt = (nn + 1.0) * (4.0 * k * nn - 3.0 * (nx * nx + ny * ny) - 2.0 * k);
    vt = vt / (45.0 * nn * nn * 4.0 * k);
    let tn = 1.0 / 6.0 + (t - et) / sqrt(45.0 * vt);
    if tn < 0.003 { return good_result(t, 1.0, nan(), nan()); }
    gauss_legendre(gx, gw, EX_BESSEL_NODES);
    let v, cv = ex_cvm_inf(tn, gx[..], gw[..]);
    if cv != OK { return fault_result(cv); }
    return good_result(t, max(1.0 - v, 0.0), nan(), nan());
}

// ---- the Epps-Singleton test ------------------------------------------

// The eigenvalues and the eigenvectors of a symmetric matrix of side p by
// the method of Jacobi. The matrix a is row by row and it ends with the
// eigenvalues on its diagonal. The columns of v are the eigenvectors. The
// sweeps are bounded at 100.
fn ex_jacobi(a: f64[>..]&, v: f64[>..]&, p: i64) {
    v.resize(0);
    for i in 0..p * p { v.push(0.0); }
    for i in 0..p { v[i * p + i] = 1.0; }
    for sweep in 0..100 {
        var off: f64 = 0.0;
        var dia: f64 = 0.0;
        for i in 0..p {
            dia += a[i * p + i] * a[i * p + i];
            for j in i + 1..p { off += a[i * p + j] * a[i * p + j]; }
        }
        if off <= 1e-34 * dia || off == 0.0 { return; }
        for i in 0..p {
            for j in i + 1..p {
                let apq = a[i * p + j];
                if apq == 0.0 { continue; }
                let theta = (a[j * p + j] - a[i * p + i]) / (2.0 * apq);
                var t: f64 = 1.0 / (abs(theta) + sqrt(theta * theta + 1.0));
                if theta < 0.0 { t = -t; }
                let c = 1.0 / sqrt(t * t + 1.0);
                let s = t * c;
                for k in 0..p {
                    let x = a[k * p + i];
                    let y = a[k * p + j];
                    a[k * p + i] = c * x - s * y;
                    a[k * p + j] = s * x + c * y;
                }
                for k in 0..p {
                    let x = a[i * p + k];
                    let y = a[j * p + k];
                    a[i * p + k] = c * x - s * y;
                    a[j * p + k] = s * x + c * y;
                }
                for k in 0..p {
                    let x = v[k * p + i];
                    let y = v[k * p + j];
                    v[k * p + i] = c * x - s * y;
                    v[k * p + j] = s * x + c * y;
                }
            }
        }
    }
}

// The biased covariance of the 2 q rows of the characteristic function of
// one sample: the rows cos(ts_j x) first, then the rows sin(ts_j x). The
// result is added into cov with the weight w. The list mean gets the mean
// of each row.
fn ex_es_cov(cov: f64[>..]&, mean: f64[>..]&, s: const f64[:], ts: const f64[:], w: f64) {
    let q = ts.len;
    let p = 2 * q;
    let n = s.len;
    var g: f64[>..] = [];
    for r in 0..p {
        for i in 0..n {
            if r < q { g.push(cos(ts[r] * s[i])); } else { g.push(sin(ts[r - q] * s[i])); }
        }
    }
    mean.resize(0);
    for r in 0..p {
        var tot = acc_new();
        for i in 0..n { acc_add(tot, g[r * n + i]); }
        mean.push(acc_total(tot) / n);
    }
    for r in 0..p {
        for c in r..p {
            var tot = acc_new();
            for i in 0..n { acc_add(tot, (g[r * n + i] - mean[r]) * (g[c * n + i] - mean[c])); }
            let v = acc_total(tot) / n * w;
            cov[r * p + c] += v;
            if c != r { cov[c * p + r] += v; }
        }
    }
}

// Signature: epps_singleton_2samp(xs, ys, t = []) -> TestResult, code. An
// empty list t means the points 0.4 and 0.8. Both samples need 5 values.
// Every t is above 0, and at most EX_MAX_POINTS points are allowed. The
// df field holds the rank of the covariance. A semi-quartile range of 0
// gives ERR_UNDEFINED, and a covariance with no rank does too.
fn epps_singleton_2samp(xs: const f64[:], ys: const f64[:], t: const f64[:] = []) -> TestResult, i64 {
    var c = check_two(xs, ys, 5);
    if c != OK { return fault_result(c); }
    if t.len > EX_MAX_POINTS { return fault_result(ERR_TOO_LONG); }
    for v in t {
        c = check_scalar(v);
        if c != OK { return fault_result(c); }
        if v <= 0.0 { return fault_result(ERR_DOMAIN); }
    }
    var tt: f64[>..] = [];
    if t.len == 0 { tt.push(0.4); tt.push(0.8); } else { for v in t { tt.push(v); } }
    let q = tt.len;
    let p = 2 * q;
    let nx = xs.len;
    let ny = ys.len;
    let n = (nx + ny) as f64;
    var pool: f64[>..] = [];
    concat2(pool, xs, ys);
    let range, ci = iqr(pool[..], METHOD_LINEAR);
    if ci != OK { return fault_result(ci); }
    let sigma = range / 2.0;
    if !(sigma > 0.0) { return fault_result(ERR_UNDEFINED); }
    var ts: f64[>..] = [];
    for v in tt { ts.push(v / sigma); }
    var cov: f64[>..] = [];
    for i in 0..p * p { cov.push(0.0); }
    var mx: f64[>..] = [];
    var my: f64[>..] = [];
    ex_es_cov(cov, mx, xs, ts[..], n / nx);
    ex_es_cov(cov, my, ys, ts[..], n / ny);
    var vec: f64[>..] = [];
    ex_jacobi(cov, vec, p);
    var top: f64 = 0.0;
    for i in 0..p { if cov[i * p + i] > top { top = cov[i * p + i]; } }
    if !(top > 0.0) { return fault_result(ERR_UNDEFINED); }
    let cut = p * EPS * top;
    var rank: i64 = 0;
    var w: f64 = 0.0;
    for e in 0..p {
        let lam = cov[e * p + e];
        if lam <= cut { continue; }
        rank += 1;
        var proj: f64 = 0.0;
        for k in 0..p { proj += vec[k * p + e] * (mx[k] - my[k]); }
        w += proj * proj / lam;
    }
    if rank == 0 { return fault_result(ERR_UNDEFINED); }
    w = w * n;
    if max(nx, ny) < 25 {
        let nxf = nx as f64;
        let nyf = ny as f64;
        w = w / (1.0 + pow(n, -0.45) + 10.1 * (pow(nxf, -1.7) + pow(nyf, -1.7)));
    }
    let pv, cp = chi2_sf(w, rank);
    if cp != OK { return fault_result(cp); }
    return good_result(w, pv, rank, nan());
}

// ---- the Baumgartner-Weiss-Schindler test -----------------------------

// The statistic B of two sorted lists of ranks. The two-sided form squares
// each numerator, and a one-sided form multiplies it by its own size.
fn ex_bws_stat(rx: const f64[:], ry: const f64[:], two: bool) -> f64 {
    let n = rx.len as f64;
    let m = ry.len as f64;
    var bx: f64 = 0.0;
    var by: f64 = 0.0;
    for i in 0..rx.len {
        let k = (i + 1) as f64;
        var num = rx[i] - (m + n) / n * k;
        if two { num = num * num; } else { num = num * abs(num); }
        let den = k / (n + 1.0) * (1.0 - k / (n + 1.0)) * m * (m + n) / n;
        bx += num / den;
    }
    for j in 0..ry.len {
        let k = (j + 1) as f64;
        var num = ry[j] - (m + n) / m * k;
        if two { num = num * num; } else { num = num * abs(num); }
        let den = k / (m + 1.0) * (1.0 - k / (m + 1.0)) * n * (m + n) / m;
        by += num / den;
    }
    bx = bx / n;
    by = by / m;
    if two { return (bx + by) / 2.0; }
    return (bx - by) / 2.0;
}

// The statistic of the split that puts the ranks z[idx[0..nx]] in the
// first sample. The list idx rises.
fn ex_bws_split(z: const f64[:], idx: const i64[:], rx: f64[>..]&, ry: f64[>..]&, two: bool) -> f64 {
    rx.resize(0);
    ry.resize(0);
    var next: i64 = 0;
    for i in 0..z.len {
        if next < idx.len && idx[next] == i {
            rx.push(z[i]);
            next += 1;
        } else {
            ry.push(z[i]);
        }
    }
    return ex_bws_stat(rx[..], ry[..], two);
}

// Signature: bws_test(xs, ys, r, alt = ALT_TWO_SIDED, resamples = 9999)
// -> TestResult, code. The p-value is a permutation p-value of the ranks.
// When the count of splits is at most MAX_EXACT_CELLS, every split counts,
// and the generator r stays unused. Above that, `resamples` random splits
// from r give the p-value, with one added to the count and to the total.
// SciPy draws random splits above 9999, so its p-value moves from run to
// run. The tolerance of 100 machine epsilons guards the comparison.
fn bws_test(xs: const f64[:], ys: const f64[:], r: rng&, alt: i64 = ALT_TWO_SIDED, resamples: i64 = 9999) -> TestResult, i64 {
    var c = check_two(xs, ys, 2);
    if c != OK { return fault_result(c); }
    if check_alt(alt) != OK || resamples < 1 { return fault_result(ERR_DOMAIN); }
    let nx = xs.len;
    let ny = ys.len;
    let total = nx + ny;
    if resamples > MAX_RESAMPLES || resamples > MAX_RESAMPLE_WORK / total { return fault_result(ERR_TOO_LONG); }
    var pool: f64[>..] = [];
    concat2(pool, xs, ys);
    var rk: f64[>..] = [];
    let cnt, cr = rankdata(rk, pool[..], TIE_AVERAGE);
    if cr != OK { return fault_result(cr); }
    let two = alt == ALT_TWO_SIDED;
    var rx: f64[>..] = [];
    var ry: f64[>..] = [];
    for i in 0..nx { rx.push(rk[i]); }
    for j in 0..ny { ry.push(rk[nx + j]); }
    sort(rx) { a, b => a < b };
    sort(ry) { a, b => a < b };
    let obs = ex_bws_stat(rx[..], ry[..], two);
    let gamma = abs(100.0 * EPS * obs);
    var z: f64[>..] = [];
    for v in rk { z.push(v); }
    sort(z) { a, b => a < b };
    let splits = ex_choose(total, min(nx, ny));
    var count: f64 = 0.0;
    var p: f64 = 0.0;
    if splits <= MAX_EXACT_CELLS {
        var idx: i64[>..] = [];
        for i in 0..nx { idx.push(i); }
        var done = false;
        var seen: f64 = 0.0;
        while !done {
            let v = ex_bws_split(z[..], idx[..], rx, ry, two);
            seen += 1.0;
            if alt == ALT_LESS {
                if v <= obs + gamma { count += 1.0; }
            } else {
                if v >= obs - gamma { count += 1.0; }
            }
            var i = nx - 1;
            while i >= 0 && idx[i] == total - nx + i { i -= 1; }
            if i < 0 { done = true; } else {
                idx[i] += 1;
                for j in i + 1..nx { idx[j] = idx[j - 1] + 1; }
            }
        }
        p = count / seen;
    } else {
        var idx: i64[>..] = [];
        var order: i64[>..] = [];
        for i in 0..total { order.push(i); }
        for k in 0..resamples {
            rs_shuffle_index(order[..], 0, total, r);
            idx.resize(0);
            for i in 0..nx { idx.push(order[i]); }
            sort(idx) { a, b => a < b };
            let v = ex_bws_split(z[..], idx[..], rx, ry, two);
            if alt == ALT_LESS {
                if v <= obs + gamma { count += 1.0; }
            } else {
                if v >= obs - gamma { count += 1.0; }
            }
        }
        p = (count + 1.0) / (resamples + 1.0);
    }
    return good_result(obs, p, nan(), nan());
}

// ---- the k-sample Anderson-Darling test -------------------------------

// The seven significance levels of the critical values.
fn ex_ksamp_levels(out: f64[>..]&) {
    out.resize(0);
    out.push(0.25); out.push(0.1); out.push(0.05); out.push(0.025);
    out.push(0.01); out.push(0.005); out.push(0.001);
}

// The p-value from a quadratic fit of the log of the levels on the
// critical values. The three normal equations are in the centered
// variable, so the system stays well conditioned.
fn ex_ksamp_pvalue(crit: const f64[:], a2: f64) -> f64 {
    var lv: f64[>..] = [];
    ex_ksamp_levels(lv);
    if a2 < crit[0] { return 0.25; }
    if a2 > crit[6] { return 0.001; }
    var mean: f64 = 0.0;
    for v in crit { mean += v; }
    mean = mean / crit.len;
    var scale: f64 = 0.0;
    for v in crit { if abs(v - mean) > scale { scale = abs(v - mean); } }
    var m: f64[12] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    // rows of the augmented matrix: 3 by 4
    for i in 0..crit.len {
        let t = (crit[i] - mean) / scale;
        let y = log(lv[i]);
        let pw: f64[5] = [1.0, t, t * t, t * t * t, t * t * t * t];
        for r in 0..3 {
            for cc in 0..3 { m[r * 4 + cc] += pw[r + cc]; }
            m[r * 4 + 3] += pw[r] * y;
        }
    }
    for col in 0..3 {
        var piv = col;
        for r in col + 1..3 { if abs(m[r * 4 + col]) > abs(m[piv * 4 + col]) { piv = r; } }
        if piv != col {
            for cc in 0..4 {
                let tmp = m[col * 4 + cc];
                m[col * 4 + cc] = m[piv * 4 + cc];
                m[piv * 4 + cc] = tmp;
            }
        }
        for r in col + 1..3 {
            let f = m[r * 4 + col] / m[col * 4 + col];
            for cc in col..4 { m[r * 4 + cc] -= f * m[col * 4 + cc]; }
        }
    }
    var co: f64[3] = [0.0, 0.0, 0.0];
    var r = 2;
    while r >= 0 {
        var s = m[r * 4 + 3];
        for cc in r + 1..3 { s -= m[r * 4 + cc] * co[cc]; }
        co[r] = s / m[r * 4 + r];
        r -= 1;
    }
    let t = (a2 - mean) / scale;
    return exp(co[0] + co[1] * t + co[2] * t * t);
}

// Signature: anderson_ksamp(out, groups, sizes, variant = KSAMP_MIDRANK)
// -> TestResult, code. The groups sit in one buffer, as in f_oneway. The
// seven critical values for the levels 0.25, 0.1, 0.05, 0.025, 0.01, 0.005
// and 0.001 go into out. The p-value comes from the quadratic fit of
// SciPy: it is capped at 0.25 below the first critical value and floored
// at 0.001 above the last. All values equal gives ERR_UNDEFINED, and so
// does a total of fewer than 4 values.
fn anderson_ksamp(out: f64[>..]&, groups: const f64[:], sizes: const i64[:], variant: i64 = KSAMP_MIDRANK) -> TestResult, i64 {
    out.resize(0);
    var c = check_groups(groups, sizes);
    if c != OK { return fault_result(c); }
    if variant < KSAMP_MIDRANK || variant > KSAMP_CONTINUOUS { return fault_result(ERR_DOMAIN); }
    let k = sizes.len;
    let total = groups.len;
    if total < 4 { return fault_result(ERR_TOO_SHORT); }
    var pairs: Vi[>..] = [];
    var g: i64 = 0;
    var left: i64 = sizes[0];
    for i in 0..total {
        while left == 0 { g += 1; left = sizes[g]; }
        pairs.push(Vi { v: groups[i], i: g });
        left -= 1;
    }
    stable_sort(pairs) { a, b => a.v < b.v };
    if pairs[0].v == pairs[total - 1].v { return fault_result(ERR_UNDEFINED); }
    var uniq: i64 = 0;
    for i in 0..total { if i == 0 || pairs[i].v != pairs[i - 1].v { uniq += 1; } }
    let steps = if variant == KSAMP_CONTINUOUS { total } else { uniq };
    if steps * k > MAX_WORK / 4 { return fault_result(ERR_TOO_LONG); }
    let nf = total as f64;
    var le: f64[>..] = [];     // the count of each group at or below the value
    var eq: f64[>..] = [];     // the count of each group at the value
    for i in 0..k { le.push(0.0); eq.push(0.0); }
    var acc: f64[>..] = [];    // the sum for each group
    for i in 0..k { acc.push(0.0); }
    var below: f64 = 0.0;      // the pooled count below the value
    var j: i64 = 0;
    var seen: f64 = 0.0;
    var run: i64 = 0;
    while j < total {
        var e = j + 1;
        while e < total && pairs[e].v == pairs[j].v { e += 1; }
        let lj = (e - j) as f64;
        for i in 0..k { eq[i] = 0.0; }
        for q in j..e { eq[pairs[q].i] += 1.0; }
        for i in 0..k { le[i] += eq[i]; }
        run += 1;
        if variant == KSAMP_MIDRANK {
            let bj = below + lj / 2.0;
            let den = bj * (nf - bj) - nf * lj / 4.0;
            for i in 0..k {
                let mij = le[i] - eq[i] / 2.0;
                let dd = nf * mij - bj * sizes[i];
                acc[i] += lj / nf * dd * dd / den;
            }
        } else if variant == KSAMP_RIGHT {
            if e < total {
                let bj = below + lj;
                for i in 0..k {
                    let dd = nf * le[i] - bj * sizes[i];
                    acc[i] += lj / nf * dd * dd / (bj * (nf - bj));
                }
            }
        } else {
            for q in j..e {
                if q + 1 < total {
                    let pos = (q + 1) as f64;
                    for i in 0..k {
                        let dd = nf * le[i] - pos * sizes[i];
                        acc[i] += dd * dd / (pos * (nf - pos));
                    }
                }
            }
        }
        below += lj;
        j = e;
    }
    var tot = acc_new();
    for i in 0..k { acc_add(tot, acc[i] / sizes[i]); }
    var a2kn = acc_total(tot);
    if variant == KSAMP_MIDRANK { a2kn = a2kn * (nf - 1.0) / nf; }
    if variant == KSAMP_CONTINUOUS { a2kn = a2kn / nf; }
    // the normalization of Scholz and Stephens
    var hh = acc_new();
    for i in 0..k { acc_add(hh, 1.0 / sizes[i]); }
    let hsum = acc_total(hh);
    var cs: f64 = 0.0;
    var gsum = acc_new();
    var cur = total - 1;
    var idx2: i64 = 2;
    while cur > 1 {
        cs += 1.0 / cur;
        acc_add(gsum, cs / idx2);
        cur -= 1;
        idx2 += 1;
    }
    let h = cs + 1.0;
    let gg = acc_total(gsum);
    let kf = k as f64;
    let aa = (4.0 * gg - 6.0) * (kf - 1.0) + (10.0 - 6.0 * gg) * hsum;
    let bb = (2.0 * gg - 4.0) * kf * kf + 8.0 * h * kf + (2.0 * gg - 14.0 * h - 4.0) * hsum - 8.0 * h + 4.0 * gg - 6.0;
    let cc = (6.0 * h + 2.0 * gg - 2.0) * kf * kf + (4.0 * h - 4.0 * gg + 6.0) * kf + (2.0 * h - 6.0) * hsum + 4.0 * h;
    let dd = (2.0 * h + 6.0) * kf * kf - 4.0 * h * kf;
    let sigsq = (aa * nf * nf * nf + bb * nf * nf + cc * nf + dd) / ((nf - 1.0) * (nf - 2.0) * (nf - 3.0));
    if !(sigsq > 0.0) { return fault_result(ERR_UNDEFINED); }
    let mm = kf - 1.0;
    let a2 = (a2kn - mm) / sqrt(sigsq);
    let b0: f64[7] = [0.675, 1.281, 1.645, 1.96, 2.326, 2.573, 3.085];
    let b1: f64[7] = [-0.245, 0.25, 0.678, 1.149, 1.822, 2.364, 3.615];
    let b2: f64[7] = [-0.105, -0.305, -0.362, -0.391, -0.396, -0.345, -0.154];
    for i in 0..7 { out.push(b0[i] + b1[i] / sqrt(mm) + b2[i] / mm); }
    let p = ex_ksamp_pvalue(out[..], a2);
    return good_result(a2, p, nan(), nan());
}

// ---- the Page trend test ----------------------------------------------

// The probability of each sum L of the products j times the value at j
// over all permutations of 1 to k. The entry l - a of out holds the
// probability of the sum l, for l from a to b. The count of the permutations
// by their partial sum runs over the subsets of the values already used.
fn ex_page_row(out: f64[>..]&, k: i64) {
    out.resize(0);
    let a = k * (k + 1) * (k + 2) / 6;
    let b = k * (k + 1) * (2 * k + 1) / 6;
    let w = b + 1;
    let full = 1 << k;
    var pos: i64[>..] = [];
    var count: i64[>..] = [];
    for t in 0..k + 1 { count.push(0); }
    for mask in 0..full {
        var bits: i64 = 0;
        for i in 0..k { if (mask >> i) & 1 == 1 { bits += 1; } }
        pos.push(count[bits]);
        count[bits] += 1;
    }
    var cur: f64[>..] = [];
    for i in 0..w { cur.push(0.0); }
    cur[0] = 1.0;
    for t in 0..k {
        var nxt: f64[>..] = [];
        for i in 0..count[t + 1] * w { nxt.push(0.0); }
        for mask in 0..full {
            var bits: i64 = 0;
            for i in 0..k { if (mask >> i) & 1 == 1 { bits += 1; } }
            if bits != t { continue; }
            let row = pos[mask] * w;
            for v in 0..k {
                if (mask >> v) & 1 == 1 { continue; }
                let add = (t + 1) * (v + 1);
                let nrow = pos[mask | (1 << v)] * w;
                for s in 0..w - add {
                    let x = cur[row + s];
                    if x != 0.0 { nxt[nrow + s + add] += x; }
                }
            }
        }
        cur.resize(0);
        for i in 0..nxt.len { cur.push(nxt[i]); }
    }
    var fact: f64 = 1.0;
    for i in 2..k + 1 { fact = fact * i; }
    for l in a..b + 1 { out.push(cur[l] / fact); }
}

// Signature: page_trend_test(data, columns, ranked = false, predicted = [],
// method = PAGE_AUTO) -> TestResult, code. The data hold one subject for
// each row, row by row, with `columns` values in a row. An empty
// predicted list means the order 1 to columns, and otherwise it holds a
// permutation of 1 to columns. With ranked set, each row holds ranks from
// 1 to columns already. The statistic is L. The df field holds PAGE_EXACT
// or PAGE_ASYMPTOTIC, the method that gave the p-value. The exact form
// needs at most MAX_PAGE_COLUMNS columns. The p-value is the upper tail,
// and the exact form cuts L to a whole number, as SciPy does.
fn page_trend_test(data: const f64[:], columns: i64, ranked: bool = false, predicted: const f64[:] = [], method: i64 = PAGE_AUTO) -> TestResult, i64 {
    var c = check_sample(data, 1);
    if c != OK { return fault_result(c); }
    if method < PAGE_AUTO || method > PAGE_ASYMPTOTIC { return fault_result(ERR_DOMAIN); }
    if columns < 3 { return fault_result(ERR_TOO_SHORT); }
    if data.len % columns != 0 { return fault_result(ERR_LENGTH); }
    let m = data.len / columns;
    if m < 2 { return fault_result(ERR_TOO_SHORT); }
    let n = columns;
    var pred: f64[>..] = [];
    if predicted.len == 0 {
        for j in 0..n { pred.push(j + 1); }
    } else {
        if predicted.len != n { return fault_result(ERR_LENGTH); }
        var seen: i64[>..] = [];
        for j in 0..n { seen.push(0); }
        for v in predicted {
            if v < 1.0 || v > n || v != floor(v) { return fault_result(ERR_DOMAIN); }
            seen[v as i64 - 1] += 1;
        }
        for j in 0..n { if seen[j] != 1 { return fault_result(ERR_DOMAIN); } }
        for v in predicted { pred.push(v); }
    }
    var sums: f64[>..] = [];
    for j in 0..n { sums.push(0.0); }
    if ranked {
        for v in data { if v < 1.0 || v > n { return fault_result(ERR_DOMAIN); } }
        for i in 0..m { for j in 0..n { sums[j] += data[i * n + j]; } }
    } else {
        var row: f64[>..] = [];
        var rk: f64[>..] = [];
        for i in 0..m {
            row.resize(0);
            for j in 0..n { row.push(data[i * n + j]); }
            let cnt, cr = rankdata(rk, row[..], TIE_AVERAGE);
            if cr != OK { return fault_result(cr); }
            for j in 0..n { sums[j] += rk[j]; }
        }
    }
    var lacc = acc_new();
    for j in 0..n { acc_add(lacc, pred[j] * sums[j]); }
    let l = acc_total(lacc);
    var exact = method == PAGE_EXACT;
    if method == PAGE_AUTO { exact = !(n > 8 || (m > 12 && n > 3) || m > 20); }
    if !exact {
        let mf = m as f64;
        let nf = n as f64;
        let e0 = mf * nf * (nf + 1.0) * (nf + 1.0) / 4.0;
        let v0 = mf * nf * nf * (nf + 1.0) * (nf * nf - 1.0) / 144.0;
        let pv, cp = norm_sf((l - e0) / sqrt(v0));
        if cp != OK { return fault_result(cp); }
        return good_result(l, pv, PAGE_ASYMPTOTIC, nan());
    }
    if n > MAX_PAGE_COLUMNS { return fault_result(ERR_TOO_LONG); }
    let a = n * (n + 1) * (n + 2) / 6;
    let b = n * (n + 1) * (2 * n + 1) / 6;
    let w1 = (b - a + 1) as f64;
    if m * (b - a) + 1 > MAX_EXACT_CELLS || m * m * w1 * w1 / 2.0 > 500000000.0 { return fault_result(ERR_TOO_LONG); }
    var p1: f64[>..] = [];
    ex_page_row(p1, n);
    var cur: f64[>..] = [];
    for v in p1 { cur.push(v); }
    for r in 1..m {
        var nxt: f64[>..] = [];
        for i in 0..cur.len + p1.len - 1 { nxt.push(0.0); }
        for i in 0..cur.len {
            for j in 0..p1.len { nxt[i + j] += cur[i] * p1[j]; }
        }
        cur.resize(0);
        for v in nxt { cur.push(v); }
    }
    let first = trunc(l) as i64 - m * a;
    var tail = acc_new();
    var from = first;
    if from < 0 { from = 0; }
    for i in from..cur.len { acc_add(tail, cur[i]); }
    return good_result(l, acc_total(tail), PAGE_EXACT, nan());
}
