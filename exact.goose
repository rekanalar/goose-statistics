// The exact tests of the package. Each call returns a TestResult and then
// a fault code. A fault gives NaN in every field of the result.
//
// The file holds the shared helpers, the nuisance parameter of the 2 by 2
// tests and the Kolmogorov-Smirnov dispatchers.
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
