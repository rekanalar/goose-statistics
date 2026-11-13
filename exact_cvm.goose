// The Poisson means test, the modified Bessel function K and the
// Cramer-von Mises tests of the package. Each call returns a TestResult
// and then a fault code. A fault gives NaN in every field of the result.
namespace statistics;

import std;
import math;

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
