// The multiple testing corrections of the package. Each call returns a
// TestResult and then a fault code. A fault gives NaN in every field of
// the result.
namespace statistics;

import std;
import math;

fn multipletests(out: f64[>..]&, ps: const f64[:], method: i64 = ADJUST_HOLM) -> i64, i64 {
    out.resize(0);
    let c = check_sample(ps, 1);
    if c != 0 { return 0, c; }
    if method > ADJUST_NONE || method < 0 { return 0, ERR_DOMAIN; }
    if method == ADJUST_HOMMEL && ps.len > HOMMEL_MAX { return 0, ERR_TOO_LONG; }
    let cp = check_probabilities(ps);
    if cp != 0 { return 0, cp; }
    let n = ps.len;
    out.resize(n, 0.0);
    let nf: f64 = n;
    if method == ADJUST_BONFERRONI {
        // The Bonferroni form keeps the order of the list, so it needs no
        // sort of the indices.
        for k in 0..n { out[k] = ps[k] * nf; }
    } else if method == ADJUST_NONE {
        for k in 0..n { out[k] = ps[k]; }
    } else {
        var order: i64[>..] = [];
        for k in 0..n { order.push(k); }
        sort(order) { a, b => ps[a] < ps[b] || (ps[a] == ps[b] && a < b) };
        if method == ADJUST_HOLM {
            var run: f64 = 0.0;
            for k in 0..n {
                let v = ps[order[k]] * (nf - k);
                if v > run { run = v; }
                out[order[k]] = run;
            }
        } else if method == ADJUST_HOCHBERG {
            var run: f64 = 1.0;
            var k = n;
            while k > 0 {
                k -= 1;
                let lim = ps[order[k]] * (nf - k);
                if lim < run { run = lim; }
                out[order[k]] = run;
            }
        } else if method == ADJUST_HOMMEL {
            hommel_adjust(out, ps, order);
        } else {
            // The Benjamini-Hochberg form. The Benjamini-Yekutieli form
            // multiplies it by the harmonic sum, added from the small
            // terms to the large ones.
            var factor: f64 = 1.0;
            if method == ADJUST_BY {
                var h: f64 = 0.0;
                var i = n;
                while i > 0 {
                    h += 1.0 / i;
                    i -= 1;
                }
                factor = h;
            }
            var run: f64 = 1.0;
            var k = n;
            while k > 0 {
                k -= 1;
                let frac = ps[order[k]] * nf / (k + 1) * factor;
                if frac < run { run = frac; }
                out[order[k]] = run;
            }
        }
    }
    for k in 0..n {
        if out[k] > 1.0 { out[k] = 1.0; }
    }
    return n, 0;
}

// The Hommel form, in the steps of Wright as R p.adjust holds them. The
// loop over m runs from n-1 down to 2 and each step costs n, so the call
// costs n squared. order lists the indices of ps from the smallest value.
// The count is at most HOMMEL_MAX. A list of one or two values gives the
// Hochberg values, because the loop over m then never runs.
fn hommel_adjust(adj: f64[>..]&, ps: const f64[:], order: const i64[:]) {
    let n = ps.len;
    var sp: f64[>..] = [];
    for k in 0..n { sp.push(ps[order[k]]); }
    var pa: f64[>..] = [];
    var lo: f64 = sp[0] * n;
    for k in 0..n {
        let v = sp[k] * n / (k + 1);
        if v < lo { lo = v; }
    }
    for k in 0..n { pa.push(lo); }
    var m = n - 1;
    while m >= 2 {
        // The indices 0 to n-m take the smaller of m times the value and
        // the smallest of m*p/(2..m) over the indices n-m+1 to n-1. The
        // later indices take no part: the floor at the p-value covers them.
        var q1: f64 = 1e300;
        for t in 0..(m - 1) {
            let v = sp[n - m + 1 + t] * m / (t + 2);
            if v < q1 { q1 = v; }
        }
        for k in 0..(n - m + 1) {
            var v = sp[k] * m;
            if q1 < v { v = q1; }
            if v > pa[k] { pa[k] = v; }
        }
        m -= 1;
    }
    for k in 0..n {
        var v = pa[k];
        if sp[k] > v { v = sp[k]; }
        adj[order[k]] = v;
    }
}

// The same values as multipletests, under the name of SciPy. The method
// is ADJUST_BH or ADJUST_BY. Another number gives the domain fault.
fn false_discovery_control(out: f64[>..]&, ps: const f64[:], method: i64 = ADJUST_BH) -> i64, i64 {
    // The sample checks come first, so a bad list gives its own fault
    // before the method fault. multipletests runs them. The number -1 is
    // the method that gives the domain fault.
    var m = method;
    if m != ADJUST_BH && m != ADJUST_BY { m = -1; }
    return multipletests(out, ps, m);
}

// Combines the p-values of independent tests of one hypothesis. The
// method is one of the COMBINE_ constants. Only the Stouffer method takes
// weights. An empty list of weights means equal weights. A p-value of 0
// gives an infinite Fisher or Stouffer statistic, and a p-value of 1 does
// the same for Stouffer, so those calls give the overflow fault.
fn combine_pvalues(ps: const f64[:], method: i64 = COMBINE_FISHER, weights: const f64[:] = []) -> TestResult, i64 {
    let c = check_sample(ps, 1);
    if c != 0 { return fault_result(c); }
    let cp = check_probabilities(ps);
    if cp != 0 { return fault_result(cp); }
    if method > COMBINE_TIPPETT || method < 0 { return fault_result(ERR_DOMAIN); }
    let n = ps.len;
    let k: f64 = n;
    // The weights are scaled by the largest one, so a tiny weight and a
    // huge weight keep their share of the evidence. The sum of the
    // squared scaled weights is the divisor of the Stouffer statistic.
    var sq: f64 = 0.0;
    var wmax: f64 = 0.0;
    if weights.len > 0 {
        if method != COMBINE_STOUFFER { return fault_result(ERR_DOMAIN); }
        if weights.len != ps.len { return fault_result(ERR_LENGTH); }
        for w in weights {
            if is_nan(w) { return fault_result(ERR_DOMAIN); }
            if is_inf(w) { return fault_result(ERR_DOMAIN); }
            if w < 0.0 { return fault_result(ERR_DOMAIN); }
            if w > wmax { wmax = w; }
        }
        if wmax <= 0.0 { return fault_result(ERR_DOMAIN); }
        for w in weights {
            let u = w / wmax;
            sq += u * u;
        }
    } else {
        sq = k;
    }
    if method == COMBINE_FISHER {
        // The exact sum of the logarithms keeps the statistic at the
        // last digits of a long list.
        var a = acc_new();
        for p in ps {
            if p == 0.0 { return fault_result(ERR_OVERFLOW); }
            acc_add(a, log(p));
        }
        let sum = acc_total(a);
        return chi_result(-2.0 * sum, 2.0 * k);
    }
    if method == COMBINE_STOUFFER {
        var top: f64 = 0.0;
        for i in 0..n {
            if ps[i] == 0.0 { return fault_result(ERR_OVERFLOW); }
            if ps[i] == 1.0 { return fault_result(ERR_OVERFLOW); }
            let z, cz = norm_ppf(ps[i]);
            if cz != 0 { return fault_result(cz); }
            var u: f64 = 1.0;
            if weights.len > 0 { u = weights[i] / wmax; }
            top -= u * z;
        }
        let zs = top / sqrt(sq);
        let p, cp = norm_sf(zs);
        if cp != 0 { return fault_result(cp); }
        return good_result(zs, p, nan(), nan());
    }
    let lo, cm = minimum(ps);
    if cm != 0 { return fault_result(cm); }
    if lo >= 1.0 { return good_result(lo, 1.0, nan(), nan()); }
    let l, cl = log1p(-lo);
    if cl != 0 { return fault_result(cl); }
    let e, ce = expm1(k * l);
    if ce != 0 { return fault_result(ce); }
    return good_result(lo, -e, nan(), nan());
}
