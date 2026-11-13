// The normality test and the correlation tests of the package. Each call
// returns a TestResult and then a fault code. A fault gives NaN in every
// field of the result.
namespace statistics;

import std;
import math;

fn poly_at(c: const f64[:], x: f64) -> f64 {
    var r: f64 = 0.0;
    var k = c.len;
    while k > 0 {
        k -= 1;
        r = r * x + c[k];
    }
    return r;
}

// The Shapiro-Wilk test of Royston (1992). It needs 3 values. A sample
// of equal values gives the undefined fault.
fn shapiro(xs: const f64[:]) -> TestResult, i64 {
    let c = check_sample(xs, 3);
    if c != 0 { return fault_result(c); }
    let s = sort_copy(xs);
    let n = s.len;
    let range = s[n - 1] - s[0];
    if is_inf(range) { return fault_result(ERR_OVERFLOW); }
    if range == 0.0 { return fault_result(ERR_UNDEFINED); }
    var y: f64[>..] = [];
    for v in s { y.push((v - s[0]) / range); }
    let half = n / 2;
    var a: f64[>..] = [];
    if n == 3 {
        a.push(sqrt(0.5));
    } else {
        let an25: f64 = n + 0.25;
        let c1: f64[6] = [0.0, 0.221157, -0.147981, -2.07119, 4.434685, -2.706056];
        let c2: f64[6] = [0.0, 0.042981, -0.293762, -1.752461, 5.682633, -3.582633];
        var m: f64[>..] = [];
        var sum2: f64 = 0.0;
        for i in 1..(half + 1) {
            let q = norm_std_ppf((i - 0.375) / an25);
            m.push(q);
            sum2 += q * q;
        }
        sum2 = 2.0 * sum2;
        let root = sqrt(sum2);
        let nf: f64 = n;
        let rsn = 1.0 / sqrt(nf);
        let a1 = poly_at(c1, rsn) - m[0] / root;
        var first = 1;
        var fac: f64 = 0.0;
        a.push(a1);
        if n > 5 {
            let a2 = poly_at(c2, rsn) - m[1] / root;
            first = 2;
            fac = sqrt((sum2 - 2.0 * m[0] * m[0] - 2.0 * m[1] * m[1]) / (1.0 - 2.0 * a1 * a1 - 2.0 * a2 * a2));
            a.push(a2);
        } else {
            fac = sqrt((sum2 - 2.0 * m[0] * m[0]) / (1.0 - 2.0 * a1 * a1));
        }
        for i in first..half { a.push(-m[i] / fac); }
    }
    var mean_y: f64 = 0.0;
    for v in y { mean_y += v; }
    let nn: f64 = n;
    mean_y = mean_y / nn;
    var ssq: f64 = 0.0;
    var b: f64 = 0.0;
    for v in y { ssq += (v - mean_y) * (v - mean_y); }
    for i in 0..half { b += a[i] * (y[n - 1 - i] - y[i]); }
    var w = b * b / ssq;
    if w > 1.0 { w = 1.0; }
    var p: f64 = 0.0;
    if n == 3 {
        p = 1.90985931710274 * (asin(sqrt(w)) - 1.04719755119660);
        if p < 0.0 { p = 0.0; }
    } else {
        let w1 = 1.0 - w;
        var ly = log(w1);
        let nf: f64 = n;
        let lx = log(nf);
        var mu: f64 = 0.0;
        var sigma: f64 = 0.0;
        if n <= 11 {
            let g: f64[2] = [-2.273, 0.459];
            let gamma = poly_at(g, nf);
            // W stays far below this limit for 4 to 11 values, so the
            // guard only keeps a NaN out of the result.
            if ly >= gamma { return fault_result(ERR_UNDEFINED); }
            ly = -log(gamma - ly);
            let c3: f64[4] = [0.544, -0.39978, 0.025054, -6.714e-4];
            let c4: f64[4] = [1.3822, -0.77857, 0.062767, -0.0020322];
            mu = poly_at(c3, nf);
            sigma = exp(poly_at(c4, nf));
        } else {
            let c5: f64[4] = [-1.5861, -0.31082, -0.083751, 0.0038915];
            let c6: f64[3] = [-0.4803, -0.082676, 0.0030302];
            mu = poly_at(c5, lx);
            sigma = exp(poly_at(c6, lx));
        }
        let q, cq = norm_sf((ly - mu) / sigma);
        if cq != 0 { return fault_result(cq); }
        p = q;
    }
    return good_result(w, p, nan(), nan());
}

// The Pearson test. The p-value uses the beta form of the distribution
// of r, so a p-value near 0 keeps its digits. Two points give p = 1.
fn pearsonr(xs: const f64[:], ys: const f64[:], alt: i64 = ALT_TWO_SIDED) -> TestResult, i64 {
    let ca = check_alt(alt);
    let c = check_pair(xs, ys, 2);
    if c != 0 { return fault_result(c); }
    if ca != 0 { return fault_result(ca); }
    let r, c2 = correlation(xs, ys);
    if c2 != 0 { return fault_result(c2); }
    let n = xs.len;
    let df: f64 = n - 2;
    if n == 2 { return good_result(r, 1.0, df, nan()); }
    let ab = n / 2.0 - 1.0;
    var x: f64 = 0.0;
    if alt == ALT_TWO_SIDED { x = 0.5 * (1.0 - abs(r)); }
    if alt == ALT_LESS { x = 0.5 * (1.0 + r); }
    if alt == ALT_GREATER { x = 0.5 * (1.0 - r); }
    if x < 0.0 { x = 0.0; }
    if x > 1.0 { x = 1.0; }
    let q, cb = betainc(ab, ab, x);
    if cb != 0 { return fault_result(cb); }
    var p = q;
    if alt == ALT_TWO_SIDED { p = 2.0 * q; }
    return good_result(r, p, df, nan());
}

// The Spearman test. The p-value uses the t form with n - 2 degrees of
// freedom. Two points give p = 1.
fn spearmanr(xs: const f64[:], ys: const f64[:], alt: i64 = ALT_TWO_SIDED) -> TestResult, i64 {
    let ca = check_alt(alt);
    let c = check_pair(xs, ys, 2);
    if c != 0 { return fault_result(c); }
    if ca != 0 { return fault_result(ca); }
    let r, c2 = spearman(xs, ys);
    if c2 != 0 { return fault_result(c2); }
    let n = xs.len;
    let df: f64 = n - 2;
    if n == 2 { return good_result(r, 1.0, df, nan()); }
    var t = (1.0 / 0.0);
    let den = (1.0 + r) * (1.0 - r);
    if den > 0.0 { t = r * sqrt(df / den); }
    if den <= 0.0 && r < 0.0 { t = -(1.0 / 0.0); }
    let p, cp = t_pvalue(t, df, alt);
    if cp != 0 { return fault_result(cp); }
    return good_result(r, p, df, nan());
}

// The sums over the runs of equal values of a sorted sample: the pairs
// t(t-1)/2, then t(t-1)(t-2) and t(t-1)(2t+5).
fn tie_sums(s: const f64[:]) -> f64, f64, f64 {
    var pairs: f64 = 0.0;
    var s0: f64 = 0.0;
    var s1: f64 = 0.0;
    var k: i64 = 0;
    while k < s.len {
        var j = k + 1;
        while j < s.len && s[j] == s[k] { j += 1; }
        let t: f64 = j - k;
        pairs += t * (t - 1.0) / 2.0;
        s0 += t * (t - 1.0) * (t - 2.0);
        s1 += t * (t - 1.0) * (2.0 * t + 5.0);
        k = j;
    }
    return pairs, s0, s1;
}

// The exact tails of the number of discordant pairs of n values with no
// ties. The first result is P(D <= c) and the second is P(D = c).
fn kendall_exact_tail(n: i64, c: i64) -> f64, f64 {
    var cur: f64[>..] = [];
    var nxt: f64[>..] = [];
    for k in 0..(c + 1) {
        cur.push(0.0);
        nxt.push(0.0);
    }
    cur[0] = 1.0;
    for j in 2..(n + 1) {
        var window: f64 = 0.0;
        for k in 0..(c + 1) {
            window += cur[k];
            if k >= j { window -= cur[k - j]; }
            nxt[k] = window / j;
        }
        for k in 0..(c + 1) { cur[k] = nxt[k]; }
    }
    var total: f64 = 0.0;
    for k in 0..(c + 1) { total += cur[k]; }
    return total, cur[c];
}

// The Kendall tau-b test. The p-value is exact for a sample with no ties
// of at most MAX_EXACT_KENDALL values, or one pair from a perfect order.
// Any other sample uses the normal form with the variance of the ties.
fn kendalltau(xs: const f64[:], ys: const f64[:], alt: i64 = ALT_TWO_SIDED) -> TestResult, i64 {
    let ca = check_alt(alt);
    let c = check_pair(xs, ys, 2);
    if c != 0 { return fault_result(c); }
    if ca != 0 { return fault_result(ca); }
    let n = xs.len;
    let cd, t = kendall_counts(xs, ys);
    // The tau-b of the pairs, from the same counts.
    let tot = n * (n - 1) / 2;
    let a = tot - t.n1;
    let b = tot - t.n2;
    if a <= 0 || b <= 0 { return fault_result(ERR_UNDEFINED); }
    // Widen before the product. Two i64 values multiply as i64 and wrap
    // above 2 to the power 63, which happens at about 78,000 values.
    let fa: f64 = a;
    let fb: f64 = b;
    let ab = fa * fb;
    let den = sqrt(ab);
    if den == 0.0 { return fault_result(ERR_UNDEFINED); }
    let tau = cd / den;
    if t.n1 == 0 && t.n2 == 0 {
        let con = (tot + cd) / 2;
        let dis = tot - con;
        let low = if dis < con { dis } else { con };
        if n <= MAX_EXACT_KENDALL || low <= 1 {
            if n == 2 { return good_result(tau, 1.0, nan(), nan()); }
            let right = con >= tot - con;
            let cc = if right { tot - con } else { con };
            let cdf, mass = kendall_exact_tail(n, cc);
            var p: f64 = 1.0;
            if alt == ALT_TWO_SIDED {
                if 4 * cc != n * (n - 1) { p = 2.0 * cdf; }
            } else {
                let agree = right == (alt == ALT_GREATER);
                if agree { p = cdf; } else { p = 1.0 - cdf + mass; }
            }
            return good_result(tau, p, nan(), nan());
        }
    }
    let px, x0, x1 = tie_sums(sort_copy(xs));
    let py, y0, y1 = tie_sums(sort_copy(ys));
    let nf: f64 = n;
    let m = nf * (nf - 1.0);
    let v = (m * (2.0 * nf + 5.0) - x1 - y1) / 18.0 + 2.0 * px * py / m + x0 * y0 / (9.0 * m * (nf - 2.0));
    if !(v > 0.0) { return fault_result(ERR_UNDEFINED); }
    let z = cd / sqrt(v);
    var p: f64 = 0.0;
    var cz: i64 = 0;
    if alt == ALT_TWO_SIDED {
        let q, cq = norm_sf(abs(z));
        p = 2.0 * q;
        cz = cq;
    } else if alt == ALT_LESS {
        let q, cq = norm_cdf(z);
        p = q;
        cz = cq;
    } else {
        let q, cq = norm_sf(z);
        p = q;
        cz = cq;
    }
    if cz != 0 { return fault_result(cz); }
    return good_result(tau, p, nan(), nan());
}
