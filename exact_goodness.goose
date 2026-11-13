// The Epps-Singleton test, the Baumgartner-Weiss-Schindler test and the
// k-sample Anderson-Darling test of the package. Each call returns a
// TestResult and then a fault code. A fault gives NaN in every field of
// the result.
namespace statistics;

import std;
import math;

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
