// The tests of variance, scale and rank: levene, bartlett, fligner,
// alexandergovern, brunnermunzel, mood, ranksums, ansari, friedmanchisquare,
// median_test, quantile_test and ttest_ind_from_stats.
//
// Every call returns a TestResult and then a fault code, like the tests of
// tests.goose. Several samples sit in one buffer and a list of sizes gives
// the length of each sample. The helpers of this file start with vr_.
namespace statistics;

import std;
import math;

// How median_test treats a value equal to the grand median.
const MEDIAN_BELOW: i64 = 0;    // the value counts as below the median
const MEDIAN_ABOVE: i64 = 1;    // the value counts as above the median
const MEDIAN_IGNORE: i64 = 2;   // the value leaves the table

// The table of the exact Ansari-Bradley count holds at most this many
// cells. Two samples of 54 values need about 160,000 cells.
const VR_ANSARI_CELLS: i64 = 400000;

// ---- shared helpers ---------------------------------------------------

// The p-value of a z statistic of the standard normal. The statistic is
// finite.
// The result of a z test: the statistic and the p-value of the normal.
fn vr_finish_z(z: f64, alt: i64) -> TestResult, i64 {
    if is_bad(z) { return fault_result(ERR_UNDEFINED); }
    let p, c = norm_z_p(z, alt);
    if c != 0 { return fault_result(c); }
    return good_result(z, p, nan(), nan());
}

// The result of a test of a chi-square statistic at df degrees of freedom.
fn vr_finish_chi(stat: f64, df: f64) -> TestResult, i64 {
    if is_nan(stat) { return fault_result(ERR_UNDEFINED); }
    if is_inf(stat) { return fault_result(ERR_OVERFLOW); }
    return chi_result(stat, df);
}

// The center of a group: the mean, the median or the trimmed mean.
fn vr_center(part: const f64[:], center: i64, trim: f64) -> f64, i64 {
    if center == CENTER_MEAN { return mean(part); }
    if center == CENTER_MEDIAN { return median(part); }
    return trim_mean(part, trim);
}

// The checks of the center and the trim of levene and fligner. They run
// after the checks of the groups.
fn vr_check_center(center: i64, trim: f64) -> i64 {
    if center == CENTER_TRIMMED {
        let c = check_scalar(trim);
        if c != 0 { return c; }
    }
    if center < CENTER_MEAN || center > CENTER_TRIMMED { return ERR_DOMAIN; }
    if center == CENTER_TRIMMED && (trim < 0.0 || trim >= 0.5) { return ERR_DOMAIN; }
    return 0;
}

// Writes the absolute deviation of every value from the center of its
// group into z, in the order of the buffer.
fn vr_deviations(z: f64[>..]&, groups: const f64[:], sizes: const i64[:], center: i64, trim: f64) -> i64 {
    z.resize(0);
    var at: i64 = 0;
    for n in sizes {
        let part = groups[at..at + n];
        at += n;
        let c, cc = vr_center(part, center, trim);
        if cc != 0 { return cc; }
        for x in part {
            let d = abs(x - c);
            if is_inf(d) { return ERR_OVERFLOW; }
            z.push(d);
        }
    }
    return 0;
}

// ---- levene, fligner and bartlett ---------------------------------------

// The Levene test of equal variances. The groups sit one after the other
// in one buffer. The center is a CENTER_ constant and trim applies to the
// trimmed mean. SciPy: levene(*samples, center='median', proportiontocut=0.05).
// A group of equal deviations only gives the undefined fault.
fn levene(groups: const f64[:], sizes: const i64[:], center: i64 = CENTER_MEDIAN, trim: f64 = 0.05) -> TestResult, i64 {
    var c = check_groups(groups, sizes);
    if c != 0 { return fault_result(c); }
    c = vr_check_center(center, trim);
    if c != 0 { return fault_result(c); }
    let k = sizes.len;
    let total = groups.len;
    if total <= k { return fault_result(ERR_TOO_SHORT); }
    var z: f64[>..] = [];
    c = vr_deviations(z, groups, sizes, center, trim);
    if c != 0 { return fault_result(c); }
    var zbars: f64[>..] = [];
    var at: i64 = 0;
    var weighted = acc_new();
    for n in sizes {
        let zb, cz = mean(z[at..at + n]);
        if cz != 0 { return fault_result(cz); }
        at += n;
        zbars.push(zb);
        acc_add(weighted, zb * n);
    }
    let zbar = acc_total(weighted) / total;
    var between = acc_new();
    var within = acc_new();
    at = 0;
    for g in 0..k {
        let n = sizes[g];
        let d = zbars[g] - zbar;
        acc_add(between, d * d * n);
        for i in at..at + n {
            let e = z[i] - zbars[g];
            acc_add(within, e * e);
        }
        at += n;
    }
    let num = acc_total(between);
    let den = acc_total(within);
    if is_inf(num) || is_inf(den) { return fault_result(ERR_OVERFLOW); }
    if den == 0.0 { return fault_result(ERR_UNDEFINED); }
    let dfn = (k - 1) as f64;
    let dfd = (total - k) as f64;
    let w = (dfd * num) / (dfn * den);
    if is_bad(w) { return fault_result(ERR_OVERFLOW); }
    let p, cp = fdist_sf(w, dfn, dfd);
    if cp != 0 { return fault_result(cp); }
    return good_result(w, p, dfn, dfd);
}

// The Fligner-Killeen test of equal variances. The arguments are those of
// levene. SciPy: fligner(*samples, center='median', proportiontocut=0.05).
fn fligner(groups: const f64[:], sizes: const i64[:], center: i64 = CENTER_MEDIAN, trim: f64 = 0.05) -> TestResult, i64 {
    var c = check_groups(groups, sizes);
    if c != 0 { return fault_result(c); }
    c = vr_check_center(center, trim);
    if c != 0 { return fault_result(c); }
    let k = sizes.len;
    let total = groups.len;
    var z: f64[>..] = [];
    c = vr_deviations(z, groups, sizes, center, trim);
    if c != 0 { return fault_result(c); }
    var rk: f64[>..] = [];
    let cnt, tt, cr = rankdata_ties(rk, z[..], TIE_AVERAGE);
    if cr != 0 { return fault_result(cr); }
    let nf = total as f64;
    var a: f64[>..] = [];
    for i in 0..total { a.push(norm_std_ppf(rk[i] / (2.0 * (nf + 1.0)) + 0.5)); }
    let abar, ca = mean(a[..]);
    if ca != 0 { return fault_result(ca); }
    let v2, cv = variance_core(a[..], 1);
    if cv != 0 { return fault_result(cv); }
    if !(v2 > 0.0) { return fault_result(ERR_UNDEFINED); }
    var s = acc_new();
    var at: i64 = 0;
    for n in sizes {
        let ai, cm = mean(a[at..at + n]);
        if cm != 0 { return fault_result(cm); }
        at += n;
        let d = ai - abar;
        acc_add(s, d * d * n);
    }
    var stat = acc_total(s) / v2;
    if stat < 0.0 { stat = 0.0; }
    return vr_finish_chi(stat, k - 1);
}

// The Bartlett test of equal variances. Each group needs two values and
// a variance above 0. SciPy: bartlett(*samples).
fn bartlett(groups: const f64[:], sizes: const i64[:]) -> TestResult, i64 {
    let c = check_groups(groups, sizes);
    if c != 0 { return fault_result(c); }
    for n in sizes { if n < 2 { return fault_result(ERR_TOO_SHORT); } }
    let k = sizes.len;
    let total = groups.len;
    var pooled = acc_new();
    var logs = acc_new();
    var recips = acc_new();
    var at: i64 = 0;
    for n in sizes {
        let v, cv = variance_core(groups[at..at + n], 1);
        if cv != 0 { return fault_result(cv); }
        at += n;
        if !(v > 0.0) { return fault_result(ERR_UNDEFINED); }
        acc_add(pooled, (n - 1) * v);
        acc_add(logs, (n - 1) * log(v));
        acc_add(recips, 1.0 / (n - 1));
    }
    let dfd = (total - k) as f64;
    let sp = acc_total(pooled) / dfd;
    if is_inf(sp) { return fault_result(ERR_OVERFLOW); }
    let numer = dfd * log(sp) - acc_total(logs);
    let denom = 1.0 + (1.0 / (3.0 * (k - 1))) * (acc_total(recips) - 1.0 / dfd);
    var t = numer / denom;
    if is_bad(t) { return fault_result(ERR_OVERFLOW); }
    if t < 0.0 { t = 0.0; }
    return vr_finish_chi(t, k - 1);
}

// The test of Alexander and Govern for equal means of groups with unequal
// variances. Each group needs two values and a variance above 0.
// SciPy: alexandergovern(*samples).
fn alexandergovern(groups: const f64[:], sizes: const i64[:]) -> TestResult, i64 {
    let c = check_groups(groups, sizes);
    if c != 0 { return fault_result(c); }
    for n in sizes { if n < 2 { return fault_result(ERR_TOO_SHORT); } }
    let k = sizes.len;
    var means: f64[>..] = [];
    var se2: f64[>..] = [];
    var at: i64 = 0;
    for n in sizes {
        let part = groups[at..at + n];
        at += n;
        let m, cm = mean(part);
        if cm != 0 { return fault_result(cm); }
        let v, cv = variance_core(part, 1);
        if cv != 0 { return fault_result(cv); }
        let s2 = v / n;
        // SciPy gives NaN when the error is no larger than eps times the mean.
        if !(sqrt(s2) > EPS * abs(m)) { return fault_result(ERR_UNDEFINED); }
        means.push(m);
        se2.push(s2);
    }
    // The weights use the inverse of the squared errors, scaled by the
    // smallest error, so a tiny error cannot overflow the sum.
    var lo = se2[0];
    for s in se2 { if s < lo { lo = s; } }
    var inv = acc_new();
    for s in se2 { acc_add(inv, lo / s); }
    let invsum = acc_total(inv);
    var vw = acc_new();
    for i in 0..k { acc_add(vw, (lo / se2[i]) / invsum * means[i]); }
    let w = acc_total(vw);
    var tot = acc_new();
    for i in 0..k {
        let n = sizes[i];
        let t = (means[i] - w) / sqrt(se2[i]);
        let v = (n - 1) as f64;
        let a = v - 0.5;
        let b = 48.0 * a * a;
        let cc = sqrt(a * log1p_core(t * t / v));
        let c2 = cc * cc;
        let c3 = c2 * cc;
        let c4 = c2 * c2;
        let c5 = c4 * cc;
        let c7 = c5 * c2;
        let z = cc + (c3 + 3.0 * cc) / b - (4.0 * c7 + 33.0 * c5 + 240.0 * c3 + 855.0 * cc) / (10.0 * b * b + 8.0 * b * c4 + 1000.0 * b);
        acc_add(tot, z * z);
    }
    return vr_finish_chi(acc_total(tot), k - 1);
}

// log(1 + x) for x at or above 0.
fn vr_log1p(x: f64) -> f64 {
    let r, c = log1p(x);
    if c != 0 { return log(1.0 + x); }
    return r;
}

// ---- two-sample rank tests ----------------------------------------------

// The rank sum test of Wilcoxon. It uses no tie correction, like SciPy.
// SciPy: ranksums(x, y, alternative='two-sided').
fn ranksums(xs: const f64[:], ys: const f64[:], alt: i64 = ALT_TWO_SIDED) -> TestResult, i64 {
    let c = check_two(xs, ys, 1);
    if c != 0 { return fault_result(c); }
    if check_alt(alt) != 0 { return fault_result(ERR_DOMAIN); }
    var pooled: f64[>..] = [];
    concat2(pooled, xs, ys);
    var rk: f64[>..] = [];
    let cnt, tt, cr = rankdata_ties(rk, pooled[..], TIE_AVERAGE);
    if cr != 0 { return fault_result(cr); }
    var s = acc_new();
    for i in 0..xs.len { acc_add(s, rk[i]); }
    let n1 = xs.len as f64;
    let n2 = ys.len as f64;
    let expected = n1 * (n1 + n2 + 1.0) / 2.0;
    let z = (acc_total(s) - expected) / sqrt(n1 * n2 * (n1 + n2 + 1.0) / 12.0);
    return vr_finish_z(z, alt);
}

// The test of Mood for equal scales. A sample with ties uses the form of
// Mielke. The two samples hold three values or more in all.
// SciPy: mood(x, y, alternative='two-sided').
fn mood(xs: const f64[:], ys: const f64[:], alt: i64 = ALT_TWO_SIDED) -> TestResult, i64 {
    let c = check_two(xs, ys, 1);
    if c != 0 { return fault_result(c); }
    if xs.len + ys.len < 3 { return fault_result(ERR_TOO_SHORT); }
    if check_alt(alt) != 0 { return fault_result(ERR_DOMAIN); }
    let m = xs.len;
    let n = ys.len;
    let tot = m + n;
    var idx: Vi[>..] = [];
    for i in 0..m { idx.push(Vi { v: xs[i], i: i }); }
    for i in 0..n { idx.push(Vi { v: ys[i], i: m + i }); }
    stable_sort(idx) { a, b => a.v < b.v };
    var tied = false;
    for k in 1..tot { if idx[k].v == idx[k - 1].v { tied = true; } }
    let mf = m as f64;
    let nf = n as f64;
    let nt = tot as f64;
    let e0 = mf * (nt * nt - 1.0) / 12.0;
    var z: f64 = 0.0;
    if !tied {
        var s = acc_new();
        for k in 0..tot {
            if idx[k].i < m {
                let d = (k + 1) - (nt + 1.0) / 2.0;
                acc_add(s, d * d);
            }
        }
        let varm = mf * nf * (nt + 1.0) * (nt + 2.0) * (nt - 2.0) / 180.0;
        z = (acc_total(s) - e0) / sqrt(varm);
        return vr_finish_z(z, alt);
    }
    // The tie blocks. The phi of a block is the mean of the squared
    // distances of its ranks from the middle rank.
    let cm = (nt + 1.0) / 2.0;
    var tsum = acc_new();
    var vsum = acc_new();
    var seen: i64 = 0;
    var k: i64 = 0;
    while k < tot {
        var j = k + 1;
        while j < tot && idx[j].v == idx[k].v { j += 1; }
        let t = (j - k) as f64;
        var a: f64 = 0.0;
        for q in k..j { if idx[q].i < m { a += 1.0; } }
        let mid_rank = seen + (t + 1.0) / 2.0;
        let phi = (mid_rank - cm) * (mid_rank - cm) + (t * t - 1.0) / 12.0;
        acc_add(tsum, a * phi);
        let s_i = (seen + (j - k)) as f64;
        let s_p = seen as f64;
        let g = nt - s_i - s_p;
        acc_add(vsum, t * (t * t - 1.0) * (t * t - 4.0 + 15.0 * g * g));
        seen += j - k;
        k = j;
    }
    let varm = mf * nf * (nt + 1.0) * (nt * nt - 4.0) / 180.0 - mf * nf / (180.0 * nt * (nt - 1.0)) * acc_total(vsum);
    if !(varm > 0.0) { return fault_result(ERR_UNDEFINED); }
    z = (acc_total(tsum) - e0) / sqrt(varm);
    return vr_finish_z(z, alt);
}

// The test of Brunner and Munzel. Each sample needs two values. The t form
// holds the Welch-style degrees of freedom in df. The normal form has no df.
// SciPy: brunnermunzel(x, y, alternative='two-sided', distribution='t').
fn brunnermunzel(xs: const f64[:], ys: const f64[:], alt: i64 = ALT_TWO_SIDED, use_t: bool = true) -> TestResult, i64 {
    let c = check_two(xs, ys, 2);
    if c != 0 { return fault_result(c); }
    if check_alt(alt) != 0 { return fault_result(ERR_DOMAIN); }
    let nx = xs.len;
    let ny = ys.len;
    var pooled: f64[>..] = [];
    concat2(pooled, xs, ys);
    var rc: f64[>..] = [];
    var rx: f64[>..] = [];
    var ry: f64[>..] = [];
    let n1, cr1 = rankdata(rc, pooled[..], TIE_AVERAGE);
    if cr1 != 0 { return fault_result(cr1); }
    let n2, cr2 = rankdata(rx, xs, TIE_AVERAGE);
    if cr2 != 0 { return fault_result(cr2); }
    let n3, cr3 = rankdata(ry, ys, TIE_AVERAGE);
    if cr3 != 0 { return fault_result(cr3); }
    let mcx, e1 = mean(rc[0..nx]);
    let mcy, e2 = mean(rc[nx..nx + ny]);
    let mrx, e3 = mean(rx[..]);
    let mry, e4 = mean(ry[..]);
    if e1 != 0 || e2 != 0 || e3 != 0 || e4 != 0 { return fault_result(ERR_UNDEFINED); }
    var sxa = acc_new();
    for i in 0..nx {
        let t = (rc[i] - rx[i]) - (mcx - mrx);
        acc_add(sxa, t * t);
    }
    var sya = acc_new();
    for i in 0..ny {
        let t = (rc[nx + i] - ry[i]) - (mcy - mry);
        acc_add(sya, t * t);
    }
    let fx = nx as f64;
    let fy = ny as f64;
    let a = acc_total(sxa) * fx / (fx - 1.0);
    let b = acc_total(sya) * fy / (fy - 1.0);
    let s = a + b;
    if !(s > 0.0) { return fault_result(ERR_UNDEFINED); }
    let w = fx * fy * (mcy - mcx) / ((fx + fy) * sqrt(s));
    if !use_t {
        let p, cp = norm_z_p(-w, alt);
        if cp != 0 { return fault_result(cp); }
        return good_result(w, p, nan(), nan());
    }
    let ra = a / s;
    let rb = b / s;
    let df = 1.0 / (ra * ra / (fx - 1.0) + rb * rb / (fy - 1.0));
    let p, cp = t_pvalue(-w, df, alt);
    if cp != 0 { return fault_result(cp); }
    return good_result(w, p, df, nan());
}

// ---- ansari -------------------------------------------------------------

// Fills out with the counts of the sums of j of the symmetric ranks of N
// values, for j from 0 to k0. The row j starts at j * (S + 1). The
// symmetric rank v = min(r, N + 1 - r) occurs twice, or once for the
// middle rank of an odd N. The call gives S, the largest sum. The counts
// reach 1e31, so a double holds them with a small relative error. No
// count is ever subtracted, so the error stays near 1e-14.
fn vr_ansari_counts(out: f64[>..]&, k0: i64, tot: i64) -> i64 {
    let vmax = (tot + 1) / 2;
    let top = k0 * vmax;
    let w = top + 1;
    out.resize(0);
    for i in 0..(k0 + 1) * w { out.push(0.0); }
    out[0] = 1.0;
    for v in 1..vmax + 1 {
        let two = !(tot % 2 == 1 && v == vmax);
        var j = k0;
        while j >= 1 {
            for s in v..w {
                var add = out[(j - 1) * w + s - v];
                if two {
                    add = add * 2.0;
                    if j >= 2 && s >= 2 * v { add += out[(j - 2) * w + s - 2 * v]; }
                }
                out[j * w + s] += add;
            }
            j -= 1;
        }
    }
    return top;
}

// The exact p-value of the Ansari-Bradley statistic ab of the first
// sample, for samples of nx and ny values with no ties.
fn vr_ansari_exact(nx: i64, ny: i64, ab: i64, alt: i64) -> f64 {
    let tot = nx + ny;
    let k0 = min(nx, ny);
    var counts: f64[>..] = [];
    let top = vr_ansari_counts(counts, k0, tot);
    let w = top + 1;
    var sum_all: i64 = 0;
    for r in 1..tot + 1 { sum_all += min(r, tot + 1 - r); }
    var le = acc_new();
    var ge = acc_new();
    var total = acc_new();
    // The first sample is the smaller one: its statistic is the sum. In the
    // other case the sum of the second sample is sum_all minus the statistic.
    let small_x = nx <= ny;
    let pivot = if small_x { ab } else { sum_all - ab };
    for s in 0..w {
        let f = counts[k0 * w + s];
        acc_add(total, f);
        if small_x {
            if s <= pivot { acc_add(le, f); }
            if s >= pivot { acc_add(ge, f); }
        } else {
            if s >= pivot { acc_add(le, f); }
            if s <= pivot { acc_add(ge, f); }
        }
    }
    let t = acc_total(total);
    let cdf = acc_total(le) / t;
    let sf = acc_total(ge) / t;
    if alt == ALT_GREATER { return cdf; }
    if alt == ALT_LESS { return sf; }
    return min(1.0, 2.0 * min(cdf, sf));
}

// The Ansari-Bradley test of equal scales. The method is MW_AUTO, MW_EXACT
// or MW_ASYMPTOTIC. The auto form is exact for two samples of at most 54
// values with no ties. The exact form refuses ties and tables above
// VR_ANSARI_CELLS cells. SciPy: ansari(x, y, alternative, method='auto').
fn ansari(xs: const f64[:], ys: const f64[:], alt: i64 = ALT_TWO_SIDED, method: i64 = MW_AUTO) -> TestResult, i64 {
    let c = check_two(xs, ys, 1);
    if c != 0 { return fault_result(c); }
    if check_alt(alt) != 0 || method < 0 || method > MW_ASYMPTOTIC { return fault_result(ERR_DOMAIN); }
    let n = xs.len;
    let m = ys.len;
    let tot = n + m;
    var pooled: f64[>..] = [];
    concat2(pooled, xs, ys);
    var rk: f64[>..] = [];
    let cnt, tt, cr = rankdata_ties(rk, pooled[..], TIE_AVERAGE);
    if cr != 0 { return fault_result(cr); }
    let nt = tot as f64;
    var sym: f64[>..] = [];
    for r in rk { sym.push(min(r, nt + 1.0 - r)); }
    var sab = acc_new();
    for i in 0..n { acc_add(sab, sym[i]); }
    let ab = acc_total(sab);
    var exact = false;
    if method == MW_EXACT {
        if tt > 0.0 { return fault_result(ERR_DOMAIN); }
        let k0 = min(n, m);
        if (k0 + 1) * (k0 * ((tot + 1) / 2) + 1) > VR_ANSARI_CELLS { return fault_result(ERR_OVERFLOW); }
        exact = true;
    }
    if method == MW_AUTO { exact = tt == 0.0 && n < 55 && m < 55; }
    if exact {
        let p = vr_ansari_exact(n, m, round(ab) as i64, alt);
        return good_result(ab, p, nan(), nan());
    }
    var mn: f64 = 0.0;
    var varab: f64 = 0.0;
    let fn_ = n as f64;
    let fm = m as f64;
    if tot % 2 == 1 { mn = fn_ * (nt + 1.0) * (nt + 1.0) / 4.0 / nt; } else { mn = fn_ * (nt + 2.0) / 4.0; }
    if tt > 0.0 {
        var fac = acc_new();
        for v in sym { acc_add(fac, v * v); }
        let f = acc_total(fac);
        if tot % 2 == 1 {
            varab = fm * fn_ * (16.0 * nt * f - (nt + 1.0) * (nt + 1.0) * (nt + 1.0) * (nt + 1.0)) / (16.0 * nt * nt * (nt - 1.0));
        } else {
            varab = fm * fn_ * (16.0 * f - nt * (nt + 2.0) * (nt + 2.0)) / (16.0 * nt * (nt - 1.0));
        }
    } else if tot % 2 == 1 {
        varab = fn_ * fm * (nt + 1.0) * (3.0 + nt * nt) / (48.0 * nt * nt);
    } else {
        varab = fm * fn_ * (nt + 2.0) * (nt - 2.0) / 48.0 / (nt - 1.0);
    }
    if !(varab > 0.0) { return fault_result(ERR_UNDEFINED); }
    let z = (mn - ab) / sqrt(varab);
    if is_bad(z) { return fault_result(ERR_UNDEFINED); }
    let p, cp = norm_z_p(z, alt);
    if cp != 0 { return fault_result(cp); }
    return good_result(ab, p, nan(), nan());
}

// ---- tests of several samples -------------------------------------------

// The Friedman test. The data hold one row for each block and one column
// for each sample, row by row: the value of block r and sample j sits at
// data[r * columns + j]. SciPy takes the samples as separate arguments;
// to convert, write sample j into the column j. It needs three columns.
// SciPy: friedmanchisquare(*samples).
fn friedmanchisquare(data: const f64[:], columns: i64) -> TestResult, i64 {
    let c = check_sample(data, 1);
    if c != 0 { return fault_result(c); }
    if columns > 0 && data.len % columns != 0 { return fault_result(ERR_LENGTH); }
    if columns < 1 { return fault_result(ERR_DOMAIN); }
    if columns < 3 { return fault_result(ERR_TOO_SHORT); }
    let k = columns;
    let rows = data.len / k;
    var colsum: f64[>..] = [];
    for j in 0..k { colsum.push(0.0); }
    var ties = acc_new();
    var rk: f64[>..] = [];
    for r in 0..rows {
        let cnt, tt, cr = rankdata_ties(rk, data[r * k..(r + 1) * k], TIE_AVERAGE);
        if cr != 0 { return fault_result(cr); }
        acc_add(ties, tt);
        for j in 0..k { colsum[j] += rk[j]; }
    }
    let kf = k as f64;
    let nf = rows as f64;
    let corr = 1.0 - acc_total(ties) / (kf * (kf * kf - 1.0) * nf);
    if !(corr > 0.0) { return fault_result(ERR_UNDEFINED); }
    var ssbn = acc_new();
    for j in 0..k { acc_add(ssbn, colsum[j] * colsum[j]); }
    var stat = (12.0 / (kf * nf * (kf + 1.0)) * acc_total(ssbn) - 3.0 * nf * (kf + 1.0)) / corr;
    if stat < 0.0 { stat = 0.0; }
    return vr_finish_chi(stat, kf - 1.0);
}

// The median test. It counts the values above and below the grand median
// in each group and tests the 2 by k table with the chi-square test. The
// table goes into out, row by row: the row 0 holds the counts above the
// median and the row 1 the counts below it. The ties argument is a MEDIAN_
// constant. The result holds the degrees of freedom in df and the grand
// median in df2. A table with an empty row or column gives the undefined
// fault, where SciPy raises. SciPy: median_test(*samples, ties='below',
// correction=True, lambda_=1).
fn median_test(out: f64[>..]&, groups: const f64[:], sizes: const i64[:], ties: i64 = MEDIAN_BELOW, correction: bool = true) -> TestResult, i64 {
    out.resize(0);
    let c = check_groups(groups, sizes);
    if c != 0 { return fault_result(c); }
    if ties < MEDIAN_BELOW || ties > MEDIAN_IGNORE { return fault_result(ERR_DOMAIN); }
    let gm, cg = median(groups);
    if cg != 0 { return fault_result(cg); }
    let k = sizes.len;
    for i in 0..2 * k { out.push(0.0); }
    var at: i64 = 0;
    for g in 0..k {
        let n = sizes[g];
        for i in at..at + n {
            let x = groups[i];
            if x > gm { out[g] += 1.0; }
            else if x < gm { out[k + g] += 1.0; }
            else if ties == MEDIAN_ABOVE { out[g] += 1.0; }
            else if ties == MEDIAN_BELOW { out[k + g] += 1.0; }
        }
        at += n;
    }
    let r, cr = chi2_contingency(out[..], 2, k, correction);
    if cr != 0 {
        out.resize(0);
        return fault_result(cr);
    }
    return good_result(r.statistic, r.pvalue, r.df, gm);
}

// The smallest count k from 0 to n with a binomial upper tail of at most
// alpha. This is the inverse survival function of SciPy.
fn vr_binom_isf(alpha: f64, n: f64, p: f64) -> f64, i64 {
    var lo: f64 = 0.0;
    var hi = n;
    var steps: i64 = 0;
    while lo < hi {
        steps += 1;
        if steps > 64 { return nan(), ERR_FIT; }
        let mid_k = lo + floor((hi - lo) / 2.0);
        let s, c = binom_sf(mid_k, n, p);
        if c != 0 { return nan(), c; }
        if s <= alpha { hi = mid_k; } else { lo = mid_k + 1.0; }
    }
    return lo, 0;
}

// The quantile test of Conover. It tests that q is the p-th quantile of
// the population. The statistic is the count T1 of values at most q or
// the count T2 of values below q. The type of the statistic, 1 or 2, goes
// into df. The confidence interval of the quantile goes into out as the
// lower and the upper bound. A missing bound is -inf or +inf, where SciPy
// gives -inf, +inf or NaN. The bounds are values of the sample.
// SciPy: quantile_test(x, q=0, p=0.5, alternative='two-sided') and the
// method confidence_interval(confidence_level).
fn quantile_test(out: f64[>..]&, xs: const f64[:], q: f64 = 0.0, p: f64 = 0.5, alt: i64 = ALT_TWO_SIDED, level: f64 = 0.95) -> TestResult, i64 {
    out.resize(0);
    var c = check_sample(xs, 1);
    if c != 0 { return fault_result(c); }
    if is_nan(q) { return fault_result(ERR_NAN); }
    c = check_scalar(p);
    if c != 0 { return fault_result(c); }
    c = check_level(level);
    if c != 0 { return fault_result(c); }
    if p <= 0.0 || p >= 1.0 || check_alt(alt) != 0 { return fault_result(ERR_DOMAIN); }
    let n = xs.len;
    let nf = n as f64;
    var t1: f64 = 0.0;
    var t2: f64 = 0.0;
    for x in xs {
        if x <= q { t1 += 1.0; }
        if x < q { t2 += 1.0; }
    }
    let pg, c1 = binom_cdf(t1, nf, p);
    let pl, c2 = tail_sf(0, t2 - 1.0, nf, p, 0.0);
    if c1 != 0 { return fault_result(c1); }
    if c2 != 0 { return fault_result(c2); }
    var stat = t1;
    var kind: f64 = 1.0;
    var pv = pg;
    if alt == ALT_LESS {
        stat = t2;
        kind = 2.0;
        pv = pl;
    } else if alt == ALT_TWO_SIDED {
        pv = 2.0 * min(pg, pl);
        if pg > pl {
            stat = t2;
            kind = 2.0;
        }
    }
    // The bounds come from the sorted sample. The index of the lower bound
    // is ppf(alpha) - 1 and the index of the upper bound is isf(alpha).
    let s = sort_copy(xs);
    var alpha = 1.0 - level;
    if alt == ALT_TWO_SIDED { alpha = (1.0 - level) / 2.0; }
    var low = -(1.0 / 0.0);
    var high = (1.0 / 0.0);
    if alt != ALT_LESS {
        let lp, cl = binom_ppf(alpha, nf, p);
        if cl != 0 { return fault_result(cl); }
        if lp - 1.0 >= 0.0 { low = s[(lp - 1.0) as i64]; }
    }
    if alt != ALT_GREATER {
        let hp, ch = vr_binom_isf(alpha, nf, p);
        if ch != 0 { return fault_result(ch); }
        if hp < nf { high = s[hp as i64]; }
    }
    let r, cr = good_result(stat, pv, kind, nan());
    if cr != 0 { return r, cr; }
    out.push(low);
    out.push(high);
    return r, 0;
}

// ---- the t test from summary values --------------------------------------

// The t test of two samples from their means, deviations and counts. The
// counts are whole numbers from 2 up. A deviation below 0 gives the domain
// fault. Equal summary values with no spread give the undefined fault, and
// so does a difference of the means with no spread, where SciPy gives an
// infinite statistic. SciPy: ttest_ind_from_stats(mean1, std1, nobs1,
// mean2, std2, nobs2, equal_var=True, alternative='two-sided').
fn ttest_ind_from_stats(mean1: f64, sd1: f64, n1: f64, mean2: f64, sd2: f64, n2: f64, equal_var: bool = true, alt: i64 = ALT_TWO_SIDED) -> TestResult, i64 {
    var c = check_scalar(mean1);
    if c == 0 { c = check_scalar(sd1); }
    if c == 0 { c = check_scalar(n1); }
    if c == 0 { c = check_scalar(mean2); }
    if c == 0 { c = check_scalar(sd2); }
    if c == 0 { c = check_scalar(n2); }
    if c != 0 { return fault_result(c); }
    if chk_count(n1) != 0 || chk_count(n2) != 0 || n1 < 2.0 || n2 < 2.0 { return fault_result(ERR_DOMAIN); }
    if sd1 < 0.0 || sd2 < 0.0 || check_alt(alt) != 0 { return fault_result(ERR_DOMAIN); }
    let d = mean1 - mean2;
    if is_inf(d) { return fault_result(ERR_OVERFLOW); }
    if equal_var {
        let df = n1 + n2 - 2.0;
        let pooled = sd1 * sd1 * ((n1 - 1.0) / df) + sd2 * sd2 * ((n2 - 1.0) / df);
        let se = sqrt(pooled * (1.0 / n1 + 1.0 / n2));
        if is_inf(se) { return fault_result(ERR_OVERFLOW); }
        if se == 0.0 { return fault_result(ERR_UNDEFINED); }
        return finish_t(d / se, df, alt);
    }
    let a = sd1 * sd1 / n1;
    let b = sd2 * sd2 / n2;
    let s = a + b;
    if is_inf(s) { return fault_result(ERR_OVERFLOW); }
    if s == 0.0 { return fault_result(ERR_UNDEFINED); }
    let ra = a / s;
    let rb = b / s;
    let df = 1.0 / (ra * ra / (n1 - 1.0) + rb * rb / (n2 - 1.0));
    return finish_t(d / sqrt(s), df, alt);
}
