// The Wilcoxon test of the package. Each call returns a TestResult and
// then a fault code. A fault gives NaN in every field of the result.
namespace statistics;

import std;
import math;

fn subset_counts(out: f64[>..]&, w: const i64[:]) {
    out.resize(0);
    var total: i64 = 0;
    for x in w { total += x; }
    for s in 0..total + 1 { out.push(0.0); }
    out[0] = 1.0;
    var top: i64 = 0;
    for x in w {
        top += x;
        var s = top;
        while s >= x {
            out[s] += out[s - x];
            s -= 1;
        }
    }
}

// The Wilcoxon test of a list of differences.
fn wilcoxon(diffs: const f64[:], alt: i64 = ALT_TWO_SIDED, zero_method: i64 = ZERO_WILCOX, correction: bool = false) -> TestResult, i64 {
    let c = check_sample(diffs, 1);
    if c != 0 { return fault_result(c); }
    if check_alt(alt) != 0 || zero_method < 0 || zero_method > ZERO_ZSPLIT { return fault_result(ERR_DOMAIN); }
    var n_zero = 0;
    for v in diffs {
        if v == 0.0 { n_zero += 1; }
    }
    if zero_method == ZERO_WILCOX && n_zero > 0 {
        // The zero differences leave the sample. The copy of the sample
        // is made only then, so the common case holds no copy.
        var kept: f64[>..] = [];
        for v in diffs {
            if v != 0.0 { kept.push(v); }
        }
        return wilcoxon_core(kept[..], diffs.len, n_zero, alt, zero_method, correction);
    }
    return wilcoxon_core(diffs, diffs.len, n_zero, alt, zero_method, correction);
}

// The Wilcoxon test of two paired samples.
fn wilcoxon_paired(xs: const f64[:], ys: const f64[:], alt: i64 = ALT_TWO_SIDED, zero_method: i64 = ZERO_WILCOX, correction: bool = false) -> TestResult, i64 {
    let c = check_pair(xs, ys, 1);
    if c != 0 { return fault_result(c); }
    if check_alt(alt) != 0 || zero_method < 0 || zero_method > ZERO_ZSPLIT { return fault_result(ERR_DOMAIN); }
    var d: f64[>..] = [];
    for i in 0..xs.len {
        let v = xs[i] - ys[i];
        if is_inf(v) { return fault_result(ERR_OVERFLOW); }
        d.push(v);
    }
    return wilcoxon(d[..], alt, zero_method, correction);
}

fn wilcoxon_core(d: const f64[:], n_diff: i64, n_zero: i64, alt: i64, zero_method: i64, correction: bool) -> TestResult, i64 {
    if d.len == 0 || n_zero == n_diff { return fault_result(ERR_UNDEFINED); }
    var absd: f64[>..] = [];
    for v in d { absd.push(abs(v)); }
    var r: f64[>..] = [];
    let cnt, t_all, cr = rankdata_ties(r, absd, TIE_AVERAGE);
    if cr != 0 { return fault_result(cr); }
    var rp: f64 = 0.0;
    var rm: f64 = 0.0;
    var rz: f64 = 0.0;
    for i in 0..d.len {
        if d[i] > 0.0 { rp += r[i]; }
        if d[i] < 0.0 { rm += r[i]; }
        if d[i] == 0.0 { rz += r[i]; }
    }
    var fixed: f64 = 0.0;
    if zero_method == ZERO_ZSPLIT {
        rp += rz / 2.0;
        rm += rz / 2.0;
        fixed = rz / 2.0;
    }
    var stat = rp;
    if alt == ALT_TWO_SIDED { stat = min(rp, rm); }
    let count = d.len;
    // The groups of equal absolute differences give the same groups as the
    // ranks. The zero differences form one group of their own. The sample
    // holds them only when the zero method keeps them.
    let keeps_zero = zero_method != ZERO_WILCOX;
    var t_nz = t_all;
    if keeps_zero && n_zero > 0 {
        let z = n_zero as f64;
        t_nz = t_all - (z * z * z - z);
    }
    var ties = t_nz > 0.0;
    var zeros = n_zero > 0;
    var exact = false;
    if n_diff <= MAX_EXACT_WILCOXON && !ties && !zeros { exact = true; }
    if n_diff <= MAX_ENUM_WILCOXON { exact = true; }
    if exact {
        var weights: i64[>..] = [];
        for i in 0..d.len {
            if d[i] != 0.0 { weights.push(round(2.0 * r[i]) as i64); }
        }
        var cts: f64[>..] = [];
        subset_counts(cts, weights);
        var total: f64 = 1.0;
        for i in 0..weights.len { total *= 2.0; }
        let at = round(2.0 * (rp - fixed)) as i64;
        var le = acc_new();
        var ge = acc_new();
        for s in 0..cts.len {
            if s <= at { acc_add(le, cts[s]); }
            if s >= at { acc_add(ge, cts[s]); }
        }
        let pl = acc_total(le) / total;
        let pg = acc_total(ge) / total;
        var p = pg;
        if alt == ALT_LESS { p = pl; }
        if alt == ALT_TWO_SIDED { p = 2.0 * min(pl, pg); }
        return good_result(stat, p, nan(), nan());
    }
    let cf = count as f64;
    var mn = cf * (cf + 1.0) * 0.25;
    var se = cf * (cf + 1.0) * (2.0 * cf + 1.0);
    if zero_method == ZERO_PRATT {
        let z = n_zero as f64;
        mn -= z * (z + 1.0) * 0.25;
        se -= z * (z + 1.0) * (2.0 * z + 1.0);
    }
    // The Pratt method corrects for the ties of the differences that are
    // not zero. The other methods correct for the ties of all differences.
    var t_corr = t_all;
    if zero_method == ZERO_PRATT && keeps_zero { t_corr = t_nz; }
    se -= 0.5 * t_corr;
    if se <= 0.0 { return fault_result(ERR_UNDEFINED); }
    se = sqrt(se / 24.0);
    var corr: f64 = 0.0;
    if correction {
        if alt == ALT_TWO_SIDED {
            if rp - mn > 0.0 { corr = 0.5; }
            if rp - mn < 0.0 { corr = -0.5; }
        }
        if alt == ALT_LESS { corr = -0.5; }
        if alt == ALT_GREATER { corr = 0.5; }
    }
    let z = (rp - mn - corr) / se;
    var p: f64 = 0.0;
    if alt == ALT_TWO_SIDED {
        let sf, cs = norm_sf(abs(z));
        if cs != 0 { return fault_result(cs); }
        p = 2.0 * sf;
    }
    if alt == ALT_GREATER {
        let sf, cs = norm_sf(z);
        if cs != 0 { return fault_result(cs); }
        p = sf;
    }
    if alt == ALT_LESS {
        let cd, cs = norm_cdf(z);
        if cs != 0 { return fault_result(cs); }
        p = cd;
    }
    return good_result(stat, p, nan(), nan());
}

fn kruskal(groups: const f64[:], sizes: const i64[:]) -> TestResult, i64 {
    let c = check_groups(groups, sizes);
    if c != 0 { return fault_result(c); }
    var rk: f64[>..] = [];
    let cnt, tt, cr = rankdata_ties(rk, groups, TIE_AVERAGE);
    if cr != 0 { return fault_result(cr); }
    let nt = groups.len as f64;
    let corr = 1.0 - tt / (nt * nt * nt - nt);
    if corr <= 0.0 { return fault_result(ERR_UNDEFINED); }
    var h: f64 = 0.0;
    var at: i64 = 0;
    for g in sizes {
        var rs: f64 = 0.0;
        for i in at..at + g { rs += rk[i]; }
        at += g;
        h += rs * rs / g;
    }
    h = 12.0 / (nt * (nt + 1.0)) * h - 3.0 * (nt + 1.0);
    h = h / corr;
    if h < 0.0 { h = 0.0; }
    let df = (sizes.len - 1) as f64;
    let p, cp = chi2_sf(h, df);
    if cp != 0 { return fault_result(cp); }
    return good_result(h, p, df, nan());
}
