// The one sided and the two sided Kolmogorov-Smirnov statistics and the
// limit of the statistic. The names and the parameters follow SciPy
// 1.18.1. The sample size n must not exceed the constant below.
namespace statistics;

import std;
import math;

const KS_MAX_N: f64 = 1000.0;       // the largest sample size of ksone and kstwo

// ---- ksone: the one sided Kolmogorov-Smirnov statistic --------------------

fn ks_check(n: f64) -> i64 {
    let c = check_scalar(n);
    if c != OK { return c; }
    if n < 1.0 || n != floor(n) || n > KS_MAX_N { return ERR_DOMAIN; }
    return OK;
}

// The survival function of the sum of Birnbaum and Tingey, and the
// density as its negative derivative. The point lies from 0 to 1.
fn ksone_both(x: f64, n: f64) -> f64, f64 {
    var sf: f64 = 0.0;
    var dd: f64 = 0.0;
    let jmax = floor(n * (1.0 - x));
    var j: f64 = 0.0;
    while j <= jmax {
        let lc, cc = lchoose(n, j);
        if cc != OK { return nan(), nan(); }
        let a = x + j / n;
        let b = 1.0 - x - j / n;
        if j == 0.0 {
            sf += pow(1.0 - x, n);
            dd += -n * pow(1.0 - x, n - 1.0);
        } else if b > 0.0 {
            let l = exp(lc + (j - 1.0) * log(a) + (n - j) * log(b));
            sf += x * l;
            dd += l * (1.0 + x * (j - 1.0) / a - x * (n - j) / b);
        } else if n - j == 0.0 {
            let l = exp(lc + (j - 1.0) * log(a));
            sf += x * l;
            dd += l * (1.0 + x * (j - 1.0) / a);
        }
        j += 1.0;
    }
    return sf, -dd;
}

fn ksone_pdf(x: f64, n: f64) -> f64, i64 {
    let k = tl_px(x, ks_check(n));
    if k != OK { return nan(), k; }
    if x <= 0.0 || x >= 1.0 { return 0.0, OK; }
    let sf, d = ksone_both(x, n);
    return tl_out(max(d, 0.0));
}

fn ksone_cdf_value(x: f64, n: f64) -> f64 {
    if x <= 0.0 { return 0.0; }
    if x >= 1.0 { return 1.0; }
    let sf, d = ksone_both(x, n);
    return min(max(1.0 - sf, 0.0), 1.0);
}

fn ksone_cdf(x: f64, n: f64) -> f64, i64 {
    let k = tl_px(x, ks_check(n));
    if k != OK { return nan(), k; }
    return tl_out(ksone_cdf_value(x, n));
}

fn ksone_ppf(q: f64, n: f64) -> f64, i64 {
    let k = tl_pq(q, ks_check(n));
    if k != OK { return nan(), k; }
    if q == 0.0 { return 0.0, OK; }
    if q == 1.0 { return 1.0, OK; }
    let x, cx = brent_root(0.0, 1.0, 1e-14, 4e-16) { v => ksone_cdf_value(v, n) - q };
    return x, cx;
}

// ---- kstwo: the two sided Kolmogorov-Smirnov statistic --------------------

// The product c = a * b of two square matrices of the side m.
fn ks_mul(a: const f64[:], b: const f64[:], c: f64[>..]&, m: i64) {
    for i in 0..m {
        for j in 0..m {
            var s: f64 = 0.0;
            for l in 0..m { s += a[i * m + l] * b[l * m + j]; }
            c[i * m + j] = s;
        }
    }
}

// Scale a matrix by 1e-140 when its largest entry grows too large. The
// call returns the added power of ten.
fn ks_rescale(v: f64[>..]&) -> i64 {
    var top: f64 = 0.0;
    for i in 0..v.len { top = max(top, abs(v[i])); }
    if top < 1e140 { return 0; }
    for i in 0..v.len { v[i] = v[i] * 1e-140; }
    return 140;
}

// The cumulative function of the statistic, by the matrix method.
fn kstwo_cdf_value(n: f64, d: f64) -> f64 {
    if d <= 0.5 / n { return 0.0; }
    if d >= 1.0 { return 1.0; }
    let s = d * d * n;
    if s > 7.24 || (s > 3.76 && n > 99.0) {
        return 1.0 - 2.0 * exp(-(2.000071 + 0.331 / sqrt(n) + 1.409 / n) * s);
    }
    let k = floor(n * d) as i64 + 1;
    let m = 2 * k - 1;
    let h = k - n * d;
    var hm: f64[>..] = [];
    for i in 0..m {
        for j in 0..m {
            if i - j + 1 < 0 { hm.push(0.0); } else { hm.push(1.0); }
        }
    }
    // The powers of h from 1 to m, one multiply each, so the two border
    // rows share one table.
    var hp: f64[>..] = [];
    var t: f64 = 1.0;
    for i in 0..m { t = t * h; hp.push(t); }
    for i in 0..m {
        hm[i * m] = hm[i * m] - hp[i];
        hm[(m - 1) * m + i] = hm[(m - 1) * m + i] - hp[m - 1 - i];
    }
    if 2.0 * h - 1.0 > 0.0 { hm[(m - 1) * m] = hm[(m - 1) * m] + pow(2.0 * h - 1.0, m); }
    for i in 0..m {
        for j in 0..m {
            if i - j + 1 > 0 {
                for g in 1..(i - j + 2) { hm[i * m + j] = hm[i * m + j] / g; }
            }
        }
    }
    // The power of the matrix by repeated squaring, with a power of ten
    // kept apart.
    var r: f64[>..] = [];
    var p: f64[>..] = [];
    var t: f64[>..] = [];
    for i in 0..m * m {
        r.push(0.0);
        t.push(0.0);
        p.push(hm[i]);
    }
    for i in 0..m { r[i * m + i] = 1.0; }
    var er: i64 = 0;
    var ep: i64 = 0;
    var e = n as i64;
    while e > 0 {
        if e % 2 == 1 {
            ks_mul(r[..], p[..], t, m);
            for i in 0..m * m { r[i] = t[i]; }
            er += ep;
            er += ks_rescale(r);
        }
        e = e / 2;
        if e > 0 {
            ks_mul(p[..], p[..], t, m);
            for i in 0..m * m { p[i] = t[i]; }
            ep = ep * 2;
            ep += ks_rescale(p);
        }
    }
    var v = r[(k - 1) * m + k - 1];
    let ni = n as i64;
    for i in 1..(ni + 1) {
        v = v * i / n;
        if v < 1e-140 {
            v = v * 1e140;
            er -= 140;
        }
    }
    return min(max(v * pow(10.0, er as f64), 0.0), 1.0);
}

fn kstwo_cdf(x: f64, n: f64) -> f64, i64 {
    let k = tl_px(x, ks_check(n));
    if k != OK { return nan(), k; }
    return tl_out(kstwo_cdf_value(n, x));
}

// The density is the derivative of the cumulative function. The call
// takes a central difference.
fn kstwo_pdf(x: f64, n: f64) -> f64, i64 {
    let k = tl_px(x, ks_check(n));
    if k != OK { return nan(), k; }
    if x <= 0.0 || x >= 1.0 { return 0.0, OK; }
    let h = 1e-6;
    let hi = min(x + h, 1.0);
    let lo = max(x - h, 0.0);
    return tl_out(max((kstwo_cdf_value(n, hi) - kstwo_cdf_value(n, lo)) / (hi - lo), 0.0));
}

fn kstwo_ppf(q: f64, n: f64) -> f64, i64 {
    let k = tl_pq(q, ks_check(n));
    if k != OK { return nan(), k; }
    if q == 0.0 { return 0.5 / n, OK; }
    if q == 1.0 { return 1.0, OK; }
    let x, cx = brent_root(0.5 / n, 1.0, 1e-14, 4e-16) { v => kstwo_cdf_value(n, v) - q };
    return x, cx;
}
