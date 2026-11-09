// The association of two samples: the covariance, the correlation and the
// weighted forms.
//
// The sum of products is exact. The product of two values is the pair that
// Dekker and Veltkamp give. Both parts go into the exact sum of v0.1.0.
// The numerator of the covariance is n * S_xy - S_x * S_y. It is the
// difference of two exact products of sums. The digits of the result
// survive a large offset.
namespace statistics;

import std;
import math;

// The split constant of Dekker and Veltkamp. It is 2^27 + 1.
const SPLIT: f64 = 134217729.0;

// The exact product of two values. The result is the pair p and e with
// p + e equal to a * b. The split needs no overflow and no underflow.
fn two_prod(a: f64, b: f64) -> f64, f64 {
    let p = a * b;
    let ca = SPLIT * a;
    let ah = ca - (ca - a);
    let al = a - ah;
    let cb = SPLIT * b;
    let bh = cb - (cb - b);
    let bl = b - bh;
    let e = ((ah * bh - p) + ah * bl + al * bh) + al * bl;
    return p, e;
}

// The exact product of two sums, added into an accumulator. The product
// of two sums is the sum of the products of all pairs of partials. Each
// pair gives one two product. The sign gives the sum or its opposite.
fn acc_add_product(a: Acc&, x: Acc&, y: Acc&, sign: f64) {
    for i in 0..x.n {
        for j in 0..y.n {
            let p, e = two_prod(x.p[i], y.p[j]);
            acc_add(a, sign * p);
            acc_add(a, sign * e);
        }
    }
}

// An accumulator that holds one value. The value must be exact.
fn acc_one(v: f64) -> Acc {
    var a = acc_new();
    acc_add(a, v);
    return a;
}

// The scale of one column. The largest magnitude becomes at most 2^400, so
// that the products and the sums stay finite. A column of small values is
// scaled up, so that the products stay normal. Two columns of very
// different magnitude need two scales, so that the smaller column does not
// fall to zero. The result is the exponent of the power of two.
fn column_scale(xs: const f64[:]) -> i64 {
    var m: f64 = 0.0;
    for x in xs {
        let a = abs(x);
        if a > m { m = a; }
    }
    if m == 0.0 { return 0; }
    let e = floor(log2(m));
    if e > 400.0 { return 400 - e as i64; }
    if e < -300.0 { return -300 - e as i64; }
    return 0;
}

// The value times a power of two. The value is split into a part and an
// exponent, so that the result stays correct when the exponent leaves the
// range of the type.
fn scale_pow2(v: f64, k: i64) -> f64 {
    if v == 0.0 { return 0.0; }
    if is_nan(v) || is_inf(v) { return v; }
    var m = v;
    var e: i64 = 0;
    while abs(m) >= 2.0 { m = m / 2.0; e += 1; }
    while abs(m) < 1.0 { m = m * 2.0; e -= 1; }
    return m * pow(2.0, (e + k) as f64);
}

// The exact sums and the exact numerators of a pair of scaled columns.
struct Cross {
    n_xy: f64,
    n_xx: f64,
    n_yy: f64,
    s_xx: f64,
    s_xy: f64,
    intercept: f64,
    exp_x: i64,
    exp_y: i64,
}

// The exact sums of a pair of columns. The sums of the products use the
// two product of Dekker and Veltkamp. The numerators are the differences
// of two exact products of sums.
fn cross_core(xs: const f64[:], ys: const f64[:], want_line: bool) -> Cross, i64 {
    let kx = column_scale(xs);
    let ky = column_scale(ys);
    let fx = pow(2.0, kx as f64);
    let fy = pow(2.0, ky as f64);
    var sx = acc_new();
    var sy = acc_new();
    var sxx = acc_new();
    var syy = acc_new();
    var sxy = acc_new();
    for k in 0..xs.len {
        let u = xs[k] * fx;
        let v = ys[k] * fy;
        acc_add(sx, u);
        acc_add(sy, v);
        let p, e = two_prod(u, u);
        acc_add(sxx, p);
        acc_add(sxx, e);
        let q, f = two_prod(u, v);
        acc_add(sxy, q);
        acc_add(sxy, f);
        let g, h = two_prod(v, v);
        acc_add(syy, g);
        acc_add(syy, h);
    }
    if sx.fault || sy.fault || sxx.fault || syy.fault || sxy.fault { return Cross { .. }, 7; }
    var nn = acc_one(xs.len);
    var xy = acc_new();
    acc_add_product(xy, nn, sxy, 1.0);
    acc_add_product(xy, sx, sy, -1.0);
    var xx = acc_new();
    acc_add_product(xx, nn, sxx, 1.0);
    acc_add_product(xx, sx, sx, -1.0);
    var yy = acc_new();
    acc_add_product(yy, nn, syy, 1.0);
    acc_add_product(yy, sy, sy, -1.0);
    if xy.fault || xx.fault || yy.fault { return Cross { .. }, 7; }
    // The exact intercept of the line. The numerator is the difference of
    // two exact products of sums. The slope of SciPy is a rounded number,
    // so its intercept loses the digits of the slope.
    let nxx = acc_total(xx);
    var intercept: f64 = 0.0;
    if want_line {
        var num = acc_new();
        acc_add_product(num, sy, xx, 1.0);
        acc_add_product(num, xy, sx, -1.0);
        if num.fault { return Cross { .. }, 7; }
        if nxx != 0.0 {
            let nn = xs.len as f64;
            intercept = acc_total(num) / (nn * nxx);
        }
    }
    let r = Cross {
        n_xy: acc_total(xy),
        n_xx: nxx,
        n_yy: acc_total(yy),
        s_xx: acc_total(sxx),
        s_xy: acc_total(sxy),
        intercept: intercept,
        exp_x: kx,
        exp_y: ky,
    };
    return r, 0;
}

// The covariance of a numerator over the divisor d. The result of the
// scaled columns is scaled back at the end. The weighted calls share this
// core.
fn covariance_value(n_xy: f64, kx: i64, ky: i64, d: f64) -> f64, i64 {
    let v = n_xy / d;
    if is_inf(v) || is_nan(v) { return nan(), 8; }
    let out = scale_pow2(v, -(kx + ky));
    if is_inf(out) || is_nan(out) { return nan(), 8; }
    if out == 0.0 && n_xy != 0.0 { return nan(), 8; }
    return out, 0;
}

// The correlation of a pair from its numerators. A column of equal values
// gives the undefined fault. The result never leaves the range from -1 to
// 1. Two columns with the same variance share the denominator, so a column
// with itself gives exactly 1. The weighted call shares this core.
fn correlation_value(n_xy: f64, n_xx: f64, n_yy: f64) -> f64, i64 {
    if is_nan(n_xy) || is_nan(n_xx) || is_nan(n_yy) { return nan(), 8; }
    if n_xx <= 0.0 || n_yy <= 0.0 { return nan(), 7; }
    var den: f64 = 0.0;
    if n_xx == n_yy { den = n_xx; } else { den = sqrt(n_xx) * sqrt(n_yy); }
    if is_inf(den) || is_nan(den) || den == 0.0 { return nan(), 8; }
    let v = n_xy / den;
    if is_nan(v) { return nan(), 8; }
    if v > 1.0 { return 1.0, 0; }
    if v < -1.0 { return -1.0, 0; }
    return v, 0;
}

// ---- public calls -----------------------------------------------------

// The population covariance of two samples. One point gives 0.
fn pcovariance(xs: const f64[:], ys: const f64[:]) -> f64, i64 {
    let c = check_pair(xs, ys, 1);
    if c != 0 { return nan(), c; }
    let r, c2 = cross_core(xs, ys, false);
    if c2 != 0 { return nan(), c2; }
    let n = xs.len as f64;
    return covariance_value(r.n_xy, r.exp_x, r.exp_y, n * n);
}

// The sample covariance of two samples. It needs two points.
fn covariance(xs: const f64[:], ys: const f64[:]) -> f64, i64 {
    let c = check_pair(xs, ys, 2);
    if c != 0 { return nan(), c; }
    let r, c2 = cross_core(xs, ys, false);
    if c2 != 0 { return nan(), c2; }
    let n = xs.len as f64;
    return covariance_value(r.n_xy, r.exp_x, r.exp_y, n * (n - 1.0));
}

// The Pearson correlation of two samples. It needs two points. A column of
// equal values gives the undefined fault.
fn correlation(xs: const f64[:], ys: const f64[:]) -> f64, i64 {
    let c = check_pair(xs, ys, 2);
    if c != 0 { return nan(), c; }
    let r, c2 = cross_core(xs, ys, false);
    if c2 != 0 { return nan(), c2; }
    return correlation_value(r.n_xy, r.n_xx, r.n_yy);
}

// The exact sums and the exact numerators of a weighted pair. The weights
// are divided by the largest weight, as in v0.1.0.
struct WCross {
    n_xy: f64,
    n_xx: f64,
    n_yy: f64,
    w: f64,
    w2: f64,
    top: f64,
    exp_x: i64,
    exp_y: i64,
}

// The divisor of the weights of a pair. The largest magnitude covers both
// columns. The policy lives in weight_top_for of weighted.goose.
fn weight_top(xs: const f64[:], ys: const f64[:], ws: const f64[:]) -> f64 {
    var big: f64 = 0.0;
    for x in xs {
        let a = abs(x);
        if a > big { big = a; }
    }
    for y in ys {
        let a = abs(y);
        if a > big { big = a; }
    }
    return weight_top_for(ws, big);
}

// The weighted sums of a pair of columns. The product of the weight, the
// value and the other value is exact.
fn weighted_core(xs: const f64[:], ys: const f64[:], ws: const f64[:]) -> WCross, i64 {
    let kx = column_scale(xs);
    let ky = column_scale(ys);
    let fx = pow(2.0, kx as f64);
    let fy = pow(2.0, ky as f64);
    let top = weight_top(xs, ys, ws);
    var sw = acc_new();
    var sw2 = acc_new();
    var swx = acc_new();
    var swy = acc_new();
    var swxy = acc_new();
    var swxx = acc_new();
    var swyy = acc_new();
    for k in 0..xs.len {
        let w = ws[k] / top;
        let u = xs[k] * fx;
        let v = ys[k] * fy;
        acc_add(sw, w);
        let ww, we = two_prod(w, w);
        acc_add(sw2, ww);
        acc_add(sw2, we);
        let p, e = two_prod(w, u);
        acc_add(swx, p);
        acc_add(swx, e);
        let q, f = two_prod(w, v);
        acc_add(swy, q);
        acc_add(swy, f);
        // The product of the weight, the value and the other value is
        // the sum of two products of two products.
        let g, h = two_prod(p, v);
        let t, u2 = two_prod(e, v);
        acc_add(swxy, g);
        acc_add(swxy, h);
        acc_add(swxy, t);
        acc_add(swxy, u2);
        let a1, a2 = two_prod(p, u);
        let a3, a4 = two_prod(e, u);
        acc_add(swxx, a1);
        acc_add(swxx, a2);
        acc_add(swxx, a3);
        acc_add(swxx, a4);
        let b1, b2 = two_prod(q, v);
        let b3, b4 = two_prod(f, v);
        acc_add(swyy, b1);
        acc_add(swyy, b2);
        acc_add(swyy, b3);
        acc_add(swyy, b4);
    }
    if sw.fault || sw2.fault || swx.fault || swy.fault || swxy.fault { return WCross { .. }, 7; }
    var xy = acc_new();
    acc_add_product(xy, sw, swxy, 1.0);
    acc_add_product(xy, swx, swy, -1.0);
    var xx = acc_new();
    acc_add_product(xx, sw, swxx, 1.0);
    acc_add_product(xx, swx, swx, -1.0);
    var yy = acc_new();
    acc_add_product(yy, sw, swyy, 1.0);
    acc_add_product(yy, swy, swy, -1.0);
    if xy.fault || xx.fault || yy.fault { return WCross { .. }, 7; }
    let r = WCross {
        n_xy: acc_total(xy),
        n_xx: acc_total(xx),
        n_yy: acc_total(yy),
        w: acc_total(sw),
        w2: acc_total(sw2),
        top: top,
        exp_x: kx,
        exp_y: ky,
    };
    return r, 0;
}

// The divisor of a weighted form. The form is one of the WEIGHT_
// constants. The divisor is the raw weight form.
fn weight_form_divisor(r: const WCross&, form: i64) -> f64 {
    if form == WEIGHT_POPULATION { return r.w * r.w; }
    if form == WEIGHT_FREQUENCY { return r.w * (r.w - 1.0 / r.top); }
    return r.w * r.w - r.w2;
}

// The weighted covariance of two samples. The form is one of the WEIGHT_
// constants. Another value gives the domain fault.
fn weighted_covariance(xs: const f64[:], ys: const f64[:], ws: const f64[:], form: i64) -> f64, i64 {
    let c = check_weighted_pair(xs, ys, ws, 1);
    if c != 0 { return nan(), c; }
    if form < 0 || form > 2 { return nan(), 5; }
    let r, c2 = weighted_core(xs, ys, ws);
    if c2 != 0 { return nan(), c2; }
    let d = weight_form_divisor(r, form);
    if d == 0.0 { return nan(), 7; }
    return covariance_value(r.n_xy, r.exp_x, r.exp_y, d);
}

// The weighted correlation of two samples. The result is the same for all
// three forms of the weight. A column of equal values gives the undefined
// fault.
fn weighted_correlation(xs: const f64[:], ys: const f64[:], ws: const f64[:]) -> f64, i64 {
    let c = check_weighted_pair(xs, ys, ws, 2);
    if c != 0 { return nan(), c; }
    let r, c2 = weighted_core(xs, ys, ws);
    if c2 != 0 { return nan(), c2; }
    return correlation_value(r.n_xy, r.n_xx, r.n_yy);
}
