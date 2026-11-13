// The gamma family of the package, the log and the exponential of pairs,
// and the beta, lbeta and lchoose.
//
// Every call returns the result first and a fault code last. A fault
// gives NaN as the result.
namespace statistics;

import std;
import math;

const BIGP: f64 = 3.273390607896142e+150;

// The gamma core on the range from 1 to 2. The fit has 24 terms.
fn gamma_core_fit(f: f64) -> f64 {
    return clenshaw(GAMMA_C, 2.0 * f - 3.0);
}

// The product gamma(1 + f) = gamma(f) * f * (f + 1) * ... for f in [1, 2)
// and k factors more. The result is m * 2^e. The pair of doubles keeps
// the digits of the product, and the exponent keeps the range.
fn gamma_core(k: f64, f: f64) -> f64, i64 {
    let n = k + 1.0;
    var ph = gamma_core_fit(f);
    var pl: f64 = 0.0;
    var e: i64 = 0;
    var j: f64 = 0.0;
    while j < n - 1.0 {
        let s, se = two_sum(f, j);
        let p, q0 = two_prod(ph, s);
        let q = q0 + (ph * se + pl * s);
        let h, l = two_sum(p, q);
        ph = h;
        pl = l;
        if ph > BIGP {
            ph *= 3.054936363499605e-151;
            pl *= 3.054936363499605e-151;
            e += 500;
        }
        j += 1.0;
    }
    return ph + pl, e;
}

// The gamma of x at or above 1. The value comes back as m * 2^e.
fn gamma_scaled(x: f64) -> f64, i64 {
    let n = floor(x);
    return gamma_core(n - 1.0, x - n + 1.0);
}

// The value m * 2^e with e a multiple of 500, scaled back up.
fn mulp(m: f64, e: i64) -> f64 {
    var r = m;
    var t = e;
    while t > 0 {
        r = r * BIGP;
        t -= 500;
    }
    return r;
}

fn divp(m: f64, e: i64) -> f64 {
    var r = m;
    var t = e;
    while t > 0 {
        r = r / BIGP;
        t -= 500;
    }
    return r;
}

// The sign of the gamma of a negative non integer x. It is +1 when the
// integer part below x is even, and -1 when it is odd.
fn gamma_sign(n: f64) -> f64 {
    let i = n as i64;
    if i % 2 == 0 { return 1.0; }
    return -1.0;
}

// The gamma function. The poles give the undefined code. A result that
// does not fit gives the overflow code. The reflection covers the
// negative side, and the far negative side runs in the log domain.
fn gamma(x: f64) -> f64, i64 {
    if is_nan(x) { return nan(), 3; }
    if is_inf(x) { return nan(), 8; }
    if x == 0.0 || (x < 0.0 && x == floor(x)) { return nan(), 7; }
    if x >= 1.0 {
        if x > 171.6243769563027 { return nan(), 8; }
        let m, e = gamma_scaled(x);
        let r = mulp(m, e);
        if is_inf(r) { return nan(), 8; }
        return r, 0;
    }
    if x > 0.0 {
        let m, e = gamma_scaled(x + 1.0);
        let r = m / x;
        if is_inf(r) || is_nan(r) { return nan(), 8; }
        return r, 0;
    }
    let n = floor(x);
    let f0 = x - n;
    var f = f0;
    if f > 0.5 { f = n + 1.0 - x; }
    let s = sin(PI * f);
    let sg = gamma_sign(n);
    if x < -201.0 {
        let lv, lc = lgamma_pos(1.0 - x);
        if lc != 0 { return nan(), lc; }
        let r = exp(log(PI / abs(s)) - lv);
        return sg * r, 0;
    }
    var m: f64 = 0.0;
    var e: i64 = 0;
    if n == -1.0 {
        let m2, e2 = gamma_scaled(1.0 - x);
        m = m2;
        e = e2;
    } else {
        let m2, e2 = gamma_core(-n - 1.0, 2.0 - f0);
        m = m2;
        e = e2;
    }
    var r = PI / (s * m);
    r = divp(r, e);
    if is_inf(r) || is_nan(r) { return nan(), 8; }
    return sg * r, 0;
}

// The log gamma of a positive x. The fit covers the range from 0.5 to 3,
// the scaled product covers the range from 3 to 15, and the Stirling
// series covers the rest.
recursive fn lgamma_pos(x: f64) -> f64, i64 {
    if x < 0.5 {
        let v, c = lgamma_pos(x + 1.0);
        if c != 0 { return nan(), c; }
        return v - log(x), 0;
    }
    if x <= 3.0 {
        if x == 1.0 || x == 2.0 { return 0.0, 0; }
        return clenshaw(LGAMMA_C, (2.0 * x - 3.5) / 2.5) * ((x - 1.0) * (x - 2.0)), 0;
    }
    if x < 15.0 {
        let m, e = gamma_scaled(x);
        if e != 0 { return nan(), 8; }
        let r = log(m);
        if is_inf(r) || is_nan(r) { return nan(), 8; }
        return r, 0;
    }
    let z = 1.0 / (x * x);
    var s: f64 = 0.0;
    var i = STIR_C.len - 1;
    while i >= 0 {
        s = s * z + STIR_C[i];
        i -= 1;
    }
    let r = (x - 0.5) * log(x) - x + 0.9189385332046727 + s / x;
    if is_inf(r) || is_nan(r) { return nan(), 8; }
    return r, 0;
}

// The log of the absolute value of the gamma, with the sign. The sign is
// +1 or -1, and 0 on a fault.
fn lgamma(x: f64) -> f64, f64, i64 {
    if is_nan(x) { return nan(), 0.0, 3; }
    if is_inf(x) { return nan(), 0.0, 8; }
    if x == 0.0 || (x < 0.0 && x == floor(x)) { return nan(), 0.0, 7; }
    if x > 0.0 {
        let v, c = lgamma_pos(x);
        if c != 0 { return nan(), 0.0, c; }
        return v, 1.0, 0;
    }
    let n = floor(x);
    // the distance to the nearest integer, exact for a tiny x
    var f = x - n;
    if f > 0.5 { f = n + 1.0 - x; }
    let s = abs(sin(PI * f));
    let sg = gamma_sign(n);
    let y = 1.0 - x;
    let v, c = lgamma_pos(y);
    if c != 0 { return nan(), 0.0, c; }
    return log(PI / s) - v, sg, 0;
}

// The digamma of a positive x. The fit around the root covers the range
// from 1 to 3, a second fit covers the range from 3 to 10, and the
// asymptotic series covers the rest.
recursive fn digamma_pos(x: f64) -> f64, i64 {
    if x >= 10.0 {
        let z = 1.0 / (x * x);
        var s: f64 = 0.0;
        var i = DIG_ASY.len - 1;
        while i >= 0 {
            s = s * z + DIG_ASY[i];
            i -= 1;
        }
        return log(x) - 0.5 / x - s * z, 0;
    }
    if x < 1.0 {
        let v, c = digamma_pos(x + 1.0);
        if c != 0 { return nan(), c; }
        return v - 1.0 / x, 0;
    }
    if x <= 3.0 {
        let d = (x - X0_HI) - X0_LO;
        return clenshaw(DIG_C, x - 2.0) * d, 0;
    }
    return clenshaw(DIG2_C, (2.0 * x - 13.0) / 7.0), 0;
}

// The digamma function. The poles give the undefined code.
fn digamma(x: f64) -> f64, i64 {
    if is_nan(x) { return nan(), 3; }
    if is_inf(x) { return nan(), 8; }
    if x == 0.0 || (x < 0.0 && x == floor(x)) { return nan(), 7; }
    if x > 0.0 {
        return digamma_pos(x);
    }
    let n = floor(x);
    let f = x - n;
    var c: f64 = 0.0;
    if f <= 0.5 {
        c = 1.0 / tan(PI * f);
    } else {
        c = -1.0 / tan(PI * (1.0 - f));
    }
    let v, cv = digamma_pos(1.0 - x);
    if cv != 0 { return nan(), cv; }
    let r = v - PI * c;
    if is_inf(r) || is_nan(r) { return nan(), 8; }
    return r, 0;
}

// ---- the log and the exponential of pairs ---------------------------

// The logarithm of a positive value, as a pair. The series of the area
// hyperbolic tangent runs on the pair z, and the exponent of two comes
// back through the log of 2.
const LOG_SERIES_TERMS: i64 = 24;

fn dd_log_d(x: f64) -> DD {
    if x <= 0.0 { return dd_from(0.0 / 0.0); }
    var m: f64 = 0.0;
    var e: i64 = 0;
    let m2, e2 = frexp2(x);
    m = m2;
    e = e2;
    if m > 1.4142135623730951 {
        m = m * 0.5;
        e += 1;
    }
    let ms, me = two_sum(m, 1.0);
    let z = dd_div(dd_from(m - 1.0), DD { hi: ms, lo: me });
    let z2 = dd_mul(z, z);
    var term = z;
    var total = z;
    var k: i64 = 1;
    while k < LOG_SERIES_TERMS {
        term = dd_mul(term, z2);
        total = dd_add(total, dd_div_d(term, 2.0 * k + 1.0));
        k += 1;
    }
    let log_m = dd_mul_d(total, 2.0);
    let ln2 = DD { hi: DD_LN2_HI, lo: DD_LN2_LO };
    return dd_add(log_m, dd_mul_d(ln2, e));
}

// The logarithm of a pair. The base comes from the one value form, and
// the second member adds a correction of the first order.
fn dd_log(a: DD) -> DD {
    if a.hi <= 0.0 { return dd_from(0.0 / 0.0); }
    let base = dd_log_d(a.hi);
    if a.lo == 0.0 { return base; }
    let t = a.lo / a.hi;
    return dd_add(base, dd_from(t - 0.5 * t * t));
}

// The exponential of a pair, as a pair and a power of two. The split
// keeps the scaling out of the series, so a subnormal result rounds once.
fn dd_exp_parts(a: DD) -> DD, i64 {
    let k = floor(a.hi / DD_LN2_HI + 0.5) as i64;
    let ln2 = DD { hi: DD_LN2_HI, lo: DD_LN2_LO };
    var r = dd_sub(a, dd_mul_d(ln2, k));
    r = dd_mul_d(r, 0.0625);
    var term = dd_from(1.0);
    var total = dd_from(1.0);
    var n: i64 = 1;
    while n < 16 {
        term = dd_div_d(dd_mul(term, r), n);
        total = dd_add(total, term);
        n += 1;
    }
    var i: i64 = 0;
    while i < 4 {
        total = dd_mul(total, total);
        i += 1;
    }
    return total, k;
}

fn dd_exp(a: DD) -> DD {
    let m, k = dd_exp_parts(a);
    return dd_scale_dd(m, k);
}

// The shift of the Stirling series. Arguments below it move up first,
// and the logs of the shifts come back off.
const STIRLING_SHIFT: f64 = 16.0;

// The log gamma of a positive pair. The Stirling series of 14 pairs
// covers the range from 16 up, and the shift covers the range below.
recursive fn dd_lngamma(z: DD) -> DD {
    if z.hi < STIRLING_SHIFT {
        let n = (floor(STIRLING_SHIFT - z.hi) + 1.0) as i64;
        var acc = dd_lngamma(dd_add_d(z, n));
        var i: i64 = 0;
        while i < n {
            acc = dd_sub(acc, dd_log(dd_add_d(z, i)));
            i += 1;
        }
        return acc;
    }
    let t1 = dd_mul(dd_log(z), dd_sub(z, dd_from(0.5)));
    var acc = dd_sub(t1, z);
    acc = dd_add(acc, dd_mul_d(DD { hi: LN2PI_H, lo: LN2PI_L }, 0.5));
    let zinv = dd_div(dd_from(1.0), z);
    let zinv2 = dd_mul(zinv, zinv);
    var term = zinv;
    var i: i64 = 0;
    while i < STIR_H.len {
        acc = dd_add(acc, dd_mul(term, DD { hi: STIR_H[i], lo: STIR_L[i] }));
        term = dd_mul(term, zinv2);
        i += 1;
    }
    return acc;
}

// The log of the beta function of two pairs. The two log gamma terms
// cancel, so the plain form loses the digits of the small shapes, and the
// pair form does not. A shape of 1e8 or more would lose its digits to the
// size of the log gamma terms, and 1e150 would overflow them. Those shapes
// take the difference of the Stirling forms, which holds only ratios.
// The product of two pairs when the first is as large as 1e308. The split
// of the two product overflows above 1e300, so the first factor goes down
// by 2^200 and the product comes back.
fn dd_mul_big(a: DD, b: DD) -> DD {
    if abs(a.hi) < 1e280 { return dd_mul(a, b); }
    return dd_scale_dd(dd_mul(dd_scale_dd(a, -200), b), 200);
}

// The ratio of two pairs of the same size, safe for a size of 1e308.
fn dd_ratio_big(a: DD, b: DD) -> DD {
    if abs(b.hi) < 1e280 { return dd_div(a, b); }
    return dd_div(dd_scale_dd(a, -200), dd_scale_dd(b, -200));
}

fn dd_lnbeta_dd(a: DD, b: DD) -> DD {
    let s = dd_add(a, b);
    var lo = a;
    var hi = b;
    if a.hi > b.hi {
        lo = b;
        hi = a;
    }
    if hi.hi < 1e8 {
        return dd_sub(dd_add(dd_lngamma(a), dd_lngamma(b)), dd_lngamma(s));
    }
    let half_ln2pi = DD { hi: HALF_LN2PI_H, lo: HALF_LN2PI_L };
    var ls = dd_log(s);
    if hi.hi >= 1e280 {
        // the sum may not fit in a double: take the log of the sum over 2^200
        let ss = dd_add(dd_scale_dd(lo, -200), dd_scale_dd(hi, -200));
        ls = dd_add(dd_log(ss), dd_mul_d(DD { hi: DD_LN2_HI, lo: DD_LN2_LO }, 200.0));
    }
    let l1 = dd_log1p_dd(dd_ratio_big(lo, hi));
    if lo.hi >= 20.0 {
        // ln(hi / s) is minus log1p(lo / hi): the plain difference of the logs loses it
        let ll = dd_log(lo);
        var t = dd_sub(dd_mul_big(lo, dd_sub(ll, ls)), dd_mul_big(hi, l1));
        t = dd_sub(t, dd_mul_d(dd_sub(ll, l1), 0.5));
        t = dd_add(t, half_ln2pi);
        let sm = dd_sub(dd_add(dd_stirling_sum(a.hi), dd_stirling_sum(b.hi)), dd_stirling_sum(a.hi + b.hi));
        return dd_add(t, sm);
    }
    // the small shape keeps its exact log gamma
    let t1 = dd_mul_big(dd_add_d(hi, -0.5), l1);
    var t = dd_sub(lo, dd_add(t1, dd_mul(lo, ls)));
    t = dd_add(t, dd_lngamma(lo));
    return dd_add(t, dd_sub(dd_stirling_sum(hi.hi), dd_stirling_sum(a.hi + b.hi)));
}

fn dd_lnbeta(a: f64, b: f64) -> DD {
    return dd_lnbeta_dd(dd_from(a), dd_from(b));
}

// ---- beta, lbeta and lchoose ----------------------------------------

fn lbeta(a: f64, b: f64) -> f64, i64 {
    let ca = check_shape(a);
    if ca != 0 { return nan(), ca; }
    let cb = check_shape(b);
    if cb != 0 { return nan(), cb; }
    let r = dd_scale(dd_lnbeta(a, b), 0);
    if is_inf(r) || is_nan(r) { return nan(), 8; }
    return r, 0;
}

fn beta(a: f64, b: f64) -> f64, i64 {
    let ca = check_shape(a);
    if ca != 0 { return nan(), ca; }
    let cb = check_shape(b);
    if cb != 0 { return nan(), cb; }
    let lnb = dd_lnbeta(a, b);
    if lnb.hi > 709.9 { return nan(), 8; }
    // below the smallest subnormal the true value is zero
    if lnb.hi < -745.3 { return 0.0, 0; }
    let m, k = dd_exp_parts(lnb);
    let r = dd_scale(m, k);
    if is_inf(r) { return nan(), 8; }
    return r, 0;
}

// The log of the binomial coefficient. The end points give the exact 0.
// The log of the binomial coefficient as a pair.
fn dd_lchoose(n: f64, k: f64) -> DD {
    if k == 0.0 || k == n { return dd_from(0.0); }
    // C(n, k) = 1 / ((n + 1) B(k + 1, n - k + 1))
    let k1 = dd_add_d(dd_from(k), 1.0);
    let nk1 = dd_add_d(dd_sub(dd_from(n), dd_from(k)), 1.0);
    return dd_neg(dd_add(dd_log(dd_add_d(dd_from(n), 1.0)), dd_lnbeta_dd(k1, nk1)));
}

fn lchoose(n: f64, k: f64) -> f64, i64 {
    let cn = check_scalar(n);
    if cn != 0 { return nan(), cn; }
    let ck = check_scalar(k);
    if ck != 0 { return nan(), ck; }
    if n < 0.0 || k < 0.0 || k > n { return nan(), 5; }
    if k == 0.0 || k == n { return 0.0, 0; }
    let r = dd_to_f64(dd_lchoose(n, k));
    if is_inf(r) || is_nan(r) { return nan(), 8; }
    return r, 0;
}
