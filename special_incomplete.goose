// The incomplete gamma and the incomplete beta of the package.
//
// Every call returns the result first and a fault code last. A fault
// gives NaN as the result. The log gamma terms of the beta cancel, so the
// calls need the layer of pairs of doubles.
namespace statistics;

import std;
import math;

fn dd_log1pmx(x: DD) -> DD {
    if abs(x.hi) > 0.5 {
        return dd_sub(dd_log(dd_add(dd_from(1.0), x)), x);
    }
    var p = DD { hi: LOG1PMX_H[LOG1PMX_H.len - 1], lo: LOG1PMX_L[LOG1PMX_L.len - 1] };
    var j = LOG1PMX_H.len - 2;
    while j >= 0 {
        p = dd_add(dd_mul(p, x), DD { hi: LOG1PMX_H[j], lo: LOG1PMX_L[j] });
        j -= 1;
    }
    return dd_mul(p, dd_mul_d(dd_mul(x, x), -0.5));
}

// psi(delta) = delta - ln(1 + delta) as a pair. It is positive on both
// sides of 0.
fn dd_psi(delta: DD) -> DD {
    if abs(delta.hi) <= 0.5 {
        return dd_neg(dd_log1pmx(delta));
    }
    return dd_sub(delta, dd_log(dd_add(dd_from(1.0), delta)));
}

// The Stirling remainder S(a) = ln Gamma(a+1) - [(a+1/2) ln a - a +
// ln(2 pi)/2]. The series of pairs covers the range from 20 up, and the
// exact log gamma covers the range below.
fn dd_stirling_sum(a: f64) -> DD {
    // the sum is below 1e-101 here
    if a >= 1e100 { return dd_from(0.0); }
    var p = dd_from(0.0);
    let inv = dd_div_d(dd_from(1.0), a);
    let inv2 = dd_mul(inv, inv);
    var i = STIR_H.len - 1;
    while i >= 0 {
        p = dd_add(dd_mul(p, inv2), DD { hi: STIR_H[i], lo: STIR_L[i] });
        i -= 1;
    }
    return dd_mul(p, inv);
}

fn dd_stirling_s(a: f64) -> DD {
    if a >= 20.0 {
        return dd_stirling_sum(a);
    }
    let lg = dd_lngamma(dd_add_d(dd_from(a), 1.0));
    let ln_a = dd_log_d(a);
    let base = dd_add(dd_sub(dd_mul(dd_add_d(dd_from(a), 0.5), ln_a), dd_from(a)), DD { hi: HALF_LN2PI_H, lo: HALF_LN2PI_L });
    return dd_sub(lg, base);
}

// The log of x^a exp(-x) / Gamma(a+1). The large terms cancel in this
// form, so nothing is lost for a shape near the argument.
fn dd_log_prefactor(a: f64, x: f64) -> DD {
    // The form with the ratio x / a loses x when x is far below a, and it
    // overflows when a is tiny or x is huge. The direct form has neither
    // fault, because no large terms cancel there.
    if a < 1e-30 || x < 0.5 * a || x > 1e6 * a {
        let lg = dd_lngamma(dd_add_d(dd_from(a), 1.0));
        return dd_sub(dd_sub(dd_mul_d(dd_log_d(x), a), dd_from(x)), lg);
    }
    let delta = dd_div_d(dd_sub(dd_from(x), dd_from(a)), a);
    var t: DD = dd_from(0.0);
    if abs(delta.hi) <= 0.5 {
        t = dd_mul_d(dd_log1pmx(delta), a);
    } else {
        t = dd_mul_d(dd_sub(dd_log(dd_add(dd_from(1.0), delta)), delta), a);
    }
    let ln_a = dd_log_d(a);
    let half = dd_mul_d(dd_add(ln_a, DD { hi: LN2PI_H, lo: LN2PI_L }), 0.5);
    return dd_add(dd_sub(t, half), dd_neg(dd_stirling_s(a)));
}

// The sum of the power series of the lower tail, as a pair. The sum
// stays a pair so the complement keeps its own digits.
fn series_p(a: f64, x: f64) -> DD, i64 {
    var term = dd_from(1.0);
    var s = dd_from(1.0);
    var k: i64 = 1;
    while k <= MAX_ITER_SERIES {
        term = dd_div(dd_mul_d(term, x), dd_add_d(dd_from(a), k));
        s = dd_add(s, term);
        if abs(term.hi) <= 1e-33 * abs(s.hi) {
            return s, 0;
        }
        k += 1;
    }
    return s, 8;
}

// The exponential of a log pair, as a pair of doubles m and k with the
// value m * 2^k. The caller multiplies before it scales.
fn scaled_exp_parts(logv: DD) -> DD, i64 {
    let k = floor(logv.hi / DD_LN2_HI) as i64;
    let ln2 = DD { hi: DD_LN2_HI, lo: DD_LN2_LO };
    let r = dd_sub(logv, dd_mul_d(ln2, k));
    let m, mk = dd_exp_parts(r);
    return m, k + mk;
}

// The Legendre continued fraction of the upper tail, in plain doubles.
// The modified Lentz form keeps the quotient stable. The bound is fixed.
fn cf_q(a: f64, x: f64) -> f64, i64 {
    var b = x + 1.0 - a;
    var c: f64 = 1e300;
    var d: f64 = 1.0 / b;
    var h = d;
    var i: i64 = 1;
    while i <= MAX_ITER_CF {
        let an = -i * (i - a);
        b += 2.0;
        d = an * d + b;
        if abs(d) < 1e-300 { d = 1e-300; }
        c = b + an / c;
        if abs(c) < 1e-300 { c = 1e-300; }
        d = 1.0 / d;
        let de = d * c;
        h = h * de;
        if abs(de - 1.0) <= 1e-17 {
            return h, 0;
        }
        i += 1;
    }
    return h, 8;
}

// The upper tail from the continued fraction, scaled by the prefactor.
fn cf_tail(a: f64, x: f64) -> f64, i64 {
    let lp = dd_log_prefactor(a, x);
    // the tail is below the smallest subnormal float
    if lp.hi < -800.0 { return 0.0, 0; }
    let m, k0 = scaled_exp_parts(lp);
    let h, c = cf_q(a, x);
    if c != 0 { return nan(), c; }
    let q = dd_scale(dd_mul_d(dd_mul_d(m, a), h), k0);
    if is_inf(q) { return nan(), 8; }
    return q, 0;
}

// The lower tail from the series, scaled by the prefactor.
fn series_tail(a: f64, x: f64) -> DD, i64 {
    let s, c = series_p(a, x);
    if c != 0 { return dd_from(nan()), c; }
    // the log of the tail is the log of the prefactor plus the log of the sum
    let lt = dd_add(dd_log_prefactor(a, x), dd_log(s));
    if lt.hi < -800.0 { return dd_from(0.0), 0; }
    return dd_exp(lt), 0;
}

// The saddle point branch of the incomplete gamma. It covers the large
// shapes near the argument, where the series and the fraction both need
// too many steps. The series of the terms b[n] gives the tail.
fn saddle_tail(a: f64, x: f64, upper: bool) -> f64, i64 {
    let delta = dd_div_d(dd_sub(dd_from(x), dd_from(a)), a);
    let big = dd_mul_d(dd_psi(delta), a);
    let bigd = dd_to_f64(big);
    let e = dd_to_f64(dd_exp(dd_neg(big)));
    // below the smallest subnormal the tail is zero
    if bigd > 745.0 || e <= 0.0 { return 0.0, 0; }
    let w = sqrt(2.0 * bigd);
    var y: f64[41];
    let sc = w / sqrt(2.0);
    if sc < 26.0 {
        let ec, cc = erfc(sc);
        if cc != 0 { return nan(), cc; }
        var ds: f64 = 0.0;
        if bigd != 0.0 {
            let st = dd_sqrt(big);
            ds = dd_to_f64(dd_sub(st, dd_from(sc)));
        }
        let hz = TWO_OVER_SQRT_PI * exp(-sc * sc) / ec;
        let i0 = sqrt(1.5707963267948966) * ec * exp(-hz * ds);
        y[0] = i0 / e;
    } else {
        let w2 = w * w;
        var pw: f64 = 1.0;
        var acc: f64 = 0.0;
        var k: i64 = 0;
        while k < ERFCA.len {
            acc += ERFCA[k] * pw;
            pw = pw / w2;
            k += 1;
        }
        y[0] = acc / w;
    }
    y[1] = 1.0;
    var n: i64 = 2;
    while n < 41 {
        var wp: f64 = 1.0;
        var j: i64 = 1;
        while j < n {
            wp = wp * w;
            j += 1;
        }
        y[n] = wp + (n - 1.0) * y[n - 2];
        n += 1;
    }
    let inv_t = 1.0 / sqrt(2.0 * a);
    var tp: f64 = 1.0;
    var ssum: f64 = 0.0;
    n = 0;
    while n < 41 {
        var bn = B_SAD[n];
        if !upper && n % 2 == 0 { bn = -bn; }
        if n == 0 {
            ssum += bn * y[0] * e;
        } else {
            ssum += bn * tp * y[n] * e;
        }
        tp = tp * inv_t;
        n += 1;
    }
    let fac = dd_to_f64(dd_exp(dd_neg(dd_stirling_s(a)))) / sqrt(PI);
    var r = fac * ssum;
    if !upper { r = -r; }
    if is_inf(r) || is_nan(r) { return nan(), 8; }
    return r, 0;
}

// The log gamma of one plus a tiny shape. The pair form above shifts to
// 16 and subtracts the logs, and the error of the shift drowns the tiny
// value. The series in a gives the value itself, so no digits are lost.
fn dd_lngamma_tiny(a: f64) -> DD {
    var s = DD { hi: LG1P_H[12], lo: LG1P_L[12] };
    var k: i64 = 11;
    while k >= 0 {
        s = dd_add(dd_mul_d(s, a), DD { hi: LG1P_H[k], lo: LG1P_L[k] });
        k -= 1;
    }
    return dd_mul_d(s, a);
}

// The upper tail for a tiny shape and a small argument. The tail 1 - P
// loses all its digits there, because P is within 1e-16 of 1. The form
// Q = (1 - w) - w a S keeps them. Here w = x^a / Gamma(a + 1) and S is the
// alternating series of (-x)^k / (k! (a + k)). The sum of 1 - w is minus
// expm1 of its log. The flag is false when the form does not apply.
fn small_upper(a: f64, x: f64) -> f64, bool {
    let u = dd_sub(dd_mul_d(dd_log_d(x), a), dd_lngamma_tiny(a));
    if abs(u.hi) >= 1e-3 { return 0.0, false; }
    // minus expm1(u) from its series; the terms fall below 1e-33 by k = 12
    var em = dd_from(0.0);
    var t = dd_from(1.0);
    var k: i64 = 1;
    while k <= 14 {
        t = dd_div_d(dd_mul(t, u), k);
        em = dd_add(em, t);
        k += 1;
    }
    let w = dd_add(dd_from(1.0), em);
    var term = dd_from(1.0);
    var ssum = dd_from(0.0);
    k = 1;
    while k <= 100 {
        term = dd_mul_d(dd_div_d(term, k), -x);
        let part = dd_div(term, dd_add_d(dd_from(a), k));
        ssum = dd_add(ssum, part);
        if abs(part.hi) <= 1e-34 { break; }
        k += 1;
    }
    let q = dd_sub(dd_neg(em), dd_mul(dd_mul_d(w, a), ssum));
    return dd_to_f64(q), true;
}

// The regularized incomplete gamma. The small tail always comes first,
// and the other tail is one minus it.
fn forward_gammainc(a: f64, x: f64, upper: bool) -> f64, i64 {
    if x == 0.0 {
        if upper { return 1.0, 0; }
        return 0.0, 0;
    }
    let delta = (x - a) / a;
    if a >= 2500.0 {
        if delta > 2.5 {
            let q, c = cf_tail(a, x);
            if c != 0 { return nan(), c; }
            if upper { return q, 0; }
            return 1.0 - q, 0;
        }
        if delta < -0.5 {
            let p, c = series_tail(a, x);
            if c != 0 { return nan(), c; }
            if upper { return dd_to_f64(dd_sub(dd_from(1.0), p)), 0; }
            return dd_to_f64(p), 0;
        }
        if delta > 0.0 {
            let q, c = saddle_tail(a, x, true);
            if c != 0 { return nan(), c; }
            if upper { return q, 0; }
            return 1.0 - q, 0;
        }
        let p, c = saddle_tail(a, x, false);
        if c != 0 { return nan(), c; }
        if upper { return 1.0 - p, 0; }
        return p, 0;
    }
    if x < a + 1.0 || x < 6.0 {
        if upper && a < 1e-3 {
            let qs, ok = small_upper(a, x);
            if ok { return qs, 0; }
        }
        let p, c = series_tail(a, x);
        if c != 0 { return nan(), c; }
        if upper { return dd_to_f64(dd_sub(dd_from(1.0), p)), 0; }
        return dd_to_f64(p), 0;
    }
    let q, c = cf_tail(a, x);
    if c != 0 { return nan(), c; }
    if upper { return q, 0; }
    return 1.0 - q, 0;
}

// The regularized lower incomplete gamma P(a, x).
fn gammainc(a: f64, x: f64) -> f64, i64 {
    let ca = check_shape(a);
    if ca != 0 { return nan(), ca; }
    if is_nan(x) { return nan(), 3; }
    if is_inf(x) {
        if x > 0.0 { return 1.0, 0; }
        return nan(), 5;
    }
    if x < 0.0 { return nan(), 5; }
    if a > MAX_SHAPE { return nan(), 5; }
    return forward_gammainc(a, x, false);
}

// The regularized upper incomplete gamma Q(a, x).
fn gammaincc(a: f64, x: f64) -> f64, i64 {
    let ca = check_shape(a);
    if ca != 0 { return nan(), ca; }
    if is_nan(x) { return nan(), 3; }
    if is_inf(x) {
        if x > 0.0 { return 0.0, 0; }
        return nan(), 5;
    }
    if x < 0.0 { return nan(), 5; }
    if a > MAX_SHAPE { return nan(), 5; }
    return forward_gammainc(a, x, true);
}

// ---- the incomplete beta ---------------------------------------------

// One step of the modified Lentz form, on pairs. The quotient aa joins
// the running values c and d.
fn lentz_step(aa: DD, c: DD, d: DD) -> DD, DD {
    let one = dd_from(1.0);
    var dn = dd_add(one, dd_mul(aa, d));
    if abs(dn.hi) < 1e-300 { dn = dd_from(1e-300); }
    var cn = dd_add(one, dd_div(aa, c));
    if abs(cn.hi) < 1e-300 { cn = dd_from(1e-300); }
    dn = dd_div(one, dn);
    return cn, dn;
}

// The continued fraction of the incomplete beta, in pairs. The count of
// steps comes back with the value. The bound is fixed.
fn betacf_dd(a: f64, b: f64, x: DD) -> DD, i64 {
    let one = dd_from(1.0);
    let sab, eab = two_sum(a, b);
    let sap, eap = two_sum(a, 1.0);
    let sam, eam = two_sum(a, -1.0);
    let qab = DD { hi: sab, lo: eab };
    let qap = DD { hi: sap, lo: eap };
    let qam = DD { hi: sam, lo: eam };
    var c = one;
    var d = dd_sub(one, dd_div(dd_mul(qab, x), qap));
    if abs(d.hi) < 1e-300 { d = dd_from(1e-300); }
    d = dd_div(one, d);
    var h = d;
    var m: i64 = 1;
    while m <= MAX_ITER_CF {
        let fm = m as f64;
        let m2 = 2.0 * fm;
        let num_e = dd_mul(dd_mul_d(dd_sub(dd_from(b), dd_from(fm)), fm), x);
        let sd, ed = two_sum(a, m2);
        let den_e = dd_mul(dd_add_d(qam, m2), DD { hi: sd, lo: ed });
        let c1, d1 = lentz_step(dd_div(num_e, den_e), c, d);
        c = c1;
        d = d1;
        h = dd_mul(h, dd_mul(d, c));
        let sn, en = two_sum(a, fm);
        let num_o = dd_neg(dd_mul(dd_mul(DD { hi: sn, lo: en }, dd_add_d(qab, fm)), x));
        let den_o = dd_mul(DD { hi: sd, lo: ed }, dd_add_d(qap, m2));
        let c2, d2 = lentz_step(dd_div(num_o, den_o), c, d);
        c = c2;
        d = d2;
        let de = dd_mul(d, c);
        h = dd_mul(h, de);
        if abs((de.hi - 1.0) + de.lo) < 1e-30 {
            return h, 0;
        }
        m += 1;
    }
    return h, 8;
}

// The small tail of the beta distribution at x in (0, 1). The result is
// the side (0 for the lower tail, 1 for the upper tail), the tail as a
// pair with its power of two, and the code. The tail that converges
// fast is the tail that the call computes.
fn beta_tail(a: f64, b: f64, x: f64) -> i64, DD, i64, i64 {
    let ys, ye = two_sum(1.0, -x);
    let y = DD { hi: ys, lo: ye };
    let lx = dd_log_d(x);
    let ly = dd_log(y);
    let lnb = dd_lnbeta(a, b);
    let logp = dd_sub(dd_add(dd_mul_d(lx, a), dd_mul_d(ly, b)), lnb);
    var side: i64 = 0;
    var cf = dd_from(0.0);
    var code: i64 = 0;
    var shape = a;
    if x < (a + 1.0) / (a + b + 2.0) {
        let h, c = betacf_dd(a, b, dd_from(x));
        cf = h;
        code = c;
    } else {
        side = 1;
        shape = b;
        let h, c = betacf_dd(b, a, y);
        cf = h;
        code = c;
    }
    if code != 0 { return side, dd_from(nan()), 0, code; }
    let mant, k = dd_exp_parts(logp);
    // A tiny shape would overflow the quotient before the power of two
    // brings it back. Scale the shape up by 2^600 and the power with it.
    if shape < 1e-100 {
        let tail = dd_div_d(dd_mul(mant, cf), shape * pow2_int(600));
        return side, tail, k + 600, 0;
    }
    let tail = dd_div_d(dd_mul(mant, cf), shape);
    return side, tail, k, 0;
}

// The tail I_x(a, b) of a tiny shape b, from the Gauss series
// x^a / B(a, b) times the sum of (1 - b)_n x^n / (n! (a + n)). For a tiny b
// the tail is as small as b, and one minus its complement holds none of its
// digits. The terms fall by the factor x, so the sum needs a bound. A sum
// that does not converge inside the bound gives the fault code 8.
fn beta_series_tail(a: f64, b: f64, xd: DD) -> f64, i64 {
    var t = dd_from(1.0);
    var ssum = dd_div_d(dd_from(1.0), a);
    var n: i64 = 1;
    var done = false;
    while n <= MAX_ITER_SERIES {
        t = dd_mul(dd_div_d(dd_mul_d(t, n - b), n), xd);
        let part = dd_div(t, dd_add_d(dd_from(a), n));
        ssum = dd_add(ssum, part);
        if abs(part.hi) <= 1e-33 * abs(ssum.hi) { done = true; break; }
        n += 1;
    }
    if !done { return nan(), 8; }
    let logp = dd_sub(dd_mul(dd_log(xd), dd_from(a)), dd_lnbeta(a, b));
    let mant, k = dd_exp_parts(logp);
    return dd_scale(dd_mul(mant, ssum), k), 0;
}

// The code of the arguments of the incomplete beta function. A shape at
// or below 0 and a probability outside the unit range are out of range.
fn check_beta_args(a: f64, b: f64, x: f64) -> i64 {
    let ca = check_shape(a);
    if ca != 0 { return ca; }
    let cb = check_shape(b);
    if cb != 0 { return cb; }
    if is_nan(x) { return 3; }
    if is_inf(x) { return 5; }
    if x < 0.0 || x > 1.0 { return 5; }
    if a > MAX_SHAPE || b > MAX_SHAPE { return 5; }
    return 0;
}

// The regularized incomplete beta function I_x(a, b), or its upper tail
// 1 - I_x(a, b). The small tail comes straight, and the other tail is
// one minus it, in pairs. A tiny shape makes the direct tail equal to
// one within 1e-15. If the wanted tail is then the other one, it comes
// from the series, because the difference holds none of its digits.
fn beta_value(a: f64, b: f64, x: f64, upper: bool) -> f64, i64 {
    if x <= 0.0 {
        if upper { return 1.0, 0; }
        return 0.0, 0;
    }
    if x >= 1.0 {
        if upper { return 0.0, 0; }
        return 1.0, 0;
    }
    let side, tail, k, code = beta_tail(a, b, x);
    if code != 0 { return nan(), code; }
    let direct = dd_scale(tail, k);
    // the direct tail is the wanted one
    if upper {
        if side == 1 { return direct, 0; }
    } else {
        if side == 0 { return direct, 0; }
    }
    // the wanted tail is one minus the direct tail
    if direct > 1.0 - 1e-15 {
        if upper {
            let ys, ye = two_sum(1.0, -x);
            return beta_series_tail(b, a, DD { hi: ys, lo: ye });
        }
        return beta_series_tail(a, b, dd_from(x));
    }
    let scaled = DD { hi: tail.hi * pow2_int(k), lo: tail.lo * pow2_int(k) };
    return dd_to_f64(dd_sub(dd_from(1.0), scaled)), 0;
}

// The regularized incomplete beta function I_x(a, b).
fn betainc(a: f64, b: f64, x: f64) -> f64, i64 {
    let c = check_beta_args(a, b, x);
    if c != 0 { return nan(), c; }
    return beta_value(a, b, x, false);
}

// The upper tail 1 - I_x(a, b) of the regularized incomplete beta
// function, with x exact. The distributions use it, because a caller
// that passes 1 - x would lose the digits of a small x.
fn betainc_upper(a: f64, b: f64, x: f64) -> f64, i64 {
    let c = check_beta_args(a, b, x);
    if c != 0 { return nan(), c; }
    return beta_value(a, b, x, true);
}
