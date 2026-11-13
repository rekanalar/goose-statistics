// The inverses of the special functions: the inverse of the error
// function, of the incomplete gamma and of the incomplete beta.
//
// Every call returns the result first and a fault code last. A fault
// gives NaN as the result.
namespace statistics;

import std;
import math;

fn norm_start(p: f64) -> f64 {
    var q = p;
    if p > 0.5 { q = 1.0 - p; }
    let r, c = erfcinv(2.0 * q);
    if c != 0 { return 0.0; }
    let z = -1.4142135623730951 * r;
    if p > 0.5 { return -z; }
    return z;
}

// The log of the tail that the call solves for. The lower tail below the
// switch comes from the series as a pair, so a tail near 1e-300 never
// underflows. The other cases take the log of the forward value.
fn log_tail_gamma(a: f64, x: f64, upper: bool) -> f64, i64 {
    if !upper && a < 2500.0 && (x < a + 1.0 || x < 6.0) {
        let s, c = series_p(a, x);
        if c != 0 { return nan(), c; }
        return dd_to_f64(dd_add(dd_log_prefactor(a, x), dd_log(s))), 0;
    }
    let v, c = forward_gammainc(a, x, upper);
    if c != 0 { return nan(), c; }
    if v <= 0.0 { return -1.0 / 0.0, 0; }
    return log(v), 0;
}

// The log of the density x^(a-1) exp(-x) / Gamma(a).
fn log_density_gamma(a: f64, x: f64) -> f64 {
    let t = dd_to_f64(dd_sub(dd_mul_d(dd_log_d(x), a - 1.0), dd_from(x)));
    let lg, c = lgamma_pos(a);
    return t - lg;
}

// Solve the tail of the incomplete gamma for a target log, in u = log x.
// The residual is monotone, so one bracket holds the root. A step that
// leaves the bracket falls back to a bisection step.
fn solve_gamma_tail(a: f64, target: f64, upper: bool) -> f64, i64 {
    let target_log = dd_to_f64(dd_log_d(target));
    var u0: f64 = 0.0;
    if upper {
        let z = -norm_start(target);
        let s = 1.0 - 1.0 / (9.0 * a) + z / (3.0 * sqrt(a));
        let x0 = a * s * s * s;
        if x0 > 0.0 && !is_inf(x0) { u0 = log(x0); } else { u0 = log(a); }
    } else {
        let lg, c = lgamma_pos(a);
        u0 = (target_log + log(a) + lg) / a;
    }
    if u0 < -740.0 { u0 = -740.0; }
    if u0 > 700.0 { u0 = 700.0; }
    var lo = u0 - 4.0;
    var hi = u0 + 4.0;
    var step: f64 = 4.0;
    var n: i64 = 0;
    // the bracket: the residual is below 0 at lo and above 0 at hi for the
    // lower tail, and the other way round for the upper tail
    while n < 60 {
        let flo, c1 = log_tail_gamma(a, exp(lo), upper);
        if c1 != 0 { return nan(), c1; }
        let fhi, c2 = log_tail_gamma(a, exp(hi), upper);
        if c2 != 0 { return nan(), c2; }
        var ok_lo = false;
        var ok_hi = false;
        if upper {
            ok_lo = flo - target_log >= 0.0;
            ok_hi = fhi - target_log <= 0.0;
        } else {
            ok_lo = flo - target_log <= 0.0;
            ok_hi = fhi - target_log >= 0.0;
        }
        if ok_lo && ok_hi { break; }
        step = step * 2.0;
        if !ok_lo { lo -= step; }
        if !ok_hi { hi += step; }
        if lo < -745.0 { lo = -745.0; }
        if hi > 709.0 { hi = 709.0; }
        if lo <= -745.0 && hi >= 709.0 { break; }
        n += 1;
    }
    var u = 0.5 * (lo + hi);
    var i: i64 = 0;
    while i < MAX_ITER_INVERSE {
        let x = exp(u);
        let lt, c = log_tail_gamma(a, x, upper);
        if c != 0 { return nan(), c; }
        let f = lt - target_log;
        if upper {
            if f > 0.0 { lo = u; } else { hi = u; }
        } else {
            if f > 0.0 { hi = u; } else { lo = u; }
        }
        let ld = log_density_gamma(a, x) - lt + u;
        var g: f64 = 0.0;
        if ld < 700.0 { g = exp(ld); }
        if upper { g = -g; }
        var un = 0.5 * (lo + hi);
        if g != 0.0 && !is_inf(g) {
            let t = u - f / g;
            if t > lo && t < hi && !is_nan(t) { un = t; }
        }
        if abs(un - u) <= 1e-16 * max(abs(un), 1.0) {
            u = un;
            break;
        }
        u = un;
        i += 1;
    }
    var xr = exp(u);
    // Polish in the value domain. The log of a tail near 1e-15 holds the
    // rounding of a number of size 35, and the value does not.
    if target > 1e-290 {
        var k: i64 = 0;
        while k < 2 {
            let v, c = forward_gammainc(a, xr, upper);
            if c != 0 { return nan(), c; }
            let dens = exp(log_density_gamma(a, xr));
            if dens <= 0.0 || is_inf(dens) || is_nan(dens) { break; }
            var xn = xr - (v - target) / dens;
            if upper { xn = xr + (v - target) / dens; }
            if !(xn > 0.0) || is_inf(xn) || abs(xn - xr) > 0.01 * xr { break; }
            xr = xn;
            k += 1;
        }
    }
    return xr, 0;
}

// The inverse of the lower regularized incomplete gamma.
fn gammaincinv(a: f64, p: f64) -> f64, i64 {
    let ca = check_shape(a);
    if ca != 0 { return nan(), ca; }
    let cp = check_prob(p);
    if cp != 0 { return nan(), cp; }
    if a > MAX_SHAPE { return nan(), 5; }
    if p == 0.0 { return 0.0, 0; }
    if p == 1.0 { return nan(), 8; }
    if p > 0.5 { return solve_gamma_tail(a, 1.0 - p, true); }
    return solve_gamma_tail(a, p, false);
}

// The inverse of the upper regularized incomplete gamma.
fn gammainccinv(a: f64, q: f64) -> f64, i64 {
    let ca = check_shape(a);
    if ca != 0 { return nan(), ca; }
    let cq = check_prob(q);
    if cq != 0 { return nan(), cq; }
    if a > MAX_SHAPE { return nan(), 5; }
    if q == 1.0 { return 0.0, 0; }
    if q == 0.0 { return nan(), 8; }
    if q > 0.5 { return solve_gamma_tail(a, 1.0 - q, false); }
    return solve_gamma_tail(a, q, true);
}

// The log of one plus a pair, for a pair above -1.
fn dd_log1p_dd(x: DD) -> DD {
    let s, e = two_sum(1.0, x.hi);
    let t = dd_add(DD { hi: s, lo: e }, dd_from(x.lo));
    if t.hi <= 0.0 { return dd_from(-1.0 / 0.0); }
    return dd_log(t);
}

// The logs of the two tails of the beta distribution and of its density
// at x. The tail that the fraction computes is exact in the log, and the
// other tail is the complement in pairs.
fn beta_log_ij(a: f64, b: f64, x: f64) -> DD, DD, DD, i64 {
    let ys, ye = two_sum(1.0, -x);
    let y = DD { hi: ys, lo: ye };
    let lx = dd_log_d(x);
    let ly = dd_log(y);
    let lnb = dd_lnbeta(a, b);
    let pref = dd_sub(dd_add(dd_mul_d(lx, a), dd_mul_d(ly, b)), lnb);
    var log_i = dd_from(0.0);
    var log_j = dd_from(0.0);
    if x < (a + 1.0) / (a + b + 2.0) {
        let cf, c = betacf_dd(a, b, dd_from(x));
        if c != 0 { return log_i, log_j, log_i, c; }
        log_i = dd_sub(dd_add(pref, dd_log(cf)), dd_log_d(a));
        log_j = dd_log1p_dd(dd_neg(dd_exp(log_i)));
    } else {
        let cf, c = betacf_dd(b, a, y);
        if c != 0 { return log_i, log_j, log_i, c; }
        log_j = dd_sub(dd_add(pref, dd_log(cf)), dd_log_d(b));
        log_i = dd_log1p_dd(dd_neg(dd_exp(log_j)));
    }
    let logpdf = dd_sub(dd_add(dd_mul_d(lx, a - 1.0), dd_mul_d(ly, b - 1.0)), lnb);
    return log_i, log_j, logpdf, 0;
}

// The start of the beta inverse. It is only a start, and the loop does
// the work.
fn beta_start(a: f64, b: f64, p: f64) -> f64 {
    let s = a + b;
    if a > 1.0 && b > 1.0 {
        let sd = sqrt(a * b / (s * s * (s + 1.0)));
        let z = norm_start(p);
        let skew = 2.0 * (b - a) * sqrt(s + 1.0) / ((s + 2.0) * sqrt(a * b));
        let x = a / s + sd * (z + skew * (z * z - 1.0) / 6.0);
        if x > 0.0 && x < 1.0 { return x; }
        return a / s;
    }
    let lb = dd_to_f64(dd_lnbeta(a, b));
    if p <= 0.5 {
        let v = (log(p * a) + lb) / a;
        if v < -1e-3 { return exp(max(v, -744.0)); }
    } else {
        let v = (log((1.0 - p) * b) + lb) / b;
        if v < -1e-3 { return 1.0 - exp(max(v, -744.0)); }
    }
    return 0.5;
}

// The inverse of the incomplete beta. The solver iterates on x with the
// log of the small tail as the residual. A geometric bisection step
// replaces every Newton step that leaves the bracket.
fn betaincinv(a: f64, b: f64, p: f64) -> f64, i64 {
    let ca = check_shape(a);
    if ca != 0 { return nan(), ca; }
    let cb = check_shape(b);
    if cb != 0 { return nan(), cb; }
    let cp = check_prob(p);
    if cp != 0 { return nan(), cp; }
    if a > MAX_SHAPE || b > MAX_SHAPE { return nan(), 5; }
    if p == 0.0 { return 0.0, 0; }
    if p == 1.0 { return 1.0, 0; }
    let upper = p > 0.5;
    var logt = dd_log_d(p);
    if upper { logt = dd_log_d(1.0 - p); }
    var x = beta_start(a, b, p);
    if !(x > 0.0 && x < 1.0) { x = 0.5; }
    var lo: f64 = 0.0;
    var hi: f64 = 1.0;
    var best_x = x;
    var best_f = 1.0 / 0.0;
    var last_f: f64 = 0.0;
    var it: i64 = 0;
    while it < MAX_ITER_INVERSE {
        let log_i, log_j, logpdf, c = beta_log_ij(a, b, x);
        if c != 0 { return nan(), c; }
        var f: f64 = 0.0;
        var lrat = dd_sub(log_i, logpdf);
        if upper {
            f = dd_to_f64(dd_sub(logt, log_j));
            lrat = dd_sub(log_j, logpdf);
        } else {
            f = dd_to_f64(dd_sub(log_i, logt));
        }
        last_f = f;
        if is_nan(f) { break; }
        if abs(f) < best_f {
            best_f = abs(f);
            best_x = x;
        }
        if abs(f) <= 4e-16 { break; }
        if f > 0.0 { hi = x; } else { lo = x; }
        var step: f64 = 0.0;
        if !is_nan(lrat.hi) && lrat.hi > -700.0 && lrat.hi < 700.0 {
            if x < 0.5 {
                let lr = dd_sub(lrat, dd_log_d(x));
                if lr.hi > -700.0 && lr.hi < 700.0 {
                    var du = -f * dd_to_f64(dd_exp(lr));
                    if !is_inf(du) && !is_nan(du) {
                        du = max(-50.0, min(50.0, du));
                        step = x * expm1_core(du);
                    }
                }
            } else {
                step = -f * dd_to_f64(dd_exp(lrat));
                if is_inf(step) || is_nan(step) { step = 0.0; }
            }
        }
        if step != 0.0 && abs(step) <= abs(x) * 2.220446049250313e-16 { break; }
        var cand: f64 = -1.0;
        if step != 0.0 {
            let t = x + step;
            if t > lo && t < hi { cand = t; }
        }
        if cand < 0.0 {
            var l = lo;
            if l <= 0.0 { l = 5e-324; }
            if hi > 8.0 * l {
                let m = sqrt(l) * sqrt(hi);
                if m > lo && m < hi { cand = m; } else { cand = 0.5 * (lo + hi); }
            } else {
                cand = 0.5 * (lo + hi);
            }
        }
        if cand == x || hi <= lo * 1.0000000000000002 { break; }
        x = cand;
        it += 1;
    }
    if best_x <= 5e-324 && last_f > 0.0 { return 0.0, 0; }
    return best_x, 0;
}
