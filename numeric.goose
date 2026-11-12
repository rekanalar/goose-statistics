// The numeric helpers of the searches and the integrals. Each search
// takes its objective as the trailing block of the call. Each loop stops
// at MAX_FIT_STEPS steps and gives ERR_FIT, so no search runs without end.
namespace statistics;

import std;
import math;

const GOLD: f64 = 1.618034;
const VERY_SMALL: f64 = 1e-21;
const BRENT_GROW: f64 = 110.0;
const BRENT_CG: f64 = 0.3819660;
const BRENT_MINTOL: f64 = 1e-11;

// A bracket (xa, xb, xc) of a minimum of F from two start points, by
// the downhill golden-ratio search of SciPy. The middle point has the
// lowest value. The code is ERR_FIT when the search does not end.
fn bracket_min<F>(xa0: f64, xb0: f64) -> f64, f64, f64, i64 {
    var xa = xa0;
    var xb = xb0;
    var fa = F(xa);
    var fb = F(xb);
    if fa < fb {
        let tx = xa; xa = xb; xb = tx;
        let tf = fa; fa = fb; fb = tf;
    }
    var xc = xb + GOLD * (xb - xa);
    var fc = F(xc);
    var iter: i64 = 0;
    while fc < fb {
        let tmp1 = (xb - xa) * (fb - fc);
        let tmp2 = (xb - xc) * (fb - fa);
        let val = tmp2 - tmp1;
        var denom = 2.0 * val;
        if abs(val) < VERY_SMALL { denom = 2.0 * VERY_SMALL; }
        var w = xb - ((xb - xc) * tmp2 - (xb - xa) * tmp1) / denom;
        let wlim = xb + BRENT_GROW * (xc - xb);
        if iter > MAX_FIT_STEPS * 5 { return nan(), nan(), nan(), ERR_FIT; }
        iter += 1;
        var fw: f64 = 0.0;
        if (w - xc) * (xb - w) > 0.0 {
            fw = F(w);
            if fw < fc {
                xa = xb; xb = w; fa = fb; fb = fw;
                return xa, xb, xc, 0;
            } else if fw > fb {
                xc = w; fc = fw;
                return xa, xb, xc, 0;
            }
            w = xc + GOLD * (xc - xb);
            fw = F(w);
        } else if (w - wlim) * (wlim - xc) >= 0.0 {
            w = wlim;
            fw = F(w);
        } else if (w - wlim) * (xc - w) > 0.0 {
            fw = F(w);
            if fw < fc {
                xb = xc; xc = w; w = xc + GOLD * (xc - xb);
                fb = fc; fc = fw; fw = F(w);
            }
        } else {
            w = xc + GOLD * (xc - xb);
            fw = F(w);
        }
        xa = xb; xb = xc; xc = w;
        fa = fb; fb = fc; fc = fw;
        if is_nan(fc) || is_nan(fb) { return nan(), nan(), nan(), ERR_FIT; }
    }
    return xa, xb, xc, 0;
}

// The minimum of F inside a bracket by the method of Brent, as SciPy
// does it. The tolerance is relative. The result is the point of the
// minimum. ERR_FIT means the search used all its steps.
fn brent_min<F>(xa: f64, xb: f64, xc: f64, tol: f64) -> f64, i64 {
    var x = xb;
    var w = xb;
    var v = xb;
    var fx = F(x);
    var fw = fx;
    var fv = fx;
    var a = xc;
    var b = xa;
    if xa < xc { a = xa; b = xc; }
    var deltax: f64 = 0.0;
    var rat: f64 = 0.0;
    var iter: i64 = 0;
    while iter < MAX_FIT_STEPS * 2 {
        let tol1 = tol * abs(x) + BRENT_MINTOL;
        let tol2 = 2.0 * tol1;
        let xmid = 0.5 * (a + b);
        if abs(x - xmid) < tol2 - 0.5 * (b - a) { return x, 0; }
        if abs(deltax) <= tol1 {
            if x >= xmid { deltax = a - x; } else { deltax = b - x; }
            rat = BRENT_CG * deltax;
        } else {
            let t1 = (x - w) * (fx - fv);
            var t2 = (x - v) * (fx - fw);
            var p = (x - v) * t2 - (x - w) * t1;
            t2 = 2.0 * (t2 - t1);
            if t2 > 0.0 { p = -p; }
            t2 = abs(t2);
            let dxtemp = deltax;
            deltax = rat;
            if p > t2 * (a - x) && p < t2 * (b - x) && abs(p) < abs(0.5 * t2 * dxtemp) {
                rat = p / t2;
                let u0 = x + rat;
                if u0 - a < tol2 || b - u0 < tol2 {
                    if xmid - x >= 0.0 { rat = tol1; } else { rat = -tol1; }
                }
            } else {
                if x >= xmid { deltax = a - x; } else { deltax = b - x; }
                rat = BRENT_CG * deltax;
            }
        }
        var u = x + rat;
        if abs(rat) < tol1 {
            if rat >= 0.0 { u = x + tol1; } else { u = x - tol1; }
        }
        let fu = F(u);
        if is_nan(fu) { return nan(), ERR_FIT; }
        if fu > fx {
            if u < x { a = u; } else { b = u; }
            if fu <= fw || w == x {
                v = w; w = u; fv = fw; fw = fu;
            } else if fu <= fv || v == x || v == w {
                v = u; fv = fu;
            }
        } else {
            if u >= x { a = x; } else { b = x; }
            v = w; w = x; x = u;
            fv = fw; fw = fx; fx = fu;
        }
        iter += 1;
    }
    return nan(), ERR_FIT;
}

// The minimum of F from two start points: a bracket, then Brent.
fn minimize_scalar<F>(xa: f64, xb: f64, tol: f64) -> f64, i64 {
    let a, b, c, cb = bracket_min(xa, xb) { x => F(x) };
    if cb != 0 { return nan(), cb; }
    let x, cx = brent_min(a, b, c, tol) { x => F(x) };
    return x, cx;
}

// A root of F between a and b, where F(a) and F(b) differ in sign, by
// the method of Brent (the zeroin form). A pair without a sign change
// gives ERR_DOMAIN. ERR_FIT means the search used all its steps.
fn brent_root<F>(a0: f64, b0: f64, xtol: f64, rtol: f64) -> f64, i64 {
    var a = a0;
    var b = b0;
    var fa = F(a);
    var fb = F(b);
    if is_nan(fa) || is_nan(fb) { return nan(), ERR_NAN; }
    if fa == 0.0 { return a, 0; }
    if fb == 0.0 { return b, 0; }
    if fa * fb > 0.0 { return nan(), ERR_DOMAIN; }
    var c = a;
    var fc = fa;
    var d = b - a;
    var e = d;
    var iter: i64 = 0;
    while iter < MAX_FIT_STEPS {
        if fb * fc > 0.0 {
            c = a; fc = fa; d = b - a; e = d;
        }
        if abs(fc) < abs(fb) {
            a = b; b = c; c = a;
            fa = fb; fb = fc; fc = fa;
        }
        let tol = 2.0 * rtol * abs(b) + xtol;
        let m = 0.5 * (c - b);
        if abs(m) <= tol || fb == 0.0 { return b, 0; }
        if abs(e) < tol || abs(fa) <= abs(fb) {
            d = m; e = m;
        } else {
            var p: f64 = 0.0;
            var q: f64 = 0.0;
            let s = fb / fa;
            if a == c {
                p = 2.0 * m * s;
                q = 1.0 - s;
            } else {
                let qq = fa / fc;
                let r = fb / fc;
                p = s * (2.0 * m * qq * (qq - r) - (b - a) * (r - 1.0));
                q = (qq - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if p > 0.0 { q = -q; } else { p = -p; }
            if 2.0 * p < min(3.0 * m * q - abs(tol * q), abs(e * q)) {
                e = d; d = p / q;
            } else {
                d = m; e = m;
            }
        }
        a = b; fa = fb;
        if abs(d) > tol { b = b + d; } else if m > 0.0 { b = b + tol; } else { b = b - tol; }
        fb = F(b);
        if is_nan(fb) { return nan(), ERR_NAN; }
        iter += 1;
    }
    return nan(), ERR_FIT;
}

// The nodes and the weights of the Gauss-Legendre rule of n points on
// the interval from -1 to 1. Newton steps on the Legendre polynomial find
// each node. n is from 2 to 256.
fn gauss_legendre(xs: f64[>..]&, ws: f64[>..]&, n: i64) -> i64 {
    xs.resize(0);
    ws.resize(0);
    if n < 2 || n > 256 { return ERR_DOMAIN; }
    for i in 0..n { xs.push(0.0); ws.push(0.0); }
    let pi = 3.14159265358979323846;
    for i in 0..(n + 1) / 2 {
        var z = cos(pi * (i + 0.75) / (n + 0.5));
        var dp: f64 = 1.0;
        for step in 0..100 {
            var p1: f64 = 1.0;
            var p2: f64 = 0.0;
            for j in 0..n {
                let p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j + 1.0) * z * p2 - j * p3) / (j + 1.0);
            }
            dp = n * (z * p1 - p2) / (z * z - 1.0);
            let z1 = z;
            z = z1 - p1 / dp;
            if abs(z - z1) < 1e-15 { break; }
        }
        xs[i] = -z;
        xs[n - 1 - i] = z;
        let w = 2.0 / ((1.0 - z * z) * dp * dp);
        ws[i] = w;
        ws[n - 1 - i] = w;
    }
    return 0;
}
