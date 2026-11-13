// The Kolmogorov-Smirnov tests of the package. Each call returns a
// TestResult and then a fault code. A fault gives NaN in every field of
// the result.
namespace statistics;

import std;
import math;

fn smirnov_sf(n: i64, x: f64) -> f64 {
    if x <= 0.0 { return 1.0; }
    if x >= 1.0 { return 0.0; }
    let nf = n as f64;
    let jmax = floor(nf * (1.0 - x)) as i64;
    var lc: f64 = 0.0;
    var sum = acc_new();
    for j in 0..jmax + 1 {
        if j > 0 { lc += log((nf - j + 1.0) / j); }
        let a = 1.0 - x - j / nf;
        if a <= 0.0 { continue; }
        let b = x + j / nf;
        acc_add(sum, exp(lc + (nf - j) * log(a) + (j - 1.0) * log(b)));
    }
    return min(1.0, x * acc_total(sum));
}

// The product of two square matrices of side m.
fn mat_mul(out: f64[>..]&, a: const f64[:], b: const f64[:], m: i64) {
    out.resize(0);
    for i in 0..m * m { out.push(0.0); }
    for i in 0..m {
        for k in 0..m {
            let v = a[i * m + k];
            if v == 0.0 { continue; }
            for j in 0..m { out[i * m + j] += v * b[k * m + j]; }
        }
    }
}

// P(D < d) of the two-sided statistic of n values, by the matrix power
// of Marsaglia, Tsang and Wang. The side of the matrix is 2 floor(n d) + 1.
fn mtw_cdf(n: i64, d: f64) -> f64 {
    let nf = n as f64;
    let k = floor(nf * d) as i64 + 1;
    let m = 2 * k - 1;
    let h = k - nf * d;
    var a: f64[>..] = [];
    for i in 0..m {
        for j in 0..m {
            if i - j + 1 < 0 { a.push(0.0); } else { a.push(1.0); }
        }
    }
    for i in 0..m {
        a[i * m] -= pow(h, i + 1.0);
        a[(m - 1) * m + i] -= pow(h, (m - i) * 1.0);
    }
    if 2.0 * h - 1.0 > 0.0 { a[(m - 1) * m] += pow(2.0 * h - 1.0, m * 1.0); }
    for i in 0..m {
        for j in 0..m {
            if i - j + 1 > 0 {
                for g in 1..i - j + 2 { a[i * m + j] /= g; }
            }
        }
    }
    var v: f64[>..] = [];
    for x in a { v.push(x); }
    var ev: i64 = 0;
    var top: i64 = 0;
    while (n >> (top + 1)) > 0 { top += 1; }
    var bit = top - 1;
    var tmp: f64[>..] = [];
    while bit >= 0 {
        mat_mul(tmp, v, v, m);
        ev = 2 * ev;
        v.resize(0);
        for x in tmp { v.push(x); }
        if ((n >> bit) & 1) == 1 {
            mat_mul(tmp, a, v, m);
            v.resize(0);
            for x in tmp { v.push(x); }
        }
        if v[(m / 2) * m + m / 2] > 1e140 {
            for i in 0..m * m { v[i] *= 1e-140; }
            ev += 140;
        }
        bit -= 1;
    }
    var s = v[(k - 1) * m + k - 1];
    for i in 1..n + 1 {
        s = s * i / nf;
        if s < 1e-140 {
            s *= 1e140;
            ev -= 140;
        }
    }
    return s * pow(10.0, ev as f64);
}

// The Pelz-Good series for the two-sided statistic: the upper tail.
fn pelz_good_sf(n: i64, x: f64) -> f64 {
    let nf = n as f64;
    let pi2 = 9.869604401089358;
    let z = sqrt(nf) * x;
    let z2 = z * z;
    let z3 = z2 * z;
    let z4 = z2 * z2;
    let z6 = z4 * z2;
    let z8 = z4 * z4;
    let qlog = -pi2 / 8.0 / z2;
    if qlog < -708.0 { return 1.0; }
    let q = exp(qlog);
    let pi4 = pi2 * pi2;
    let pi6 = pi4 * pi2;
    let k1a = -z2;
    let k1b = pi2 / 4.0;
    let k2a = 6.0 * z6 + 2.0 * z4;
    let k2b = (2.0 * z4 - 5.0 * z2) * pi2 / 4.0;
    let k2c = pi4 * (1.0 - 2.0 * z2) / 16.0;
    let k3d = pi6 * (5.0 - 30.0 * z2) / 64.0;
    let k3c = pi4 * (-60.0 * z2 + 212.0 * z4) / 16.0;
    let k3b = pi2 * (135.0 * z4 - 96.0 * z6) / 4.0;
    let k3a = -30.0 * z6 - 90.0 * z8;
    let maxk = ceil(16.0 * z / 3.141592653589793) as i64;
    var s0: f64 = 0.0;
    var s1: f64 = 0.0;
    var s2: f64 = 0.0;
    var s3: f64 = 0.0;
    var k = maxk;
    while k >= 1 {
        let m = (2 * k - 1) as f64;
        let ms = m * m;
        let m4 = ms * ms;
        let m6 = m4 * ms;
        let qp = pow(q, 8.0 * k);
        s0 = s0 * qp + 1.0;
        s1 = s1 * qp + k1a + k1b * ms;
        s2 = s2 * qp + k2a + k2b * ms + k2c * m4;
        s3 = s3 * qp + k3a + k3b * ms + k3c * m4 + k3d * m6;
        k -= 1;
    }
    let root = sqrt(2.0 * 3.141592653589793);
    let c0 = s0 * q * root / z;
    var c1 = s1 * q * root / (6.0 * z4);
    var c2 = s2 * q * root / (72.0 * z6 * z);
    var c3 = s3 * q * root / (6480.0 * z8 * z2);
    let q2 = exp(-pi2 / 2.0 / z2);
    let sq3 = sqrt(3.0) * z;
    var e2: f64 = 0.0;
    var e3: f64 = 0.0;
    for j in 1..maxk + 1 {
        let jf = j as f64;
        let w = jf * jf * pow(q2, jf * jf);
        e2 += w;
        e3 += (sq3 + 3.141592653589793 * jf) * (sq3 - 3.141592653589793 * jf) * w;
    }
    c2 += e2 * pi2 * root / (-36.0 * z3);
    c3 += e3 * pi2 * root / (216.0 * z6);
    c1 = c1 / sqrt(nf);
    c2 = c2 / nf;
    c3 = c3 / (nf * sqrt(nf));
    return (1.0 - c0) - c1 - c2 - c3;
}

fn clip01(p: f64) -> f64 {
    if p < 0.0 { return 0.0; }
    if p > 1.0 { return 1.0; }
    return p;
}

// The upper tail of the two-sided statistic of n values: P(D >= x).
fn kolmogn_sf(n: i64, x: f64) -> f64 {
    if x >= 1.0 { return 0.0; }
    if x <= 0.0 { return 1.0; }
    let nf = n as f64;
    let t = nf * x;
    if t <= 1.0 {
        if t <= 0.5 { return 1.0; }
        var prob: f64 = 1.0;
        if n <= 140 {
            for i in 1..n + 1 { prob = prob * i / nf * (2.0 * t - 1.0); }
        } else {
            let lf, sg, cl = lgamma(nf + 1.0);
            prob = exp(lf - nf * log(nf) + nf * log(2.0 * t - 1.0));
        }
        return clip01(1.0 - prob);
    }
    if t >= nf - 1.0 { return clip01(2.0 * pow(1.0 - x, nf)); }
    if x >= 0.5 { return clip01(2.0 * smirnov_sf(n, x)); }
    let nx2 = t * x;
    if n <= 140 {
        if nx2 <= 4.0 { return clip01(1.0 - mtw_cdf(n, x)); }
        return clip01(2.0 * smirnov_sf(n, x));
    }
    if nx2 >= 370.0 { return 0.0; }
    if nx2 >= 2.2 { return clip01(2.0 * smirnov_sf(n, x)); }
    if n <= 100000 && nf * pow(x, 1.5) <= 1.4 { return clip01(1.0 - mtw_cdf(n, x)); }
    return clip01(pelz_good_sf(n, x));
}

// The greatest common divisor of two positive values.
// The probability that a lattice path stays inside the band of width h
// around the diagonal from (0, 0) to (m, n). The count of the paths is
// scaled at each column, so the values stay in 0 to 1. This is the
// method of Hodges.
fn ks_inside(m: i64, n: i64, g: i64, h: i64, two_sided: bool) -> f64 {
    let mg = m / g;
    let ng = n / g;
    var minj: i64 = 0;
    var maxj = min(ceil(h as f64 / mg) as i64, n + 1);
    var curlen = maxj - minj;
    var len_a = min(2 * maxj + 2, n + 1);
    var a: f64[>..] = [];
    for i in 0..len_a { a.push(1.0); }
    for i in minj..maxj { a[i] = 0.0; }
    for i in 1..m + 1 {
        let last_minj = minj;
        let last_len = curlen;
        minj = max(floor((ng * i - h) as f64 / mg) as i64 + 1, 0);
        if minj > n { minj = n; }
        if two_sided {
            maxj = min(ceil((ng * i + h) as f64 / mg) as i64, n + 1);
        } else {
            maxj = n + 1;
        }
        if maxj <= minj { return 1.0; }
        var val: f64 = 0.0;
        if minj != 0 { val = 1.0; }
        for jj in minj..maxj {
            let j = jj;
            val = (a[j - last_minj] * i + val * j) / (i + j);
            a[jj - minj] = val;
        }
        curlen = maxj - minj;
        if last_len > curlen {
            for i in curlen..last_len { a[i] = 1.0; }
        }
    }
    return a[maxj - minj - 1];
}

// The two-sided exact p-value of the two-sample test. The two forms
// follow the source of SciPy: the square for equal sizes, the band for
// the others.
fn ks_two_sided_exact(d: f64, n1: i64, n2: i64, g: i64) -> f64 {
    let lcm = (n1 / g) * n2;
    let h = round(d * lcm) as i64;
    if h == 0 { return 1.0; }
    if n1 == n2 { return ks_outside_square(n1, h); }
    return ks_inside(n1, n2, g, h, true);
}

// The probability that a path leaves the band of width h around the
// diagonal of two equal samples.
fn ks_outside_square(n: i64, h: i64) -> f64 {
    var p: f64 = 0.0;
    var k = n / h;
    while k >= 0 {
        var p1: f64 = 1.0;
        for j in 0..h {
            p1 = (n - k * h - j) * p1 / (n + k * h + j + 1);
        }
        p = p1 * (1.0 - p);
        k -= 1;
    }
    return 2.0 * p;
}

// The binomial coefficient of a and b as a float. The loop runs b times.
fn binom_f(a: i64, b: i64) -> f64 {
    if b < 0 || b > a { return 0.0; }
    let k = if b > a - b { a - b } else { b };
    var r: f64 = 1.0;
    for i in 1..k + 1 { r = r * (a - k + i) / i; }
    return r;
}

// The number of lattice paths that touch the one-sided boundary at height
// h, for samples of m >= n values. This is the count of SciPy.
fn ks_paths_outside(m: i64, n: i64, g: i64, h: i64) -> f64 {
    let mg = m / g;
    let ng = n / g;
    // The floor of (mg - h) / mg, as the quotient of a negative number.
    var q = (mg - h) / mg;
    if (mg - h) % mg != 0 && mg - h < 0 { q -= 1; }
    let lxj = n + q;
    if lxj <= 0 { return binom_f(m + n, n); }
    var xj: i64[>..] = [];
    for j in 0..lxj { xj.push((h + mg * j + ng - 1) / ng); }
    var bs: f64[>..] = [];
    for j in 0..lxj { bs.push(0.0); }
    bs[0] = 1.0;
    for j in 1..lxj {
        var bj = binom_f(xj[j] + j, j);
        for i in 0..j { bj -= binom_f(xj[j] - xj[i] + j - i, j - i) * bs[i]; }
        bs[j] = bj;
    }
    var num: f64 = 0.0;
    for j in 0..lxj { num += bs[j] * binom_f(m - xj[j] + n - j, n - j); }
    return num;
}

// The one-sided exact p-value of the two-sample test.
fn ks_one_sided_exact(d: f64, n1: i64, n2: i64, g: i64) -> f64 {
    let lcm = (n1 / g) * n2;
    let h = round(d * lcm) as i64;
    if h == 0 { return 1.0; }
    if n1 == n2 {
        var p: f64 = 1.0;
        for j in 0..h { p = p * (n1 - j) / (n1 + j + 1); }
        return p;
    }
    return ks_paths_outside(n1, n2, g, h) / binom_f(n1 + n2, n2);
}

// The two-sample Kolmogorov-Smirnov test. The statistic is the largest
// gap between the two step functions.
fn ks_2samp(xs: const f64[:], ys: const f64[:], alt: i64 = ALT_TWO_SIDED) -> TestResult, i64 {
    var c = check_sample(xs, 1);
    if c != 0 { return fault_result(c); }
    c = check_sample(ys, 1);
    if c != 0 { return fault_result(c); }
    if check_alt(alt) != 0 { return fault_result(ERR_DOMAIN); }
    let n1 = xs.len;
    let n2 = ys.len;
    let sx = sort_copy(xs);
    let sy = sort_copy(ys);
    var dplus: f64 = 0.0;
    var dminus: f64 = 0.0;
    var i: i64 = 0;
    var j: i64 = 0;
    while i < n1 || j < n2 {
        var v: f64 = 0.0;
        if i < n1 && j < n2 {
            v = if sx[i] <= sy[j] { sx[i] } else { sy[j] };
        } else if i < n1 {
            v = sx[i];
        } else {
            v = sy[j];
        }
        while i < n1 && sx[i] == v { i += 1; }
        while j < n2 && sy[j] == v { j += 1; }
        let gap: f64 = i / (n1 as f64) - j / (n2 as f64);
        if gap > dplus { dplus = gap; }
        if -gap > dminus { dminus = -gap; }
    }
    var d = max(dplus, dminus);
    if alt == ALT_GREATER { d = dplus; }
    if alt == ALT_LESS { d = dminus; }
    if d == 0.0 { return good_result(0.0, 1.0, nan(), nan()); }
    let g = gcd_counts(n1, n2);
    let big = max(n1, n2);
    let small = min(n1, n2);
    let product = n1 * n2;
    let tx, _sx0, _sx1 = tie_sums(sx);
    let ty, _sy0, _sy1 = tie_sums(sy);
    let tied = tx > 0.0 || ty > 0.0;
    var p: f64 = 1.0;
    if product <= MAX_EXACT_KS && !tied {
        if alt == ALT_TWO_SIDED {
            p = ks_two_sided_exact(d, big, small, g);
        } else {
            p = ks_one_sided_exact(d, big, small, g);
        }
        if alt == ALT_TWO_SIDED { p = min(1.0, p); }
        return good_result(d, p, nan(), nan());
    }
    let m = max(n1, n2) as f64;
    let n = min(n1, n2) as f64;
    let en = m * n / (m + n);
    if alt == ALT_TWO_SIDED {
        let q = kolmogn_sf(round(en) as i64, d);
        return good_result(d, q, nan(), nan());
    }
    let z = sqrt(en) * d;
    let expt = -2.0 * z * z - 2.0 * z * (m + 2.0 * n) / (3.0 * sqrt(m * n * (m + n)));
    p = exp(expt);
    return good_result(d, p, nan(), nan());
}

// The sweep of one sample against a family. It gives the two deviations,
// the location of the larger one and its sign. The call checks the
// sample, the family and the parameters first.
fn ks_1samp_sweep(xs: const f64[:], family: i64, params: const f64[:], alt: i64) -> f64, f64, f64, f64, i64 {
    var c = check_sample(xs, 1);
    if c != 0 { return nan(), nan(), nan(), nan(), c; }
    if check_alt(alt) != 0 { return nan(), nan(), nan(), nan(), ERR_DOMAIN; }
    if family < 0 || family > FAMILY_F { return nan(), nan(), nan(), nan(), ERR_DOMAIN; }
    let need = family_arity(family);
    if params.len != need { return nan(), nan(), nan(), nan(), ERR_LENGTH; }
    for p in params {
        let cp = check_scalar(p);
        if cp != 0 { return nan(), nan(), nan(), nan(), cp; }
    }
    let s = sort_copy(xs);
    let n = s.len;
    var dplus: f64 = 0.0;
    var dminus: f64 = 0.0;
    var lplus: f64 = s[0];
    var lminus: f64 = s[0];
    var k: i64 = 0;
    while k < n {
        let f, cf = family_cdf(family, s[k], params);
        if cf != 0 { return nan(), nan(), nan(), nan(), cf; }
        let up = (k + 1) as f64 / n - f;
        let down = f - k as f64 / n;
        if up > dplus { dplus = up; lplus = s[k]; }
        if down > dminus { dminus = down; lminus = s[k]; }
        k += 1;
    }
    return dplus, dminus, lplus, lminus, 0;
}

// The Kolmogorov-Smirnov test of one sample against a family. The family
// is one of the FAMILY_ constants. The buffer holds the parameters of
// that family.
fn ks_1samp(xs: const f64[:], family: i64, params: const f64[:], alt: i64 = ALT_TWO_SIDED) -> TestResult, i64 {
    let dplus, dminus, _lp, _lm, c = ks_1samp_sweep(xs, family, params, alt);
    if c != 0 { return fault_result(c); }
    let n = xs.len;
    var d = max(dplus, dminus);
    if alt == ALT_GREATER { d = dplus; }
    if alt == ALT_LESS { d = dminus; }
    if d == 0.0 { return good_result(0.0, 1.0, nan(), nan()); }
    var p: f64 = 1.0;
    if alt == ALT_TWO_SIDED {
        p = kolmogn_sf(n, d);
    } else {
        p = smirnov_sf(n, d);
    }
    return good_result(d, p, nan(), nan());
}

// The count of the parameters of a family.
fn family_arity(family: i64) -> i64 {
    if family == FAMILY_EXPONENTIAL { return 1; }
    if family == FAMILY_CHI2 { return 1; }
    if family == FAMILY_T { return 1; }
    return 2;
}

// The cumulative function of a family at one value. The value is the
// result of the call and the code is the fault of the call.
fn family_cdf(family: i64, x: f64, params: const f64[:]) -> f64, i64 {
    if family == FAMILY_NORMAL { let v, c = norm_cdf(x, params[0], params[1]); return v, c; }
    if family == FAMILY_UNIFORM { let v, c = unif_cdf(x, params[0], params[1]); return v, c; }
    if family == FAMILY_EXPONENTIAL { let v, c = expon_cdf(x, params[0]); return v, c; }
    if family == FAMILY_LOGNORMAL { let v, c = lognorm_cdf(x, params[0], params[1]); return v, c; }
    if family == FAMILY_WEIBULL { let v, c = weibull_cdf(x, params[0], params[1]); return v, c; }
    if family == FAMILY_GAMMA { let v, c = gammadist_cdf(x, params[0], params[1]); return v, c; }
    if family == FAMILY_BETA { let v, c = betadist_cdf(x, params[0], params[1]); return v, c; }
    if family == FAMILY_CHI2 { let v, c = chi2_cdf(x, params[0]); return v, c; }
    if family == FAMILY_T { let v, c = tdist_cdf(x, params[0]); return v, c; }
    let v, c = fdist_cdf(x, params[0], params[1]);
    return v, c;
}
