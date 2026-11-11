// The normality tests. Each call asks whether one sample departs from
// the normal shape, in the skewness, in the kurtosis or as a whole.
// Each call returns a TestResult and then a fault code. A fault gives
// NaN in every field of the result.
//
// The names follow SciPy. The component tests carry the three
// alternatives; the omnibus tests and the Anderson-Darling test give
// only the two-sided form, as the references do.
namespace statistics;

import std;
import math;

// The p-value of a standard normal z in the alternative alt. The
// two-sided form covers both tails. A fault gives NaN.
fn norm_z_p(z: f64, alt: i64) -> f64, i64 {
    if alt == ALT_TWO_SIDED {
        let upper, c1 = norm_sf(abs(z));
        if c1 != 0 { return nan(), c1; }
        return 2.0 * upper, 0;
    }
    if alt == ALT_LESS {
        return norm_cdf(z);
    }
    return norm_sf(z);
}

// The z value of the skewness test, Anscombe and Glynn 1983. A constant
// sample faults. The call needs 8 values, as SciPy and statrs do.
fn skew_z(xs: const f64[:]) -> f64, i64 {
    let b1, c = skew(xs, true);
    if c != 0 { return nan(), c; }
    let n: f64 = xs.len;
    var y = b1 * sqrt(((n + 1.0) * (n + 3.0)) / (6.0 * (n - 2.0)));
    // A symmetric sample gives y = 0. SciPy and statrs both replace it by 1
    // before the transform, and so does this call.
    if y == 0.0 { y = 1.0; }
    let beta2 = (3.0 * (n * n + 27.0 * n - 70.0) * (n + 1.0) * (n + 3.0))
        / ((n - 2.0) * (n + 5.0) * (n + 7.0) * (n + 9.0));
    let w2 = -1.0 + sqrt(2.0 * (beta2 - 1.0));
    // For 8 values or more, w2 stays above 1, so the logarithm is positive.
    let delta = 1.0 / sqrt(0.5 * log(w2));
    let alpha = sqrt(2.0 / (w2 - 1.0));
    let t = y / alpha;
    return delta * log(t + sqrt(t * t + 1.0)), 0;
}

// The z value of the kurtosis test, Anscombe and Glynn 1983. The call
// needs 5 values, as SciPy does. A zero denominator faults.
fn kurt_z(xs: const f64[:]) -> f64, i64 {
    let b2, c = kurtosis(xs, false, true);
    if c != 0 { return nan(), c; }
    let n: f64 = xs.len;
    let e = 3.0 * (n - 1.0) / (n + 1.0);
    let varb2 = 24.0 * n * (n - 2.0) * (n - 3.0)
        / ((n + 1.0) * (n + 1.0) * (n + 3.0) * (n + 5.0));
    let x = (b2 - e) / sqrt(varb2);
    let sqrtbeta1 = 6.0 * (n * n - 5.0 * n + 2.0) / ((n + 7.0) * (n + 9.0))
        * sqrt(6.0 * (n + 3.0) * (n + 5.0) / (n * (n - 2.0) * (n - 3.0)));
    let a = 6.0 + 8.0 / sqrtbeta1
        * (2.0 / sqrtbeta1 + sqrt(1.0 + 4.0 / (sqrtbeta1 * sqrtbeta1)));
    let denom = 1.0 + x * sqrt(2.0 / (a - 4.0));
    if denom == 0.0 { return nan(), ERR_UNDEFINED; }
    // The cube root of the signed quotient is the sign of the denominator
    // times the cube root of the magnitude, as SciPy computes it. The call
    // uses pow of the magnitude, because pow of a negative base gives NaN
    // and the older Goose releases hold no cbrt.
    let q = (1.0 - 2.0 / a) / denom;
    var term2 = pow(abs(q), 1.0 / 3.0);
    if q < 0.0 { term2 = -term2; }
    let term1 = 1.0 - 2.0 / (9.0 * a);
    return (term1 - term2) / sqrt(2.0 / (9.0 * a)), 0;
}

// The test of the skewness against the normal shape. The sample needs
// at least 8 values.
fn skewtest(xs: const f64[:], alt: i64 = ALT_TWO_SIDED) -> TestResult, i64 {
    let c = check_sample(xs, 8);
    if c != 0 { return fault_result(c); }
    if check_alt(alt) != 0 { return fault_result(ERR_DOMAIN); }
    let z, cz = skew_z(xs);
    if cz != 0 { return fault_result(cz); }
    let p, cp = norm_z_p(z, alt);
    if cp != 0 { return fault_result(cp); }
    return good_result(z, p, nan(), nan());
}

// The test of the kurtosis against the normal shape. The sample needs
// at least 5 values.
fn kurtosistest(xs: const f64[:], alt: i64 = ALT_TWO_SIDED) -> TestResult, i64 {
    let c = check_sample(xs, 5);
    if c != 0 { return fault_result(c); }
    if check_alt(alt) != 0 { return fault_result(ERR_DOMAIN); }
    let z, cz = kurt_z(xs);
    if cz != 0 { return fault_result(cz); }
    let p, cp = norm_z_p(z, alt);
    if cp != 0 { return fault_result(cp); }
    return good_result(z, p, nan(), nan());
}

// The omnibus test of D'Agostino and Pearson. The statistic is the sum
// of the squares of the two z values, and the p-value comes from the
// chi-square upper tail with 2 degrees of freedom. The sample needs at
// least 8 values, as SciPy does.
fn normaltest(xs: const f64[:]) -> TestResult, i64 {
    let c = check_sample(xs, 8);
    if c != 0 { return fault_result(c); }
    let zs, c1 = skew_z(xs);
    if c1 != 0 { return fault_result(c1); }
    let zk, c2 = kurt_z(xs);
    if c2 != 0 { return fault_result(c2); }
    let k2 = zs * zs + zk * zk;
    return chi_result(k2, 2.0);
}

// The test of Jarque and Bera. The statistic is n over 6, times the
// square of the biased skewness plus a quarter of the square of the
// biased excess kurtosis. The p-value uses the same chi-square tail.
// The sample needs at least 2 values, as SciPy allows.
fn jarque_bera(xs: const f64[:]) -> TestResult, i64 {
    let c = check_sample(xs, 2);
    if c != 0 { return fault_result(c); }
    let b1, c1 = skew(xs, true);
    if c1 != 0 { return fault_result(c1); }
    let b2, c2 = kurtosis(xs, true, true);
    if c2 != 0 { return fault_result(c2); }
    let n: f64 = xs.len;
    let jb = n / 6.0 * (b1 * b1 + b2 * b2 / 4.0);
    return chi_result(jb, 2.0);
}

// The p-value of the Anderson-Darling statistic, the asymptotic rule of
// D'Agostino and Stephens adjusted by the sample size. The rule holds a
// value in [0,1] over the whole range of the statistic.
fn ad_p_value(a2: f64, n: f64) -> f64 {
    let adj = a2 * (1.0 + 0.75 / n + 2.25 / (n * n));
    // The quadratic of the top branch turns at adj 153.4. Above the turn
    // the rule rises back toward 1 and then above it, so statrs and R
    // report a value above 1 there. A statistic that large means a
    // p-value of 0, and this rule returns 0.
    if adj >= 153.4 { return 0.0; }
    if adj >= 0.6 {
        return exp(1.2937 - 5.709 * adj + 0.0186 * adj * adj);
    }
    if adj >= 0.34 {
        return exp(0.9177 - 4.279 * adj - 1.38 * adj * adj);
    }
    if adj >= 0.2 {
        return 1.0 - exp(-8.318 + 42.796 * adj - 59.938 * adj * adj);
    }
    return 1.0 - exp(-13.436 + 101.14 * adj - 223.73 * adj * adj);
}

// The logarithm of one standard normal tail. Beyond 30 standard units
// the tail underflows, so the logarithm splits the scaled tail from the
// exponential: log(0.5 q) - z^2 / 2. The other tail rounds to 1, and its
// logarithm is 0.
fn log_tail(z: f64, upper: bool) -> f64 {
    let a = abs(z);
    if a > 30.0 {
        if upper != (z > 0.0) { return 0.0; }
        return log(0.5 * scaled_q(a * INV_SQRT2)) - 0.5 * a * a;
    }
    if upper {
        let s, c = norm_sf(z);
        if c != 0 { return nan(); }
        return log(s);
    }
    let v, c = norm_cdf(z);
    if c != 0 { return nan(); }
    return log(v);
}

// The test of Anderson and Darling against the normal distribution with
// the mean and the unbiased scale of the sample. The sample needs at
// least 5 values, the documented range of the asymptotic p-value rule.
fn anderson_darling(xs: const f64[:]) -> TestResult, i64 {
    let c = check_sample(xs, 5);
    if c != 0 { return fault_result(c); }
    let n: f64 = xs.len;
    // The mean and the unbiased scale, from the shared helpers. They scale
    // the sample before they square it, so a spread that does not fit
    // faults with ERR_OVERFLOW.
    let mu, cm = exact_mean(xs);
    if cm != 0 { return fault_result(cm); }
    let scale2, cs = variance_core(xs, 1);
    if cs != 0 { return fault_result(cs); }
    if scale2 <= 0.0 { return fault_result(ERR_UNDEFINED); }
    let sigma = sqrt(scale2);
    // The sum runs once over the sorted copy, one term at a time, in the
    // order that matches the reference value. Each tail enters through
    // its logarithm. Beyond 30 standard units the tail underflows, so the
    // term uses the scaled tail.
    let sorted = sort_copy(xs);
    var sum = 0.0;
    for i in 0..sorted.len {
        let z1 = (sorted[i] - mu) / sigma;
        let z2 = (sorted[sorted.len - 1 - i] - mu) / sigma;
        let w: f64 = 2 * (i + 1) - 1;
        sum += w * (log_tail(z1, false) + log_tail(z2, true));
    }
    let a2 = -n - sum / n;
    return good_result(a2, ad_p_value(a2, n), nan(), nan());
}
