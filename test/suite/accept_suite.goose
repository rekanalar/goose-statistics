// One test for each documented case of the package. Each case states the
// sample and the answer that the reference libraries give.
namespace suite;

import std;
import statistics from "../..";

// Counts a failure and prints the case.
fn check(bad: i64&, ok: bool, what: u8[:]) {
    if ok { return; }
    bad += 1;
    print("FAIL ", what, "\n");
}

// A sample from a literal list.
fn list(vals: const f64[:]) -> f64[>..] {
    var out: f64[>..] = [];
    for v in vals { out.push(v); }
    return out;
}

fn run_accept() -> i64 {
    var bad: i64 = 0;

    // Location and spread.
    let a = list([1.0, 2.0, 3.0, 4.0]);
    let m, c = statistics::mean(a[..]);
    check(bad, m == 2.5 && c == 0, "mean of 1 2 3 4 is 2.5");
    let md, c2 = statistics::median(a[..]);
    check(bad, md == 2.5 && c2 == 0, "median of 1 2 3 4 is 2.5");
    let b = list([1.0, 3.0, 5.0, 7.0]);
    let lo, c3 = statistics::median_low(b[..]);
    let hi, c4 = statistics::median_high(b[..]);
    check(bad, lo == 3.0 && hi == 5.0 && c3 == 0 && c4 == 0, "low and high median of 1 3 5 7 are 3 and 5");
    let sd, c5 = statistics::stdev(a[..]);
    let pd, c6 = statistics::pstdev(a[..]);
    check(bad, within(sd, 1.2909944487358056, 2.0) && within(pd, 1.118033988749895, 2.0) && c5 == 0 && c6 == 0, "deviations of 1 2 3 4");
    let t = list([1.0, 1.0, 2.0, 2.0]);
    let mo, c7 = statistics::mode(t[..]);
    var all: f64[>..] = [];
    let cnt, c8 = statistics::multimode(all, t[..]);
    check(bad, mo == 1.0 && c7 == 0 && cnt == 2 && all[0] == 1.0 && all[1] == 2.0 && c8 == 0, "mode of 1 1 2 2 is 1, all modes 1 and 2");
    let g = list([1.0, 2.0, 4.0]);
    let gm, c9 = statistics::geometric_mean(g[..]);
    let hm, c10 = statistics::harmonic_mean(g[..]);
    check(bad, within(gm, 2.0, 8.0) && within(hm, 12.0 / 7.0, 4.0) && c9 == 0 && c10 == 0, "geometric 2 and harmonic 12/7 of 1 2 4");
    let off = list([1e9 + 4.0, 1e9 + 7.0, 1e9 + 13.0, 1e9 + 16.0]);
    let v, c11 = statistics::variance(off[..]);
    check(bad, v == 30.0 && c11 == 0, "variance of the large offset sample is exactly 30");
    let cancel = list([1e16, 1.0, -1e16]);
    let s, c12 = statistics::sum(cancel[..]);
    check(bad, s == 1.0 && c12 == 0, "sum of 1e16 1 -1e16 is 1");
    let grp = list([1.0, 2.0, 2.0, 3.0, 4.0, 4.0, 4.0, 4.0, 4.0, 5.0]);
    let gmed, c13 = statistics::median_grouped(grp[..], 1.0);
    check(bad, within(gmed, 3.7, 2.0) && c13 == 0, "grouped median is 3.7");
    let _w, c14 = statistics::median_grouped(grp[..], 0.0);
    check(bad, c14 == statistics::ERR_DOMAIN, "a class width of 0 is a domain fault");

    // Quantiles on 1 to 10.
    let ten = list([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]);
    let want: f64[9] = [3.0, 3.5, 3.0, 3.0, 3.5, 3.3, 3.6999999999999997, 3.4333333333333336, 3.45];
    for k in 0..9 {
        let q, cq = statistics::quantile(ten[..], 0.3, k);
        check(bad, within(q, want[k], 2.0) && cq == 0, "quantile 0.3 of 1 to 10 by one of the nine methods");
    }
    let ex, c15 = statistics::quantile(ten[..], 0.25, statistics::METHOD_EXCLUSIVE);
    let inc, c16 = statistics::quantile(ten[..], 0.25, statistics::METHOD_INCLUSIVE);
    check(bad, within(ex, 2.75, 2.0) && within(inc, 3.25, 2.0) && c15 == 0 && c16 == 0, "exclusive 2.75 and inclusive 3.25");
    let ew, c17 = statistics::quantile(ten[..], 0.25, statistics::METHOD_WEIBULL);
    let el, c18 = statistics::quantile(ten[..], 0.25, statistics::METHOD_LINEAR);
    check(bad, ew == ex && el == inc && c17 == 0 && c18 == 0, "exclusive equals weibull and inclusive equals linear");
    var cuts: f64[>..] = [];
    let nc, c19 = statistics::cut_points(cuts, ten[..], 4, statistics::METHOD_EXCLUSIVE);
    check(bad, nc == 3 && within(cuts[0], 2.75, 2.0) && within(cuts[1], 5.5, 2.0) && within(cuts[2], 8.25, 2.0) && c19 == 0, "exclusive cut points 2.75 5.5 8.25");
    let _q1, c20 = statistics::quantile(ten[..], -0.1, 6);
    let _q2, c21 = statistics::quantile(ten[..], 1.1, 6);
    let _q3, c22 = statistics::quantile(ten[..], 0.5, 9);
    check(bad, c20 == statistics::ERR_DOMAIN && c21 == statistics::ERR_DOMAIN && c22 == statistics::ERR_DOMAIN, "a fraction or a method out of range is a domain fault");

    // Shape.
    let sh = list([2.0, 8.0, 0.0, 4.0, 1.0, 9.0, 9.0, 0.0]);
    let sb, c23 = statistics::skew(sh[..], true);
    let su, c24 = statistics::skew(sh[..], false);
    check(bad, within(sb, 0.2650554122698573, 8.0) && within(su, 0.3305821804079746, 8.0) && c23 == 0 && c24 == 0, "skewness 0.265 and 0.331");
    let kf, c25 = statistics::kurtosis(sh[..], true, true);
    let kp, c26 = statistics::kurtosis(sh[..], false, true);
    check(bad, within(kf, -1.6660010752838508, 8.0) && within(kp, 1.3339989247161492, 8.0) && c25 == 0 && c26 == 0, "kurtosis Fisher and Pearson");
    let cm, c27 = statistics::central_moment(sh[..], 3);
    check(bad, within(cm, 13.67578125, 8.0) && c27 == 0, "central moment 3 is 13.67578125");
    let five = list([1.0, 2.0, 3.0, 4.0, 5.0]);
    var zs: f64[>..] = [];
    let nz, c28 = statistics::zscore(zs, five[..], 0);
    check(bad, nz == 5 && within(zs[0], -1.414213562373095, 4.0) && zs[2] == 0.0 && within(zs[4], 1.414213562373095, 4.0) && c28 == 0, "z-scores of 1 to 5");
    let same = list([3.0, 3.0, 3.0, 3.0]);
    let _s1, c29 = statistics::skew(same[..], true);
    let _s2, c30 = statistics::zscore(zs, same[..], 0);
    check(bad, c29 == statistics::ERR_UNDEFINED && c30 == statistics::ERR_UNDEFINED, "equal values give an undefined fault");

    // Trimmed and weighted forms.
    let zero_nine = list([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
    let tr, c31 = statistics::trim_mean(zero_nine[..], 0.1);
    check(bad, tr == 4.5 && c31 == 0, "trimmed mean of 0 to 9 is 4.5");
    let _t2, c32 = statistics::trim_mean(zero_nine[..], 0.5);
    check(bad, c32 == statistics::ERR_DOMAIN, "a cut of 0.5 is a domain fault");
    let tm, c33 = statistics::tmean(five[..], 2.0, 4.0);
    check(bad, tm == 3.0 && c33 == 0, "mean within limits 2 and 4 is 3");
    let x3 = list([1.0, 2.0, 3.0]);
    let w3 = list([3.0, 2.0, 1.0]);
    let wm, c34 = statistics::weighted_mean(x3[..], w3[..]);
    check(bad, within(wm, 1.6666666666666667, 2.0) && c34 == 0, "weighted mean is 1.6666666666666667");
    let w4 = list([1.0, 2.0, 3.0, 4.0]);
    let wf, c35 = statistics::weighted_variance(a[..], w4[..], statistics::WEIGHT_FREQUENCY);
    let wr, c36 = statistics::weighted_variance(a[..], w4[..], statistics::WEIGHT_RELIABILITY);
    let wp, c37 = statistics::weighted_variance(a[..], w4[..], statistics::WEIGHT_POPULATION);
    check(bad, within(wf, 1.1111111111111112, 4.0) && within(wr, 1.4285714285714284, 4.0) && within(wp, 1.0, 4.0) && c35 == 0 && c36 == 0 && c37 == 0, "weighted variance forms");
    let neg = list([1.0, -1.0, 1.0]);
    let _n1, c38 = statistics::weighted_mean(x3[..], neg[..]);
    check(bad, c38 == statistics::ERR_DOMAIN, "a negative weight is a domain fault");
    let zw = list([0.0, 0.0, 0.0]);
    let _n2, c39 = statistics::weighted_mean(x3[..], zw[..]);
    check(bad, c39 == statistics::ERR_UNDEFINED, "weights that add to zero give an undefined fault");
    let w2 = list([1.0, 2.0]);
    let _n3, c40 = statistics::weighted_mean(x3[..], w2[..]);
    check(bad, c40 == statistics::ERR_LENGTH, "a weight list of another length is a length fault");
    let nanw = list([1.0, 0.0 / 0.0, 1.0]);
    let _n4, c41 = statistics::weighted_mean(x3[..], nanw[..]);
    check(bad, c41 == statistics::ERR_NAN, "a NaN weight is a NaN fault");

    // Safe numerics.
    var empty: f64[>..] = [];
    let e1, c42 = statistics::mean(empty[..]);
    check(bad, is_nan(e1) && c42 == statistics::ERR_EMPTY, "an empty sample gives an empty fault and NaN");
    let one = list([5.0]);
    let _o1, c43 = statistics::variance(one[..]);
    let o2, c44 = statistics::pvariance(one[..]);
    check(bad, c43 == statistics::ERR_TOO_SHORT && o2 == 0.0 && c44 == 0, "one point: sample variance too short, population variance 0");
    let hasnan = list([1.0, 0.0 / 0.0]);
    let hasinf = list([1.0, 1.0 / 0.0]);
    let _f1, c45 = statistics::median(hasnan[..]);
    let _f2, c46 = statistics::median(hasinf[..]);
    check(bad, c45 == statistics::ERR_NAN && c46 == statistics::ERR_INFINITY, "a NaN and an infinity name their fault");
    let big = list([1e308, 1e308]);
    let bm, c47 = statistics::mean(big[..]);
    check(bad, bm == 1e308 && c47 == 0, "the mean of two 1e308 values is 1e308");
    let opp = list([1e308, -1e308]);
    let _b2, c48 = statistics::variance(opp[..]);
    check(bad, c48 == statistics::ERR_OVERFLOW, "the variance of 1e308 and -1e308 does not fit");
    var long: f64[>..] = [];
    for k in 0..statistics::MAX_INPUT + 1 { long.push(1.0); }
    let _l1, c49 = statistics::mean(long[..]);
    let _l2, c50 = statistics::median(long[..]);
    check(bad, c49 == statistics::ERR_TOO_LONG && c50 == statistics::ERR_TOO_LONG, "a sample above the bound faults before it reads data");
    // A sum that passes the largest double must still give the mean.
    let near = list([9.9e307, 9.9e307]);
    let nm, c51 = statistics::mean(near[..]);
    check(bad, nm == 9.9e307 && c51 == 0, "the mean of two 9.9e307 values is 9.9e307");
    let pair = list([9.9e307, 9e307]);
    let pm, c57b = statistics::mean(pair[..]);
    check(bad, pm == 9.45e307 && c57b == 0, "the mean of 9.9e307 and 9e307 is 9.45e307");
    let nz = list([2.0, 0.0, 4.0]);
    let hz, c52 = statistics::harmonic_mean(nz[..]);
    check(bad, hz == 0.0 && c52 == 0, "the harmonic mean with a zero is 0");
    let hn = list([2.0, -0.5, 4.0]);
    let _hn1, c53 = statistics::harmonic_mean(hn[..]);
    let _hn2, c54 = statistics::geometric_mean(hn[..]);
    check(bad, c53 == statistics::ERR_DOMAIN && c54 == statistics::ERR_DOMAIN, "a negative value gives a domain fault in both ratio means");
    let cm, c55 = statistics::central_moment(a[..], 1);
    check(bad, cm == 0.0 && c55 == 0, "the first central moment is 0");
    let tenth = list([0.1, 0.1, 0.1]);
    let tm, c56 = statistics::mean(tenth[..]);
    check(bad, tm == 0.1 && c56 == 0, "the mean of equal values is that value");
    var fifty: f64[>..] = [];
    for k in 0..49 { fifty.push(0.1); }
    let fm, c57 = statistics::mean(fifty[..]);
    check(bad, fm == 0.1 && c57 == 0, "the mean of 49 equal values is that value");
    // A tie in the last bit breaks toward the digits that lie below it.
    let tie1 = list([1.0, 1.1102230246251565e-16, 1e-30]);
    let tie2 = list([1e-30, 1.1102230246251565e-16, 1.0]);
    let tie3 = list([-1.0, -1.1102230246251565e-16, -1e-30]);
    let s1, d1 = statistics::sum(tie1[..]);
    let s2, d2 = statistics::sum(tie2[..]);
    let s3, d3 = statistics::sum(tie3[..]);
    check(bad, s1 == 1.0000000000000002 && s2 == s1 && s3 == -s1 && d1 + d2 + d3 == 0, "a tie in the sum rounds by the lower partials");
    // A sum whose last bit is a tie keeps the digit of the larger part.
    let tied = list([16777216.0, 3.725290298461914e-09, 9.313225746154785e-10, 33554432.0]);
    let ts, d4 = statistics::sum(tied[..]);
    check(bad, ts == 50331648.00000001 && d4 == 0, "a tie in the sum of four parts keeps the larger digit");
    // Many partial sums: values spaced across the whole exponent range.
    // Four copies of a block sum to four times the sum of the block, and the
    // exact sum does not depend on the order of the values.
    var part: f64[>..] = [];
    var spaced: f64[>..] = [];
    var backward: f64[>..] = [];
    for e in 0..40 {
        part.push(pow(2.0, -1070.0 + 52.0 * (e as f64)));
    }
    for r in 0..4 {
        for v in part {
            spaced.push(v);
        }
    }
    for k in 0..spaced.len {
        backward.push(spaced[spaced.len - 1 - k]);
    }
    let sb, c58a = statistics::sum(part[..]);
    let sp, c58 = statistics::sum(spaced[..]);
    let sr, c58b = statistics::sum(backward[..]);
    check(bad, c58a == 0 && c58 == 0 && c58b == 0 && sp == 4.0 * sb && sr == sp, "a sample of 160 spaced powers of two still sums");
    // Signed zeros: the mean of zeros of both signs is positive zero.
    let zz = list([-0.0, 0.0]);
    let zm, c59 = statistics::mean(zz[..]);
    check(bad, c59 == 0 && zm == 0.0 && 1.0 / zm > 0.0, "the mean of -0 and 0 is +0");
    bad += run_accept_assoc();
    return bad;
}

// One test for each acceptance scenario of the association calls.
fn run_accept_assoc() -> i64 {
    var bad: i64 = 0;

    // The covariance and the correlation of one pair.
    let x = list([1.0, 2.0, 3.0, 4.0, 5.0]);
    let y = list([2.0, 4.0, 5.0, 4.0, 5.0]);
    let cov, c1 = statistics::covariance(x[..], y[..]);
    let pcov, c2 = statistics::pcovariance(x[..], y[..]);
    check(bad, within(cov, 1.5, 4.0) && within(pcov, 1.2, 4.0) && c1 == 0 && c2 == 0, "the covariance of 1 2 3 4 5 and 2 4 5 4 5 is 1.5 and 1.2");
    let r, c3 = statistics::correlation(x[..], y[..]);
    check(bad, within(r, 0.7745966692414834, 4.0) && c3 == 0, "the correlation of the pair is 0.7745966692414834");
    let off = list([1e9 + 1.0, 1e9 + 2.0, 1e9 + 3.0]);
    let lin = list([1.0, 2.0, 3.0]);
    let ro, c4 = statistics::correlation(off[..], lin[..]);
    check(bad, within(ro, 1.0, 4.0) && c4 == 0, "the correlation of a large offset and a line is 1");
    // A perfectly correlated pair: the rounding of the sums pushes the
    // correlation past one, and the clamp holds it.
    let perfect_x = list([23.0, 5.0, 27.0]);
    let perfect_y = list([69.0, 15.0, 81.0]);
    let rp, c4b = statistics::correlation(perfect_x[..], perfect_y[..]);
    check(bad, rp == 1.0 && c4b == 0, "the correlation of a perfectly correlated pair is exactly 1");
    // A sample with a large magnitude: the products of the values do not
    // fit in a float, but the correlation does.
    let big_x = list([1e200, 2e200, 3e200, 4e200]);
    let big_y = list([2e200, 4e200, 5e200, 9e200]);
    let rb, c4c = statistics::correlation(big_x[..], big_y[..]);
    check(bad, within(rb, 0.9647638212377322, 4.0) && c4c == 0, "the correlation of a large sample keeps its value");
    let cb, c4d = statistics::covariance(big_x[..], big_y[..]);
    check(bad, is_nan(cb) && c4d == statistics::ERR_OVERFLOW, "the covariance of a large sample gives the overflow fault");
    // Two columns of a very different magnitude. The small column must not
    // fall to zero. The covariance is 1 and the correlation is 1.
    let mix_x = list([1e300, 2e300, 3e300]);
    let mix_y = list([1e-300, 2e-300, 3e-300]);
    let cm, c4e = statistics::covariance(mix_x[..], mix_y[..]);
    check(bad, within(cm, 1.0, 4.0) && c4e == 0, "the covariance of columns of a different magnitude keeps its value");
    let rm, c4f = statistics::correlation(mix_x[..], mix_y[..]);
    check(bad, rm == 1.0 && c4f == 0, "the correlation of columns of a different magnitude is exactly 1");
    let dec = list([1.0, 2.0, 3.0, 4.0]);
    let rev = list([4.0, 3.0, 2.0, 1.0]);
    let rn, c5 = statistics::correlation(dec[..], rev[..]);
    check(bad, within(rn, -1.0, 4.0) && c5 == 0, "the correlation of a line and its reverse is -1");
    let flat = list([2.0, 2.0, 2.0, 2.0]);
    let rc, c6 = statistics::correlation(dec[..], flat[..]);
    check(bad, is_nan(rc) && c6 == statistics::ERR_UNDEFINED, "the correlation of a constant column is undefined");
    let w1 = list([1.0, 2.0, 3.0]);
    let w2 = list([1.0, 1.0, 2.0]);
    let w3 = list([2.0, 3.0, 4.0]);
    // The weights 2, 3 and 4 repeat the points of the sample.
    let rep_x = list([1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0]);
    let rep_y = list([1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0]);
    let wc, c7 = statistics::weighted_covariance(w1[..], w2[..], w3[..], statistics::WEIGHT_FREQUENCY);
    let rc2, c8 = statistics::covariance(rep_x[..], rep_y[..]);
    check(bad, within(wc, rc2, 4.0) && c7 == 0 && c8 == 0, "the weighted covariance in the frequency form equals the repeated sample covariance");

    // The ranks and the rank correlations.
    let rs = list([3.0, 1.0, 3.0, 2.0]);
    var avg: f64[>..] = [];
    var mn: f64[>..] = [];
    var mx: f64[>..] = [];
    var dense: f64[>..] = [];
    var ord: f64[>..] = [];
    let n1, c9 = statistics::rankdata(avg, rs[..], statistics::TIE_AVERAGE);
    let n2, c10 = statistics::rankdata(mn, rs[..], statistics::TIE_MIN);
    let n3, c11 = statistics::rankdata(mx, rs[..], statistics::TIE_MAX);
    let n4, c12 = statistics::rankdata(dense, rs[..], statistics::TIE_DENSE);
    let n5, c13 = statistics::rankdata(ord, rs[..], statistics::TIE_ORDINAL);
    check(bad, c9 == 0 && c10 == 0 && c11 == 0 && c12 == 0 && c13 == 0
        && avg[0] == 3.5 && avg[1] == 1.0 && avg[2] == 3.5 && avg[3] == 2.0
        && mn[0] == 3.0 && mn[1] == 1.0 && mn[2] == 3.0 && mn[3] == 2.0
        && mx[0] == 4.0 && mx[1] == 1.0 && mx[2] == 4.0 && mx[3] == 2.0
        && dense[0] == 3.0 && dense[1] == 1.0 && dense[2] == 3.0 && dense[3] == 2.0
        && ord[0] == 3.0 && ord[1] == 1.0 && ord[2] == 4.0 && ord[3] == 2.0, "the five rank methods of 3 1 3 2");
    let rx = list([1.0, 2.0, 2.0, 3.0]);
    let ry = list([1.0, 3.0, 2.0, 2.0]);
    let rho, c14 = statistics::spearman(rx[..], ry[..]);
    check(bad, within(rho, 0.5, 4.0) && c14 == 0, "the Spearman coefficient of the pair is 0.5");
    let tb, c15 = statistics::kendall(rx[..], ry[..], statistics::KENDALL_B);
    let tc, c16 = statistics::kendall(rx[..], ry[..], statistics::KENDALL_C);
    check(bad, within(tb, 0.4, 4.0) && within(tc, 0.375, 4.0) && c15 == 0 && c16 == 0, "the Kendall tau-b is 0.4 and tau-c is 0.375");
    let rc3, c17 = statistics::spearman(dec[..], flat[..]);
    let tc2, c18 = statistics::kendall(dec[..], flat[..], statistics::KENDALL_B);
    check(bad, is_nan(rc3) && c17 == statistics::ERR_UNDEFINED && is_nan(tc2) && c18 == statistics::ERR_UNDEFINED, "the Spearman and the Kendall coefficient of a constant column are undefined");

    // The line and the Theil-Sen estimate.
    let l, c19 = statistics::linear_regression(x[..], y[..], false);
    check(bad, within(l.slope, 0.6, 4.0) && within(l.intercept, 2.2, 4.0) && c19 == 0, "the line of the pair has the slope 0.6 and the intercept 2.2");
    let se, ie, c20 = statistics::regression_errors(x[..], y[..]);
    check(bad, within(se, 0.28284271247461906, 4.0) && within(ie, 0.9380831519646861, 4.0) && c20 == 0, "the standard errors of the line are 0.28284271247461906 and 0.9380831519646861");
    let ox = list([1.0, 2.0, 3.0]);
    let oy = list([2.0, 4.0, 7.0]);
    let l2, c21 = statistics::linear_regression(ox[..], oy[..], true);
    check(bad, within(l2.slope, 2.2142857142857144, 4.0) && l2.intercept == 0.0 && c21 == 0, "the line through the origin has the slope 2.2142857142857144");
    let ts, ti, c22 = statistics::theil_sen(x[..], y[..]);
    check(bad, within(ts, 0.7083333333333333, 4.0) && within(ti, 1.875, 4.0) && c22 == 0, "the Theil-Sen estimate of the pair is 0.7083333333333333 and 1.875");
    let lc, c23 = statistics::linear_regression(dec[..], flat[..], false);
    check(bad, is_nan(lc.slope) && c23 == statistics::ERR_UNDEFINED, "the line of a constant column is undefined");
    var big: f64[>..] = [];
    for k in 0..statistics::MAX_THEIL + 1 {
        let v: f64 = k;
        big.push(v);
    }
    let ts2, ti2, c24 = statistics::theil_sen(big[..], big[..]);
    check(bad, is_nan(ts2) && is_nan(ti2) && c24 == statistics::ERR_TOO_LONG, "the Theil-Sen call above its bound gives the too long fault");

    // The matrices.
    var data: f64[>..] = [1.0, 2.0, 3.0, 4.0, 2.0, 4.0, 5.0, 9.0, 4.0, 3.0, 2.0, 1.0];
    var cells: f64[>..] = [];
    let nc, c25 = statistics::covariance_matrix(cells, data[..], 3, 1);
    check(bad, nc == 9 && c25 == 0 && within(cells[1], 3.6666666666666665, 4.0) && within(cells[0], 1.6666666666666667, 4.0), "the cell (0,1) of the covariance matrix is 3.6666666666666665");
    var ccells: f64[>..] = [];
    let nc2, c26 = statistics::correlation_matrix(ccells, data[..], 3);
    check(bad, nc2 == 9 && c26 == 0 && ccells[0] == 1.0 && ccells[4] == 1.0 && ccells[8] == 1.0 && within(ccells[2], -1.0, 4.0), "the diagonal of the correlation matrix is 1 and the cell (0,2) is -1");
    var one: f64[>..] = [1.0, 2.0, 3.0, 4.0];
    var onecells: f64[>..] = [];
    let nc3, c27 = statistics::covariance_matrix(onecells, one[..], 1, 1);
    let v, c28 = statistics::variance(one[..]);
    check(bad, nc3 == 1 && c27 == 0 && within(onecells[0], v, 4.0) && c28 == 0, "the covariance matrix of one column is the variance of the column");
    var cdata: f64[>..] = [1.0, 2.0, 3.0, 4.0, 5.0, 5.0, 5.0, 5.0];
    var cccells: f64[>..] = [];
    let nc4, c29 = statistics::correlation_matrix(ccells, cdata[..], 2);
    check(bad, nc4 == 0 && c29 == statistics::ERR_UNDEFINED && cccells.len == 0, "the correlation matrix of a constant column gives no partial matrix");
    var onecol: f64[>..] = [1.0, 2.0];
    var rowcells: f64[>..] = [];
    let nc5, c31 = statistics::covariance_matrix(rowcells, onecol[..], 2, 1);
    check(bad, nc5 == 0 && c31 == statistics::ERR_TOO_SHORT && rowcells.len == 0, "a matrix of one row gives the too short fault");
    let nc6, c32 = statistics::covariance_matrix(cells, data[..], statistics::MAX_COLUMNS + 1, 1);
    check(bad, nc6 == 0 && c32 == statistics::ERR_TOO_LONG && cells.len == 0, "a matrix call above the column bound gives the too long fault");

    // The faults of the pair calls.
    let short = list([1.0]);
    let cs, c33 = statistics::covariance(short[..], short[..]);
    check(bad, is_nan(cs) && c33 == statistics::ERR_TOO_SHORT, "the sample covariance of one point gives the too short fault");
    let other = list([1.0, 2.0]);
    let cl, c34 = statistics::covariance(short[..], other[..]);
    check(bad, is_nan(cl) && c34 == statistics::ERR_LENGTH, "a pair of columns of another length gives the length fault");
    let nan_col = list([1.0, 0.0 / 0.0, 3.0, 4.0]);
    let cn, c35 = statistics::correlation(nan_col[..], dec[..]);
    check(bad, is_nan(cn) && c35 == statistics::ERR_NAN, "a NaN in a column gives the NaN fault");
    let inf_col = list([1.0, 1.0 / 0.0, 3.0, 4.0]);
    let ci, c36 = statistics::correlation(inf_col[..], dec[..]);
    check(bad, is_nan(ci) && c36 == statistics::ERR_INFINITY, "an infinity in a column gives the infinity fault");
    let huge = list([1e308, -1e308, 1.0]);
    let ch, c37 = statistics::covariance(huge[..], huge[..]);
    check(bad, is_nan(ch) && c37 == statistics::ERR_OVERFLOW, "the covariance of 1e308 -1e308 1 gives the overflow fault");
    return bad;
}
