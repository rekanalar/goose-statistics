// One test for each documented case of the normality tests. The values
// come from the measurements against SciPy and statrs. The random
// vectors live in normality_vectors.goose, which the generator
// test/make_normality_vectors.py writes.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_normality_tests() -> i64 {
    var bad: i64 = 0;
    let nan_v = statistics::nan();

    // The acceptance sample.
    let acc: f64[16] = [1.1, 2.3, 3.7, 4.2, 5.9, 6.1, 7.4, 8.8,
                        2.0, 2.9, 3.1, 5.5, 6.6, 7.2, 9.1, 9.9];

    // The component tests carry the three alternatives.
    let s1, sc1 = statistics::skewtest(acc[..], statistics::ALT_TWO_SIDED);
    check(bad, sc1 == 0 && near_tol(s1.statistic, 0.18964323837234598, 1e-14)
        && near_tol(s1.pvalue, 0.8495887026717367, 1e-14), "skewtest, two-sided");
    let s2, sc2 = statistics::skewtest(acc[..], statistics::ALT_LESS);
    check(bad, sc2 == 0 && near_tol(s2.pvalue, 0.5752056486641317, 1e-14), "skewtest, less");
    let s3, sc3 = statistics::skewtest(acc[..], statistics::ALT_GREATER);
    check(bad, sc3 == 0 && near_tol(s3.pvalue, 0.42479435133586835, 1e-14), "skewtest, greater");

    let k1, kc1 = statistics::kurtosistest(acc[..], statistics::ALT_TWO_SIDED);
    check(bad, kc1 == 0 && near_tol(k1.statistic, -1.2823982121365567, 1e-14)
        && near_tol(k1.pvalue, 0.19970299106588685, 1e-14), "kurtosistest, two-sided");
    let k2, kc2 = statistics::kurtosistest(acc[..], statistics::ALT_LESS);
    check(bad, kc2 == 0 && near_tol(k2.pvalue, 0.09985149553294342, 1e-14), "kurtosistest, less");
    let k3, kc3 = statistics::kurtosistest(acc[..], statistics::ALT_GREATER);
    check(bad, kc3 == 0 && near_tol(k3.pvalue, 0.9001485044670565, 1e-14), "kurtosistest, greater");

    // The tails of one z add to one plus half of the two-sided p-value.
    check(bad, near(s2.pvalue + s3.pvalue, 1.0), "skewtest tails add to one");
    check(bad, near(k2.pvalue + k3.pvalue, 1.0), "kurtosistest tails add to one");
    check(bad, near(k1.pvalue, 2.0 * k2.pvalue), "kurtosistest two-sided is twice the smaller tail");

    // The omnibus tests.
    let n1, nc1 = statistics::normaltest(acc[..]);
    check(bad, nc1 == 0 && near_tol(n1.statistic, 1.6805097323513876, 1e-14)
        && near_tol(n1.pvalue, 0.43160050903903557, 1e-14), "normaltest, the acceptance sample");
    check(bad, n1.df == 2.0, "normaltest holds 2 degrees of freedom");
    let j1, jc1 = statistics::jarque_bera(acc[..]);
    check(bad, jc1 == 0 && near_tol(j1.statistic, 0.90763058673989838, 1e-14)
        && near_tol(j1.pvalue, 0.6352000480714689, 1e-14), "jarque-bera, the acceptance sample");

    // The Anderson-Darling test against the fitted normal.
    let a1, ac1 = statistics::anderson_darling(acc[..]);
    check(bad, ac1 == 0 && near_tol(a1.statistic, 0.23247860960977107, 1e-14)
        && near_tol(a1.pvalue, 0.7595249955902043, 1e-14), "anderson-darling, the acceptance sample");

    // A sample of 1490 values with one outlier at 1. The fitted normal
    // puts the outlier 38.6 standard units out, where the package tail
    // rounds to 0. The statistic uses the asymptotic tail form there.
    var big: f64[>..] = [];
    big.resize(1490, 0.0);
    big[1489] = 1.0;
    let bx, bcx = statistics::anderson_darling(big[..]);
    check(bad, bcx == 0 && near_tol(bx.statistic, 575.2895733499386, 1e-8)
        && bx.pvalue == 0.0, "anderson-darling, one outlier past the tail cutoff");
    big[0] = -1.0;
    let lo, lco = statistics::anderson_darling(big[..]);
    check(bad, lco == 0 && near_tol(lo.statistic, 573.3113471097017, 1e-8)
        && lo.pvalue == 0.0, "anderson-darling, one low outlier past the cutoff");

    // A statistic of 309 sits above the turn of the p-value quadratic at
    // 153.4. The rule would rise above 1 there, so the call gives 0.
    var mid: f64[>..] = [];
    mid.resize(800, 0.0);
    mid[799] = 1.0;
    let mc, mcc = statistics::anderson_darling(mid[..]);
    check(bad, mcc == 0 && near_tol(mc.statistic, 308.73273228529456, 1e-8)
        && mc.pvalue == 0.0, "anderson-darling, the p-value clamps at the quadratic turn");

    // The smallest samples. The kurtosis test runs from 5 values, the
    // skewness test from 8.
    let five: f64[5] = [1.1, 2.3, 3.7, 4.2, 5.9];
    let f5, fc5 = statistics::kurtosistest(five[..]);
    check(bad, fc5 == 0 && near_tol(f5.statistic, -0.1743311402893598, 1e-14), "kurtosistest, five values");
    let seven: f64[7] = [1.1, 2.3, 3.7, 4.2, 5.9, 6.1, 7.4];
    let s7, sc7 = statistics::skewtest(seven[..]);
    check(bad, sc7 == 2 && is_nan(s7.statistic), "skewtest, seven values is too short");
    let four: f64[4] = [1.1, 2.3, 3.7, 4.2];
    let f4, fc4 = statistics::kurtosistest(four[..]);
    check(bad, fc4 == 2, "kurtosistest, four values is too short");
    let n7, nc7 = statistics::normaltest(seven[..]);
    check(bad, nc7 == 2, "normaltest, seven values is too short");
    let two: f64[2] = [1.5, 2.5];
    let j2, jc2 = statistics::jarque_bera(two[..]);
    check(bad, jc2 == 0 && near(j2.statistic, 0.3333333333333333), "jarque-bera, two values");
    let one: f64[1] = [4.0];
    let j1s, jc1s = statistics::jarque_bera(one[..]);
    check(bad, jc1s == 2 && is_nan(j1s.statistic), "jarque-bera, one value is too short");
    let a4, ac4 = statistics::anderson_darling(four[..]);
    check(bad, ac4 == 2, "anderson-darling, four values is too short");

    // The hostile inputs.
    let none: f64[0] = [];
    let e1, ec1 = statistics::skewtest(none[..]);
    check(bad, ec1 == 1, "skewtest, an empty sample");
    let e2, ec2 = statistics::normaltest(none[..]);
    check(bad, ec2 == 1, "normaltest, an empty sample");
    let e3, ec3 = statistics::jarque_bera(none[..]);
    check(bad, ec3 == 1, "jarque-bera, an empty sample");
    let e4, ec4 = statistics::anderson_darling(none[..]);
    check(bad, ec4 == 1, "anderson-darling, an empty sample");

    var na: f64[>..] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
    na.push(nan_v);
    let nn, nnc = statistics::normaltest(na[..]);
    check(bad, nnc == 3, "normaltest, a NaN in the sample");
    var ia: f64[>..] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
    ia.push(1.0 / 0.0);
    let ii, iic = statistics::jarque_bera(ia[..]);
    check(bad, iic == 4, "jarque-bera, an infinity in the sample");
    let iad, iadc = statistics::anderson_darling(ia[..]);
    check(bad, iadc == 4, "anderson-darling, an infinity in the sample");

    let con: f64[8] = [2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0];
    let q1, qc1 = statistics::skewtest(con[..]);
    check(bad, qc1 == 7, "skewtest, a constant sample");
    let q2, qc2 = statistics::kurtosistest(con[..]);
    check(bad, qc2 == 7, "kurtosistest, a constant sample");
    let q3, qc3 = statistics::normaltest(con[..]);
    check(bad, qc3 == 7, "normaltest, a constant sample");
    let q4, qc4 = statistics::jarque_bera(con[..]);
    check(bad, qc4 == 7, "jarque-bera, a constant sample");
    let q5, qc5 = statistics::anderson_darling(con[..]);
    check(bad, qc5 == 7, "anderson-darling, a constant sample has no scale");

    let d3, dc3 = statistics::skewtest(acc[..], 3);
    check(bad, dc3 == 5, "skewtest, an unknown alternative");

    // The round trips. The K-squared statistic is the sum of the squares
    // of the two z values, and both omnibus p-values are the upper
    // chi-square tail with 2 degrees of freedom, which is exp(-x/2).
    let k2r = s1.statistic * s1.statistic + k1.statistic * k1.statistic;
    check(bad, near_tol(n1.statistic, k2r, 1e-14), "the K-squared round trip");
    check(bad, near_tol(n1.pvalue, exp(-0.5 * n1.statistic), 1e-14), "the K-squared tail");
    check(bad, near_tol(j1.pvalue, exp(-0.5 * j1.statistic), 1e-14), "the Jarque-Bera tail");

    // No call writes into the sample of the caller.
    var keep: f64[>..] = [3.1, 1.4, 1.5, 9.2, 5.2, 6.1, 2.7, 4.4, 8.8, 2.2];
    let before: f64[10] = [3.1, 1.4, 1.5, 9.2, 5.2, 6.1, 2.7, 4.4, 8.8, 2.2];
    let _u1, _u1c = statistics::skewtest(keep[..]);
    let _u2, _u2c = statistics::kurtosistest(keep[..]);
    let _u3, _u3c = statistics::normaltest(keep[..]);
    let _u4, _u4c = statistics::jarque_bera(keep[..]);
    let _u5, _u5c = statistics::anderson_darling(keep[..]);
    var same = true;
    for i in 0..10 {
        if keep[i] != before[i] { same = false; }
    }
    check(bad, same, "the calls leave the sample unchanged");

    // A sample above the input bound gives the too long fault in every
    // call and no work starts.
    var long: f64[>..] = [];
    for i in 0..statistics::MAX_INPUT + 1 { long.push(1.0); }
    let g1, gc1 = statistics::skewtest(long[..]);
    let g2, gc2 = statistics::kurtosistest(long[..]);
    let g3, gc3 = statistics::normaltest(long[..]);
    let g4, gc4 = statistics::jarque_bera(long[..]);
    let g5, gc5 = statistics::anderson_darling(long[..]);
    check(bad, gc1 == 9 && gc2 == 9 && gc3 == 9 && gc4 == 9 && gc5 == 9,
        "the calls refuse a sample above the input bound");

    return bad;
}
