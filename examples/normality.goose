// The normality tests: the skewness and kurtosis tests, the omnibus
// tests and the Anderson-Darling test, each on one sample.
//
// Run it with: goose run examples/normality.goose
import std;
import statistics from "..";

fn main() {
    let sample: f64[12] = [5.1, 4.9, 5.3, 4.8, 5.0, 5.4, 4.7, 5.2,
                           4.6, 5.5, 4.4, 5.6];

    // The test of the skewness. The alternatives cover one tail or both.
    let s, c1 = statistics::skewtest(sample[..], statistics::ALT_TWO_SIDED);
    print("skewtest z = ", s.statistic, ", p = ", s.pvalue, ", code ", c1, "\n");

    // The test of the kurtosis.
    let k, c2 = statistics::kurtosistest(sample[..], statistics::ALT_TWO_SIDED);
    print("kurtosistest z = ", k.statistic, ", p = ", k.pvalue, ", code ", c2, "\n");

    // The omnibus test of D'Agostino and Pearson.
    let n, c3 = statistics::normaltest(sample[..]);
    print("normaltest K-squared = ", n.statistic, ", p = ", n.pvalue, ", code ", c3, "\n");

    // The test of Jarque and Bera.
    let j, c4 = statistics::jarque_bera(sample[..]);
    print("jarque-bera = ", j.statistic, ", p = ", j.pvalue, ", code ", c4, "\n");

    // The Anderson-Darling test against the fitted normal.
    let a, c5 = statistics::anderson_darling(sample[..]);
    print("anderson-darling A-squared = ", a.statistic, ", p = ", a.pvalue, ", code ", c5, "\n");

    // A constant sample has no shape and every call faults.
    let flat: f64[8] = [2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0];
    let q, c6 = statistics::normaltest(flat[..]);
    print("constant sample gives code ", c6, "\n");
}
