// The probability distributions: the density, the cumulative function, the
// quantile and the moments of a few families.
//
// Run it with: goose run examples/distributions.goose
import std;
import statistics from "..";

fn main() {
    // The normal family. The survival function keeps the digits of the
    // far tail, where one minus the cumulative function gives zero.
    let d, c1 = statistics::norm_pdf(1.5);
    let f, c2 = statistics::norm_cdf(1.5);
    let s, c3 = statistics::norm_sf(1.5);
    print("normal pdf(1.5)   ", d, "\n");
    print("normal cdf(1.5)   ", f, "\n");
    print("normal sf(1.5)    ", s, "\n");

    // The quantile of the normal is the inverse of the cumulative function.
    let q, c4 = statistics::norm_ppf(0.975);
    print("normal ppf(0.975) ", q, "\n");

    // The Student t family. The p value of a test with seven degrees of
    // freedom uses the survival function of the family.
    let t, c5 = statistics::tdist_sf(2.364624251592784, 7.0);
    print("t sf(2.3646, 7)   ", t, "\n");
    let tp, c6 = statistics::tdist_ppf(0.975, 7.0);
    print("t ppf(0.975, 7)   ", tp, "\n");

    // The F family. The quantile gives the critical value of a test.
    let f, c7 = statistics::fdist_ppf(0.95, 3.0, 6.0);
    print("F ppf(0.95, 3, 6) ", f, "\n");

    // A discrete family. The mass and the cumulative function of the
    // binomial family.
    let m, c8 = statistics::binom_pmf(3.0, 10.0, 0.5);
    let cd, c9 = statistics::binom_cdf(3.0, 10.0, 0.5);
    print("binom pmf(3,10)   ", m, "\n");
    print("binom cdf(3,10)   ", cd, "\n");

    // The quantile of a discrete family is a count.
    let bq, c10 = statistics::binom_ppf(0.5, 10.0, 0.5);
    print("binom ppf(0.5)    ", bq, "\n");

    // The moments of a family.
    let mean, c11 = statistics::gammadist_mean(3.0, 2.0);
    let variance, c12 = statistics::gammadist_var(3.0, 2.0);
    let skew, c13 = statistics::gammadist_skew(3.0, 2.0);
    let kurt, c14 = statistics::gammadist_kurt(3.0, 2.0);
    print("gamma mean        ", mean, "\n");
    print("gamma var         ", variance, "\n");
    print("gamma skew        ", skew, "\n");
    print("gamma kurt        ", kurt, "\n");

    // A fault gives NaN and a code. A shape at or below 0 is out of range.
    let bad, c15 = statistics::gammadist_pdf(1.0, 0.0, 1.0);
    print("gamma pdf shape 0 ", bad, " code ", c15, "\n");
    print("codes             ", c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, "\n");
    print("codes             ", c11, c12, c13, c14, "\n");
    return 0;
}
