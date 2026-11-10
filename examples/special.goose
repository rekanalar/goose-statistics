// The special functions: the error function, the gamma family and the
// incomplete gamma and beta functions with their inverses.
//
// Run it with: goose run examples/special.goose
import std;
import statistics from "..";

fn main() {
    // erfc keeps its digits in the far tail, where 1 - erf gives zero.
    let tail, c1 = statistics::erfc(10.0);
    print("erfc(10)          ", tail, "\n");

    // gamma is exact on small integers. lgamma gives the log and the sign.
    let g5, c2 = statistics::gamma(5.0);
    let lg, sign, c3 = statistics::lgamma(-2.5);
    print("gamma(5)          ", g5, "\n");
    print("lgamma(-2.5)      ", lg, " sign ", sign, "\n");

    // The chi-square distribution with one degree of freedom is the
    // incomplete gamma with shape 1/2 and argument x/2.
    let q, c4 = statistics::gammainccinv(0.5, 0.05);
    print("chi-square 95%    ", 2.0 * q, "\n");
    let p, c5 = statistics::gammaincc(0.5, 3.841458820694124 / 2.0);
    print("p value           ", p, "\n");

    // The beta distribution: the probability and its inverse.
    let cdf, c6 = statistics::betainc(2.0, 5.0, 0.3);
    let x, c7 = statistics::betaincinv(2.0, 5.0, cdf);
    print("betainc(2,5,0.3)  ", cdf, " and back ", x, "\n");

    // A fault gives NaN and a code. A pole of the gamma is a fault.
    let pole, c8 = statistics::gamma(-3.0);
    print("gamma(-3)         ", pole, " code ", c8, "\n");
    print("codes             ", c1, c2, c3, c4, c5, c6, c7, "\n");
    return 0;
}
