// One test for each documented case of the distributions. Each case
// states the arguments and the value that the reference libraries give.
namespace suite;

import std;
import statistics from "../..";

fn run_accept_dist() -> i64 {
    var bad: i64 = 0;

    // The densities and the cumulative functions of the continuous
    // families.
    let a1, c1 = statistics::norm_pdf(0.5);
    check(bad, within(a1, 0.3520653267642995, 4.0) && c1 == 0, "the normal density at 0.5");
    let a2, c2 = statistics::norm_cdf(1.5);
    check(bad, within(a2, 0.9331927987311419, 4.0) && c2 == 0, "the normal cumulative function at 1.5");
    let a3, c3 = statistics::norm_cdf(-0.75);
    check(bad, within(a3, 0.2266273523768682, 4.0) && c3 == 0, "the normal cumulative function at -0.75");
    let a4, c4 = statistics::norm_sf(6.0);
    check(bad, within(a4, 9.86587645037698e-10, 4.0) && c4 == 0, "the normal survival function at 6");
    let a5, c5 = statistics::unif_pdf(0.5);
    check(bad, within(a5, 1.0, 4.0) && c5 == 0, "the uniform density at 0.5");
    let a6, c6 = statistics::expon_pdf(2.0, 3.0);
    check(bad, within(a6, 0.0074362565299990755, 8.0) && c6 == 0, "the exponential density at 2 with a rate of 3");
    let a7, c7 = statistics::gammadist_cdf(3.0, 2.0);
    check(bad, within(a7, 0.8008517265285442, 16.0) && c7 == 0, "the gamma cumulative function at 3 with a shape of 2");
    let a8, c8 = statistics::betadist_pdf(0.4, 2.0, 3.0);
    check(bad, within(a8, 1.7279999999999998, 8.0) && c8 == 0, "the beta density at 0.4");
    let a9, c9 = statistics::chi2_cdf(3.0, 4.0);
    check(bad, within(a9, 0.4421745996289252, 16.0) && c9 == 0, "the chi-square cumulative function at 3");
    let a10, c10 = statistics::tdist_cdf(1.5, 7.0);
    check(bad, within(a10, 0.911350756505015, 16.0) && c10 == 0, "the Student t cumulative function at 1.5");
    let a11, c11 = statistics::fdist_cdf(2.0, 3.0, 6.0);
    check(bad, within(a11, 0.7844465853788261, 16.0) && c11 == 0, "the F cumulative function at 2");

    // The quantiles of the continuous families.
    let a12, c12 = statistics::norm_ppf(0.975);
    check(bad, within(a12, 1.959963984540054, 4.0) && c12 == 0, "the normal quantile at 0.975");
    let a13, c13 = statistics::norm_ppf(0.001);
    check(bad, within(a13, -3.090232306167813, 8.0) && c13 == 0, "the normal quantile at 0.001");
    let a14, c14 = statistics::tdist_ppf(0.975, 7.0);
    check(bad, within(a14, 2.364624251592784, 16.0) && c14 == 0, "the Student t quantile at 0.975");
    let a15, c15 = statistics::fdist_ppf(0.95, 3.0, 6.0);
    check(bad, within(a15, 4.757062663089412, 16.0) && c15 == 0, "the F quantile at 0.95");

    // The masses and the cumulative functions of the discrete families.
    let a16, c16 = statistics::binom_pmf(3.0, 10.0, 0.5);
    check(bad, within(a16, 0.1171875, 4.0) && c16 == 0, "the binomial mass at 3");
    let a17, c17 = statistics::binom_cdf(3.0, 10.0, 0.5);
    check(bad, within(a17, 0.171875, 4.0) && c17 == 0, "the binomial cumulative function at 3");
    let a18, c18 = statistics::poisson_pmf(2.0, 3.0);
    check(bad, within(a18, 0.22404180765538775, 4.0) && c18 == 0, "the Poisson mass at 2");
    let a19, c19 = statistics::poisson_cdf(2.0, 3.0);
    check(bad, within(a19, 0.42319008112684364, 8.0) && c19 == 0, "the Poisson cumulative function at 2");
    let a20, c20 = statistics::geom_pmf(2.0, 0.5);
    check(bad, within(a20, 0.25, 4.0) && c20 == 0, "the geometric mass at 2");
    let a21, c21 = statistics::nbinom_pmf(2.0, 3.0, 0.5);
    check(bad, within(a21, 0.1875, 4.0) && c21 == 0, "the negative binomial mass at 2");
    let a22, c22 = statistics::hypergeom_pmf(1.0, 20.0, 10.0, 4.0);
    check(bad, within(a22, 0.24767801857585137, 8.0) && c22 == 0, "the hypergeometric mass at 1");

    // The quantiles of the discrete families.
    let a23, c23 = statistics::binom_ppf(0.5, 10.0, 0.5);
    check(bad, a23 == 5.0 && c23 == 0, "the binomial quantile at 0.5");
    let a24, c24 = statistics::poisson_ppf(0.5, 3.0);
    check(bad, a24 == 3.0 && c24 == 0, "the Poisson quantile at 0.5");
    let a25, c25 = statistics::hypergeom_ppf(0.5, 20.0, 10.0, 4.0);
    check(bad, a25 == 2.0 && c25 == 0, "the hypergeometric quantile at 0.5");
    let a26, c26 = statistics::binom_ppf(0.0, 10.0, 0.5);
    check(bad, a26 == 0.0 && c26 == 0, "the binomial quantile at 0");
    let a27, c27 = statistics::binom_ppf(1.0, 10.0, 0.5);
    check(bad, a27 == 10.0 && c27 == 0, "the binomial quantile at 1");

    // The moments.
    let a28, c28 = statistics::norm_mean();
    let a29, c29 = statistics::norm_var();
    check(bad, a28 == 0.0 && a29 == 1.0 && c28 == 0 && c29 == 0, "the mean and the variance of the normal family");
    let a30, c30 = statistics::binom_mean(10.0, 0.5);
    let a31, c31 = statistics::binom_var(10.0, 0.5);
    check(bad, a30 == 5.0 && a31 == 2.5 && c30 == 0 && c31 == 0, "the mean and the variance of the binomial family");
    let a32, c32 = statistics::poisson_mean(3.0);
    let a33, c33 = statistics::poisson_var(3.0);
    check(bad, a32 == 3.0 && a33 == 3.0 && c32 == 0 && c33 == 0, "the mean and the variance of the Poisson family");

    // The faults.
    let nan = 0.0 / 0.0;
    let a34, c34 = statistics::norm_pdf(nan);
    check(bad, is_nan(a34) && c34 == statistics::ERR_NAN, "the normal density of a NaN");
    let a35, c35 = statistics::norm_pdf(1.0, 0.0, 0.0);
    check(bad, is_nan(a35) && c35 == statistics::ERR_DOMAIN, "the normal density of a scale of 0");
    let a36, c36 = statistics::gammadist_pdf(1.0, 0.0);
    check(bad, is_nan(a36) && c36 == statistics::ERR_DOMAIN, "the gamma density of a shape of 0");
    let a37, c37 = statistics::norm_ppf(0.0);
    check(bad, is_nan(a37) && c37 == statistics::ERR_DOMAIN, "the normal quantile of 0");
    let a38, c38 = statistics::norm_ppf(1.5);
    check(bad, is_nan(a38) && c38 == statistics::ERR_DOMAIN, "the normal quantile of 1.5");
    let a39, c39 = statistics::binom_pmf(-1.0, 10.0, 0.5);
    check(bad, is_nan(a39) && c39 == statistics::ERR_DOMAIN, "the binomial mass of a count of -1");

    return bad;
}
