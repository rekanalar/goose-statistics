// The families of the truncation, the gamma, the generalized, the
// symmetric and the Johnson groups: the density, the cumulative function,
// the quantile and the moments of a few of them.
//
// Run it with: goose run examples/families.goose
import std;
import statistics from "..";

fn main() {
    // A family of the elementary group. Its density is the derivative of
    // its cumulative function.
    let d, c1 = statistics::arcsine_pdf(0.5);
    let f, c2 = statistics::arcsine_cdf(0.5);
    print("arcsine pdf(0.5)  ", d, "\n");
    print("arcsine cdf(0.5)  ", f, "\n");

    // A family of the bounded group. The quantile of a probability gives a
    // point whose cumulative function is that probability.
    let q, c3 = statistics::triang_ppf(0.75, 0.4);
    let f2, c4 = statistics::triang_cdf(q, 0.4);
    print("triang ppf(0.75)  ", q, "\n");
    print("triang cdf of it  ", f2, "\n");

    // A family of the Weibull group. Its mean is 0, because its law is
    // symmetric around 0.
    let m, c5 = statistics::dweibull_mean(1.7);
    let v, c6 = statistics::dweibull_var(1.7);
    print("dweibull mean     ", m, "\n");
    print("dweibull var      ", v, "\n");

    // A family of the power group. Its moments come from a rule of
    // Gauss-Legendre on the density of the law.
    let pm, c7 = statistics::powerlognorm_mean(1.4, 2.0);
    let pk, c8 = statistics::powerlognorm_kurt(1.4, 2.0);
    print("powerlognorm mean ", pm, "\n");
    print("powerlognorm kurt ", pk, "\n");

    // A family of the power group. Its moments have a closed form.
    let fm, c9 = statistics::fisk_mean(5.0);
    let fk, c10 = statistics::genpareto_kurt(-0.5);
    print("fisk mean         ", fm, "\n");
    print("genpareto kurt    ", fk, "\n");

    // A family of the growth group. Its moments come from the Hurwitz zeta
    // function or from a rule of Gauss-Legendre.
    let gm, c11 = statistics::genlogistic_mean(2.5);
    let gk, c12 = statistics::gompertz_kurt(1.5);
    print("genlogistic mean  ", gm, "\n");
    print("gompertz kurt     ", gk, "\n");

    // A family of the truncation group. Its support is bounded.
    let tm, c13 = statistics::truncpareto_mean(2.0, 5.0);
    let tq, c14 = statistics::truncnorm_ppf(0.75, -1.0, 2.0);
    print("truncpareto mean  ", tm, "\n");
    print("truncnorm ppf     ", tq, "\n");

    // A family that comes from the gamma law. Its moments have a closed
    // form in the gamma function.
    let cm, c15 = statistics::chi_mean(3.0);
    let cv, c16 = statistics::invgamma_var(3.0);
    print("chi mean          ", cm, "\n");
    print("invgamma var      ", cv, "\n");

    // A family of the generalized group. Its moments have a closed form.
    let gm2, c17 = statistics::gengamma_mean(3.0, 2.0);
    let gk2, c18 = statistics::mielke_mean(2.0, 3.0);
    print("gengamma mean     ", gm2, "\n");
    print("mielke mean       ", gk2, "\n");

    // A family of the symmetric group. Its law is symmetric around 0.
    let sm, c19 = statistics::gennorm_mean(1.5);
    let sk, c20 = statistics::halfgennorm_kurt(1.5);
    print("gennorm mean      ", sm, "\n");
    print("halfgennorm kurt  ", sk, "\n");

    // A family of the Johnson group. Its quantile is a logistic function.
    let jq, c21 = statistics::johnsonsb_ppf(0.75, 1.0, 2.0);
    let ju, c22 = statistics::johnsonsu_mean(1.0, 2.0);
    print("johnsonsb ppf     ", jq, "\n");
    print("johnsonsu mean    ", ju, "\n");

    // A family of the special group. Its quantile searches the root of
    // its survival function.
    let sp, c23 = statistics::invgauss_ppf(0.3, 0.5);
    let sc, c24 = statistics::invgauss_cdf(sp, 0.5);
    print("invgauss ppf     ", sp, "\n");
    print("invgauss cdf     ", sc, "\n");

    // A family of the special group. Its moments have a closed form.
    let km, c25 = statistics::nakagami_mean(2.5);
    let kk, c26 = statistics::kappa3_mean(3.5);
    print("nakagami mean    ", km, "\n");
    print("kappa3 mean      ", kk, "\n");

    // A family of the skew group. Its cumulative function uses the T
    // function of Owen.
    let sn, c27 = statistics::skewnorm_cdf(0.5, 2.0);
    let la, c28 = statistics::laplace_asymmetric_mean(0.5);
    print("skewnorm cdf     ", sn, "\n");
    print("laplace asym mean", la, "\n");

    // A fault gives NaN and a code. A shape at or below 0 is out of range.
    let bad, c29 = statistics::exponpow_pdf(1.0, 0.0);
    print("exponpow shape 0  ", bad, " code ", c29, "\n");
    print("codes             ", c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, "\n");
    print("codes             ", c11, c12, c13, c14, c15, c16, c17, c18, c19, c20, "\n");
    print("codes             ", c21, c22, c23, c24, c25, c26, c27, c28, c29, "\n");
    return 0;
}
