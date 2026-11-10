// Tests of the distributions against the reference vectors. Each row holds
// the family number, the operation number, the argument, up to three
// parameters, the exact value and the tolerance in units in the last place.
namespace suite;

import std;
import statistics from "../..";
import distribution_vectors;

// The result of the call of one family and one operation.
fn dist_call(fam: i64, op: i64, x: f64, q: f64, a: f64, b: f64, c: f64) -> f64, i64 {
    if fam == 0 {
        if op == 0 { return statistics::norm_pdf(x, a, b); }
        if op == 1 { return statistics::norm_cdf(x, a, b); }
        if op == 2 { return statistics::norm_sf(x, a, b); }
        if op == 3 { return statistics::norm_ppf(q, a, b); }
        if op == 4 { return statistics::norm_mean(a, b); }
        if op == 5 { return statistics::norm_var(a, b); }
        if op == 6 { return statistics::norm_skew(a, b); }
        if op == 7 { return statistics::norm_kurt(a, b); }
    }
    if fam == 1 {
        if op == 0 { return statistics::unif_pdf(x, a, b); }
        if op == 1 { return statistics::unif_cdf(x, a, b); }
        if op == 2 { return statistics::unif_sf(x, a, b); }
        if op == 3 { return statistics::unif_ppf(q, a, b); }
        if op == 4 { return statistics::unif_mean(a, b); }
        if op == 5 { return statistics::unif_var(a, b); }
        if op == 6 { return statistics::unif_skew(a, b); }
        if op == 7 { return statistics::unif_kurt(a, b); }
    }
    if fam == 2 {
        if op == 0 { return statistics::expon_pdf(x, a); }
        if op == 1 { return statistics::expon_cdf(x, a); }
        if op == 2 { return statistics::expon_sf(x, a); }
        if op == 3 { return statistics::expon_ppf(q, a); }
        if op == 4 { return statistics::expon_mean(a); }
        if op == 5 { return statistics::expon_var(a); }
        if op == 6 { return statistics::expon_skew(a); }
        if op == 7 { return statistics::expon_kurt(a); }
    }
    if fam == 3 {
        if op == 0 { return statistics::lognorm_pdf(x, a, b); }
        if op == 1 { return statistics::lognorm_cdf(x, a, b); }
        if op == 2 { return statistics::lognorm_sf(x, a, b); }
        if op == 3 { return statistics::lognorm_ppf(q, a, b); }
        if op == 4 { return statistics::lognorm_mean(a, b); }
        if op == 5 { return statistics::lognorm_var(a, b); }
        if op == 6 { return statistics::lognorm_skew(a, b); }
        if op == 7 { return statistics::lognorm_kurt(a, b); }
    }
    if fam == 4 {
        if op == 0 { return statistics::weibull_pdf(x, a, b); }
        if op == 1 { return statistics::weibull_cdf(x, a, b); }
        if op == 2 { return statistics::weibull_sf(x, a, b); }
        if op == 3 { return statistics::weibull_ppf(q, a, b); }
        if op == 4 { return statistics::weibull_mean(a, b); }
        if op == 5 { return statistics::weibull_var(a, b); }
        if op == 6 { return statistics::weibull_skew(a, b); }
        if op == 7 { return statistics::weibull_kurt(a, b); }
    }
    if fam == 5 {
        if op == 0 { return statistics::gammadist_pdf(x, a, b); }
        if op == 1 { return statistics::gammadist_cdf(x, a, b); }
        if op == 2 { return statistics::gammadist_sf(x, a, b); }
        if op == 3 { return statistics::gammadist_ppf(q, a, b); }
        if op == 4 { return statistics::gammadist_mean(a, b); }
        if op == 5 { return statistics::gammadist_var(a, b); }
        if op == 6 { return statistics::gammadist_skew(a, b); }
        if op == 7 { return statistics::gammadist_kurt(a, b); }
    }
    if fam == 6 {
        if op == 0 { return statistics::betadist_pdf(x, a, b); }
        if op == 1 { return statistics::betadist_cdf(x, a, b); }
        if op == 2 { return statistics::betadist_sf(x, a, b); }
        if op == 3 { return statistics::betadist_ppf(q, a, b); }
        if op == 4 { return statistics::betadist_mean(a, b); }
        if op == 5 { return statistics::betadist_var(a, b); }
        if op == 6 { return statistics::betadist_skew(a, b); }
        if op == 7 { return statistics::betadist_kurt(a, b); }
    }
    if fam == 7 {
        if op == 0 { return statistics::chi2_pdf(x, a); }
        if op == 1 { return statistics::chi2_cdf(x, a); }
        if op == 2 { return statistics::chi2_sf(x, a); }
        if op == 3 { return statistics::chi2_ppf(q, a); }
        if op == 4 { return statistics::chi2_mean(a); }
        if op == 5 { return statistics::chi2_var(a); }
        if op == 6 { return statistics::chi2_skew(a); }
        if op == 7 { return statistics::chi2_kurt(a); }
    }
    if fam == 8 {
        if op == 0 { return statistics::tdist_pdf(x, a); }
        if op == 1 { return statistics::tdist_cdf(x, a); }
        if op == 2 { return statistics::tdist_sf(x, a); }
        if op == 3 { return statistics::tdist_ppf(q, a); }
        if op == 4 { return statistics::tdist_mean(a); }
        if op == 5 { return statistics::tdist_var(a); }
        if op == 6 { return statistics::tdist_skew(a); }
        if op == 7 { return statistics::tdist_kurt(a); }
    }
    if fam == 9 {
        if op == 0 { return statistics::fdist_pdf(x, a, b); }
        if op == 1 { return statistics::fdist_cdf(x, a, b); }
        if op == 2 { return statistics::fdist_sf(x, a, b); }
        if op == 3 { return statistics::fdist_ppf(q, a, b); }
        if op == 4 { return statistics::fdist_mean(a, b); }
        if op == 5 { return statistics::fdist_var(a, b); }
        if op == 6 { return statistics::fdist_skew(a, b); }
        if op == 7 { return statistics::fdist_kurt(a, b); }
    }
    if fam == 10 {
        if op == 0 { return statistics::binom_pmf(x, a, b); }
        if op == 1 { return statistics::binom_cdf(x, a, b); }
        if op == 2 { return statistics::binom_sf(x, a, b); }
        if op == 3 { return statistics::binom_ppf(q, a, b); }
        if op == 4 { return statistics::binom_mean(a, b); }
        if op == 5 { return statistics::binom_var(a, b); }
        if op == 6 { return statistics::binom_skew(a, b); }
        if op == 7 { return statistics::binom_kurt(a, b); }
    }
    if fam == 11 {
        if op == 0 { return statistics::poisson_pmf(x, a); }
        if op == 1 { return statistics::poisson_cdf(x, a); }
        if op == 2 { return statistics::poisson_sf(x, a); }
        if op == 3 { return statistics::poisson_ppf(q, a); }
        if op == 4 { return statistics::poisson_mean(a); }
        if op == 5 { return statistics::poisson_var(a); }
        if op == 6 { return statistics::poisson_skew(a); }
        if op == 7 { return statistics::poisson_kurt(a); }
    }
    if fam == 12 {
        if op == 0 { return statistics::geom_pmf(x, a); }
        if op == 1 { return statistics::geom_cdf(x, a); }
        if op == 2 { return statistics::geom_sf(x, a); }
        if op == 3 { return statistics::geom_ppf(q, a); }
        if op == 4 { return statistics::geom_mean(a); }
        if op == 5 { return statistics::geom_var(a); }
        if op == 6 { return statistics::geom_skew(a); }
        if op == 7 { return statistics::geom_kurt(a); }
    }
    if fam == 13 {
        if op == 0 { return statistics::nbinom_pmf(x, a, b); }
        if op == 1 { return statistics::nbinom_cdf(x, a, b); }
        if op == 2 { return statistics::nbinom_sf(x, a, b); }
        if op == 3 { return statistics::nbinom_ppf(q, a, b); }
        if op == 4 { return statistics::nbinom_mean(a, b); }
        if op == 5 { return statistics::nbinom_var(a, b); }
        if op == 6 { return statistics::nbinom_skew(a, b); }
        if op == 7 { return statistics::nbinom_kurt(a, b); }
    }
    if fam == 14 {
        if op == 0 { return statistics::hypergeom_pmf(x, a, b, c); }
        if op == 1 { return statistics::hypergeom_cdf(x, a, b, c); }
        if op == 2 { return statistics::hypergeom_sf(x, a, b, c); }
        if op == 3 { return statistics::hypergeom_ppf(q, a, b, c); }
        if op == 4 { return statistics::hypergeom_mean(a, b, c); }
        if op == 5 { return statistics::hypergeom_var(a, b, c); }
        if op == 6 { return statistics::hypergeom_skew(a, b, c); }
        if op == 7 { return statistics::hypergeom_kurt(a, b, c); }
    }
    return 0.0 / 0.0, 99;
}

fn dist_family_name(fam: i64) -> u8[:] {
    if fam == 0 { return "norm"; }
    if fam == 1 { return "unif"; }
    if fam == 2 { return "expon"; }
    if fam == 3 { return "lognorm"; }
    if fam == 4 { return "weibull"; }
    if fam == 5 { return "gammadist"; }
    if fam == 6 { return "betadist"; }
    if fam == 7 { return "chi2"; }
    if fam == 8 { return "tdist"; }
    if fam == 9 { return "fdist"; }
    if fam == 10 { return "binom"; }
    if fam == 11 { return "poisson"; }
    if fam == 12 { return "geom"; }
    if fam == 13 { return "nbinom"; }
    if fam == 14 { return "hypergeom"; }
    return "?";
}

// The names of the operations.
fn dist_op_name(op: i64) -> u8[:] {
    if op == 0 { return "pdf"; }
    if op == 1 { return "cdf"; }
    if op == 2 { return "sf"; }
    if op == 3 { return "ppf"; }
    if op == 4 { return "mean"; }
    if op == 5 { return "var"; }
    if op == 6 { return "skew"; }
    return "kurt";
}

// Run all the rows. A row passes when the call gives no fault and the
// value lies within the tolerance of the exact value.
fn run_distributions(report: bool) -> i64 {
    var bad: i64 = 0;
    var worst: f64[120];
    for k in 0..120 { worst[k] = 0.0; }
    for r in 0..D_ROWS {
        let fam = D_FAM[r];
        let op = D_OP[r];
        let got, code = dist_call(fam, op, D_X[r], D_X[r], D_A[r], D_B[r], D_C[r]);
        let want = D_E[r];
        let err = ulp_error(got, want);
        if code == 0 && err > worst[fam * 8 + op] { worst[fam * 8 + op] = err; }
        let ok = code == 0 && err <= D_TOL[r];
        if !ok {
            bad += 1;
            if bad <= 30 {
                print("FAIL ", dist_family_name(fam), " ", dist_op_name(op), " row ", r, " x ", D_X[r], " a ", D_A[r], " b ", D_B[r], " c ", D_C[r], " got ", got, " want ", want, " code ", code, " ulp ", err, "\n");
            }
        }
    }
    if report {
        for k in 0..120 {
            print(dist_family_name(k / 8), " ", dist_op_name(k % 8), " worst ", worst[k], "\n");
        }
    }
    return bad;
}
