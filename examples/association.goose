// The association of two samples: the covariance, the correlation, the
// rank correlations and the line fit.
//
// Run it with: goose run examples/association.goose
import std;
import statistics from "..";

fn main() {
    let xs: f64[>..] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
    let ys: f64[>..] = [3.0, 5.0, 4.0, 7.0, 6.0, 9.0, 8.0, 11.0];
    let cov, c1 = statistics::covariance(xs[..], ys[..]);
    let pcov, c2 = statistics::pcovariance(xs[..], ys[..]);
    let r, c3 = statistics::correlation(xs[..], ys[..]);
    let rho, c4 = statistics::spearman(xs[..], ys[..]);
    let tau, c5 = statistics::kendall(xs[..], ys[..], statistics::KENDALL_B);
    let l, c6 = statistics::linear_regression(xs[..], ys[..], false);
    print("covariance  ", cov, "\n");
    print("pcovariance ", pcov, "\n");
    print("correlation ", r, "\n");
    print("spearman    ", rho, "\n");
    print("kendall b   ", tau, "\n");
    print("slope       ", l.slope, "\n");
    print("intercept   ", l.intercept, "\n");
    print("r2          ", l.r2, "\n");
    print("codes       ", c1, " ", c2, " ", c3, " ", c4, " ", c5, " ", c6, "\n");
    return 0;
}
