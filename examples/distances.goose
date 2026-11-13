// The distances of two samples, the multiscale graph correlation of a
// curve, the rank correlation and a closed-form law.
//
// Run it with: goose run examples/distances.goose
import std;
import statistics from "..";

fn main() {
    // The distance of two sets of points in two dimensions.
    let u: f64[8] = [0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 3.0, 3.0];
    let v: f64[6] = [1.0, 1.0, 2.0, 0.0, 4.0, 4.0];
    let w, cw = statistics::wasserstein_distance_nd(u[..], v[..], 2);
    print("distance of two sets of points: ", w, ", code ", cw, "\n");

    // The dependence of a curve, which a linear correlation misses.
    var xs: f64[>..] = [];
    var ys: f64[>..] = [];
    for k in 0..30 {
        let t = (k as f64 - 15.0) / 5.0;
        xs.push(t);
        ys.push(t * t);
    }
    var map: f64[>..] = [];
    var gen = rng { s: 21 };
    let m, cm = statistics::multiscale_graphcorr(map, xs[..], ys[..], 1, 1, 99, gen);
    print("graph correlation ", m.statistic, ", p ", m.pvalue, ", scale ", m.kx, " ", m.ky, ", code ", cm, "\n");
    let rho, cr = statistics::spearmanrho(xs[..], ys[..], statistics::ALT_TWO_SIDED);
    print("rank correlation ", rho.statistic, ", code ", cr, "\n");

    // The Cauchy law at its median and at its quartile.
    let half, c1 = statistics::cauchy_cdf(0.0, 0.0, 1.0);
    let q, c2 = statistics::cauchy_ppf(0.75, 0.0, 1.0);
    print("cdf at 0 ", half, ", upper quartile ", q, ", codes ", c1, " ", c2, "\n");
}
