// The round trip of the quantiles. Each quantile of a continuous family
// composes with its cumulative function. The quantile of a discrete
// family gives the smallest count whose cumulative function reaches the
// probability.
namespace suite;

import std;
import statistics from "../..";

// The parameters of one family. The discrete families hold their own
// counts.
fn round_args(fam: i64, p: f64) -> f64, f64, f64 {
    if fam == 0 { return 0.0, 1.0, 0.0; }
    if fam == 1 { return 0.0, 1.0, 0.0; }
    if fam == 2 { return 1.0, 0.0, 0.0; }
    if fam == 3 { return 0.0, 1.0, 0.0; }
    if fam == 4 { return 1.0, 1.0, 0.0; }
    if fam == 5 { return 2.0, 1.0, 0.0; }
    if fam == 6 { return 2.0, 3.0, 0.0; }
    if fam == 7 { return 4.0, 0.0, 0.0; }
    if fam == 8 { return 7.0, 0.0, 0.0; }
    if fam == 9 { return 3.0, 6.0, 0.0; }
    if fam == 10 { return 10.0, 0.5, 0.0; }
    if fam == 11 { return 3.0, 0.0, 0.0; }
    if fam == 12 { return 0.5, 0.0, 0.0; }
    if fam == 13 { return 3.0, 0.5, 0.0; }
    return 20.0, 10.0, 4.0;
}

// The count of the points of the round trip. The default gives 2,000
// points for each family, so 30,000 points in all. The mutation script
// sets it lower, because it runs the suite many times.
fn round_points() -> i64 {
    var buf: u8[>..] = [];
    if !env("STATISTICS_ROUND_N", buf) { return 2000; }
    var n: i64 = 0;
    for ch in buf {
        if ch < 48 || ch > 57 { return 2000; }
        n = n * 10 + (ch - 48);
    }
    if n < 2 { return 2000; }
    return n;
}

fn run_accept_round() -> i64 {
    var bad: i64 = 0;
    let points = round_points();
    for fam in 0..15 {
        var worst: f64 = 0.0;
        var worst_p: f64 = 0.0;
        for i in 0..points {
            // The probability runs from 1e-4 to 1 - 1e-4.
            let p = 1e-4 + (1.0 - 2e-4) * (i as f64) / ((points - 1) as f64);
            let a, b, c = round_args(fam, p);
            let x, cx = dist_call(fam, 3, p, p, a, b, c);
            if cx != 0 {
                bad += 1;
                print("FAIL the quantile of family ", fam, " at ", p, " gives the code ", cx, "\n");
                continue;
            }
            let v, cv = dist_call(fam, 1, x, p, a, b, c);
            if cv != 0 {
                bad += 1;
                print("FAIL the cumulative function of family ", fam, " gives the code ", cv, "\n");
                continue;
            }
            if fam < 10 {
                // The cumulative function of a continuous family returns
                // the probability of the quantile.
                let err = ulp_error(v, p);
                if err > worst { worst = err; worst_p = p; }
            } else {
                // The cumulative function of a discrete family reaches
                // the probability at the count, and stays below it at the
                // count below.
                if v < p {
                    bad += 1;
                    print("FAIL the cumulative function of family ", fam, " at the count ", x, " is ", v, " below ", p, "\n");
                }
                let w, cw = dist_call(fam, 1, x - 1.0, p, a, b, c);
                if cw == 0 && w >= p && x > 0.0 {
                    bad += 1;
                    print("FAIL the quantile of family ", fam, " is not the smallest count at ", p, "\n");
                }
            }
        }
        print("round family ", fam, " worst ", worst, " at ", worst_p, "\n");
    }
    return bad;
}
