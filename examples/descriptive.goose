// The descriptive helpers: the summary of a sample, the percentiles, the
// power mean, the L-moments, the expectiles and the geometric deviation.
//
// Run it with: goose run examples/descriptive.goose
import std;
import statistics from "..";

fn main() {
    let xs: f64[8] = [2.5, -1.0, 4.0, 3.25, 10.5, 0.5, 7.75, 3.0];
    let pos: f64[8] = [0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0];

    // The summary of the sample.
    let s, c1 = statistics::describe(xs[..], 1, false);
    print("describe: n ", s.n, ", mean ", s.mean, ", variance ", s.variance,
        ", skewness ", s.skewness, ", kurtosis ", s.kurtosis, ", code ", c1, "\n");

    // The percentiles, with the three methods of the interpolation.
    let q1, c2 = statistics::scoreatpercentile(xs[..], 25.0);
    let q2, c3 = statistics::scoreatpercentile_method(xs[..], 25.0, statistics::INTERP_LOWER);
    let q3, c4 = statistics::scoreatpercentile_method(xs[..], 25.0, statistics::INTERP_HIGHER);
    print("scoreatpercentile 25: ", q1, " lower ", q2, " higher ", q3,
        ", codes ", c2, " ", c3, " ", c4, "\n");

    // The share of the sample at or below a score.
    let p1, c5 = statistics::percentileofscore(xs[..], 3.0, statistics::PCT_WEAK);
    let p2, c6 = statistics::percentileofscore(xs[..], 3.0, statistics::PCT_STRICT);
    print("percentileofscore 3: weak ", p1, ", strict ", p2, ", codes ", c5, " ", c6, "\n");

    // The median and the mean of the absolute distance to the median.
    let mad, c7 = statistics::median_abs_deviation(xs[..], 1.0, false);
    let madn, c8 = statistics::median_abs_deviation_normal(xs[..]);
    print("median absolute deviation: ", mad, ", scaled ", madn, ", codes ", c7, " ", c8, "\n");

    // The power means of a sample of positive values.
    let g, c9 = statistics::pmean(pos[..], 0.0);
    let a1, c10 = statistics::pmean(pos[..], 1.0);
    let a2, c11 = statistics::pmean(pos[..], 2.0);
    print("pmean: geometric ", g, ", arithmetic ", a1, ", quadratic ", a2,
        ", codes ", c9, " ", c10, " ", c11, "\n");

    // The geometric deviation and its z-scores.
    let gs, c12 = statistics::gstd(pos[..], 1);
    print("geometric deviation: ", gs, ", code ", c12, "\n");
    var scores: f64[>..] = [];
    let _n, c13 = statistics::gzscore(scores, pos[..], 1);
    print("geometric z-scores: ", scores, ", code ", c13, "\n");

    // The L-moments, with and without the standardization.
    var lm: f64[>..] = [];
    let _n2, c14 = statistics::lmoment(lm, xs[..], 4, true);
    print("L-moments: ", lm, ", code ", c14, "\n");

    // The expectiles of a sample.
    let e1, c15 = statistics::expectile(xs[..], 0.5);
    let e2, c16 = statistics::expectile(xs[..], 0.8);
    print("expectiles: 0.5 gives ", e1, ", 0.8 gives ", e2, ", codes ", c15, " ", c16, "\n");
}
