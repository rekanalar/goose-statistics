// The frequencies of a sample: the empirical distribution function, the
// estimated distribution function, the cumulative and the relative
// frequencies, the entropy and the k-statistics.
//
// Run it with: goose run examples/frequencies.goose
import std;
import statistics from "..";

fn main() {
    let xs: f64[8] = [2.5, -1.0, 4.0, 3.25, 10.5, 0.5, 7.75, 3.0];
    let ranks: f64[4] = [1.0, 2.5, 2.5, 4.0];

    // The empirical distribution function.
    var quantiles: f64[>..] = [];
    var probs: f64[>..] = [];
    let _n, c1 = statistics::ecdf(quantiles, probs, xs[..]);
    print("ecdf quantiles: ", quantiles, "\n");
    print("ecdf probabilities: ", probs, ", code ", c1, "\n");

    // The Greenwood band of the empirical distribution function.
    var low: f64[>..] = [];
    var high: f64[>..] = [];
    let _n2, c2 = statistics::ecdf_interval(low, high, xs[..], 0.95);
    print("ecdf band low: ", low, "\n");
    print("ecdf band high: ", high, ", code ", c2, "\n");

    // The Kaplan-Meier curve of a censored sample.
    let uncensored: f64[3] = [1.0, 3.0, 5.0];
    let right: f64[2] = [2.0, 4.0];
    var times: f64[>..] = [];
    var survival: f64[>..] = [];
    let _n3, c3 = statistics::ecdf_censored(times, survival, uncensored[..], right[..]);
    print("kaplan-meier times: ", times, "\n");
    print("kaplan-meier survival: ", survival, ", code ", c3, "\n");

    // The estimated distribution function, with two methods.
    let e1, c4 = statistics::estimated_cdf(xs[..], 3.0, statistics::CDF_LINEAR);
    let e2, c5 = statistics::estimated_cdf(xs[..], 3.0, statistics::CDF_HAZEN);
    print("estimated cdf at 3: linear ", e1, ", hazen ", e2, ", codes ", c4, " ", c5, "\n");

    // The cumulative and the relative frequencies.
    var counts: f64[>..] = [];
    let f1, c6 = statistics::cumfreq_limits(counts, xs[..], 4, 0.0, 11.0);
    print("cumfreq counts: ", counts, ", lower limit ", f1.lowerlimit,
        ", bin size ", f1.binsize, ", code ", c6, "\n");
    var freqs: f64[>..] = [];
    let f2, c7 = statistics::relfreq_limits(freqs, xs[..], 4, 0.0, 11.0);
    print("relfreq shares: ", freqs, ", code ", c7, "\n");

    // The tie correction of a rank test.
    let t, c8 = statistics::tiecorrect(ranks[..]);
    print("tie correction: ", t, ", code ", c8, "\n");

    // The entropy of a set of probabilities.
    let pk: f64[2] = [0.5, 0.5];
    let qk: f64[2] = [0.9, 0.1];
    let e, c9 = statistics::entropy(pk[..], 2.0);
    let r, c10 = statistics::relative_entropy(pk[..], qk[..], 0.0);
    print("entropy in bits: ", e, ", relative entropy: ", r, ", codes ", c9, " ", c10, "\n");

    // The k-statistics of a sample.
    let k2, c11 = statistics::kstat(xs[..], 2);
    let k3, c12 = statistics::kstat(xs[..], 3);
    let kv, c13 = statistics::kstatvar(xs[..], 2);
    print("k-statistics: k2 ", k2, ", k3 ", k3, ", variance of k2 ", kv,
        ", codes ", c11, " ", c12, " ", c13, "\n");
}
