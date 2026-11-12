// The resampling calls: the bootstrap interval, the permutation test,
// the Monte Carlo test, the power of a test and the Sobol indices.
//
// Run it with: goose run examples/resampling.goose
import std;
import statistics from "..";

// The mean of a list.
fn mean_of(s: const f64[:]) -> f64 {
    var t: f64 = 0.0;
    for v in s { t += v; }
    return t / s.len;
}

// The difference of the means of the two halves of a list.
fn diff_of(s: const f64[:]) -> f64 {
    let h = s.len / 2;
    return mean_of(s[0..h]) - mean_of(s[h..s.len]);
}

// The p-value of the one-sample t test of a list.
fn ttest_p(s: const f64[:]) -> f64 {
    let r, c = statistics::ttest_1samp(s, 0.0);
    return r.pvalue;
}

// The Ishigami function of a point of the unit cube.
fn ishigami(p: const f64[:]) -> f64 {
    let pi = 3.141592653589793;
    let x1 = (2.0 * p[0] - 1.0) * pi;
    let x2 = (2.0 * p[1] - 1.0) * pi;
    let x3 = (2.0 * p[2] - 1.0) * pi;
    return sin(x1) + 7.0 * sin(x2) * sin(x2) + 0.1 * pow(x3, 4.0) * sin(x1);
}

fn main() {
    let a: f64[8] = [2.1, 3.4, 1.9, 5.6, 4.4, 3.8, 2.7, 4.9];
    let b: f64[8] = [3.9, 5.2, 4.8, 6.1, 5.5, 4.1, 6.6, 5.0];
    var ab: f64[16];
    for i in 0..8 { ab[i] = a[i]; ab[8 + i] = b[i]; }
    let one: i64[1] = [8];
    let two: i64[2] = [8, 8];
    var nul: f64[>..] = [];

    // The caller owns the generator. One seed gives one answer.
    var r = rng { s: 3 };

    // The bootstrap interval of the mean. The method holds the kind, the
    // resample count and the form of the interval.
    let bm = statistics::ResampleMethod { resamples: 999, interval: statistics::INTERVAL_PERCENTILE, .. };
    let br, bc = statistics::bootstrap(nul, a[..], one[..], bm, 0.95, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    print("bootstrap percentile: ", br.interval.low, " to ", br.interval.high, ", error ", br.standard_error, ", code ", bc, "\n");
    let bam = statistics::ResampleMethod { resamples: 999, interval: statistics::INTERVAL_BCA, .. };
    let ba, bac = statistics::bootstrap(nul, a[..], one[..], bam, 0.95, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    print("bootstrap bca: ", ba.interval.low, " to ", ba.interval.high, ", code ", bac, "\n");

    // The permutation test of the difference of the means. The resample
    // count of 1000000 is above the count of the splits, so the test lists
    // every split and the p-value is exact.
    let pm = statistics::ResampleMethod { kind: statistics::RESAMPLE_PERMUTATION, resamples: 1000000, .. };
    let pt, pc = statistics::permutation_test(nul, ab[..], two[..], pm, statistics::ALT_TWO_SIDED, r) { s => diff_of(s) };
    print("permutation: statistic ", pt.statistic, ", p ", pt.pvalue, ", splits ", nul.len, ", code ", pc, "\n");

    // The Monte Carlo test against a normal null of mean 3.
    let mm = statistics::ResampleMethod { kind: statistics::RESAMPLE_MONTE_CARLO, resamples: 9999, .. };
    let np: f64[2] = [3.0, 1.5];
    let mt, mc = statistics::monte_carlo_test(nul, a[..], statistics::FAMILY_NORMAL, np[..], mm, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    print("monte carlo: statistic ", mt.statistic, ", p ", mt.pvalue, ", code ", mc, "\n");

    // The power of the one-sample t test for 20 values. The alternative
    // has a mean of 0.5 and the test checks against a mean of 0.
    let pn: f64[2] = [0.5, 1.0];
    let pw, pwc = statistics::power(statistics::FAMILY_NORMAL, pn[..], 20, 0.05, 2000, r) { s => ttest_p(s) };
    print("power: ", pw, ", code ", pwc, "\n");

    // The Sobol indices of the Ishigami function. The first-order indices
    // lie near 0.3139, 0.4424 and 0. The total indices lie near 0.5576,
    // 0.4424 and 0.2437.
    var first: f64[>..] = [];
    var total: f64[>..] = [];
    var rs = rng { s: 5 };
    let sn, sk = statistics::sobol_indices(first, total, 3, 4096, rs) { p => ishigami(p) };
    print("sobol: ", sn, " indices, code ", sk, "\n");
    for i in 0..first.len {
        print("  first ", first[i], ", total ", total[i], "\n");
    }
}
