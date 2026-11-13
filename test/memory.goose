// The memory probe of the hypothesis tests. It runs one call on the largest
// input of the package. The variable OP picks the call, from 1 to 27. The
// call 0 runs no test and gives the baseline. The two samples hold the
// input of the call, so the peak above the baseline is the memory of the
// call alone. Build the program and read the peak size:
//   goose build -o memory memory.goose
//   OP=1 /usr/bin/time -f "%M KB" ./memory
// The peak above the baseline stays under 8 times the input of the call.
import std;
import os;
import statistics from "..";

// The mean of a list. The resampling calls take this as their statistic.
fn mem_mean(xs: const f64[:]) -> f64 {
    if xs.len == 0 { return statistics::nan(); }
    var t: f64 = 0.0;
    for v in xs { t += v; }
    return t / xs.len;
}

fn main() {
    var arg: u8[>..] = [];
    var op: i64 = 0;
    if env("OP", arg) { op = arg[0] - 48; if arg.len > 1 { op = op * 10 + (arg[1] - 48); } }
    // The calls that read two samples hold half of the input in each one.
    // The variable HALF sets the size for the baseline run.
    var half_flag: u8[>..] = [];
    let half = statistics::MAX_INPUT / 2;
    var n1 = statistics::MAX_INPUT;
    if env("HALF", half_flag) { n1 = half; }
    var size: u8[>..] = [];
    if env("N", size) {
        n1 = 0;
        for ch in size { n1 = n1 * 10 + (ch - 48); }
    }
    var n2 = n1;
    var xs: f64[>..] = [];
    var ys: f64[>..] = [];
    var s: u64 = 12345;
    for i in 0..n1 {
        s = s * 6364136223846793005 + 1442695040888963407;
        xs.push((s >> 33) % 100003);
    }
    for i in 0..n2 {
        s = s * 6364136223846793005 + 1442695040888963407;
        ys.push((s >> 33) % 100003 + 0.5);
    }
    var c: i64 = 0;
    var sizes: i64[>..] = [n1 / 2, n1 - n1 / 2];
    var out: f64[>..] = [];
    var out2: f64[>..] = [];
    var out3: f64[>..] = [];
    var out4: f64[>..] = [];
    var out5: f64[>..] = [];
    if op == 1 { let r, cc = statistics::wilcoxon(xs[..]); c = cc; }
    if op == 2 { let r, cc = statistics::wilcoxon_paired(xs[..], ys[..]); c = cc; }
    if op == 3 { let r, cc = statistics::kruskal(xs[..], sizes[..]); c = cc; }
    if op == 4 { let r, cc = statistics::shapiro(xs[..]); c = cc; }
    if op == 5 { let r, cc = statistics::spearmanr(xs[..], ys[..]); c = cc; }
    if op == 6 { let r, cc = statistics::kendalltau(xs[..], ys[..]); c = cc; }
    if op == 7 { let r, cc = statistics::ks_2samp(xs[..], ys[..]); c = cc; }
    if op == 8 { let r, cc = statistics::mannwhitneyu(xs[..], ys[..]); c = cc; }
    if op == 9 { let r, cc = statistics::ks_1samp(xs[..], statistics::FAMILY_NORMAL, [50000.0, 30000.0][..]); c = cc; }
    if op == 10 { let r, cc = statistics::ttest_rel(xs[..], ys[..]); c = cc; }
    if op == 11 {
        // The probabilities for the correction lie in 0 to 1.
        var ps: f64[>..] = [];
        for y in ys { ps.push(1.0 / (1.0 + y)); }
        let k, cc = statistics::multipletests(out, ps[..], statistics::ADJUST_BH);
        c = cc;
    }
    if op == 12 { let r, cc = statistics::f_oneway(xs[..], sizes[..]); c = cc; }
    if op == 13 {
        var ps: f64[>..] = [];
        for y in ys { ps.push(1.0 / (1.0 + y)); }
        let r, cc = statistics::combine_pvalues(ps[..], statistics::COMBINE_STOUFFER, xs[..]);
        c = cc;
    }
    if op == 14 {
        var ps: f64[>..] = [];
        for y in ys { ps.push(1.0 / (1.0 + y)); }
        let k, cc = statistics::false_discovery_control(out, ps[..], statistics::ADJUST_BY);
        c = cc;
    }
    // The new calls of the variance tests, the transforms, the resampling
    // calls and the survival test.
    if op == 15 { let r, cc = statistics::levene(xs[..], sizes[..]); c = cc; }
    if op == 16 { let r, cc = statistics::fligner(xs[..], sizes[..]); c = cc; }
    if op == 17 { let r, cc = statistics::ranksums(xs[..], ys[..]); c = cc; }
    if op == 18 { let r, cc = statistics::mood(xs[..], ys[..]); c = cc; }
    if op == 19 { let r, cc = statistics::brunnermunzel(xs[..], ys[..]); c = cc; }
    if op == 20 { let r, cc = statistics::ansari(xs[..], ys[..]); c = cc; }
    if op == 21 { let r, cc = statistics::friedmanchisquare(xs[..], 2); c = cc; }
    if op == 22 {
        // The transform needs positive values for the Box-Cox form.
        var pos: f64[>..] = [];
        for v in xs { pos.push(v + 1.0); }
        let r, cc = statistics::boxcox(out, pos[..], 0.5);
        c = cc;
    }
    if op == 23 {
        var pos: f64[>..] = [];
        for v in xs { pos.push(v + 1.0); }
        let r, cc = statistics::boxcox_fit(out, pos[..]);
        c = cc;
    }
    if op == 24 {
        // The bootstrap draws one resample at a time, so the null list
        // holds the statistic of each resample.
        var gen = rng { s: 7 };
        let m = statistics::ResampleMethod { kind: statistics::RESAMPLE_BOOTSTRAP, resamples: 200, .. };
        let r, cc = statistics::bootstrap(out, xs[..], sizes[..], m, 0.95, statistics::ALT_TWO_SIDED, gen) { s => mem_mean(s) };
        c = cc;
    }
    if op == 25 {
        var gen = rng { s: 7 };
        let m = statistics::ResampleMethod { kind: statistics::RESAMPLE_PERMUTATION, resamples: 200, .. };
        let r, cc = statistics::permutation_test(out, xs[..], sizes[..], m, statistics::ALT_TWO_SIDED, gen) { s => mem_mean(s) };
        c = cc;
    }
    if op == 26 {
        var gen = rng { s: 7 };
        let r, cc = statistics::goodness_of_fit(out, xs[0..1000000], statistics::FAMILY_NORMAL, statistics::GOF_KOLMOGOROV, 1, gen);
        c = cc;
    }
    if op == 27 { let r, cc = statistics::epps_singleton_2samp(xs[..], ys[..]); c = cc; }
    if op == 28 { let r, cc = statistics::describe(xs[..], 1, false); c = cc; }
    if op == 29 { let r, cc = statistics::lmoment(out, xs[..], 4, true); c = cc; }
    if op == 30 { let r, cc = statistics::expectile(xs[..], 0.3); c = cc; }
    if op == 31 { let r, cc = statistics::kstat(xs[..], 4); c = cc; }
    if op == 32 { let r, cc = statistics::differential_entropy(xs[..], statistics::ENTROPY_VASICEK, 0, 0.0); c = cc; }
    if op == 33 { let r, cc = statistics::estimated_cdf(xs[..], 0.1, statistics::CDF_LINEAR); c = cc; }
    if op == 34 { let r, cc = statistics::ecdf(out, out2, xs[..]); c = cc; }
    if op == 35 { let r, cc = statistics::ecdf_interval(out, out2, xs[..], 0.95); c = cc; }
    if op == 36 { let r, cc = statistics::ecdf_censored(out, out2, xs[..], ys[..]); c = cc; }
    if op == 37 { let r, cc = statistics::cumfreq(out, xs[..], 1000); c = cc; }
    if op == 38 { let r, cc = statistics::trimboth(out, xs[..], 0.1); c = cc; }
    if op == 39 { let r, cc = statistics::circmean(xs[..], 0.0, 6.283185307179586); c = cc; }
    if op == 40 { let r, cc = statistics::binned_statistic(out, out2, out3, out4, out5, xs[..], ys[..], statistics::STAT_MEDIAN, 1000); c = cc; }
    print(c, "\n");
}
