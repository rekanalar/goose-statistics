// The memory probe of the hypothesis tests. It runs one call on the largest
// input of the package. The variable OP picks the call, from 1 to 14. The
// call 0 runs no test and gives the baseline. The two samples hold the
// input of the call, so the peak above the baseline is the memory of the
// call alone. Build the program and read the peak size:
//   goose build -o memory memory.goose
//   OP=1 /usr/bin/time -f "%M KB" ./memory
// The peak above the baseline stays under 8 times the input of the call.
import std;
import os;
import statistics from "..";

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
    print(c, "\n");
}
