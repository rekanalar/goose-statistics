// The comparisons of many group means and the Bayesian intervals:
// tukey_hsd, dunnett, bayes_mvs and mvsdist.
//
// Run it with: goose run examples/comparisons.goose
import std;
import statistics from "..";

fn main() {
    // Three groups of eight values, one after the other in one list.
    let groups: f64[24] = [
        2.1, 3.4, 1.9, 5.6, 4.4, 3.8, 2.7, 4.9,
        3.9, 5.2, 4.8, 6.1, 5.5, 4.1, 6.6, 5.0,
        1.2, 2.5, 2.2, 3.1, 1.8, 2.9, 2.4, 3.3];
    let sizes: i64[3] = [8, 8, 8];

    // Each row of the table is one pair: the difference of the means, the
    // p-value, and the low and the high end of the interval.
    var table: f64[>..] = [];
    let pairs, c = statistics::tukey_hsd(table, groups[..], sizes[..], 0.95, true);
    if c != 0 { print("tukey_hsd fault ", c, "\n"); return; }
    var at: i64 = 0;
    for i in 0..sizes.len - 1 {
        for j in i + 1..sizes.len {
            print("group ", i, " and ", j, ": diff ", table[at], " p ", table[at + 1], " interval ", table[at + 2], " to ", table[at + 3], "\n");
            at += 4;
        }
    }

    // Each other group against the control, which is the group 0. The call
    // is deterministic: the same input gives the same p-values.
    var versus: f64[>..] = [];
    let rows, cd = statistics::dunnett(versus, groups[..], sizes[..], 0, statistics::ALT_TWO_SIDED, 0.95);
    if cd != 0 { print("dunnett fault ", cd, "\n"); return; }
    for r in 0..rows {
        print("control against group ", r + 1, ": t ", versus[4 * r], " p ", versus[4 * r + 1], "\n");
    }

    // The posterior of the mean, the variance and the deviation of group 0.
    let first = groups[0..8];
    let m, cm = statistics::bayes_mvs(first, 0.9);
    if cm != 0 { print("bayes_mvs fault ", cm, "\n"); return; }
    print("mean ", m.mean, " in ", m.mean_interval.low, " to ", m.mean_interval.high, "\n");
    print("variance ", m.variance, " in ", m.variance_interval.low, " to ", m.variance_interval.high, "\n");
    print("deviation ", m.deviation, " in ", m.deviation_interval.low, " to ", m.deviation_interval.high, "\n");

    // The parameters of the three posteriors.
    let d, cp = statistics::mvsdist(first);
    print("the mean follows a t law with ", d.mean_df, " degrees of freedom and scale ", d.mean_scale, "\n");
}
