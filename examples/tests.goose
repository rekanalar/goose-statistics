// The hypothesis tests: a t test, a rank test, an exact count test, a
// correlation test and a correction of p-values for many tests.
//
// Run it with: goose run examples/tests.goose
import std;
import statistics from "..";

fn main() {
    let control: f64[8] = [4.1, 5.2, 4.8, 5.9, 5.1, 4.4, 5.6, 4.9];
    let treated: f64[8] = [5.8, 6.1, 5.5, 6.9, 6.4, 5.7, 7.2, 6.0];

    // The Welch t test. Each result holds the statistic, the p-value and
    // the degrees of freedom. A fault gives NaN and a code.
    let t, c1 = statistics::ttest_ind(control[..], treated[..], statistics::ALT_TWO_SIDED, false);
    print("welch t = ", t.statistic, ", p = ", t.pvalue, ", df = ", t.df, ", code ", c1, "\n");

    // The Mann-Whitney test needs no normal shape. With no ties and small
    // samples it counts the exact distribution.
    let u, c2 = statistics::mannwhitneyu(control[..], treated[..]);
    print("mann-whitney U = ", u.statistic, ", p = ", u.pvalue, ", code ", c2, "\n");

    // The exact binomial test: 7 successes in 10 trials against p = 0.5.
    let b, c3 = statistics::binomtest(7.0, 10.0, 0.5);
    print("binomial p = ", b.pvalue, ", code ", c3, "\n");

    // The Fisher exact test of a 2x2 table, row by row.
    let table: f64[4] = [8.0, 2.0, 1.0, 5.0];
    let f, c4 = statistics::fisher_exact(table[..]);
    print("fisher odds ratio = ", f.statistic, ", p = ", f.pvalue, ", code ", c4, "\n");

    // A correlation with its p-value. Kendall's tau is exact for small
    // samples without ties.
    let k, c5 = statistics::kendalltau(control[..], treated[..]);
    print("kendall tau = ", k.statistic, ", p = ", k.pvalue, ", code ", c5, "\n");

    // Three p-values of many tests, adjusted by the Holm step-down rule.
    let ps: f64[3] = [0.01, 0.04, 0.03];
    var adjusted: f64[>..] = [];
    let n, c6 = statistics::multipletests(adjusted, ps[..], statistics::ADJUST_HOLM);
    print("holm adjusted:");
    for p in adjusted { print(" ", p); }
    print(" (", n, " values, code ", c6, ")\n");

    // A bad input gives NaN and a code, never a wrong p-value.
    let bad, c7 = statistics::shapiro([1.0, 2.0][..]);
    print("shapiro of two values: p = ", bad.pvalue, ", code ", c7, "\n");
}
