// The tests of variance, scale and rank: equal variances of three groups,
// two rank tests of two samples, a test of several samples and a median test.
//
// Run it with: goose run examples/variance_tests.goose
import std;
import statistics from "..";

fn main() {
    // Three groups in one buffer, with a list of their sizes.
    let groups: f64[18] = [
        2.1, 3.4, 1.9, 5.6, 4.4, 3.8,
        3.9, 5.2, 4.8, 6.1, 5.5, 4.1,
        1.2, 2.5, 2.2, 3.1, 1.8, 2.9];
    let sizes: i64[3] = [6, 6, 6];

    // Equal variances. The center is the median unless the call says more.
    let lv, c1 = statistics::levene(groups[..], sizes[..]);
    print("levene W = ", lv.statistic, ", p = ", lv.pvalue, ", code ", c1, "\n");
    let fl, c2 = statistics::fligner(groups[..], sizes[..]);
    print("fligner chi2 = ", fl.statistic, ", p = ", fl.pvalue, ", code ", c2, "\n");
    let bt, c3 = statistics::bartlett(groups[..], sizes[..]);
    print("bartlett T = ", bt.statistic, ", p = ", bt.pvalue, ", code ", c3, "\n");

    // Equal means of groups with unequal variances.
    let ag, c4 = statistics::alexandergovern(groups[..], sizes[..]);
    print("alexander-govern A = ", ag.statistic, ", p = ", ag.pvalue, ", code ", c4, "\n");

    // Two samples: a test of scale and a test of order.
    let xs: f64[6] = [1.0, 2.0, 2.5, 3.0, 4.0, 4.5];
    let ys: f64[7] = [0.5, 2.2, 5.0, 6.0, 7.5, 8.0, 9.5];
    let an, c5 = statistics::ansari(xs[..], ys[..]);
    print("ansari AB = ", an.statistic, ", p = ", an.pvalue, ", code ", c5, "\n");
    let bm, c6 = statistics::brunnermunzel(xs[..], ys[..]);
    print("brunner-munzel W = ", bm.statistic, ", p = ", bm.pvalue, ", df = ", bm.df, ", code ", c6, "\n");

    // Three samples measured on the same four blocks, row by row.
    let blocks: f64[12] = [1.0, 2.0, 3.0, 2.0, 3.0, 1.0, 1.5, 3.5, 2.5, 4.0, 2.0, 3.0];
    let fr, c7 = statistics::friedmanchisquare(blocks[..], 3);
    print("friedman chi2 = ", fr.statistic, ", p = ", fr.pvalue, ", code ", c7, "\n");

    // The median test writes its 2 by k table into out.
    var table: f64[>..] = [];
    let md, c8 = statistics::median_test(table, groups[..], sizes[..]);
    print("median chi2 = ", md.statistic, ", p = ", md.pvalue, ", code ", c8, "\n");

    // The t test from summary values.
    let tt, c9 = statistics::ttest_ind_from_stats(5.0, 1.2, 20.0, 5.8, 1.5, 25.0, false);
    print("t = ", tt.statistic, ", p = ", tt.pvalue, ", df = ", tt.df, ", code ", c9, "\n");
}
