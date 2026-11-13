// The kernel density of a sample, the odds ratio of a table and the
// Anderson-Darling test of a Gumbel sample.
//
// Run it with: goose run examples/smoothing.goose
import std;
import statistics from "..";

fn main() {
    // A density of one dimension, with the rule of Scott.
    let xs: f64[6] = [1.0, 2.0, 2.5, 3.0, 4.5, 6.0];
    let k, ck = statistics::kde_fit(xs[..], 1, [], statistics::KDE_SCOTT, 0.0);
    let at: f64[3] = [1.0, 3.0, 5.0];
    var dens: f64[>..] = [];
    let cp = statistics::kde_pdf(k, dens, at[..]);
    print("density: ", dens, ", codes ", ck, " ", cp, "\n");
    let share, cs = statistics::kde_integrate_box_1d(k, 2.0, 4.0);
    print("share between 2 and 4: ", share, ", code ", cs, "\n");

    // The odds ratio of a 2 by 2 table, with its interval.
    let table: f64[4] = [12.0, 5.0, 3.0, 9.0];
    let r, co = statistics::odds_ratio(table[..], statistics::ODDS_SAMPLE, statistics::ALT_TWO_SIDED, 0.95);
    print("odds ratio ", r.statistic, " in [", r.low, ", ", r.high, "], code ", co, "\n");

    // The test of a Gumbel sample.
    let ys: f64[8] = [0.4, 1.1, 1.9, 0.9, 2.6, 1.4, 3.8, 0.2];
    var crit: f64[>..] = [];
    var level: f64[>..] = [];
    var params: f64[>..] = [];
    let t, ct = statistics::anderson_gumbel_r(crit, level, params, ys[..]);
    print("statistic ", t.statistic, ", p ", t.pvalue, ", code ", ct, "\n");
}
