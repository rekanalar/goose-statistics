// The corrections of a list of p-values: the adjusted values, the
// combined value of independent tests and the false discovery control.
//
// Run it with: goose run examples/multipletests.goose
import std;
import statistics from "..";

fn main() {
    let ps: f64[8] = [0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205];

    // The adjustment of a list of p-values. The method is one of the
    // ADJUST_ constants.
    var adj: f64[>..] = [];
    let n1, c1 = statistics::multipletests(adj, ps[..], statistics::ADJUST_BH);
    print("benjamini-hochberg: ", n1, " values, code ", c1, "\n");
    for i in 0..n1 { print("  ", adj[i], "\n"); }

    let n2, c2 = statistics::multipletests(adj, ps[..], statistics::ADJUST_HOMMEL);
    print("hommel: ", n2, " values, code ", c2, "\n");
    for i in 0..n2 { print("  ", adj[i], "\n"); }

    let n3, c3 = statistics::multipletests(adj, ps[..], statistics::ADJUST_BY);
    print("benjamini-yekutieli: ", n3, " values, code ", c3, "\n");
    for i in 0..n3 { print("  ", adj[i], "\n"); }

    // The same values under the name of SciPy.
    let n4, c4 = statistics::false_discovery_control(adj, ps[..], statistics::ADJUST_BY);
    print("false discovery control: ", n4, " values, code ", c4, "\n");

    // The combined p-value of independent tests of one hypothesis.
    let f, cf = statistics::combine_pvalues(ps[..], statistics::COMBINE_FISHER);
    print("fisher: ", f.statistic, ", p = ", f.pvalue, ", code ", cf, "\n");
    let s, cs = statistics::combine_pvalues(ps[..], statistics::COMBINE_STOUFFER);
    print("stouffer: ", s.statistic, ", p = ", s.pvalue, ", code ", cs, "\n");
    let t, ct = statistics::combine_pvalues(ps[..], statistics::COMBINE_TIPPETT);
    print("tippett: ", t.statistic, ", p = ", t.pvalue, ", code ", ct, "\n");

    // The Stouffer method takes weights. Each weight says how much its
    // test counts.
    let weights: f64[8] = [1.0, 1.0, 2.0, 2.0, 1.0, 1.0, 3.0, 3.0];
    let w, cw = statistics::combine_pvalues(ps[..], statistics::COMBINE_STOUFFER, weights[..]);
    print("stouffer with weights: ", w.statistic, ", p = ", w.pvalue, ", code ", cw, "\n");

    // A p-value outside 0 to 1 gives the domain fault.
    let bad: f64[2] = [0.5, 1.5];
    let e, ce = statistics::multipletests(adj, bad[..], statistics::ADJUST_BH);
    print("a p-value above 1 gives code ", ce, "\n");
}
