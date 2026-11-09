// The nine quantile methods on one sample.
//
// Run it with: goose run examples/quantiles.goose
import std;
import statistics from "..";

fn main() {
    let xs: f64[>..] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
    let methods: f64[>..] = [];
    let names: const u8[:][9] = [
        "inverted_cdf          ",
        "averaged_inverted_cdf ",
        "closest_observation   ",
        "interpolated_inverted ",
        "hazen                 ",
        "weibull               ",
        "linear                ",
        "median_unbiased       ",
        "normal_unbiased       "
    ];
    for m in 0..9 {
        let q, c = statistics::quantile(xs[..], 0.3, m);
        print(names[m], " ", q, "  code ", c, "\n");
    }
    let iq, c2 = statistics::iqr(xs[..], statistics::METHOD_LINEAR);
    print("iqr ", iq, "  code ", c2, "\n");
    var cuts: f64[>..] = [];
    let n, c3 = statistics::cut_points(cuts, xs[..], 4, statistics::METHOD_LINEAR);
    print("cut points ", n, " code ", c3, "\n");
    for v in cuts { print("  ", v, "\n"); }
    return 0;
}
