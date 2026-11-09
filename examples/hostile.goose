// Each fault code of the package, one case for each.
//
// Run it with: goose run examples/hostile.goose
import std;
import statistics from "..";

fn show(what: u8[:], v: f64, code: i64) {
    print(what, " -> ", v, "  code ", code, "\n");
}

fn main() {
    let empty: f64[>..] = [];
    let v1, c1 = statistics::mean(empty[..]);
    show("empty sample      ", v1, c1);

    let one: f64[>..] = [1.0];
    let v2, c2 = statistics::variance(one[..]);
    show("one point         ", v2, c2);

    let nan: f64[>..] = [1.0, 0.0 / 0.0, 3.0];
    let v3, c3 = statistics::mean(nan[..]);
    show("a NaN in the sample", v3, c3);

    let inf: f64[>..] = [1.0, 1.0 / 0.0, 3.0];
    let v4, c4 = statistics::mean(inf[..]);
    show("an infinity       ", v4, c4);

    let xs: f64[>..] = [1.0, 2.0, 3.0];
    let v5, c5 = statistics::quantile(xs[..], 1.5, statistics::METHOD_LINEAR);
    show("a fraction above 1", v5, c5);

    let v6, c6 = statistics::quantile(xs[..], 0.5, 99);
    show("an unknown method ", v6, c6);

    let ws: f64[>..] = [1.0, -2.0, 3.0];
    let v7, c7 = statistics::weighted_mean(xs[..], ws[..]);
    show("a negative weight ", v7, c7);

    let v8, c8 = statistics::harmonic_mean([-1.0, 2.0][..]);
    show("a negative value  ", v8, c8);

    let big: f64[>..] = [1e308, -1e308];
    let v9, c9 = statistics::variance(big[..]);
    show("a variance too big", v9, c9);

    let v10, c10 = statistics::trim_mean(xs[..], 0.5);
    show("a trim of 0.5     ", v10, c10);

    return 0;
}
