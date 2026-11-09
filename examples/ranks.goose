// The ranks of one sample and the Theil-Sen estimate of two samples.
//
// Run it with: goose run examples/ranks.goose
import std;
import statistics from "..";

fn main() {
    let xs: f64[>..] = [3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0];
    var out: f64[>..] = [];
    let n, c1 = statistics::rankdata(out, xs[..], statistics::TIE_AVERAGE);
    var out2: f64[>..] = [];
    let n2, c2 = statistics::rankdata(out2, xs[..], statistics::TIE_DENSE);
    let sx: f64[>..] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
    let sy: f64[>..] = [2.1, 3.9, 6.2, 8.1, 9.8, 12.2, 14.1];
    let ts, ti, c3 = statistics::theil_sen(sx[..], sy[..]);
    print("average ranks ", n, ": ");
    for r in out {
        print(r, " ");
    }
    print("\n");
    print("dense ranks   ", n2, ": ");
    for r in out2 {
        print(r, " ");
    }
    print("\n");
    print("theil slope    ", ts, "\n");
    print("theil intercept", ti, "\n");
    print("codes          ", c1, " ", c2, " ", c3, "\n");
    return 0;
}
