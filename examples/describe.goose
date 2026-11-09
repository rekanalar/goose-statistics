// The location and the spread of one sample.
//
// Run it with: goose run examples/describe.goose
import std;
import statistics from "..";

fn main() {
    let xs: f64[>..] = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
    let m, c1 = statistics::mean(xs[..]);
    let md, c2 = statistics::median(xs[..]);
    let sd, c3 = statistics::stdev(xs[..]);
    let mo, c4 = statistics::mode(xs[..]);
    let sk, c5 = statistics::skew(xs[..], true);
    print("mean   ", m, "\n");
    print("median ", md, "\n");
    print("stdev  ", sd, "\n");
    print("mode   ", mo, "\n");
    print("skew   ", sk, "\n");
    print("codes  ", c1, " ", c2, " ", c3, " ", c4, " ", c5, "\n");
    return 0;
}
