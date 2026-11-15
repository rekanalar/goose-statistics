// The timing probe of the package. Run from the package directory with:
//   goose run test/speed_probe.goose
// One row prints per target: name, baseline milliseconds, current time.
// The program asserts nothing about wall-clock time.
import std;
import statistics from "..";
import suite from "./suite";

fn main() {
    let bad = suite::run_speed();
    if bad != 0 {
        print(bad, " checks failed\n");
        assert(false);
    }
}
