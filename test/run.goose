// Tests for the statistics package. Run with: goose run test/run.goose
// The program must start in the package directory.

import std;
import statistics from "..";
import suite from "./suite";

fn main() {
    var bad: i64 = 0;
    bad += suite::run_vectors();
    bad += suite::run_weighted();
    bad += suite::run_accept();
    bad += suite::run_hostile();
    if bad > 0 {
        print(bad, " checks failed\n");
        assert(false);
    }
    print("ok\n");
    return 0;
}
