// Tests for the statistics package. Run with: goose run test/run.goose
// The program must start in the package directory.

import std;
import os;
import statistics from "..";
import suite from "./suite";

fn main() {
    var bad: i64 = 0;
    // The variable STATISTICS_SUITES set to "tests" runs the suites of the
    // hypothesis tests only. The mutation probes of those calls use it.
    var which: u8[>..] = [];
    let only_tests = env("STATISTICS_SUITES", which) && is_word_tests(which[..]);
    if !only_tests {
        bad += run_older();
    }
    bad += suite::run_accept_tests();
    bad += suite::run_test_vectors();
    bad += suite::run_normality_vectors();
    bad += suite::run_normality_tests();
    bad += suite::run_hostile_tests();
    bad += suite::run_time_tests();
    if bad > 0 {
        print(bad, " checks failed\n");
        assert(false);
    }
    print("ok\n");
    return 0;
}

// Whether the bytes hold the word tests.
fn is_word_tests(b: const u8[:]) -> bool {
    if b.len != 5 { return false; }
    return b[0] == 116 && b[1] == 101 && b[2] == 115 && b[3] == 116 && b[4] == 115;
}

fn run_older() -> i64 {
    var bad: i64 = 0;
    bad += suite::run_vectors();
    bad += suite::run_assoc();
    bad += suite::run_weighted();
    bad += suite::run_accept();
    bad += suite::run_hostile();
    bad += suite::run_special(false);
    bad += suite::run_special_extra();
    bad += suite::run_distributions(false);
    bad += suite::run_distribution_extra();
    bad += suite::run_accept_dist();
    bad += suite::run_accept_round();
    return bad;
}
