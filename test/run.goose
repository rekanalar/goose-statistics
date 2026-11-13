// Tests for the statistics package. Run with: goose run test/run.goose
// The program must start in the package directory.

import std;
import os;
import statistics from "..";
import suite from "./suite";

fn main() {
    var bad: i64 = 0;
    // The variable STATISTICS_SUITES holds the name of a group of suites.
    // The empty name runs every suite. The name "tests" runs the suites of
    // the hypothesis tests. A group name runs the older suites of one file.
    // The mutation probes use a group, so a probe costs only the suite that
    // covers the rule it changes.
    var which: u8[>..] = [];
    let has = env("STATISTICS_SUITES", which);
    if !has || which.len == 0 {
        bad += run_older();
        bad += run_hypothesis();
    } else if is_word(which[..], "tests") {
        bad += run_hypothesis();
    } else {
        bad += run_group(which[..]);
    }
    if bad > 0 {
        print(bad, " checks failed\n");
        assert(false);
    }
    print("ok\n");
    return 0;
}

// Whether the bytes hold the word.
fn is_word(b: const u8[:], word: const u8[:]) -> bool {
    if b.len != word.len { return false; }
    for i in 0..b.len {
        if b[i] != word[i] { return false; }
    }
    return true;
}

// The suites of the hypothesis tests.
fn run_hypothesis() -> i64 {
    var bad: i64 = 0;
    bad += run_group("hypothesis");
    bad += run_group("normality");
    bad += run_group("multipletests");
    bad += run_group("variance");
    bad += run_group("transform");
    bad += run_group("exact");
    bad += run_group("resample");
    bad += run_group("multicomp");
    bad += run_group("fitting");
    bad += run_group("descriptive");
    bad += run_group("frequency");
    bad += run_group("entropy");
    bad += run_group("trimmed");
    bad += run_group("circular");
    bad += run_group("binned");
    bad += run_group("time");
    return bad;
}

// The groups of the older suites. A group holds the suites that share the
// calls of one file, so a probe of that file runs only this group.
fn run_group(name: u8[:]) -> i64 {
    var bad: i64 = 0;
    if is_word(name, "vectors") {
        bad += suite::run_vectors();
        bad += suite::run_weighted();
    }
    if is_word(name, "assoc") {
        bad += suite::run_assoc();
    }
    if is_word(name, "accept") {
        bad += suite::run_accept();
    }
    if is_word(name, "hostile") {
        bad += suite::run_hostile();
    }
    if is_word(name, "special") {
        bad += suite::run_special(false);
        bad += suite::run_special_extra();
    }
    if is_word(name, "distributions") {
        bad += suite::run_distributions(false);
        bad += suite::run_distribution_extra();
        bad += suite::run_accept_dist();
        bad += suite::run_accept_round();
    }
    if is_word(name, "hypothesis") {
        bad += suite::run_accept_tests();
        bad += suite::run_test_vectors();
        bad += suite::run_hostile_tests();
    }
    if is_word(name, "normality") {
        bad += suite::run_normality_vectors();
        bad += suite::run_normality_tests();
    }
    if is_word(name, "multipletests") {
        bad += suite::run_multipletests_vectors();
        bad += suite::run_multipletests_tests();
    }
    if is_word(name, "variance") {
        bad += suite::run_variance_vectors();
        bad += suite::run_variance_tests();
    }
    if is_word(name, "transform") {
        bad += suite::run_transform_vectors();
        bad += suite::run_transform_tests();
    }
    if is_word(name, "exact") {
        bad += suite::run_exact_vectors();
        bad += suite::run_exact_tests();
    }
    if is_word(name, "resample") {
        bad += suite::run_resample_vectors();
        bad += suite::run_resample_tests();
    }
    if is_word(name, "multicomp") {
        bad += suite::run_multicomp_vectors();
        bad += suite::run_multicomp_tests();
    }
    if is_word(name, "fitting") {
        bad += suite::run_fitting_vectors();
        bad += suite::run_fitting_tests();
    }
    if is_word(name, "descriptive") {
        bad += suite::run_descriptive_vectors();
        bad += suite::run_descriptive_tests();
        bad += suite::run_descriptive_hostile();
    }
    if is_word(name, "frequency") {
        bad += suite::run_frequency_vectors();
        bad += suite::run_frequency_tests();
        bad += suite::run_frequency_hostile();
    }
    if is_word(name, "entropy") {
        bad += suite::run_entropy_vectors();
        bad += suite::run_entropy_tests();
        bad += suite::run_entropy_hostile();
    }
    if is_word(name, "trimmed") {
        bad += suite::run_trimmed_vectors();
        bad += suite::run_trimmed_tests();
        bad += suite::run_trimmed_hostile();
    }
    if is_word(name, "circular") {
        bad += suite::run_circular_vectors();
        bad += suite::run_circular_tests();
        bad += suite::run_circular_hostile();
    }
    if is_word(name, "binned") {
        bad += suite::run_binned_vectors();
        bad += suite::run_binned_tests();
        bad += suite::run_binned_hostile();
    }
    if is_word(name, "time") {
        bad += suite::run_time_tests();
        bad += suite::run_time_descriptive(200000);
        bad += suite::run_time_entropy(200000);
        bad += suite::run_time_frequency(200000);
        bad += suite::run_time_trimmed(200000);
        bad += suite::run_time_circular(200000);
        bad += suite::run_time_binned(200000);
    }
    return bad;
}

fn run_older() -> i64 {
    var bad: i64 = 0;
    bad += run_group("vectors");
    bad += run_group("assoc");
    bad += run_group("accept");
    bad += run_group("hostile");
    bad += run_group("special");
    bad += run_group("distributions");
    return bad;
}
