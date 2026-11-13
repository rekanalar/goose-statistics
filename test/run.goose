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
    bad += run_group("distributions_location");
    bad += run_group("contingency");
    bad += run_group("divergence");
    bad += run_group("clip");
    bad += run_group("assoc_rank");
    bad += run_group("distance");
    bad += run_group("distance_nd");
    bad += run_group("assoc_mgc");
    bad += run_group("line_siegel");
    bad += run_group("laws");
    bad += run_group("probplot");
    bad += run_group("anderson_families");
    bad += run_group("fitting_mle");
    bad += run_group("assoc_somers");
    bad += run_group("contingency_risk");
    bad += run_group("kde");
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
    if is_word(name, "distributions_location") {
        bad += suite::run_distributions_location_vectors();
        bad += suite::run_distributions_location_tests();
        bad += suite::run_distributions_location_hostile();
    }
    if is_word(name, "contingency") {
        bad += suite::run_contingency_vectors();
        bad += suite::run_contingency_tests();
        bad += suite::run_contingency_hostile();
    }
    if is_word(name, "divergence") {
        bad += suite::run_divergence_vectors();
        bad += suite::run_divergence_tests();
        bad += suite::run_divergence_hostile();
    }
    if is_word(name, "clip") {
        bad += suite::run_clip_vectors();
        bad += suite::run_clip_tests();
        bad += suite::run_clip_hostile();
    }
    if is_word(name, "assoc_rank") {
        bad += suite::run_assoc_rank_vectors();
        bad += suite::run_assoc_rank_tests();
        bad += suite::run_assoc_rank_hostile();
    }
    if is_word(name, "distance") {
        bad += suite::run_distance_vectors();
        bad += suite::run_distance_tests();
        bad += suite::run_distance_hostile();
    }
    if is_word(name, "distance_nd") {
        bad += suite::run_distance_nd_vectors();
        bad += suite::run_distance_nd_tests();
        bad += suite::run_distance_nd_hostile();
    }
    if is_word(name, "assoc_mgc") {
        bad += suite::run_assoc_mgc_vectors();
        bad += suite::run_assoc_mgc_tests();
        bad += suite::run_assoc_mgc_hostile();
    }
    if is_word(name, "line_siegel") {
        bad += suite::run_line_siegel_vectors();
        bad += suite::run_line_siegel_tests();
        bad += suite::run_line_siegel_hostile();
    }
    if is_word(name, "laws") {
        bad += suite::run_laws_vectors();
        bad += suite::run_laws_tests();
        bad += suite::run_laws_hostile();
    }
    if is_word(name, "probplot") {
        bad += suite::run_probplot_vectors();
        bad += suite::run_probplot_tests();
        bad += suite::run_probplot_hostile();
    }
    if is_word(name, "anderson_families") {
        bad += suite::run_anderson_families_vectors();
        bad += suite::run_anderson_families_tests();
        bad += suite::run_anderson_families_hostile();
    }
    if is_word(name, "fitting_mle") {
        bad += suite::run_fitting_mle_vectors();
        bad += suite::run_fitting_mle_tests();
        bad += suite::run_fitting_mle_hostile();
    }
    if is_word(name, "assoc_somers") {
        bad += suite::run_assoc_somers_vectors();
        bad += suite::run_assoc_somers_tests();
        bad += suite::run_assoc_somers_hostile();
    }
    if is_word(name, "contingency_risk") {
        bad += suite::run_contingency_risk_vectors();
        bad += suite::run_contingency_risk_tests();
        bad += suite::run_contingency_risk_hostile();
    }
    if is_word(name, "kde") {
        bad += suite::run_kde_vectors();
        bad += suite::run_kde_tests();
        bad += suite::run_kde_hostile();
    }
    if is_word(name, "time") {
        bad += suite::run_time_tests();
        bad += suite::run_time_descriptive(200000);
        bad += suite::run_time_entropy(200000);
        bad += suite::run_time_frequency(200000);
        bad += suite::run_time_trimmed(200000);
        bad += suite::run_time_circular(200000);
        bad += suite::run_time_binned(200000);
        bad += suite::run_time_distributions_location(200000);
        bad += suite::run_time_contingency(0);
        bad += suite::run_time_divergence(200000);
        bad += suite::run_time_clip(200000);
        bad += suite::run_time_assoc_rank(20000);
        bad += suite::run_time_distance(200000);
        bad += suite::run_time_distance_nd(256);
        bad += suite::run_time_assoc_mgc(1024);
        bad += suite::run_time_line_siegel(2000);
        bad += suite::run_time_laws(4096);
        bad += suite::run_time_probplot(200000);
        bad += suite::run_time_anderson_families(20000);
        bad += suite::run_time_fitting_mle(200);
        bad += suite::run_time_assoc_somers(200000);
        bad += suite::run_time_contingency_risk();
        bad += suite::run_time_kde();
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
