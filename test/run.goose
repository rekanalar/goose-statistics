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
        // The variable may hold several names, separated by commas.
        var start: i64 = 0;
        var i: i64 = 0;
        while i < which.len {
            if which[i] as i64 == 44 {
                bad += run_group(which[start..i]);
                start = i + 1;
            }
            i += 1;
        }
        if start < which.len { bad += run_group(which[start..]); }
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
    bad += run_group("fitting_shapes");
    bad += run_group("assoc_somers");
    bad += run_group("contingency_risk");
    bad += run_group("kde");
    bad += run_group("distributions_elementary");
    bad += run_group("distributions_bounded");
    bad += run_group("distributions_weibull");
    bad += run_group("distributions_exponent");
    bad += run_group("distributions_discrete_a");
    bad += run_group("distributions_discrete_b");
    bad += run_group("distributions_rice");
    bad += run_group("distributions_search_a");
    bad += run_group("distributions_search_b");
    bad += run_group("distributions_search_c");
    bad += run_group("distributions_search_d");
    bad += run_group("distributions_closed");
    bad += run_group("distributions_noncentral");
    bad += run_group("distributions_nch");
    bad += run_group("distributions_landau");
    bad += run_group("distributions_multivariate");
    bad += run_group("distributions_random_matrix");
    bad += run_group("distributions_objects");
    bad += run_group("distributions_mixture");
    bad += run_group("distributions_covariance");
    bad += run_group("distributions_histogram");
    bad += run_group("distributions_power");
    bad += run_group("distributions_growth");
    bad += run_group("distributions_truncated");
    bad += run_group("distributions_gamma");
    bad += run_group("distributions_generalized");
    bad += run_group("distributions_symmetric");
    bad += run_group("distributions_johnson");
    bad += run_group("distributions_special");
    bad += run_group("distributions_extreme");
    bad += run_group("distributions_skew");
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
    if is_word(name, "fitting_shapes") {
        bad += suite::run_fitting_shapes_vectors();
        bad += suite::run_fitting_shapes_tests();
        bad += suite::run_fitting_shapes_hostile();
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
    if is_word(name, "distributions_elementary") {
        bad += suite::run_distributions_elementary_vectors();
        bad += suite::run_distributions_elementary_tests();
        bad += suite::run_distributions_elementary_hostile();
        bad += suite::run_distributions_elementary_time();
    }
    if is_word(name, "distributions_bounded") {
        bad += suite::run_distributions_bounded_vectors();
        bad += suite::run_distributions_bounded_tests();
        bad += suite::run_distributions_bounded_hostile();
        bad += suite::run_distributions_bounded_time();
    }
    if is_word(name, "distributions_weibull") {
        bad += suite::run_distributions_weibull_vectors();
        bad += suite::run_distributions_weibull_tests();
        bad += suite::run_distributions_weibull_hostile();
        bad += suite::run_distributions_weibull_time();
    }
    if is_word(name, "distributions_exponent") {
        bad += suite::run_distributions_exponent_vectors();
        bad += suite::run_distributions_exponent_tests();
        bad += suite::run_distributions_exponent_hostile();
        bad += suite::run_distributions_exponent_time();
    }
    if is_word(name, "distributions_discrete_a") {
        bad += suite::run_distributions_discrete_a_vectors();
        bad += suite::run_distributions_discrete_a_tests();
        bad += suite::run_distributions_discrete_a_hostile();
        bad += suite::run_distributions_discrete_a_time();
    }
    if is_word(name, "distributions_discrete_b") {
        bad += suite::run_distributions_discrete_b_vectors();
        bad += suite::run_distributions_discrete_b_tests();
        bad += suite::run_distributions_discrete_b_hostile();
        bad += suite::run_distributions_discrete_b_time();
    }
    if is_word(name, "distributions_rice") {
        bad += suite::run_distributions_rice_vectors();
        bad += suite::run_distributions_rice_tests();
        bad += suite::run_distributions_rice_hostile();
        bad += suite::run_distributions_rice_time();
    }
    if is_word(name, "distributions_search_a") {
        bad += suite::run_distributions_search_a_vectors();
        bad += suite::run_distributions_search_a_tests();
        bad += suite::run_distributions_search_a_hostile();
        bad += suite::run_distributions_search_a_time();
    }
    if is_word(name, "distributions_search_b") {
        bad += suite::run_distributions_search_b_vectors();
        bad += suite::run_distributions_search_b_tests();
        bad += suite::run_distributions_search_b_hostile();
        bad += suite::run_distributions_search_b_time();
    }
    if is_word(name, "distributions_search_c") {
        bad += suite::run_distributions_search_c_vectors();
        bad += suite::run_distributions_search_c_tests();
        bad += suite::run_distributions_search_c_hostile();
        bad += suite::run_distributions_search_c_time();
    }
    if is_word(name, "distributions_search_d") {
        bad += suite::run_distributions_search_d_vectors();
        bad += suite::run_distributions_search_d_tests();
        bad += suite::run_distributions_search_d_hostile();
        bad += suite::run_distributions_search_d_time();
    }
    if is_word(name, "distributions_closed") {
        bad += suite::run_distributions_closed_vectors();
        bad += suite::run_distributions_closed_tests();
        bad += suite::run_distributions_closed_hostile();
        bad += suite::run_distributions_closed_time();
    }
    if is_word(name, "distributions_noncentral") {
        bad += suite::run_distributions_noncentral_vectors();
        bad += suite::run_distributions_noncentral_tests();
        bad += suite::run_distributions_noncentral_hostile();
        bad += suite::run_distributions_noncentral_time();
    }
    if is_word(name, "distributions_nch") {
        bad += suite::run_distributions_nch_vectors();
        bad += suite::run_distributions_nch_tests();
        bad += suite::run_distributions_nch_hostile();
        bad += suite::run_distributions_nch_time();
    }
    if is_word(name, "distributions_landau") {
        bad += suite::run_distributions_landau_vectors();
        bad += suite::run_distributions_landau_tests();
        bad += suite::run_distributions_landau_hostile();
        bad += suite::run_distributions_landau_time();
    }
    if is_word(name, "distributions_multivariate") {
        bad += suite::run_distributions_multivariate_tests();
        bad += suite::run_distributions_multivariate_hostile();
        bad += suite::run_distributions_multivariate_time();
    }
    if is_word(name, "distributions_random_matrix") {
        bad += suite::run_distributions_random_matrix_tests();
        bad += suite::run_distributions_random_matrix_hostile();
        bad += suite::run_distributions_random_matrix_time();
    }
    if is_word(name, "distributions_objects") {
        bad += suite::run_distributions_objects_vectors();
        bad += suite::run_distributions_objects_tests();
        bad += suite::run_distributions_objects_hostile();
        bad += suite::run_distributions_objects_time();
    }
    if is_word(name, "distributions_mixture") {
        bad += suite::run_distributions_mixture_vectors();
        bad += suite::run_distributions_mixture_tests();
        bad += suite::run_distributions_mixture_time();
    }
    if is_word(name, "distributions_covariance") {
        bad += suite::run_distributions_covariance_vectors();
        bad += suite::run_distributions_covariance_tests();
        bad += suite::run_distributions_covariance_time();
    }
    if is_word(name, "distributions_histogram") {
        bad += suite::run_distributions_histogram_vectors();
        bad += suite::run_distributions_histogram_tests();
    }
    if is_word(name, "distributions_power") {
        bad += suite::run_distributions_power_vectors();
        bad += suite::run_distributions_power_tests();
        bad += suite::run_distributions_power_hostile();
        bad += suite::run_distributions_power_time();
    }
    if is_word(name, "distributions_growth") {
        bad += suite::run_distributions_growth_vectors();
        bad += suite::run_distributions_growth_tests();
        bad += suite::run_distributions_growth_hostile();
        bad += suite::run_distributions_growth_time();
    }
    if is_word(name, "distributions_truncated") {
        bad += suite::run_distributions_truncated_vectors();
        bad += suite::run_distributions_truncated_tests();
        bad += suite::run_distributions_truncated_hostile();
        bad += suite::run_distributions_truncated_time();
    }
    if is_word(name, "distributions_gamma") {
        bad += suite::run_distributions_gamma_vectors();
        bad += suite::run_distributions_gamma_tests();
        bad += suite::run_distributions_gamma_hostile();
        bad += suite::run_distributions_gamma_time();
    }
    if is_word(name, "distributions_generalized") {
        bad += suite::run_distributions_generalized_vectors();
        bad += suite::run_distributions_generalized_tests();
        bad += suite::run_distributions_generalized_hostile();
        bad += suite::run_distributions_generalized_time();
    }
    if is_word(name, "distributions_symmetric") {
        bad += suite::run_distributions_symmetric_vectors();
        bad += suite::run_distributions_symmetric_tests();
        bad += suite::run_distributions_symmetric_hostile();
        bad += suite::run_distributions_symmetric_time();
    }
    if is_word(name, "distributions_johnson") {
        bad += suite::run_distributions_johnson_vectors();
        bad += suite::run_distributions_johnson_tests();
        bad += suite::run_distributions_johnson_hostile();
        bad += suite::run_distributions_johnson_time();
    }
    if is_word(name, "distributions_special") {
        bad += suite::run_distributions_special_vectors();
        bad += suite::run_distributions_special_tests();
        bad += suite::run_distributions_special_hostile();
        bad += suite::run_distributions_special_time();
    }
    if is_word(name, "distributions_extreme") {
        bad += suite::run_distributions_extreme_vectors();
        bad += suite::run_distributions_extreme_tests();
        bad += suite::run_distributions_extreme_hostile();
        bad += suite::run_distributions_extreme_time();
    }
    if is_word(name, "distributions_skew") {
        bad += suite::run_distributions_skew_vectors();
        bad += suite::run_distributions_skew_tests();
        bad += suite::run_distributions_skew_hostile();
        bad += suite::run_distributions_skew_time();
    }
    if is_word(name, "time") {
        bad += suite::run_time_tests();
        bad += suite::run_time_descriptive(20000);
        bad += suite::run_time_entropy(20000);
        bad += suite::run_time_frequency(20000);
        bad += suite::run_time_trimmed(20000);
        bad += suite::run_time_circular(20000);
        bad += suite::run_time_binned(20000);
        bad += suite::run_time_distributions_location(20000);
        bad += suite::run_time_contingency(0);
        bad += suite::run_time_divergence(20000);
        bad += suite::run_time_clip(20000);
        bad += suite::run_time_assoc_rank(2000);
        bad += suite::run_time_distance(20000);
        bad += suite::run_time_distance_nd(128);
        bad += suite::run_time_assoc_mgc(512);
        bad += suite::run_time_line_siegel(1000);
        bad += suite::run_time_laws(2048);
        bad += suite::run_time_probplot(20000);
        bad += suite::run_time_anderson_families(5000);
        bad += suite::run_time_fitting_mle(100);
        bad += suite::run_time_fitting_shapes(50);
        bad += suite::run_time_assoc_somers(20000);
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
