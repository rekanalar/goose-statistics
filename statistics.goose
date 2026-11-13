// statistics - descriptive statistics for Goose.
//
// The package gives the location, the spread, the quantiles, the shape,
// the trimmed forms and the weighted forms of a sample of 64 bit floats.
// The sum is exact, so the mean and the variance keep their digits.
//
// Every call returns the result first and a fault code second. A fault
// gives NaN as the result. A call checks the length of its input before
// any other work. A call never changes the data of the caller, and a fault
// never holds data of the caller.
namespace statistics;

import check;
import exactsum;
import spread;
import order;
import shape;
import weighted;
import cross;
import ranks;
import line;
import matrix;
import special;
import distributions;
import tests;
import normality;
import numeric;
import variance;
import transform;
import exact;
import resample;
import multicomp;
import fitting;
import descriptive;
import frequency;
import entropy;
import trimmed;
import binned;
import circular;
import transform_yeo;
import variance_multi;
import resample_sobol;
import resample_permutation;
import exact_page;
import exact_goodness;
import exact_cvm;
import tests_multiplicity;
import tests_correlation;
import tests_ks;
import tests_wilcoxon;
import special_inverse;
import special_incomplete;
import special_gamma;
import special_erf;
import distributions_discrete;
import distributions_t;
import distributions_shape;
import distributions_location;
import contingency;
import divergence;
import clip;
import assoc_rank;
import distance;
import line_siegel;
import laws;
import probplot;
import anderson_families;
import fitting_mle;
import assoc_somers;
import contingency_risk;
import kde;
import distance_nd;
import assoc_mgc;

// ---- fault codes ----------------------------------------------------
const OK: i64 = 0;
const ERR_EMPTY: i64 = 1;       // the sample has no value
const ERR_TOO_SHORT: i64 = 2;   // the sample has fewer values than the call needs
const ERR_NAN: i64 = 3;         // the sample, a weight or an argument holds a NaN
const ERR_INFINITY: i64 = 4;    // the sample or a weight holds an infinity
const ERR_DOMAIN: i64 = 5;      // an argument or a value is out of range for the call
const ERR_LENGTH: i64 = 6;      // the list of weights has another length than the sample
const ERR_UNDEFINED: i64 = 7;   // the result has no defined value
const ERR_OVERFLOW: i64 = 8;    // the true result does not fit in a 64 bit float
const ERR_TOO_LONG: i64 = 9;    // the sample has more than MAX_INPUT values
const ERR_FIT: i64 = 10;        // a search or a fit does not end within its step bound

// ---- bounds ---------------------------------------------------------
const MAX_INPUT: i64 = 4194304;  // values in a sample

// The difference between 1 and the next float above it. The searches
// and the tests compare against it.
const EPS: f64 = 2.220446049250313e-16;

// The relative bound of a root search.
const BRENT_RTOL: f64 = 8.881784197001252e-16;
const MAX_ORDER: i64 = 8;        // the largest moment order
const MAX_THEIL: i64 = 4096;     // points in a Theil-Sen estimate
const MAX_COLUMNS: i64 = 1024;   // columns in a matrix
const MAX_WORK: i64 = 268435456; // cell computations in a matrix

// The bounds of the special functions. The measurement of the prototypes
// set each bound. A loop that meets its bound gives the overflow code.
const MAX_SHAPE: f64 = 10000000000.0;   // the largest shape of an incomplete call
const MAX_ITER_SERIES: i64 = 2048;      // terms of a power series
const MAX_ITER_CF: i64 = 40000;         // steps of a continued fraction
const MAX_ITER_INVERSE: i64 = 100;      // steps of an inverse loop

// The bounds of the distributions.
const MAX_COUNT: f64 = 9007199254740992.0;  // 2^53, the largest exact count
const MAX_DIST_SUPPORT: i64 = 100000;   // counts in the support of a hypergeometric sum
const MAX_ITER_SEARCH: i64 = 200;       // steps of the search of a discrete quantile

// ---- tie methods -------------------------------------------------------
const TIE_AVERAGE: i64 = 0;   // the tied values share the mean of their ranks
const TIE_MIN: i64 = 1;       // the tied values share the lowest rank of the run
const TIE_MAX: i64 = 2;       // the tied values share the highest rank of the run
const TIE_DENSE: i64 = 3;     // the tied values share the count of the distinct smaller values plus one
const TIE_ORDINAL: i64 = 4;   // the tied values take the next ranks in the order of the input

// ---- Kendall variants --------------------------------------------------
const KENDALL_B: i64 = 0;   // the tau-b coefficient
const KENDALL_C: i64 = 1;   // the tau-c coefficient

// ---- quantile methods -------------------------------------------------
const METHOD_INVERTED_CDF: i64 = 0;
const METHOD_AVERAGED_INVERTED_CDF: i64 = 1;
const METHOD_CLOSEST_OBSERVATION: i64 = 2;
const METHOD_INTERPOLATED_INVERTED_CDF: i64 = 3;
const METHOD_HAZEN: i64 = 4;
const METHOD_WEIBULL: i64 = 5;
const METHOD_LINEAR: i64 = 6;
const METHOD_MEDIAN_UNBIASED: i64 = 7;
const METHOD_NORMAL_UNBIASED: i64 = 8;
const METHOD_EXCLUSIVE: i64 = 5;   // the exclusive method of Python, the same as WEIBULL
const METHOD_INCLUSIVE: i64 = 6;   // the inclusive method of Python, the same as LINEAR

// ---- weighted variance forms ------------------------------------------
const WEIGHT_POPULATION: i64 = 0;   // divide by the sum of the weights
const WEIGHT_FREQUENCY: i64 = 1;    // divide by the sum minus 1
const WEIGHT_RELIABILITY: i64 = 2;  // divide by the sum minus the sum of squares over the sum

// The alternative of a test.
const ALT_TWO_SIDED: i64 = 0;   // the p-value covers both tails
const ALT_LESS: i64 = 1;        // the first sample lies below the second
const ALT_GREATER: i64 = 2;     // the first sample lies above the second

// The method of the Mann-Whitney test.
const MW_AUTO: i64 = 0;         // the exact form where the rules allow it
const MW_EXACT: i64 = 1;        // the exact form, or a fault
const MW_ASYMPTOTIC: i64 = 2;   // the normal form

// The treatment of a zero difference in the Wilcoxon test.
const ZERO_WILCOX: i64 = 0;     // drop the zero differences
const ZERO_PRATT: i64 = 1;      // rank the zeros, then adjust the mean and the variance
const ZERO_ZSPLIT: i64 = 2;     // give half of the zero ranks to each sum

// The correction of a list of p-values.
const ADJUST_BONFERRONI: i64 = 0;   // multiply by the count
const ADJUST_HOLM: i64 = 1;         // the step-down form of Holm
const ADJUST_BH: i64 = 2;           // the step-up form of Benjamini and Hochberg
const ADJUST_HOCHBERG: i64 = 3;     // the step-up form of Holm
const ADJUST_HOMMEL: i64 = 4;       // the closed form of Hommel, in the square of the count
const ADJUST_BY: i64 = 5;           // Benjamini and Yekutieli, for any dependence
const ADJUST_NONE: i64 = 6;         // the values unchanged
const HOMMEL_MAX: i64 = 10000;      // the longest list of the Hommel form

// The method that combines p-values of independent tests.
const COMBINE_FISHER: i64 = 0;      // minus twice the sum of the logarithms
const COMBINE_STOUFFER: i64 = 1;    // the sum of the normal quantiles, with weights
const COMBINE_TIPPETT: i64 = 2;     // the smallest p-value

// The family of the one-sample Kolmogorov-Smirnov test.
const FAMILY_NORMAL: i64 = 0;
const FAMILY_UNIFORM: i64 = 1;
const FAMILY_EXPONENTIAL: i64 = 2;
const FAMILY_LOGNORMAL: i64 = 3;
const FAMILY_WEIBULL: i64 = 4;
const FAMILY_GAMMA: i64 = 5;
const FAMILY_BETA: i64 = 6;
const FAMILY_CHI2: i64 = 7;
const FAMILY_T: i64 = 8;
const FAMILY_F: i64 = 9;

// The bounds of the exact forms of the tests.
const MAX_EXACT_MW: i64 = 8;          // the smaller sample of an exact Mann-Whitney test
const MAX_MW_CELLS: i64 = 100000;     // cells in the table of the exact Mann-Whitney test
const MAX_EXACT_WILCOXON: i64 = 50;   // pairs in an exact Wilcoxon test
const MAX_ENUM_WILCOXON: i64 = 13;    // pairs whose sign patterns the call counts
const MAX_EXACT_KENDALL: i64 = 33;    // pairs in an exact Kendall test
const MAX_EXACT_KS: i64 = 10000;      // the size bound of an exact Kolmogorov-Smirnov test

// The center of the Levene and Fligner-Killeen tests.
const CENTER_MEAN: i64 = 0;       // the mean of each group
const CENTER_MEDIAN: i64 = 1;     // the median of each group
const CENTER_TRIMMED: i64 = 2;    // the trimmed mean of each group

// The kind of a resampling method.
const RESAMPLE_BOOTSTRAP: i64 = 0;
const RESAMPLE_PERMUTATION: i64 = 1;
const RESAMPLE_MONTE_CARLO: i64 = 2;

// The interval of a bootstrap.
const INTERVAL_PERCENTILE: i64 = 0;   // the quantiles of the bootstrap values
const INTERVAL_BASIC: i64 = 1;        // the quantiles reflected about the estimate
const INTERVAL_BCA: i64 = 2;          // the bias-corrected and accelerated form

// The permutation type of a permutation test.
const PERMUTE_INDEPENDENT: i64 = 0;   // the samples swap their values
const PERMUTE_SAMPLES: i64 = 1;       // the values swap within each pair
const PERMUTE_PAIRINGS: i64 = 2;      // each sample is shuffled on its own

// The statistic of a goodness-of-fit test.
const GOF_ANDERSON: i64 = 0;
const GOF_KOLMOGOROV: i64 = 1;
const GOF_CRAMER: i64 = 2;

// The bounds of resampling, of the exact forms and of the fits.
const MAX_RESAMPLES: i64 = 1000000;          // resamples in one call
const MAX_RESAMPLE_WORK: i64 = 500000000;    // resamples times the sample size
const MAX_FIT_WORK: i64 = 2000000;            // resamples times the sample size of a goodness-of-fit test
const MAX_GOF_REPS: i64 = 100000;         // Monte Carlo draws of a goodness-of-fit test
const MAX_SOBOL_POINTS: i64 = 1048576;       // base points of a Sobol estimate
const MAX_SOBOL_DIMS: i64 = 32;              // input variables of a Sobol estimate
const MAX_EXACT_TABLE: i64 = 400;            // the total of a table in an exact 2 by 2 test
const MAX_NUISANCE_GRID: i64 = 1024;         // grid points of the nuisance parameter
const MAX_PAGE_COLUMNS: i64 = 12;            // columns of an exact Page test
const MAX_EXACT_CELLS: i64 = 100000;         // cells of an exact count table
const MAX_FIT_STEPS: i64 = 200;              // steps of a search or a fit

// ---- bounds of the descriptive and shape calls ----------------------
const MAX_BINS: i64 = 1048576;         // bins of a binned statistic
const MAX_CELLS: i64 = 4194304;        // cells of a binned table
const MAX_DIMS: i64 = 16;              // dimensions of a binned table
const MAX_LORDER: i64 = 16;            // the largest L-moment order
const MAX_EXPECTILE_STEPS: i64 = 200;  // steps of the search of an expectile

// ---- kinds of the percentile calls ---------------------------------
const PCT_RANK: i64 = 0;      // the mean rank of the tied values
const PCT_WEAK: i64 = 1;      // the share of the values at or below the score
const PCT_STRICT: i64 = 2;    // the share of the values below the score
const PCT_MEAN: i64 = 3;      // the mean of the weak and the strict share

// The methods of the percentile of a score.
const INTERP_FRACTION: i64 = 0;  // the linear interpolation
const INTERP_LOWER: i64 = 1;     // the lower value
const INTERP_HIGHER: i64 = 2;    // the higher value

// The methods of the estimated distribution function.
const CDF_INTERPOLATED: i64 = 0;      // the interpolated inverted function
const CDF_HAZEN: i64 = 1;             // the Hazen function
const CDF_WEIBULL: i64 = 2;           // the Weibull function
const CDF_LINEAR: i64 = 3;            // the linear function
const CDF_MEDIAN_UNBIASED: i64 = 4;   // the median unbiased function
const CDF_NORMAL_UNBIASED: i64 = 5;   // the normal unbiased function
const CDF_INVERTED: i64 = 6;          // the inverted function
const CDF_AVERAGED_INVERTED: i64 = 7; // the averaged inverted function
const CDF_CLOSEST: i64 = 8;           // the closest observation

// The estimators of the differential entropy.
const ENTROPY_AUTO: i64 = 0;       // the estimator of the size of the sample
const ENTROPY_VASICEK: i64 = 1;    // the Vasicek estimator
const ENTROPY_VAN_ES: i64 = 2;     // the van Es estimator
const ENTROPY_EBRAHIMI: i64 = 3;   // the Ebrahimi estimator
const ENTROPY_CORREA: i64 = 4;     // the Correa estimator
const ENTROPY_SMALL: i64 = 10;     // the largest sample that ENTROPY_AUTO gives to van Es
const ENTROPY_MEDIUM: i64 = 1000;  // the largest sample that ENTROPY_AUTO gives to Ebrahimi

// The statistics of a binned table.
const STAT_MEAN: i64 = 0;      // the mean of the values of the bin
const STAT_STD: i64 = 1;       // the deviation of the values of the bin
const STAT_MEDIAN: i64 = 2;    // the median of the values of the bin
const STAT_COUNT: i64 = 3;     // the number of the values of the bin
const STAT_SUM: i64 = 4;       // the sum of the values of the bin
const STAT_MIN: i64 = 5;       // the lowest value of the bin
const STAT_MAX: i64 = 6;       // the highest value of the bin

// The tail that a trim call cuts.
const TAIL_LEFT: i64 = 0;      // the lower tail
const TAIL_RIGHT: i64 = 1;     // the upper tail

// ---- structs of the descriptive and shape calls ---------------------
// The summary of a sample. The count is at least 1. The variance is
// for the given ddof. The skewness and the kurtosis are for the given
// bias. The kurtosis is the Fisher form.
struct Summary {
    n: i64,
    min: f64,
    max: f64,
    mean: f64,
    variance: f64,
    skewness: f64,
    kurtosis: f64,
}

// The counts of a histogram. The lower limit is the low edge of the
// first bin. The bin size is above 0. The extra points count the
// values outside the real limits.
struct Frequencies {
    lowerlimit: f64,
    binsize: f64,
    extrapoints: i64,
}

// The sizes of a binned table. The dims is the number of dimensions.
// The cells is the number of bins in all. The points is the number of
// the values of the caller.
struct BinnedTable {
    dims: i64,
    cells: i64,
    points: i64,
}

// ---- bounds of the association, distance and distribution calls -------
const MAX_TABLE_CELLS: i64 = 1048576;   // cells of a contingency table
const MAX_TABLE_SIDE: i64 = 1024;       // rows or columns of a contingency table
const MAX_MGC_PAIRS: i64 = 1024;        // pairs of a multiscale graph correlation
const MAX_MGC_REPS: i64 = 10000;        // permutations of a multiscale graph correlation
const MAX_MGC_WORK: i64 = 8000000;      // cells of all the maps of a multiscale graph correlation
const MAX_NDW_POINTS: i64 = 256;        // points of a set of an n-dimensional distance
const MAX_NDW_DIMS: i64 = 16;           // dimensions of an n-dimensional distance
const MAX_KDE_POINTS: i64 = 65536;      // sample points of a kernel density
const MAX_KDE_DIMS: i64 = 16;           // dimensions of a kernel density
const MAX_KDE_WORK: i64 = 67108864;     // sample points times evaluation points
const MAX_PPCC_GRID: i64 = 4096;        // grid points of a probability plot correlation
const MAX_PPCC_POINTS: i64 = 65536;     // sample values of a probability plot search
const MAX_DIRMULT_K: i64 = 4096;        // categories of a Dirichlet multinomial
const MAX_ORDER_N: i64 = 1048576;       // size of an order statistic
const MAX_SIEGEL_POINTS: i64 = 2048;      // points of a Siegel estimate

// The methods of the association of a contingency table.
const ASSOC_CRAMER: i64 = 0;
const ASSOC_TSCHUPROW: i64 = 1;
const ASSOC_PEARSON: i64 = 2;

// The kinds of an odds ratio.
const ODDS_CONDITIONAL: i64 = 0;   // the maximum likelihood estimate of the conditional law
const ODDS_SAMPLE: i64 = 1;        // the ratio of the cross products

// The rules of the bandwidth of a kernel density.
const KDE_SCOTT: i64 = 0;
const KDE_SILVERMAN: i64 = 1;
const KDE_SCALAR: i64 = 2;          // the caller gives the factor

// The families that truncate and order_statistic take. Every family reads
// a location and a scale, and the uniform family reads its two ends.
const LAW_NORMAL: i64 = 0;
const LAW_UNIFORM: i64 = 1;
const LAW_LOGISTIC: i64 = 2;
const LAW_CAUCHY: i64 = 3;
const LAW_LAPLACE: i64 = 4;
const LAW_RAYLEIGH: i64 = 5;
const LAW_GUMBEL_R: i64 = 6;
const LAW_GUMBEL_L: i64 = 7;

// The families of a probability plot. The Tukey lambda family and the
// Weibull family read a shape.
const PLOT_NORMAL: i64 = 0;
const PLOT_EXPON: i64 = 1;
const PLOT_UNIFORM: i64 = 2;
const PLOT_TUKEY: i64 = 3;
const PLOT_LOGISTIC: i64 = 4;
const PLOT_LAPLACE: i64 = 5;
const PLOT_CAUCHY: i64 = 6;
const PLOT_GUMBEL_R: i64 = 7;
const PLOT_GUMBEL_L: i64 = 8;
const PLOT_RAYLEIGH: i64 = 9;
const PLOT_WEIBULL: i64 = 10;

// The methods of the Siegel slopes.
const SIEGEL_HIERARCHICAL: i64 = 0;
const SIEGEL_SEPARATE: i64 = 1;

// The families of the Anderson-Darling test.
const ANDERSON_EXPON: i64 = 1;
const ANDERSON_LOGISTIC: i64 = 2;
const ANDERSON_GUMBEL: i64 = 3;
const ANDERSON_GUMBEL_L: i64 = 4;
const ANDERSON_GUMBEL_R: i64 = 5;
const ANDERSON_EXTREME1: i64 = 6;
const ANDERSON_WEIBULL: i64 = 7;

// The rank rule of the weighted tau.

// ---- structs of the association, distance and distribution calls -----
// The statistic and the p-value of a correlation or a dependence test.
struct CorrelationResult {
    statistic: f64,
    pvalue: f64,
}

// A correlation result with a fault code. Every field is NaN.
fn corr_fault(code: i64) -> CorrelationResult {
    return CorrelationResult { statistic: nan(), pvalue: nan() };
}

// The estimate and the confidence interval of a risk measure.
struct IntervalResult {
    statistic: f64,
    low: f64,
    high: f64,
}

// An interval result with a fault code. Every field is NaN.
fn interval_fault(code: i64) -> IntervalResult {
    return IntervalResult { statistic: nan(), low: nan(), high: nan() };
}

// A correlation result with a statistic and a p-value.
fn corr_good(statistic: f64, pvalue: f64, df: f64) -> CorrelationResult {
    return CorrelationResult { statistic: statistic, pvalue: pvalue };
}

// The p-value of a correlation of n points. The p-value comes from the
// beta law of the distribution of the correlation.
fn corr_beta_p(r: f64, n: i64, alt: i64) -> f64, i64 {
    let ab = n / 2.0 - 1.0;
    var x: f64 = 0.0;
    if alt == ALT_TWO_SIDED { x = 0.5 * (1.0 - abs(r)); }
    if alt == ALT_LESS { x = 0.5 * (1.0 + r); }
    if alt == ALT_GREATER { x = 0.5 * (1.0 - r); }
    if x < 0.0 { x = 0.0; }
    if x > 1.0 { x = 1.0; }
    let q, cb = betainc(ab, ab, x);
    if cb != OK { return nan(), cb; }
    var p = q;
    if alt == ALT_TWO_SIDED { p = 2.0 * q; }
    return p, OK;
}

// The result of a probability plot: the line of the ordered data on the
// ordered medians.
struct ProbabilityPlot {
    slope: f64,
    intercept: f64,
    r: f64,
}

// The bounds of the last pass of a sigma clip, and the kept count.
struct ClipResult {
    low: f64,
    high: f64,
    kept: i64,
}

