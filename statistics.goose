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

// ---- bounds ---------------------------------------------------------
const MAX_INPUT: i64 = 4194304;  // values in a sample
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
