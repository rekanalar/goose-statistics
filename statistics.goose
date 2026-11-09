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
