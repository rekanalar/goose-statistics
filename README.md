# goose-statistics

[![ci](https://github.com/rekanalar/goose-statistics/actions/workflows/ci.yml/badge.svg)](https://github.com/rekanalar/goose-statistics/actions/workflows/ci.yml)

Descriptive statistics for the [Goose](https://github.com/aardappel/goose)
programming language.

The package gives the location, the spread, the quantiles and the shape of
a sample of 64 bit floats. It gives the trimmed forms and the weighted
forms too. The sum is exact, so the mean and the variance keep their
digits.

Every call returns the result first and a fault code second. A fault gives
NaN as the result. A call never changes the data of the caller.

Pure Goose, no C.

## Install

Add this to the root file of the program that uses it:

```goose
import statistics from "rekanalar/goose-statistics@v0.2.0";
```

## API

| Function | What it does |
|---|---|
| `sum(xs)` | the exact sum of the sample, rounded once |
| `mean(xs)` | the mean of the sample |
| `fmean(xs)` | the same as `mean`; the name is for a reader who knows Python |
| `geometric_mean(xs)` | the geometric mean; a zero gives 0 and a negative value gives `ERR_DOMAIN` |
| `harmonic_mean(xs)` | the harmonic mean; a zero gives 0 and a negative value gives `ERR_DOMAIN` |
| `minimum(xs)` | the smallest value of the sample |
| `maximum(xs)` | the largest value of the sample |
| `ptp(xs)` | the largest value minus the smallest value |
| `median(xs)` | the median of the sample |
| `median_low(xs)` | the lower middle value of an even sample |
| `median_high(xs)` | the upper middle value of an even sample |
| `median_grouped(xs, width)` | the grouped median; `width` above 0, or `ERR_DOMAIN` |
| `mode(xs)` | the most common value; a tie gives the first value of the sample |
| `multimode(out, xs)` | writes every most common value into `out`; gives the count and the code |
| `pvariance(xs)` | the population variance |
| `variance(xs)` | the sample variance |
| `pstdev(xs)` | the population standard deviation |
| `stdev(xs)` | the sample standard deviation |
| `sem(xs)` | the standard error of the mean |
| `variation(xs)` | the coefficient of variation; a mean of 0 gives `ERR_UNDEFINED` |
| `iqr(xs, method)` | the interquartile range |
| `quantile(xs, q, method)` | the quantile at `q`; `q` below 0 or above 1 gives `ERR_DOMAIN` |
| `cut_points(out, xs, groups, method)` | writes `groups - 1` cut points into `out`; gives the count and the code |
| `moment(xs, order)` | the raw moment; `order` from 1 to `MAX_ORDER`, or `ERR_DOMAIN` |
| `central_moment(xs, order)` | the central moment; the same range of order |
| `skew(xs, bias)` | the skewness; a constant sample gives `ERR_UNDEFINED` |
| `kurtosis(xs, fisher, bias)` | the kurtosis; a constant sample gives `ERR_UNDEFINED` |
| `zscore(out, xs, ddof)` | writes each z-score into `out`; gives the count and the code |
| `tmean(xs, lower, upper)` | the mean of the values from `lower` to `upper`, both included |
| `trim_mean(xs, proportion)` | the mean after the cut; `proportion` at or above 0.5 gives `ERR_DOMAIN` |
| `weighted_mean(xs, ws)` | the weighted mean |
| `weighted_variance(xs, ws, form)` | the weighted variance in one of the three forms |
| `weighted_median(xs, ws)` | the weighted median |
| `pcovariance(xs, ys)` | the population covariance of two samples |
| `covariance(xs, ys)` | the sample covariance of two samples |
| `correlation(xs, ys)` | the Pearson correlation of two samples |
| `weighted_covariance(xs, ys, ws, form)` | the weighted covariance of two samples in one of the three forms |
| `weighted_correlation(xs, ys, ws)` | the weighted correlation of two samples |
| `rankdata(out, xs, method)` | writes each rank into `out`; gives the count and the code |
| `spearman(xs, ys)` | the Spearman rank correlation of two samples |
| `kendall(xs, ys, variant)` | the Kendall rank correlation, tau-b or tau-c |
| `linear_regression(xs, ys, proportional)` | the least squares line of two samples; gives the line and the code |
| `regression_errors(xs, ys)` | the standard error of the slope and the standard error of the intercept |
| `theil_sen(xs, ys)` | the Theil-Sen slope and the Theil-Sen intercept |
| `covariance_matrix(out, data, ncols, ddof)` | writes the covariance matrix into `out`; gives the count and the code |
| `correlation_matrix(out, data, ncols)` | writes the correlation matrix into `out`; gives the count and the code |

The argument `method` is one of the nine `METHOD_` constants. The
constants `METHOD_EXCLUSIVE` and `METHOD_INCLUSIVE` stand for the two
methods of Python.

The argument `form` is one of `WEIGHT_POPULATION`, `WEIGHT_FREQUENCY` and
`WEIGHT_RELIABILITY`.

The argument `method` of `rankdata` is one of the five `TIE_` constants.
The constants `TIE_AVERAGE`, `TIE_MIN`, `TIE_MAX`, `TIE_DENSE` and
`TIE_ORDINAL` name the five tie methods of SciPy.

The argument `variant` of `kendall` is one of the two `KENDALL_`
constants. The constant `KENDALL_B` gives tau-b and the constant
`KENDALL_C` gives tau-c.

The argument `bias` and the argument `fisher` are booleans. The argument
`ddof` is 0 or 1, and any other value gives `ERR_DOMAIN`.

The call `linear_regression` gives a line. The line holds the fields
`slope`, `intercept`, `r` and `r2`. The call `regression_errors` gives the
standard error of the slope first and the standard error of the intercept
second. The call `theil_sen` gives the slope first and the intercept
second.

The buffer of a matrix call holds the columns one after the other. The
result buffer holds the cells row by row. The argument `ncols` gives the
count of the columns. The argument `ddof` of `covariance_matrix` is 0 for
the population form and 1 for the sample form.

The constant `MAX_THEIL` of 4,096 points gives the largest sample of the
Theil-Sen call. The constant `MAX_COLUMNS` of 1,024 gives the largest
count of the columns of a matrix call. The constant `MAX_WORK` of
268,435,456 gives the largest work of a matrix call.

Every call returns the result first and the fault code second. A fault
gives NaN as the result and one of these codes.

| Code | Meaning |
|---|---|
| `OK` | the call succeeded |
| `ERR_EMPTY` | the sample has no value |
| `ERR_TOO_SHORT` | the sample has fewer values than the call needs |
| `ERR_NAN` | the sample, a weight or an argument holds a NaN |
| `ERR_INFINITY` | the sample or a weight holds an infinity |
| `ERR_DOMAIN` | an argument or a value is out of range for the call |
| `ERR_LENGTH` | the list of weights has another length than the sample |
| `ERR_UNDEFINED` | the result has no defined value |
| `ERR_OVERFLOW` | the true result does not fit in a 64 bit float |
| `ERR_TOO_LONG` | the sample has more than `MAX_INPUT` values |

The order of the checks never changes. The order is the length bound, the
empty sample, a NaN and an infinity. The arguments and the length of the
sample come last. The same input gives the same code on every run.

```goose
let xs: f64[>..] = [1.0, 2.0, 3.0, 4.0];
let m, code = statistics::mean(xs[..]);
if code == statistics::OK {
    print(m);   // 2.5
}
```

## Numerics

The sum of the package is the exact sum of Shewchuk. It keeps a list of
partial sums, and each addition uses the two sum of Knuth. The result is
rounded once, at the end. The error of the sum stays below one unit in the
last place.

The variance uses the corrected two pass method. The first pass gets the
mean with the exact sum. The second pass sums the squares of the
deviations. The result subtracts the correction term of the rounding of
the mean. A large offset therefore keeps its spread.

The quantiles use the nine methods of NumPy. The method is a required
argument, because NumPy uses `linear` and Python uses `exclusive`. A
silent default would give a different answer to half of the callers.

The calls that take two columns scale each column by its own power of
two. The largest value of a column becomes at most 2 to the power 400.
It becomes at least 2 to the power -300. The products and the sums then
stay finite and stay normal. The scale of a column is an exponent, so
the result comes back without a false overflow.

## Limits

- The package holds no p-value and no confidence interval. SciPy holds
  `pearsonr`, `spearmanr` and `kendalltau` with these values.
- The package holds no confidence slope of the Theil-Sen estimate. SciPy
  holds `theilslopes` with the low and the high slope.
- The package holds no partial correlation and no partial rank
  correlation. SciPy holds `partialcorr`.
- The package holds no complex sample. NumPy holds the complex forms.
- The package holds no significance test for the rank correlations.
- The package holds no special function. SciPy holds `special`. The crate
  `statrs` of Rust holds `function`.
- The package holds no distribution. Python holds `NormalDist`. SciPy holds
  the distributions of `stats`. gonum holds `distuv`.
- The package holds no hypothesis test. SciPy holds `ttest_ind`,
  `chisquare` and `mannwhitneyu`.
- The package holds no resampling call. SciPy holds `bootstrap` and
  `permutation_test`.
- The package holds no histogram and no density estimate. NumPy holds
  `histogram`. SciPy holds `gaussian_kde`.
- The package holds no entropy. gonum holds `Entropy`. SciPy holds
  `entropy`.
- The package holds no 32 bit float and no complex value. NumPy holds
  both.
- The package holds no streaming form. Every call reads the whole sample.
- The weights exist only in the three named forms. The quantile methods
  other than `inverted_cdf` have no weighted form. NumPy refuses them too.
- The variation is the sample form. The reference libraries hold no
  coefficient of variation.
- The geometric mean and the harmonic mean use `log` and `exp` of the
  standard library. Their error grows near the limits of the type. The
  tests allow 8 units in the last place.
- The exact sum scales a sample by a power of two when a value reaches
  2 to the power 1000. The scaling drops the bits of a value below 2 to
  the power -974. The error stays below one unit of the largest value.
- The tolerance of the tests is relative. A value near zero has no
  absolute tolerance.
- The Theil-Sen call reaches its bound at 4,096 points. It then gives
  `ERR_TOO_LONG`. A larger sample gives no slope and no intercept.
- A matrix call takes at most 1,024 columns. The work of the call is the
  count of the cells times the count of the rows. It stays at most
  `MAX_WORK`. A larger work gives `ERR_TOO_LONG` and writes no cell.
- The intercept of a line with a large offset keeps the digits of the
  exact sums. SciPy rounds the slope first, so its intercept differs from
  the exact value.
- A slope below 2 to the power -1074 fits in no float. Two columns whose
  magnitudes differ by more than about 300 orders of magnitude give such
  a slope. The call gives 0 for it and no fault.
- The weighted forms divide the weights by the largest weight. The
  weights of a sample that reach 1e150 lose a few units in the last
  place.

## Security

The package reads a sample from people we do not know. These protections
go beyond the reference libraries.

- Every call checks the length of its input before any other work. The
  limit is `MAX_INPUT` of 4,194,304 values. A longer sample gives
  `ERR_TOO_LONG` and no work starts.
- Every call rejects a NaN and an infinity in the sample and in the
  weights. The reference libraries let a NaN reach the result, and some
  of them let an infinity reach it too.
- Every call rejects a negative weight. NumPy and gonum accept a negative
  weight and give a result that no reader expects.
- Every call checks its result for an infinity and for a NaN. A result
  that does not fit gives `ERR_OVERFLOW`. The reference libraries return
  an infinity without a warning in many cases.
- A call that sorts copies the sample first. The copy holds at most
  `MAX_INPUT` values, which is 32 MiB. The data of the caller never
  changes.
- The sort comes from the standard library. It is a pattern defeating
  quicksort with a heapsort fallback, so its worst case is n log n. The
  test suite times the sort on five crafted samples. It fails when a
  crafted sample takes more than five times as long as a random one.
- A fault holds a code and NaN. It never holds a value of the sample.
- The exact sum keeps a growable list of partial sums. The list holds at
  most 2098 entries, one for each bit position of a double. A guard of
  2200 entries gives `ERR_UNDEFINED` and writes nothing past the list.
- A call that takes two columns checks that both have the same length.
  A different length gives `ERR_LENGTH` and no work starts. The
  reference libraries raise or give a wrong answer.
- A call that takes two columns rejects a NaN and an infinity in both
  columns. The reference libraries let them reach the result.
- The correlation stays between -1 and 1. The clamp holds a result that
  the rounding pushes past the range.
- A call that takes two columns scales each column on its own. The
  reference libraries use one scale for both columns. A small column
  next to a huge column falls to zero there. The call gives a wrong
  result with no warning. This package scales each column, so the value
  of the small column survives.
- The Theil-Sen call stops at 4,096 points. It then gives
  `ERR_TOO_LONG` and computes no slope. The call of SciPy has no such
  bound and uses the memory of the caller without a limit.
- A matrix call checks the work before it starts. A work above
  `MAX_WORK` gives `ERR_TOO_LONG` and writes no cell. The count of the
  cells stays at most 1,048,576.
- The Kendall count sorts the pairs and counts the swaps. It runs in the
  time of a sort. The test suite times the count on five crafted samples
  of the maximum length.

These risks stay.

- The package does not check the meaning of the values. A value that
  holds a price, a weight or a count is a number. The package leaves the
  choice to the caller.
- The package does not check the units. A sample of lengths and a sample
  of weights give the same mean.
- The tests compare the values with the reference libraries within a
  tolerance in units in the last place. A tolerance of 8 units in the
  last place is not an exact match.

## Test

```bash
goose run test/run.goose
goose run examples/describe.goose
goose run examples/quantiles.goose
goose run examples/hostile.goose
goose run examples/association.goose
goose run examples/ranks.goose
```

The scripts `test/make_vectors.py` and `test/make_assoc_vectors.py` build
the vectors of the test suite. They write numbers only into the package. Each vector is compared with a second
source, and the script fails on a difference beyond the tolerance.

The second sources are in four languages. Python gives the exact value
through the modules `fractions` and `decimal`. NumPy, SciPy and Python
`statistics` give the reference values. Three more libraries give their
own opinion. They are gonum of Go, the crate `statrs` of Rust and the
package `simple-statistics` of Node.

The oracle programs are under `test/`. The program `test/oracle` reads a
list of samples and gives the values of Go with gonum. The program
`test/oracle-rust` gives the values of the crate `statrs`. The program
`test/oracle-node` gives the values of the package `simple-statistics`.
Each program needs its own toolchain and its own packages.

The script `test/mutate.sh` flips one rule of the code at a time in a copy
of the package. Each flip must make the test program fail. The last line
reads `N probes, 0 failures`.

## License

MIT
