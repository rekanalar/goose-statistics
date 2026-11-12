# goose-statistics

[![ci](https://github.com/rekanalar/goose-statistics/actions/workflows/ci.yml/badge.svg)](https://github.com/rekanalar/goose-statistics/actions/workflows/ci.yml)

Descriptive statistics for the [Goose](https://github.com/aardappel/goose)
programming language.

The package gives the location, the spread, the quantiles and the shape of
a sample of 64 bit floats. It gives the trimmed forms and the weighted
forms too. The sum is exact, so the mean and the variance keep their
digits.

The package gives fifteen probability families too. Ten are continuous
and five are discrete. Each family has its density or mass, its
cumulative function, its survival function, its quantile and its four
moments.

Every call returns the result first and a fault code second. A fault gives
NaN as the result. A call never changes the data of the caller.

Pure Goose, no C.

## Install

Add this to the root file of the program that uses it:

```goose
import statistics from "rekanalar/goose-statistics@v0.6.0";
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
| `erf(x)`, `erfc(x)` | the error function and its complement; `erfc` keeps its digits far in the tail |
| `erfinv(y)`, `erfcinv(p)` | the inverses; `erfinv` takes -1 to 1 and `erfcinv` takes 0 to 2 |
| `gamma(x)` | the gamma function; a pole gives `ERR_UNDEFINED` and a result that does not fit gives `ERR_OVERFLOW` |
| `lgamma(x)` | gives the log of the absolute gamma, the sign and the code |
| `digamma(x)` | the logarithmic derivative of the gamma function |
| `beta(a, b)`, `lbeta(a, b)` | the beta function and its log, for `a` and `b` above 0 |
| `lchoose(n, k)` | the log of the binomial coefficient, for `k` from 0 to `n` |
| `gammainc(a, x)`, `gammaincc(a, x)` | the regularized lower and upper incomplete gamma, for `a` up to `MAX_SHAPE` |
| `gammaincinv(a, p)`, `gammainccinv(a, q)` | the inverses of the two incomplete gamma functions |
| `betainc(a, b, x)` | the regularized incomplete beta, for `x` from 0 to 1 |
| `betaincinv(a, b, p)` | the inverse of the regularized incomplete beta |
| `log1p(x)`, `expm1(x)` | the log of 1 + x and the exponential minus 1, exact near 0 |
| `norm_pdf(x, mu, sigma)`, `norm_cdf`, `norm_sf`, `norm_ppf` | the normal family; `mu` = 0 and `sigma` = 1 |
| `norm_mean(mu, sigma)`, `norm_var`, `norm_skew`, `norm_kurt` | the moments of the normal family |
| `unif_pdf(x, lo, hi)`, and the same seven calls | the uniform family; `lo` = 0 and `hi` = 1 |
| `expon_pdf(x, rate)`, and the same seven calls | the exponential family; `rate` = 1 |
| `lognorm_pdf(x, mu, sigma)`, and the same seven calls | the lognormal family; `mu` = 0 and `sigma` = 1 |
| `weibull_pdf(x, shape, scale)`, and the same seven calls | the Weibull family; `shape` = 1 and `scale` = 1 |
| `gammadist_pdf(x, shape, scale)`, and the same seven calls | the gamma family; `shape` = 1 and `scale` = 1 |
| `betadist_pdf(x, a, b)`, and the same seven calls | the beta family; `a` = 1 and `b` = 1 |
| `chi2_pdf(x, df)`, and the same seven calls | the chi-square family; `df` = 1 |
| `tdist_pdf(x, df)`, and the same seven calls | the Student t family; `df` = 1 |
| `fdist_pdf(x, df1, df2)`, and the same seven calls | the F family; `df1` = 1 and `df2` = 1 |
| `binom_pmf(k, n, p)`, and the same seven calls | the binomial family; `p` = 0.5 |
| `poisson_pmf(k, rate)`, and the same seven calls | the Poisson family; `rate` = 1 |
| `geom_pmf(k, p)`, and the same seven calls | the geometric family in the trials form; `p` = 0.5 |
| `nbinom_pmf(k, r, p)`, and the same seven calls | the negative binomial family; `r` = 1 and `p` = 0.5 |
| `hypergeom_pmf(k, M, n, N)`, and the same seven calls | the hypergeometric family, no default |
| `ttest_1samp(xs, popmean, alt)` | the one-sample t test; gives a `TestResult` |
| `ttest_ind(xs, ys, alt, equal_var)` | the two-sample t test; `equal_var` false gives the Welch form |
| `ttest_rel(xs, ys, alt)` | the paired t test |
| `f_oneway(groups, sizes)` | the one-way analysis of variance of the groups in one flat buffer |
| `mean_confidence_interval(xs, confidence)` | the t interval of the mean; gives the lower bound, the upper bound and the code |
| `chisquare(obs, expected, ddof)` | the chi-square test of a fit; an empty `expected` gives the uniform fit; a zero expected count gives `ERR_UNDEFINED`, or `ERR_DOMAIN` when the sums differ |
| `chi2_contingency(table, rows, cols, correction)` | the chi-square test of a table; `correction` applies the Yates form to a 2x2 table |
| `binomtest(k, n, p, alt)` | the exact binomial test |
| `fisher_exact(table, alt)` | the exact test of a 2x2 table; the statistic is the odds ratio |
| `mannwhitneyu(xs, ys, continuity, alt, method)` | the Mann-Whitney U test |
| `wilcoxon(diffs, alt, zero_method, correction)` | the Wilcoxon signed-rank test of the differences |
| `wilcoxon_paired(xs, ys, alt, zero_method, correction)` | the same test of two paired samples |
| `kruskal(groups, sizes)` | the Kruskal-Wallis test with the correction for ties |
| `ks_1samp(xs, family, params, alt)` | the Kolmogorov-Smirnov test against one of the ten families |
| `ks_2samp(xs, ys, alt)` | the two-sample Kolmogorov-Smirnov test |
| `shapiro(xs)` | the Shapiro-Wilk test of normality; it needs 3 values |
| `skewtest(xs, alt)` | the test of the skewness against the normal shape; it needs 8 values |
| `kurtosistest(xs, alt)` | the test of the kurtosis against the normal shape; it needs 5 values |
| `normaltest(xs)` | the omnibus K-squared test of D'Agostino and Pearson; it needs 8 values |
| `jarque_bera(xs)` | the test of Jarque and Bera; it needs 2 values |
| `anderson_darling(xs)` | the Anderson-Darling test against the fitted normal; it needs 5 values |
| `pearsonr(xs, ys, alt)` | the Pearson correlation with its p-value |
| `spearmanr(xs, ys, alt)` | the Spearman correlation with its p-value |
| `kendalltau(xs, ys, alt)` | the Kendall tau-b with its p-value |
| `multipletests(out, ps, method)` | writes the adjusted p-values into `out`; gives the count and the code |

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

The special functions take one to three numbers and no sample. The call
`lgamma` returns three values: the log of the absolute gamma, the sign
and the code. A pole gives `ERR_UNDEFINED`. A NaN gives `ERR_NAN`.

Each family name stands for eight calls. The endings are `pdf` or `pmf`,
`cdf`, `sf`, `ppf`, `mean`, `var`, `skew` and `kurt`. The names `gamma`
and `beta` belong to the special functions. The family calls carry the
prefixes `gammadist_` and `betadist_`. A continuous family takes a
density, and a discrete family takes a mass. The call `kurt` gives the
excess kurtosis, as gonum gives it.

A test returns a `TestResult` and then the fault code. The fields are
`statistic`, `pvalue`, `df` and `df2`. A field that the test does not have
holds NaN. A fault gives NaN in every field. A p-value with the code `OK`
lies from 0 to 1 and is never NaN.

The argument `alt` is one of `ALT_TWO_SIDED`, `ALT_LESS` and `ALT_GREATER`.
Another value gives `ERR_DOMAIN`. The argument `method` of `mannwhitneyu`
is one of `MW_AUTO`, `MW_EXACT` and `MW_ASYMPTOTIC`. The argument
`zero_method` is one of `ZERO_WILCOX`, `ZERO_PRATT` and `ZERO_ZSPLIT`. The
argument `method` of `multipletests` is one of `ADJUST_BONFERRONI`,
`ADJUST_HOLM` and `ADJUST_BH`.

A group call takes one flat buffer of all the values and an array of the
group sizes. A table call takes the counts row by row. The argument
`family` of `ks_1samp` is one of the ten `FAMILY_` constants, from
`FAMILY_NORMAL` to `FAMILY_F`. The buffer `params` holds the parameters in
the order of the family call. Most families take two. The exponential,
the chi-square and the Student t family take one.

The exact forms follow the rules of SciPy. The Mann-Whitney test counts the
exact distribution when one sample holds at most 8 values and no value
repeats. The Wilcoxon test counts exactly for at most 50 pairs without
ties and zeros. It also counts exactly for at most 13 pairs with them. The Kendall test is
exact for at most 33 pairs without ties. The two-sample Kolmogorov-Smirnov
test is exact when the product of the sizes is at most 10,000 and no value
repeats.

The other cases use the normal form with the correction for ties.

The quantile calls take the probability first.

A discrete quantile gives the smallest count whose cumulative function
reaches the probability. The probability 0 gives the smallest count of
the support. The probability 1 gives the largest count of the support.

The support of the geometric family starts at 1, because the call counts
the trials. The negative binomial family counts the failures before the
r-th success. The support of the hypergeometric family holds at most
`MAX_DIST_SUPPORT` of 100,000 counts. A larger support gives
`ERR_TOO_LONG`. The mass of a count outside the support is 0.

An infinity gives `ERR_INFINITY`, except for `erf` and `erfc`, where it is
a valid limit. A shape or a probability out of range gives `ERR_DOMAIN`.
A result that does not fit in a float gives `ERR_OVERFLOW`. The constant
`MAX_SHAPE` of 1e10 gives the largest shape of an incomplete function.

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

The special functions use pairs of doubles where the plain form loses
digits. The log of the beta function subtracts log gamma terms of equal
size. The incomplete functions form the small tail directly and give the
large tail as one minus the small tail. So neither tail loses its digits.

The tests allow 4 units in the last place for the error functions and 8
for the gamma family. They allow 16 for the incomplete functions up to a
shape of 1e4. They allow 64 up to a shape of 1e10.

The forward call of an inverse returns the input within 20 units. This
holds for a shape up to 10. For a larger shape the root itself stays
within 64 units.

The distributions use the same pairs. The normal cumulative function is
the tail itself, and so is the survival function. Neither loses the
digits of a far tail. The tails of the Student t and of the F family
take a first order correction on the ratio argument. The discrete
quantiles run a bounded search. A search that meets its bound gives
`ERR_OVERFLOW`, never a silent wrong count.

The tests allow 4 units for the density of the normal family. They allow
8 for the exponential and the lognormal forms. They allow 16 for the
gamma, the beta and the tailed families, and for every moment. The
counts of the discrete quantiles are exact.

## Limits

- The package holds no confidence interval of a correlation. SciPy holds
  `pearsonr`, `spearmanr` and `kendalltau` with these values.
- The package holds no confidence slope of the Theil-Sen estimate. SciPy
  holds `theilslopes` with the low and the high slope.
- The package holds no partial correlation and no partial rank
  correlation. SciPy holds `partialcorr`.
- The package holds no complex sample. NumPy holds the complex forms.
- The package holds no significance test for the rank correlations.
- The special functions hold no unregularized incomplete gamma. The crate
  `statrs` holds `gamma_li` and `gamma_ui`.
- The package holds no integer factorial and no exact binomial
  coefficient. Python holds `math.factorial` and `math.comb`.
- The package holds no complex gamma, no Bessel function, no zeta
  function and no elliptic integral. SciPy holds them in `special`.
- The package holds no `ndtr` and no `ndtri`. SciPy holds both.
- The package holds no inverse of the digamma. The crate `statrs` holds
  `inv_digamma`.
- A shape above `MAX_SHAPE` of 1e10 gives `ERR_DOMAIN` in the incomplete
  functions and their inverses. SciPy gives a number there.
- A tail below the smallest subnormal float gives 0 and no fault. The
  call `erfc` gives 0 from 27.25. The call `beta` gives 0 when the true
  value is below the smallest subnormal.
- The log gamma of an argument near a root keeps a relative error only
  away from the zeros. Near 1 and 2 the error is absolute, near 1e-15.
- An inverse of the incomplete functions loses digits when the input
  probability is near 1. Pass the small tail instead, with `gammainccinv`
  or the complement of the shape.
- A tail below the smallest subnormal gives 0 and no fault. This happens
  for the survival function of the normal family from 38.6 standard
  deviations. It happens for the mass of a Poisson count far from the
  rate. The call `beta` behaves the same.
- The survival function of the reference libraries loses the digits of a
  small tail. At 6 standard deviations the value of SciPy differs from
  the exact tail by 8 units in the last place. The value of gonum
  differs by 630 units. The value of the crate `statrs` differs by 242
  units. This package gives the exact tail there.
- SciPy and gonum return NaN or an infinity for a domain error, without a
  signal. The crate `statrs` panics for some arguments. This package
  gives a fault code.
- The kurtosis call of SciPy gives the kurtosis. This call `kurt` gives
  the excess kurtosis, which is 3 less.
- gonum holds no geometric, negative binomial or hypergeometric family.
  It holds no quantile of the binomial and the Poisson families. The
  package `jstat` holds no geometric family and no skewness or kurtosis.
- The mean of the Student t family needs `df` above 1. The variance
  needs it above 2, the skewness above 3 and the kurtosis above 4. The
  moments of the F family need `df2` above 2, 4, 6 and 8 in the same
  order. A pole gives `ERR_UNDEFINED`.
- The quantile of a discrete family is the smallest count whose
  cumulative function reaches the probability. The references agree on
  this rule.
- The package holds no resampling test. SciPy holds `permutation_test`
  and the random permutation form of `wilcoxon`. The Wilcoxon test counts
  all sign patterns for at most 13 pairs and uses the normal form above.
- The Shapiro-Wilk p-value is not accurate
  above 5,000 values, as in SciPy.
- The normality tests beyond Shapiro-Wilk follow SciPy. The minimum
  sizes are the sizes of SciPy:
  `skewtest` needs 8 values.
  `kurtosistest` needs 5.
  `normaltest` needs 8.
  `jarque_bera` needs 2.
  `anderson_darling` needs 5.

  SciPy returns NaN below its sizes for the first three, and it has no
  floor for `anderson`. This package faults below the sizes instead. The
  asymptotic p-value rule holds for 5 values or more.
- The Anderson-Darling p-value rule is a reference disagreement. SciPy
  gives no p-value by default and warns. Its `method=interpolate`
  reads a table of five points and floors the p-value at 0.15. This
  package follows the asymptotic rule of D'Agostino and Stephens, the
  rule of Rust `statrs`.

  On the 16 values of the acceptance sample, SciPy gives 0.15.
  This package gives 0.7595.

  The statistic follows the SciPy formula with the lower and the
  upper tail logarithms.
  statrs takes the logarithm of one minus the cumulative function.
  It loses digits when a tail is small.

  The fit of this package uses exact sums.
  It differs from the float fit of SciPy by about 1e-11.
  The cause is that the statistic cancels against the sample size.

  R `nortest::ad.test`
  documents the same asymptotic family of rules. R was not run; the R
  rules come from its documentation.
- The smallest size of `jarque_bera` is a reference disagreement too.
  SciPy accepts 2 values. The R documentation of
  `tseries::jarque.bera.test` asks for more. This package follows
  SciPy.
- The multiple-testing call holds three methods. R `p.adjust` also holds
  Hochberg, Hommel and Benjamini-Yekutieli.
- The Mann-Whitney call limits each of its two samples to `MAX_INPUT`
  values. The call copies both samples into one pool, so the pool can
  hold twice the limit. A sample above the limit gives `ERR_TOO_LONG`.
- The Mann-Whitney table has at most 100,000 cells. A larger exact case
  uses the normal form. SciPy keeps the exact form.
- The two-sample Kolmogorov-Smirnov test uses the exact form for a product
  of the sizes up to 10,000. SciPy attempts it for sizes up to 10,000
  each, so it differs for large samples. A tie in either sample gives the
  normal form, because the exact tables assume distinct values. SciPy
  computes the exact value on tied data.
- The one-sample Kolmogorov-Smirnov test gives the exact p-value for any
  sample, tied or not, up to 10,000 values. It uses the series above that.
- `chi2_contingency` accepts a count with a fraction, as SciPy does.
  `fisher_exact` rejects it.
- The Wilcoxon test of a sample of all zeros gives NaN with
  `ERR_UNDEFINED`. SciPy gives a p-value of 1. The statistic of such a
  sample has no meaning, so this package refuses it.
- R and SciPy differ on a constant sample. R stops with an error and SciPy
  gives NaN. This package gives NaN and `ERR_UNDEFINED`.
- The package returns the statistic k/n for `binomtest`, as SciPy does.
- The `simple-statistics` function `tTest` divides by the population
  deviation, so its statistic differs from SciPy by the factor
  `sqrt(n / (n - 1))`. This package follows SciPy and R.
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
- The review of this release did not check some parts of the code.
  It did not check the coefficients of the Shapiro-Wilk test. It did
  not check the bounds of the series of the exact Kolmogorov-Smirnov
  test. It did not check the asymptotic branches of the one-sample
  Kolmogorov-Smirnov test. A wrong value in those places stays.

## Security

The package reads a sample from people we do not know. These protections
go beyond the reference libraries.

- Every call checks the length of its input before any other work. The
  limit is `MAX_INPUT` of 4,194,304 values. A longer sample gives
  `ERR_TOO_LONG` and no work starts.
- A test that reads a sample checks the length and the empty sample
  first. It checks NaN and infinity before it reads a table or a group.
  A test that reads counts checks each count for NaN, for an infinity
  and for a negative value.
  Every exact count has a bound: 8 values, 50 pairs, 13 pairs, 33 pairs
  and a product of 10,000.
  The Mann-Whitney table stops at 100,000 cells. A p-value with the code
  `OK` is a number from 0 to 1, never NaN.
- An exact count stays below 2 to the power 53, or it is a probability.
  An exact value never overflows in silence.
- `chi2_contingency` and `fisher_exact` reject a negative count and a
  count above 2 to the power 53. `fisher_exact` rejects a fraction too.
- `multipletests` rejects a p-value outside 0 to 1 and clamps each
  adjusted value to 1.
- `mean_confidence_interval` rejects a confidence outside 0 to 1.
- The normality tests fault below their minimum sizes. SciPy returns
  NaN for those inputs. A caller who ignores codes sees no NaN here.
- The p-value rule of `anderson_darling` holds only to the turn of its
  top quadratic, at an adjusted statistic of 153.4. Above the turn,
  the rule of statrs and of R `nortest` rises above 1.
  A caller reads a reversed verdict. This call gives 0 there.
- `anderson_darling` rejects a fitted scale of 0. A value past 30
  standard units of the fit has a tail near the limit of a float.
  The statistic takes the logarithm of the scaled tail there.
  The result stays finite and keeps the digits of the tail.
- statrs takes one minus the cumulative function in the
  Anderson-Darling test. Its statistic is infinite there.
  SciPy takes the logarithm of the tail and stays exact.
- Every call rejects a NaN and an infinity in the sample and in the
  weights. The reference libraries let a NaN reach the result, and some
  of them let an infinity reach it too.
- Every call rejects a negative weight. NumPy and gonum accept a negative
  weight and give a result that no reader expects.
- Every call checks its result for an infinity and for a NaN. A result
  that does not fit gives `ERR_OVERFLOW`. The reference libraries return
  an infinity without a warning in many cases.
- A call that sorts copies the sample first. The copy of one sample
  holds at most `MAX_INPUT` values, which is 32 MiB. The calls that
  read two samples copy both, so their copy can reach 64 MiB. The data
  of the caller never changes.
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
  next to a column of large values falls to zero there. The call gives a wrong
  result with no warning. This package scales each column, so the value
  of the small column survives.
- The Theil-Sen call stops at 4,096 points. It then gives
  `ERR_TOO_LONG` and computes no slope. The call of SciPy has no such
  bound and uses the memory of the caller without a limit.
- A matrix call checks the work before it starts. A work above
  `MAX_WORK` gives `ERR_TOO_LONG` and writes no cell. The count of the
  cells stays at most 1,048,576.
- The special functions return a fault code and NaN for a NaN, a domain
  error, a pole and an overflow. SciPy and Go return NaN or an infinity
  without a signal. gonum and the crate `statrs` panic on a shape out of
  range. Python raises an exception.
- Every iteration of the special functions has a fixed bound. A series
  has 2048 terms, a fraction has 40,000 steps and an inverse has 100
  steps. A loop that meets its bound gives `ERR_OVERFLOW`. It never gives
  a wrong value without a signal.
- A shape above 1e10 gives `ERR_DOMAIN`. The bounds above hold for every
  shape up to that value.
- No special function allocates memory that grows with its argument, and
  none keeps state between calls.
- The test suite runs 100,000 seeded inputs on every special function.
  It times crafted shapes against a random call.
- The Kendall count sorts the pairs and counts the swaps. It runs in the
  time of a sort. The test suite times the count on five crafted samples
  of the maximum length.
- The support of the hypergeometric family is bounded. A support above
  `MAX_DIST_SUPPORT` gives `ERR_TOO_LONG` and no sum starts. The sum
  holds pairs, so a mass far below the largest mass keeps its digits.
- The search of a discrete quantile has two fixed bounds. The doubling
  of an unbounded support has 200 steps, and the bisection has 200 steps.
  A search that meets a bound gives a fault.
- A count of a discrete family has to be a whole number from 0. The
  geometric family starts at 1. A fraction gives `ERR_DOMAIN`. gonum and
  SciPy accept a fraction and truncate it.
- The moments check their poles before any division. A moment with no
  defined value gives `ERR_UNDEFINED`, never an infinity.
- A density at a pole gives `ERR_OVERFLOW`. An example is the gamma
  density at 0 with a shape below 1. SciPy returns an infinity there. A density that does not fit also gives
  `ERR_OVERFLOW`.
- A call with a fault code of 0 never returns an infinity or a NaN. The
  test suite checks 40,960 calls with extreme arguments, from 1e-300 to
  1e300, on every continuous family. A product of two large values keeps
  a finite result. A ratio that underflows uses the difference of two logs.
- The checks run in one order for every call: NaN, infinity, domain. An
  infinite value of x gives `ERR_INFINITY` in every density and
  cumulative call, and an infinite scale gives `ERR_INFINITY` too.

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
goose run examples/special.goose
goose run examples/distributions.goose
goose run examples/tests.goose
```

The scripts `test/make_vectors.py`, `test/make_assoc_vectors.py`,
`test/make_special_vectors.py`, `test/make_distribution_vectors.py` and
`test/make_test_vectors.py` build the vectors of the test suite. They write numbers only into the
package. Each vector is compared with a second source, and the script
fails on a difference beyond the tolerance.

The second sources are in four languages. Python gives the exact value
through the modules `fractions` and `decimal`. NumPy, SciPy and Python
`statistics` give the reference values. Three more libraries give their
own opinion. They are gonum of Go, the crate `statrs` of Rust and the
package `simple-statistics` of Node.

The oracle programs are under `test/`. The program `test/oracle` reads a
list of samples and gives the values of Go with gonum. The program
`test/oracle-rust` gives the values of the crate `statrs`. The program
`test/oracle-node` gives the values of the package `simple-statistics`,
and of the package `jstat` in the distribution mode. Each program needs
its own toolchain and its own packages.

The distribution vectors carry 3,190 rows. The exact value of each row
comes from `mpmath` at 50 digits. Run the vector script with `--oracles`
to compare the three libraries with the rows. The package `jstat` holds
a beta solver with an error near 1e-8, so its rows take a wider bound.
gonum forms a survival function as one minus a value. Its small tails
take the bound that this loss gives.

For the special functions each program runs in the mode `special`. The script
`test/make_special_vectors.py` needs `mpmath` for the exact values. Run it
with `--oracles` to compare gonum, `statrs` and `simple-statistics` with
them. The `gamma` of `simple-statistics` errs by up to half a
percent. The script leaves it out of the gamma comparison.

The script `test/mutate.sh` flips one rule of the code at a time in a copy
of the package. Each flip must make the test program fail. The last line
reads `N probes, 0 failures`.

The environment variable `STATISTICS_HOSTILE_N` sets the count of the
seeded inputs of the hostile suite. The default is 100000. The script
`test/mutate.sh` sets it to 2000, because it runs the suite many times.
The three inverse calls take a fiftieth of that count. Each one costs
milliseconds, so a full count would take too long.

The distributions take a fiftieth of that count too. The cumulative
function of every family costs a continued fraction. The quantile costs
a search. The moments of the hypergeometric cost a sum of the masses of
its support.

The script `test/make_test_vectors.py` writes the vectors of the hypothesis
tests. It checks each exact p-value of the Mann-Whitney and the Wilcoxon
test against a count in Python integers. It also compares the t tests, the
analysis of variance and the Pearson correlation with `jstat` and
`simple-statistics`. Each call has its own bound.

The worst relative difference from SciPy is 2.4e-8 for a p-value. The
statistics agree to 1e-13. Those libraries are less exact than this
package, so a difference there is no defect here.

The same script compares the crate `statrs` with SciPy. The crate holds
seven tests. The statistics and the p-values agree to 1e-9 relative, with
these differences:

- The exact Mann-Whitney p-value of the crate is wrong. For 1, 2, 3, 4
  against 5, 6, 7, 8.5 it gives 1.0. SciPy and this package give
  0.02857142857142857. The script compares the asymptotic form only.
- The crate computes the Kolmogorov-Smirnov p-values with the asymptotic
  series by default. The p-value of its two-sample series differs from
  SciPy, so the script compares the statistic only.
- The crate loses the relative digits of a small chi-square p-value.
  One case differs by 5% at an absolute difference below 1e-12.
- The crate panics on a constant sample with the mean it tests. This package
  gives `ERR_UNDEFINED`.
- The crate gives an infinite statistic and p-value 0 for a zero expected
  count in the chi-square test. This package gives a fault.

gonum holds the chi-square and Kolmogorov-Smirnov statistics but no
p-value. R is not installed, so the rules of R come from its
documentation only.

The hypothesis tests allocate at most 8 times the bytes of their input. The
script `test/memory.sh` measures the peak of every call at the largest
input of the package. The worst case is `wilcoxon` at 7.1 times its input.
The Mann-Whitney test takes 6.1 times its two samples, and the rest take
less. The exact tables hold at most 100,000 cells.

A zero expected count in `chisquare` is a fault. SciPy gives an infinite
statistic and p-value 0 when the sums are equal. It raises an error when the
sums differ. This package gives `ERR_UNDEFINED` for equal sums and
`ERR_DOMAIN` for different sums, so no infinity reaches the caller.

The environment variable `STATISTICS_ROUND_N` sets the count of the
points of the round trip of the quantiles. The default is 2000 points
for each family, so 30,000 points in all. The script `test/mutate.sh`
sets it to 200, because it runs the suite many times.

## License

MIT
