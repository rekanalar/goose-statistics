# goose-statistics

[![ci](https://github.com/rekanalar/goose-statistics/actions/workflows/ci.yml/badge.svg)](https://github.com/rekanalar/goose-statistics/actions/workflows/ci.yml)

Descriptive statistics for the [Goose](https://github.com/aardappel/goose)
programming language.

The package gives the location, the spread, the quantiles and the shape of
a sample of 64 bit floats. It gives the trimmed forms and the weighted
forms too. The sum is exact, so the mean and the variance keep their
digits.

The package gives sixty-seven probability families too. Sixty-two are
continuous and five are discrete. Each family has its density or mass,
its cumulative function, its survival function, its quantile and its
four moments.

Every call returns the result first and a fault code second. A fault gives
NaN as the result. A call never changes the data of the caller.

Pure Goose, no C.

## Install

Add this to the root file of the program that uses it:

```goose
import statistics from "rekanalar/goose-statistics@v0.10.0";
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
| `anglit_pdf(x)`, and the same seven calls | the anglit family, no default |
| `arcsine_pdf(x)`, and the same seven calls | the arcsine family, no default |
| `cosine_pdf(x)`, and the same seven calls | the cosine family, no default |
| `gibrat_pdf(x)`, and the same seven calls | the Gibrat family, no default |
| `halfcauchy_pdf(x)`, and the same seven calls | the half Cauchy family, no default |
| `halflogistic_pdf(x)`, and the same seven calls | the half logistic family, no default |
| `halfnorm_pdf(x)`, and the same seven calls | the half normal family, no default |
| `hypsecant_pdf(x)`, and the same seven calls | the hyperbolic secant family, no default |
| `semicircular_pdf(x)`, and the same seven calls | the semicircular family, no default |
| `triang_pdf(x, c)`, and the same seven calls | the triangular family; `c` is the mode |
| `trapezoid_pdf(x, c, d)`, and the same seven calls | the trapezoid family; `c` and `d` are the ends of the plateau |
| `bradford_pdf(x, c)`, and the same seven calls | the Bradford family; `c` is the shape |
| `wald_pdf(x)`, and the same seven calls | the inverse Gaussian family with a mean of 1, no default |
| `maxwell_pdf(x)`, and the same seven calls | the Maxwell family, no default |
| `moyal_pdf(x)`, and the same seven calls | the Moyal family, no default |
| `dweibull_pdf(x, c)`, and the same seven calls | the double Weibull family of the whole line; `c` is the shape |
| `invweibull_pdf(x, c)`, and the same seven calls | the inverse Weibull family; `c` is the shape |
| `weibull_max_pdf(x, c)`, and the same seven calls | the Weibull family of the negative line; `c` is the shape |
| `truncweibull_min_pdf(x, c, a, b)`, and the same seven calls | the Weibull family on `a` to `b`; `c` is the shape |
| `powernorm_pdf(x, c)`, and the same seven calls | the power of the normal family; `c` is the power |
| `powerlognorm_pdf(x, c, s)`, and the same seven calls | the power of the log-normal family; `c` is the power and `s` the scale |
| `exponpow_pdf(x, b)`, and the same seven calls | the generalized exponential family; `b` is the shape |
| `exponweib_pdf(x, a, c)`, and the same seven calls | the Weibull family of the exponential kind; `a` and `c` are the shapes |
| `fisk_pdf(x, c)`, and the same seven calls | the log-logistic family; `c` is the shape |
| `burr_pdf(x, c, d)`, and the same seven calls | the Burr family of the third kind; `c` and `d` are the shapes |
| `burr12_pdf(x, c, d)`, and the same seven calls | the Burr family of the twelfth kind; `c` and `d` are the shapes |
| `betaprime_pdf(x, a, b)`, and the same seven calls | the beta prime family; `a` and `b` are the shapes |
| `genpareto_pdf(x, c)`, and the same seven calls | the generalized Pareto family; `c` is the shape |
| `genlogistic_pdf(x, c)`, and the same seven calls | the generalized logistic family; `c` is the shape |
| `genhalflogistic_pdf(x, c)`, and the same seven calls | the generalized half logistic family; `c` is the shape |
| `gompertz_pdf(x, c)`, and the same seven calls | the Gompertz family; `c` is the shape |
| `truncexpon_pdf(x, b)`, and the same seven calls | the truncated exponential family; `b` is the upper bound |
| `truncpareto_pdf(x, b, c)`, and the same seven calls | the truncated Pareto family; `b` is the shape and `c` the upper bound |
| `truncnorm_pdf(x, a, b)`, and the same seven calls | the truncated normal family; `a` and `b` are the bounds |
| `chi_pdf(x, df)`, and the same seven calls | the chi family; `df` is the degrees of freedom |
| `dgamma_pdf(x, a)`, and the same seven calls | the double gamma family; `a` is the shape |
| `invgamma_pdf(x, a)`, and the same seven calls | the inverse gamma family; `a` is the shape |
| `loggamma_pdf(x, c)`, and the same seven calls | the log gamma family; `c` is the shape |
| `gengamma_pdf(x, a, c)`, and the same seven calls | the generalized gamma family; `a` and `c` are the shapes |
| `genexpon_pdf(x, a, b, c)`, and the same seven calls | the generalized exponential family; `a`, `b` and `c` are the shapes |
| `mielke_pdf(x, k, s)`, and the same seven calls | the Mielke family; `k` and `s` are the shapes |
| `gennorm_pdf(x, beta)`, and the same seven calls | the generalized normal family; `beta` is the shape |
| `halfgennorm_pdf(x, beta)`, and the same seven calls | the generalized half normal family; `beta` is the shape |
| `tukeylambda_pdf(x, lam)`, and the same seven calls | the law of Tukey and of lambda; `lam` is the shape |
| `jf_skew_t_pdf(x, a, b)`, and the same seven calls | the skew t family of the Johnson kind; `a` and `b` are the shapes |
| `johnsonsb_pdf(x, a, b)`, and the same seven calls | the Johnson family of the B kind; `a` and `b` are the shapes |
| `johnsonsu_pdf(x, a, b)`, and the same seven calls | the Johnson family of the U kind; `a` and `b` are the shapes |
| `invgauss_pdf(x, mu)`, and the same seven calls | the inverse Gaussian family; `mu` is the mean |
| `nakagami_pdf(x, nu)`, and the same seven calls | the Nakagami family; `nu` is the shape |
| `kappa3_pdf(x, a)`, and the same seven calls | the Kappa family of the third kind; `a` is the shape |
| `pearson3_pdf(x, skew)`, and the same seven calls | the Pearson family of the third kind; `skew` is the skewness |
| `genextreme_pdf(x, c)`, and the same seven calls | the generalized extreme value family; `c` is the shape |
| `laplace_asymmetric_pdf(x, kappa)`, and the same seven calls | the asymmetric Laplace family; `kappa` is the shape |
| `skewnorm_pdf(x, a)`, and the same seven calls | the skew normal family; `a` is the shape |
| `bernoulli_pmf(k, p)`, `cdf` and `ppf` | the Bernoulli family; `p` is the chance of 1 |
| `betabinom_pmf(k, n, a, b)`, `cdf` and `ppf` | the beta binomial family; `n` is the count, `a` and `b` the shapes |
| `betanbinom_pmf(k, n, a, b)`, `cdf` and `ppf` | the beta negative binomial family; `n` is the count of successes |
| `boltzmann_pmf(k, lambda, n)`, `cdf` and `ppf` | the truncated geometric family of Boltzmann; `n` is the count of values |
| `dlaplace_pmf(k, a)`, `cdf` and `ppf` | the discrete Laplace family; `a` is the shape |
| `logser_pmf(k, p)`, `cdf` and `ppf` | the logarithmic series family; `p` is the shape |
| `nhypergeom_pmf(k, M, n, r)`, `cdf` and `ppf` | the negative hypergeometric family; `r` is the count of successes |
| `planck_pmf(k, lambda)`, `cdf` and `ppf` | the Planck family, a geometric law from 0; `lambda` is the shape |
| `poisson_binom_pmf(k, ps)`, `cdf` and `ppf` | the Poisson binomial family; `ps` is the list of the chances |
| `randint_pmf(k, lo, hi)`, `cdf` and `ppf` | the uniform family of whole numbers from `lo` up to `hi` less 1 |
| `skellam_pmf(k, mu1, mu2)`, `cdf` and `ppf` | the Skellam family, the difference of two Poisson counts |
| `yulesimon_pmf(k, alpha)`, `cdf` and `ppf` | the Yule-Simon family; `alpha` is the shape |
| `zipf_pmf(k, a)`, `cdf` and `ppf` | the Zipf family without an end; `a` is the shape above 1 |
| `zipfian_pmf(k, a, n)`, `cdf` and `ppf` | the Zipf family on 1 to `n`; `a` is the shape |
| `rice_pdf(x, b)`, and the same seven calls | the Rice family; `b` is the shape, 0 or more and at most 3000 |
| `argus_pdf(x, chi)`, `cdf` and `ppf` | the Argus family; `chi` is the shape |
| `dpareto_lognorm_pdf(x, s, f, scale)`, `cdf` and `ppf` | the Pareto of the log-normal kind; `s` is the shape, `f` the scale and `scale` the shift |
| `erlang_pdf(x, a, loc)`, `cdf` and `ppf` | the Erlang family; `a` is the shape and `loc` the shift |
| `exponnorm_pdf(x, k, loc)`, `cdf` and `ppf` | the exponential-normal family; `k` is the shape and `loc` the shift |
| `foldcauchy_pdf(x, c)`, `cdf` and `ppf` | the folded Cauchy family; `c` is the shift |
| `foldnorm_pdf(x, c)`, `cdf` and `ppf` | the folded normal family; `c` is the shift |
| `loguniform_pdf(x, a, b)`, `cdf` and `ppf` | the log uniform family; `a` and `b` are the ends |
| `irwinhall_pdf(x, n)`, `cdf` and `ppf` | the Irwin-Hall family; `n` is the count of the summed values |
| `gausshyper_pdf(x, a, b, c)`, `cdf` and `ppf` | the Gauss hypergeometric family; `a`, `b` and `c` are the shapes |
| `geninvgauss_pdf(x, p, b)`, `cdf` and `ppf` | the generalized inverse Gaussian family; `p` is the shape and `b` the scale |
| `genhyperbolic_pdf(x, p, a, b)`, `cdf` and `ppf` | the generalized hyperbolic family; `p` is the shape, `a` the location and `b` the scale |
| `norminvgauss_pdf(x, a, b)`, `cdf` and `ppf` | the normal inverse Gaussian family; `a` is the shape and `b` the skewness |
| `recipinvgauss_pdf(x, mu)`, `cdf` and `ppf` | the reciprocal inverse Gaussian family; `mu` is the mean |
| `rel_breitwigner_pdf(x, rho)`, `cdf` and `ppf` | the relative Breit-Wigner family; `rho` is the shape |
| `studentized_range_pdf(x, k, df)`, `cdf` and `ppf` | the studentized range family; `k` is the count of the samples and `df` the degrees of freedom |
| `vonmises_pdf(x, kappa)`, `cdf` and `ppf` | the von Mises family of an angle; `kappa` is the shape |
| `levy_stable_pdf(x, alpha, beta)`, `cdf` and `ppf` | the stable family of Levy; `alpha` is the shape and `beta` the skewness |
| `alpha_pdf(x, a)`, `cdf` and `ppf` | the alpha family; `a` is the shape |
| `fatiguelife_pdf(x, c)`, `cdf` and `ppf` | the fatigue life family; `c` is the shape |
| `levy_pdf(x)`, `cdf` and `ppf` | the Levy family, a stable law of shape one half |
| `levy_l_pdf(x)`, `cdf` and `ppf` | the mirror of the Levy family, on the negative line |
| `logistic_pdf(x)`, `cdf` and `ppf` | the logistic family |
| `loglaplace_pdf(x, c)`, `cdf` and `ppf` | the log Laplace family; `c` is the shape |
| `lomax_pdf(x, c)`, `cdf` and `ppf` | the Lomax family; `c` is the shape |
| `pareto_pdf(x, b)`, `cdf` and `ppf` | the Pareto family from 1; `b` is the shape |
| `powerlaw_pdf(x, a)`, `cdf` and `ppf` | the power law on 0 to 1; `a` is the shape |
| `reciprocal_pdf(x, a, b)`, `cdf` and `ppf` | the same law as `loguniform`; `a` and `b` are the ends |
| `rdist_pdf(x, c)`, `cdf` and `ppf` | the R family on -1 to 1; `c` is the shape |
| `skewcauchy_pdf(x, a)`, `cdf` and `ppf` | the skew Cauchy family; `a` is the skew, below 1 in size |
| `wrapcauchy_pdf(x, c)`, `cdf` and `ppf` | the wrapped Cauchy family on 0 to 2 pi; `c` is between 0 and 1 |
| `kappa4_pdf(x, h, k)`, `cdf` and `ppf` | the four parameter Kappa family; `h` and `k` are the shapes |
| `crystalball_pdf(x, beta, m)`, `cdf` and `ppf` | the crystal ball family; `m` is above 1 |
| `kstwobign_pdf(x)`, `cdf` and `ppf` | the limit of the Kolmogorov-Smirnov statistic |
| `vonmises_line_pdf(x, kappa)`, `cdf` and `ppf` | the von Mises law on the line from -pi to pi |
| `ncx2_pdf(x, df, nc)`, `cdf` and `ppf` | the non-central chi-square family; `nc` is up to 100,000 |
| `ncf_pdf(x, dfn, dfd, nc)`, `cdf` and `ppf` | the non-central F family; `nc` is up to 100,000 |
| `nct_pdf(x, df, nc)`, `cdf` and `ppf` | the non-central t family |
| `ksone_pdf(x, n)`, `cdf` and `ppf` | the one sided Kolmogorov-Smirnov statistic; `n` is up to 1000 |
| `kstwo_pdf(x, n)`, `cdf` and `ppf` | the two sided Kolmogorov-Smirnov statistic; `n` is up to 1000 |
| `landau_pdf(x)`, `cdf` and `ppf` | the Landau family |
| `nchypergeom_fisher_pmf(k, M, n, N, odds)`, `cdf` and `ppf` | the non-central hypergeometric family of Fisher; `M` is up to 100,000 |
| `nchypergeom_wallenius_pmf(k, M, n, N, odds)`, `cdf` and `ppf` | the non-central hypergeometric family of Wallenius; `M` is up to 500 |
| `multivariate_normal_pdf(x, mean, cov, d)`, `mean` and `cov` | the multivariate normal family; the matrices hold their values in row order |
| `multivariate_t_rvs(out, mean, shape, df, d, size, r)` | the draws of the multivariate t family |
| `multivariate_t_mean(out, mean, shape, d)`, `cov` and `df` | the moments of the multivariate t family |
| `matrix_normal_rvs(out, mean, cov, d, size, r)` | the draws of the matrix normal family |
| `matrix_normal_mean(out, mean, cov, d)` and `cov` | the moments of the matrix normal family |
| `matrix_t_rvs(out, mean, shape, df, d, size, r)` | the draws of the matrix t family |
| `matrix_t_mean(out, mean, shape, d)` | the mean of the matrix t family |
| `wishart_rvs(out, shape, df, d, size, r)` and `wishart_mean(out, shape, df, d)` | the Wishart family and its mean |
| `invwishart_rvs(out, shape, df, d, size, r)` | the inverse Wishart family |
| `normal_inverse_gamma_rvs(out, mu, lam, a, b, d, size, r)` | the normal inverse gamma family |
| `multinomial_pmf(x, n, p)`, `mean` and `cov` | the multinomial family |
| `multinomial_rvs(out, n, p, size, r)` | the draws of the multinomial family |
| `multivariate_hypergeom_pmf(x, n, N)` and `mean` | the multivariate hypergeometric family |
| `dirichlet_pdf(x, alpha)`, `logpdf`, `mean`, `var` and `cov` | the Dirichlet family |
| `dirichlet_rvs(out, alpha, size, r)` | the draws of the Dirichlet family |
| `uniform_direction_rvs(out, dim, size, r)` | the draws of the uniform direction on a sphere |
| `ortho_group_rvs(out, dim, size, r)` | the draws of the orthogonal group |
| `special_ortho_group_rvs(out, dim, size, r)` | the draws of the special orthogonal group; each determinant is 1 |
| `unitary_group_rvs(out, dim, size, r)` | the draws of the unitary group; the real parts come first, then the imaginary parts |
| `random_correlation_rvs(out, eigs, r)` | the draw of a correlation matrix with the given eigenvalues |
| `vonmises_fisher_pdf(x, mu, kappa)`, `logpdf` and `rvs` | the von Mises-Fisher family on a sphere |
| `random_table_rvs(out, row, col, size, r)` | the draw of a table with the given row sums and column sums |
| `Dist { kind, a, b }` | the object of a distribution; the kind is 1 for the normal, 2 for the uniform, 3 for the logistic and 4 for the binomial law |
| `dist_pdf(d, x)`, `dist_cdf`, `dist_ccdf`, `dist_icdf` and `dist_iccdf` | the density, the cumulative function, the survival function and their inverses of an object |
| `dist_logpdf(d, x)`, `dist_logcdf`, `dist_logccdf`, `dist_ilogcdf` and `dist_ilogccdf` | the logarithms of the calls of an object |
| `dist_mean(d)`, `dist_median`, `dist_mode`, `dist_variance`, `dist_standard_deviation`, `dist_skewness`, `dist_kurtosis`, `dist_entropy` and `dist_logentropy` | the moments of an object |
| `dist_moment(d, order)` | the raw moment of an object; the order is up to 64 |
| `dist_support(d)` | the lower bound, the upper bound and the code of an object |
| `mixture_pdf(comps, w, x)`, `mixture_cdf`, `mixture_ccdf` and `mixture_icdf` | the mixture of continuous components |
| `mixture_mean(comps, w)`, `mixture_variance`, `mixture_standard_deviation`, `mixture_skewness`, `mixture_kurtosis` and `mixture_moment` | the moments of a mixture |
| `mixture_support(comps, w)` | the lower bound, the upper bound and the code of a mixture |
| `covariance_matrix(out, kind, data, d)` | the covariance matrix of one of the four forms |
| `covariance_log_pdet(kind, data, d)` and `covariance_rank(kind, data, d)` | the log of the pseudo-determinant and the rank of a form |
| `covariance_whiten(out, kind, data, d, x)` and `covariance_colorize(out, kind, data, d, x)` | the whitening of a point and its inverse |
| `histogram_pdf(content, edges, density, x)`, `cdf`, `sf` and `ppf` | the law of a histogram |
| `histogram_mean(content, edges, density)`, `histogram_var`, `histogram_entropy` and `histogram_moment` | the moments of the law of a histogram |
| `histogram_rvs(out, content, edges, density, size, r)` | the draws of the law of a histogram |
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
| `false_discovery_control(out, ps, method)` | the adjusted p-values under the name of SciPy; gives the count and the code |
| `combine_pvalues(ps, method, weights)` | combines the p-values of independent tests; gives the statistic, the p-value and the code |
| `levene(groups, sizes, center, trim)` | the Levene test of equal variance; `center` is a `CENTER_` constant |
| `fligner(groups, sizes, center, trim)` | the Fligner-Killeen test of equal variance |
| `bartlett(groups, sizes)` | the Bartlett test of equal variance |
| `alexandergovern(groups, sizes)` | the Alexander-Govern test of equal variance |
| `ranksums(xs, ys, alt)` | the Wilcoxon rank-sum test of two samples |
| `mood(xs, ys, alt)` | the Mood test of median equality |
| `brunnermunzel(xs, ys, alt, use_t)` | the Brunner-Munzel test; `use_t` picks the t form |
| `ansari(xs, ys, alt, method)` | the Ansari-Bradley test of equal dispersion |
| `friedmanchisquare(data, columns)` | the Friedman test of k samples; `data` holds the blocks row by row |
| `median_test(out, groups, sizes, ties, correction)` | the median test; writes the 2 by k table into `out` |
| `quantile_test(out, xs, q, p, alt, level)` | the test of a quantile; writes the bounds into `out` |
| `ttest_ind_from_stats(mean1, sd1, n1, mean2, sd2, n2, equal_var, alt)` | the two-sample t-test from the summary values |
| `boxcox_llf(lambda, xs)` | the log-likelihood of the Box-Cox transform |
| `yeojohnson_llf(lambda, xs)` | the log-likelihood of the Yeo-Johnson transform |
| `boxcox(out, xs, lambda)` | the Box-Cox transform; writes the values into `out` |
| `yeojohnson(out, xs, lambda)` | the Yeo-Johnson transform; writes the values into `out` |
| `boxcox_fit(out, xs)` | the Box-Cox transform with the lambda of the maximum likelihood |
| `yeojohnson_fit(out, xs)` | the Yeo-Johnson transform with the lambda of the maximum likelihood |
| `boxcox_normmax(xs, method)` | the lambda of the best fit; `method` is a `NORMMAX_` constant |
| `yeojohnson_normmax(xs)` | the lambda of the best fit |
| `boxcox_interval(xs, alpha)` | the low and the high lambda of the confidence interval |
| `boxcox_normplot(out, xs, la, lb, n)` | writes n rows of lambda and r into `out` |
| `yeojohnson_normplot(out, xs, la, lb, n)` | writes n rows of lambda and r into `out` |
| `bootstrap(out, data, sizes, method, level, alt, rng)` | the bootstrap interval and its standard error; writes the null list into `out` |
| `permutation_test(out, data, sizes, method, alt, rng)` | the permutation test; writes the null list into `out` |
| `monte_carlo_test(out, xs, family, params, method, alt, rng)` | the Monte Carlo test against a family |
| `power(family, params, n, level, resamples, rng)` | the power of the test at a sample size |
| `sobol_indices(first, total, d, n, rng)` | the first-order and total Sobol indices; writes the two lists |
| `sobol_from_evaluations(first, total, f_a, f_b, f_ab, d)` | the same indices from three lists of evaluations that the caller made |
| `ResampleMethod` | the struct that holds the method of the resampling calls |
| `Interval` | the struct of a low and a high bound |
| `BootstrapResult` | the struct of the interval and the standard error that `bootstrap` returns |
| `CensoredData` | the struct of four lists: uncensored, left, right and interval times |
| `MvsResult` | the struct of the three intervals that `bayes_mvs` returns |
| `MvsDist` | the struct of the three posterior distributions that `mvsdist` returns |
| `logrank(x, y, alt)` | the log-rank test of two groups of `CensoredData` |
| `goodness_of_fit(out, xs, family, statistic, resamples, rng)` | the goodness-of-fit test over the ten families with `GOF_ANDERSON`, `GOF_KOLMOGOROV`, `GOF_CRAMER` or `GOF_FILLIBEN`; writes the fitted parameters into `out` |
| `barnard_exact(a, b, c, d, alt, pooled, grid)` | the Barnard exact test of a 2 by 2 table |
| `boschloo_exact(a, b, c, d, alt, grid)` | the Boschloo exact test of a 2 by 2 table |
| `kstest(out, xs, family, params, alt)` | the one-sample Kolmogorov-Smirnov test; writes the location and the sign into `out` |
| `kstest_2samp(out, xs, ys, alt)` | the two-sample Kolmogorov-Smirnov test; writes the location and the sign into `out` |
| `cramervonmises(xs, family, params)` | the Cramer-von Mises test against a family |
| `cramervonmises_2samp(xs, ys, method)` | the two-sample Cramer-von Mises test |
| `epps_singleton_2samp(xs, ys, t)` | the Epps-Singleton test of two samples |
| `bws_test(xs, ys, rng, alt, resamples)` | the BWS test of two samples |
| `anderson_ksamp(out, groups, sizes, variant)` | the Anderson-Darling test of k samples; writes the critical values into `out` |
| `poisson_means_test(k1, n1, k2, n2, diff, alt)` | the test of two Poisson means |
| `page_trend_test(data, columns, ranked, predicted, method)` | the Page trend test of k samples |
| `tukey_hsd(out, groups, sizes, level, equal_var)` | the Tukey intervals of every pair; writes the tables into `out` |
| `dunnett(out, groups, sizes, control, alt, level)` | the Dunnett intervals of every group against the control |
| `bayes_mvs(xs, level)` | the Bayesian intervals of the mean, the variance and the deviation |
| `mvsdist(xs)` | the parameters of the three posterior distributions |
| `describe(xs, ddof, bias)` | the summary of a sample; gives a `Summary` and the code |
| `scoreatpercentile(xs, per)` | the value at a percentile; `per` outside 0 to 100 gives `ERR_DOMAIN` |
| `scoreatpercentile_method(xs, per, method)` | the same value with a method of the interpolation |
| `scoreatpercentile_limit(xs, per, low, high)` | the value at a percentile of the values inside the limits |
| `scoreatpercentile_limit_method(xs, per, low, high, method)` | the same value with a method of the interpolation |
| `percentileofscore(xs, score, kind)` | the share of the sample at or below a score |
| `median_abs_deviation(xs, scale, use_mean)` | the median of the absolute distance to the median or the mean |
| `median_abs_deviation_normal(xs)` | the same value with the scale of the normal distribution |
| `pmean(xs, p)` | the power mean of a sample of positive values |
| `pmean_weighted(xs, ws, p)` | the power mean of a sample with weights |
| `gstd(xs, ddof)` | the geometric deviation of a sample of positive values |
| `gzscore(out, xs, ddof)` | the z-scores of the logarithms of a sample; writes them into `out` |
| `zmap(out, scores, compare, ddof)` | the z-scores of a sample against a reference sample |
| `lmoment(out, xs, max_order, standardize)` | the L-moments of a sample; writes them into `out` |
| `expectile(xs, alpha)` | the expectile of a sample at a level from 0 to 1 |
| `expectile_weighted(xs, ws, alpha)` | the expectile of a sample with weights |
| `entropy(pk, base)` | the Shannon entropy of a set of probabilities |
| `relative_entropy(pk, qk, base)` | the entropy of one set of probabilities against another |
| `kstat(xs, order)` | the k-statistic of a sample, for the orders 1 to 4 |
| `kstatvar(xs, order)` | the variance of the k-statistic, for the orders 1 and 2 |
| `differential_entropy(xs, method, window, base)` | the differential entropy of a sample |
| `tiecorrect(ranks)` | the tie correction of a rank test |
| `ecdf(quantiles, probs, xs)` | the empirical distribution function of a sample |
| `ecdf_value(xs, x)` | the share of the values of a sample at or below a point |
| `ecdf_censored(quantiles, survival, uncensored, right)` | the Kaplan-Meier curve of a censored sample |
| `ecdf_censored_interval(low, high, uncensored, right, level)` | the Greenwood band of that curve |
| `ecdf_interval(low, high, xs, level)` | the Greenwood band of the empirical distribution function |
| `estimated_cdf(xs, y, method)` | the estimated distribution function at a point |
| `cumfreq(counts, xs, numbins)` | the cumulative frequencies of a sample; gives a `Frequencies` and the code |
| `cumfreq_limits(counts, xs, numbins, low, high)` | the same counts over a range of the values |
| `cumfreq_weighted(counts, xs, ws, numbins)` | the same counts with weights |
| `relfreq(freqs, xs, numbins)` | the relative frequencies of a sample |
| `relfreq_limits(freqs, xs, numbins, low, high)` | the same shares over a range of the values |
| `relfreq_weighted(freqs, xs, ws, numbins)` | the same shares with weights |
| `tmin(xs, low, inclusive)` | the smallest value of a sample at or above a limit |
| `tmax(xs, high, inclusive)` | the largest value of a sample at or below a limit |
| `tvar(xs, low, high, use_low, use_high, lo_inc, hi_inc, ddof)` | the variance of the values inside the limits |
| `tstd(xs, low, high, use_low, use_high, lo_inc, hi_inc, ddof)` | the deviation of those values |
| `tsem(xs, low, high, use_low, use_high, lo_inc, hi_inc, ddof)` | the standard error of the mean of those values |
| `trimboth(out, xs, proportion)` | the sample without a share of the values at each end |
| `trim1(out, xs, proportion, tail)` | the sample without a share of the values at one end |
| `circmean(xs, low, high)` | the circular mean of a sample of angles |
| `circvar(xs, low, high)` | the circular variance of a sample of angles |
| `circstd(xs, low, high, normalize)` | the circular deviation of a sample of angles |
| `directional_stats(dir, xs, dims, normalize)` | the mean direction of a set of vectors and its resultant length |
| `binned_statistic(stats, edges, shape, binnumber, empty, xs, values, stat, nbins)` | the statistic of the values of each bin of one dimension |
| `binned_statistic_edges(stats, edges, shape, binnumber, empty, xs, values, stat, given)` | the same statistic with the edges of the caller |
| `binned_statistic_with(stats, edges, shape, binnumber, empty, xs, values, nbins)` | the same statistic, where the caller gives the statistic |
| `binned_statistic_2d(stats, edges, shape, binnumber, empty, xs, ys, values, stat, nx, ny)` | the statistic of the values of each bin of two dimensions |
| `binned_statistic_2d_with(stats, edges, shape, binnumber, empty, xs, ys, values, nx, ny)` | the same statistic, where the caller gives the statistic |
| `binned_statistic_dd(stats, edges, shape, binnumber, empty, sample, dims, values, stat, nbins)` | the statistic of the values of each bin of many dimensions |
| `binned_statistic_dd_with(stats, edges, shape, binnumber, empty, sample, dims, values, nbins)` | the same statistic, where the caller gives the statistic |
| `spearmanrho(xs, ys, alt)` | the Spearman rank correlation, with its p-value |
| `pointbiserialr(xs, ys, alt)` | the point-biserial correlation; `xs` holds two values |
| `chatterjeexi(xs, ys, y_continuous)` | the xi correlation of Chatterjee, with its asymptotic p-value |
| `weightedtau(xs, ys, additive)` | the weighted tau of Vigna, with the hyperbolic weigher |
| `weightedtau_ranked(xs, ys, ranks, additive)` | the weighted tau with the ranks of the caller |
| `somersd(xs, ys, alt)` | Somers' D of two coded samples, with its p-value |
| `somersd_table(table, rows, cols, alt)` | Somers' D of a contingency table |
| `sigmaclip(out, xs, low, high)` | the iterative clip by standard deviations; writes the kept values into `out` |
| `power_divergence(obs, expected, lam, ddof)` | the power divergence statistic; `lam` 1 is Pearson and 0 is the G test |
| `power_divergence_uniform(obs, lam, ddof)` | the same, against equal expected counts |
| `margins(out, table, shape, axis)` | the sums of a table over the axes that `axis` names |
| `expected_freq(out, table, shape)` | the expected counts of a table under independence |
| `association(table, rows, cols, method, correction, lam)` | Cramer's V, Tschuprow's T or the contingency coefficient |
| `odds_ratio(table, kind, alt, level)` | the sample or conditional odds ratio of a 2 by 2 table, with its interval |
| `relative_risk(exposed_cases, exposed_total, control_cases, control_total, level)` | the relative risk with the Katz interval |
| `wasserstein_distance_nd(u, v, dims)` | the exact first Wasserstein distance of two sets of points, with the Euclidean cost |
| `wasserstein_distance_nd_weighted(u, v, dims, uw, vw)` | the same, with a weight for each point |
| `multiscale_graphcorr(map, xs, ys, xdims, ydims, reps, r)` | the multiscale graph correlation of two samples of points; writes the map of local correlations into `map` |
| `multiscale_graphcorr_dist(map, dx, dy, reps, r)` | the same, from two distance tables of n by n values |
| `multiscale_graphcorr_twosample(map, u, v, dims, reps, r)` | the test of two samples of points of the same width |
| `wasserstein_distance(u, v)` | the first Wasserstein distance of two samples |
| `wasserstein_distance_weighted(u, v, uw, vw)` | the same, with weights |
| `energy_distance(u, v)` | the energy distance of two samples |
| `energy_distance_weighted(u, v, uw, vw)` | the same, with weights |
| `obrientransform(out, groups, sizes)` | the transform of O'Brien for a test of equal variances |
| `siegelslopes(ys, xs, method)` | the repeated median line; gives the slope and the intercept |
| `probplot(osm, osr, xs, family, shape)` | the quantile plot of a sample; the fit line and the correlation |
| `ppcc_plot(out, xs, a, b, n, family)` | the correlation of a quantile plot over a range of shapes |
| `ppcc_max(xs, family, lo, hi)` | the shape with the best correlation |
| `anderson_exponential(crit, level, params, xs)` | the Anderson-Darling test of the exponential family |
| `anderson_logistic(crit, level, params, xs)` | the Anderson-Darling test of the logistic family |
| `anderson_gumbel_r(crit, level, params, xs)` | the Anderson-Darling test of the right Gumbel family |
| `anderson_gumbel_l(crit, level, params, xs)` | the Anderson-Darling test of the left Gumbel family |
| `anderson_weibull(crit, level, params, xs)` | the Anderson-Darling test of the Weibull family |
| `cauchy_pdf`, `cauchy_cdf`, `cauchy_sf`, `cauchy_ppf` | the Cauchy law with `loc` and `scale`; also `cauchy_mean`, `cauchy_var`, `cauchy_skew`, `cauchy_kurt` |
| `laplace_pdf`, `laplace_cdf`, `laplace_sf`, `laplace_ppf` | the Laplace law; the same moment calls |
| `rayleigh_pdf`, `rayleigh_cdf`, `rayleigh_sf`, `rayleigh_ppf` | the Rayleigh law; the same moment calls |
| `gumbel_r_pdf`, `gumbel_r_cdf`, `gumbel_r_sf`, `gumbel_r_ppf` | the right Gumbel law; the same moment calls |
| `gumbel_l_pdf`, `gumbel_l_cdf`, `gumbel_l_sf`, `gumbel_l_ppf` | the left Gumbel law; the same moment calls |
| `truncate_pdf`, `truncate_cdf`, `truncate_ppf`, `truncate_mean`, `truncate_var` | a law truncated to `lb` and `ub` |
| `order_statistic_pdf`, `order_statistic_cdf`, `order_statistic_ppf`, `order_statistic_mean`, `order_statistic_var` | the r-th of n draws from a law |
| `dirichlet_multinomial_logpmf(x, alpha, n)` | the log mass; also `dirichlet_multinomial_pmf`, `_mean`, `_var`, `_cov` |
| `kde_fit(data, dims, weights, bw, factor)` | a Gaussian kernel estimate; `bw` is `KDE_SCOTT`, `KDE_SILVERMAN` or `KDE_SCALAR` |
| `kde_pdf(k, out, xs)` | the density of the estimate; `kde_logpdf` gives the logarithm |
| `kde_covariance(k, out)` | the covariance of the estimate, with the bandwidth in it |
| `kde_integrate_box_1d(k, low, high)` | the share of a one-dimensional estimate between two bounds |
| `kde_resample(k, out, size, r)` | draws points from the estimate with a seeded stream |

The argument `method` is one of the nine `METHOD_` constants. The
constants `METHOD_EXCLUSIVE` and `METHOD_INCLUSIVE` stand for the two
methods of Python.

The struct `Summary` holds the summary of a sample. Its fields are `n`,
`min`, `max`, `mean`, `variance`, `skewness` and `kurtosis`. The struct
`Frequencies` holds the limits of a histogram. Its fields are
`lowerlimit`, `binsize` and `extrapoints`. The struct `BinnedTable` holds
the sizes of a table of bins. Its fields are `dims`, `cells` and
`points`.

The argument `kind` of `percentileofscore` is one of the four `PCT_`
constants. They are `PCT_RANK`, `PCT_WEAK`, `PCT_STRICT` and `PCT_MEAN`.

The argument `method` of `scoreatpercentile_method` is one of the three
`INTERP_` constants. They are `INTERP_FRACTION`, `INTERP_LOWER` and
`INTERP_HIGHER`.

The argument `method` of `estimated_cdf` is one of the nine `CDF_`
constants, from `CDF_INTERPOLATED` to `CDF_CLOSEST`.

The argument `method` of `differential_entropy` is one of the five
`ENTROPY_` constants, from `ENTROPY_AUTO` to `ENTROPY_CORREA`.

The argument `stat` of the binned calls is one of the seven `STAT_`
constants. They are `STAT_MEAN`, `STAT_STD`, `STAT_MEDIAN`, `STAT_COUNT`,
`STAT_SUM`, `STAT_MIN` and `STAT_MAX`. The argument `tail` of `trim1` is
`TAIL_LEFT` or `TAIL_RIGHT`.

The binned calls fill five output lists. The caller passes them in this
order: `stats`, `edges`, `shape`, `binnumber` and `empty`. The list
`stats` holds one value for each cell, in the row-major order. The list
`edges` holds the edges of each dimension, one dimension after the other.
The list `shape` holds the count of the edges of each dimension.

The list `binnumber` holds one bin number for each point of the caller.
The list `empty` holds 1 for each empty cell and 0 for the others. An
empty cell holds the value 0 in the list `stats`.

The struct `ResampleMethod` holds the method of a resampling call. Its
fields are `kind`, `resamples`, `interval` and `permutation`. The field
`kind` is one of `RESAMPLE_BOOTSTRAP`, `RESAMPLE_PERMUTATION` and
`RESAMPLE_MONTE_CARLO`. The field `interval` is one of
`INTERVAL_PERCENTILE`, `INTERVAL_BASIC` and `INTERVAL_BCA`. The field
`permutation` is one of `PERMUTE_INDEPENDENT`, `PERMUTE_SAMPLES` and
`PERMUTE_PAIRINGS`. The seed comes from the `rng` argument.

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
argument `method` of `multipletests` is one of the seven `ADJUST_`
constants. They are `ADJUST_BONFERRONI`, `ADJUST_HOLM`, `ADJUST_BH`,
`ADJUST_HOCHBERG`, `ADJUST_HOMMEL`, `ADJUST_BY` and `ADJUST_NONE`.

The argument `method` of `false_discovery_control` is one of
`ADJUST_BH` and `ADJUST_BY`. The argument `method` of `combine_pvalues`
is one of the three `COMBINE_` constants. They are `COMBINE_FISHER`,
`COMBINE_STOUFFER` and `COMBINE_TIPPETT`. Only the Stouffer method takes
the list of weights.

The calls `truncate_pdf`, `truncate_cdf`, `truncate_ppf`,
`truncate_mean`, `truncate_var` and the calls `order_statistic_pdf`,
`order_statistic_cdf`, `order_statistic_ppf`, `order_statistic_mean`,
`order_statistic_var` take a family number. The numbers are the
`LAW_` constants. They are `LAW_NORMAL`, `LAW_UNIFORM`, `LAW_LOGISTIC`,
`LAW_CAUCHY`, `LAW_LAPLACE`, `LAW_RAYLEIGH`, `LAW_GUMBEL_R` and
`LAW_GUMBEL_L`. The first two take a location and a scale. The logistic
law takes them too. The Cauchy, the Laplace, the Rayleigh and the two
Gumbel laws take them.

The calls `probplot`, `ppcc_plot` and `ppcc_max` take a family number.
The numbers are the `PLOT_` constants. They are `PLOT_NORMAL`,
`PLOT_EXPON`, `PLOT_UNIFORM`, `PLOT_TUKEY`, `PLOT_LOGISTIC`,
`PLOT_LAPLACE`, `PLOT_CAUCHY`, `PLOT_GUMBEL_R`, `PLOT_GUMBEL_L`,
`PLOT_RAYLEIGH` and `PLOT_WEIBULL`. The Tukey and the Weibull forms take
a shape. The other forms take no shape.

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
| `ERR_LENGTH` | a list has another length than the call needs |
| `ERR_UNDEFINED` | the result has no defined value |
| `ERR_OVERFLOW` | the true result does not fit in a 64 bit float |
| `ERR_TOO_LONG` | the sample has more than `MAX_INPUT` values |
| `ERR_FIT` | a fit or a search did not end, or it gave a scale that is not usable |

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

The moments of the families with no closed form come from a rule of
Gauss-Legendre on panels. The rule finds the peak of the weighted
density first. It cuts the range where the weighted density falls below
1e-15 of that peak. A law of a nonnegative support is integrated in the
logarithm of the point. That keeps the accuracy of a law with a long
tail. The rule uses 16 nodes on 100 panels.

## Limits

- The `levene` and `fligner` calls take a center of the mean, the median
  or the trimmed mean. SciPy takes the same three. The trimmed mean of
  SciPy uses the partition of NumPy, which can end one rounding step
  off the exact mean. The two centers then break a tie of the deviations
  differently. The package keeps the exact mean.
- The `ansari` call counts the exact distribution with exact sums. The
  table of SciPy holds single-precision floats. Its exact p-value loses
  about 1e-9 of its digits on a large sample.
- The `boxcox` and `yeojohnson` calls refuse a value that the form does
  not take. SciPy gives NaN or an infinity there. A sample of fewer than
  two values gives `ERR_TOO_SHORT`.
- The `boxcox_normmax` and `yeojohnson_normmax` calls give a lambda that
  differs from SciPy by about 1e-7 relative. The search of SciPy follows
  another path through its Brent minimizer. The difference of the
  mle lambda of SciPy is about 1.6e-7. The difference of the pearsonr
  lambda is about 2.5e-7.
- The `boxcox_normmax` call refuses data that needs the Lambert W guard
  of SciPy. It gives `ERR_OVERFLOW` there. SciPy swaps in a smaller
  lambda.
- The `tukey_hsd` call computes the upper tail of the studentized range
  directly. SciPy computes one minus the cumulative function, so its
  p-value loses its digits below about 1e-12.
- The `dunnett` call computes its probability with a deterministic
  integral. SciPy uses a randomized quasi Monte Carlo with no fixed
  seed. Its result errs by about 1.5e-4 and it changes between calls.
  The one-sided critical value of SciPy searches a noisy objective with
  a tolerance of 1e-3.
- The `bayes_mvs` and `mvsdist` calls need at least four values. The
  variance posterior of SciPy has no mean below that. SciPy gives an
  infinity or a NaN there.
- The `sobol_indices` call uses pseudo-random matrices from the given
  generator. SciPy uses a scrambled Sobol sequence. SciPy converges
  faster for the same count of points.
- The `goodness_of_fit` call fits the normal, exponential, uniform,
  lognormal, Weibull, gamma and beta families. The fit of the Weibull,
  gamma and beta families uses a Newton search with a fixed step bound.
  SciPy calls its optimizer for them. The fitted parameters differ from
  the ones of SciPy. The gamma family differs by about 1e-6.
  The Weibull family differs by about 1e-4.
- The beta family of `goodness_of_fit` differs from SciPy by about
  2e-3. Its fit searches two parameters at once.
- The chi-square, student and F families of `goodness_of_fit` use a
  Nelder-Mead search. SciPy calls its own optimizer with the location
  free. The likelihood of the chi-square family has no unique optimum
  there. This package fixes the location at 0 for the chi-square and F
  families. It does the same for the gamma and Weibull families. The
  parameters agree with SciPy to about 1e-4.
- The student family of `goodness_of_fit` keeps its location free.
- The `GOF_FILLIBEN` statistic of `goodness_of_fit` is a correlation,
  so a small value is extreme. Its p-value counts the null values at
  most as large as the observed one. SciPy counts the same way for
  this statistic.
- The `cramervonmises` call stops its series at an absolute term size of
  1e-7, as SciPy does. Its p-value then differs from SciPy by about
  1e-5 on a random sample. The statistic differs by about 1e-9.
- The `epps_singleton_2samp` call gives a p-value that differs from
  SciPy by about 1e-6 on a random sample. Its pseudo-inverse of the
  covariance matrix loses digits.
- The `poisson_means_test` call differs from SciPy by about 2e-9 of its
  digits. The boundary of the Poisson quantile differs slightly.
- The `kstest` and `kstest_2samp` calls give the dispatcher result of
  the package `ks_1samp` and `ks_2samp`. A tie inside a sample uses the
  asymptotic form, where SciPy would try the exact form.
- The `bws_test` call above 9999 splits draws from the given generator.
  SciPy draws without a seed, so the two cannot match. The call below
  that bound counts exactly.
- The `anderson_ksamp` call gives the statistic of SciPy to 1e-9. Its
  p-value differs by about 1e-9.
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
- The Wilcoxon test counts all sign patterns for at most 13 pairs and
  uses the normal form above. The random permutation form of SciPy is
  not built. The call `permutation_test` covers it.
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
- The multiple-testing call holds seven methods. They are the
  Bonferroni, Holm, Benjamini-Hochberg, Hochberg, Hommel,
  Benjamini-Yekutieli and the identity. R `p.adjust` holds the same
  seven. SciPy 1.18.1 does not hold the call. The values of Hochberg,
  Hommel and Benjamini-Yekutieli come from the rules of R `p.adjust`.
  R was not run, so statsmodels 0.15.0 and SciPy 1.18.1 gave the
  measured values.
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
- The `sobol_indices` call uses pseudo-random matrices from the given
  generator. SciPy uses a scrambled Sobol sequence. SciPy converges
  faster for the same count of points.
- The `bootstrap` call with the BCa interval refuses more than one
  sample. The jackknife of several samples is not built.
- The draws of the resampling calls come from the given generator. The
  seed of SciPy uses another stream. A seeded run agrees with SciPy only
  within the error of the Monte Carlo.
- The `rand_int` of the standard library gives another stream on v2026.3
  and on v2026.4. The package draws its own counts, so one seed gives
  one answer on both releases.
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

- The `scoreatpercentile` call has the three methods of the
  interpolation of SciPy: `fraction`, `lower` and `higher`. SciPy has
  five more, which use a midpoint and an inverted function.
- The `percentileofscore` call has the four kinds of SciPy: `rank`,
  `weak`, `strict` and `mean`.
- The `kstat` call has the orders 1 to 4. The `kstatvar` call has the
  orders 1 and 2. SciPy has the same orders.
- The `lmoment` call writes the orders 1 to `max_order`. SciPy gives
  one order at a time. The constant `MAX_LORDER` of 16 is the largest
  order of the call.
- The `expectile` call solves the equation of SciPy by bisection. It
  ends with a tolerance of 1e-14 relative. The value of SciPy differs
  by about 1e-9.
- The `pmean` call uses the linearized form of SciPy for an exponent
  with an absolute size of at most 2e-6. The value of SciPy uses the
  same form. The other exponents use the power mean directly. A zero
  value gives `ERR_UNDEFINED` in that band, as SciPy gives NaN. The
  exponent 0 gives the geometric mean, as in SciPy.
- The `differential_entropy` call has the four estimators of SciPy and
  the choice of the size of the sample. The window must satisfy
  `2 * window < n`. A tie gives `ERR_UNDEFINED`, as SciPy gives NaN.
- The `ecdf_interval` and the `ecdf_censored_interval` calls use the
  Greenwood formula of SciPy. The formula has no value at the last
  observation, so the band holds NaN there. The bounds are clipped to
  0 and 1.
- The `ecdf_censored` call takes uncensored and right-censored values.
  SciPy also takes left-censored and interval-censored values. This
  package gives `ERR_DOMAIN` for them.
- The `cumfreq` and the `relfreq` calls take the range of the sample
  when the caller gives no limits. SciPy widens that range by a half
  step of a bin on each side. The package widens it the same way.
- The `trimboth` and the `trim1` calls write the kept values in
  ascending order. SciPy leaves them in the order of the sample. The
  package sorts every result that it gives as a list.
- The `binned_statistic` calls number the bins in the row-major order
  over `nbins + 2` slots per dimension. SciPy numbers them the same
  way. The value that rounds to the last edge stays in the last bin.
- The `directional_stats` call gives the mean resultant length as its
  result. SciPy gives it as a field of a struct.
- `wasserstein_distance_nd` solves the transport problem exactly with
  the method of the shortest augmenting path. SciPy solves it with a
  linear program. Each set holds at most `MAX_NDW_POINTS` of 256 points
  and `MAX_NDW_DIMS` of 16 dimensions. A mass below 1e-14 of the total
  counts as used up.
- `multiscale_graphcorr` gives the statistic and the scale of SciPy to
  1e-12. The p-value uses the seeded stream of the caller, so it matches
  SciPy within the sampling error. A `reps` of 0 gives the p-value 1, as
  SciPy does. A sample holds at most `MAX_MGC_PAIRS` of 1,024 points and
  a call at most `MAX_MGC_REPS` of 10,000 repetitions. Its pairs squared
  times one plus the repetitions stay under `MAX_MGC_WORK` of 8,000,000.
- When every distance of a sample is equal, SciPy gives the scale as one
  number. This package gives the pair of counts.
- `weightedtau` holds the hyperbolic weigher and the rank list of the
  caller. A weigher that is a function of the caller is not possible,
  because Goose cannot store a block.
- `odds_ratio` gives the conditional estimate as the root of the mean of
  the noncentral hypergeometric law. A table with a zero row sum or a
  zero column sum gives `ERR_UNDEFINED`.
- The `kde` calls keep one point per row. SciPy keeps one point per
  column. The evaluation points use the same layout as the data.
- The `anderson` calls of the new families read the p-value from the
  critical values by linear interpolation. They clamp it to the first
  and the last level.
- The logistic fit of `anderson_logistic` agrees with SciPy to about
  1e-5, because SciPy stops its own solver at that tolerance.
- The `powernorm`, `powerlognorm`, `exponpow` and `exponweib` moments
  come from an integration of the density of the law. SciPy integrates
  the quantile function instead. The two routes agree on the mean of
  `powernorm`, `exponpow` and `exponweib` to about 1e-11. They differ
  by about 1e-8 on the skewness and on the kurtosis of those three
  families.
- The `powerlognorm` moments differ from SciPy by more. SciPy gives a
  kurtosis of 21062 for the shape 1.4 and the scale 2. The package
  gives 148180.
- The rule of the moments finds a peak that is narrow. It also finds a
  peak that lies far from the origin. The moments of `powerlognorm` with
  the shape 1 and the scale 0.001 exist. The moments of `exponweib` with
  the shape 1 and the shape 0.03 exist.
- The moments of `exponweib` with the shape 1 and the shape 0.03 agree
  with the exact law. That law is the Weibull family. SciPy gives a
  mean of 6.81e32. The exact law gives 2.80e37. The package gives
  2.80e37.
- The variance, the skewness and the kurtosis of a law with a narrow
  peak lose digits. They come from differences of the raw moments. The
  `powerlognorm` law of the scale 0.001 has a variance of 1.0000015e-6.
  The package gives it to about 3e-9 relative. Its kurtosis is 1.6e-5.
  The package gives -0.0037 and SciPy gives -0.0219.
- Three independent routes agree on the value of the package. They are
  the density, the survival function, and the exact law of the
  log-normal family at the shape 1.
- A sample of four million draws gives a fourth raw moment near 3.6e8.
  That value is below the value of the package. It matches the known
  bias of the estimator. The tail of SciPy loses its digits here.
- The `gompertz` moments come from an integration of the density of the
  law. SciPy integrates the quantile function instead. The two routes
  differ by about 1e-8 on the kurtosis of the shape 1.5. They agree to
  about 1e-11 on the other three moments.
- The `genhalflogistic` moments come from an integration of the density
  of the law too. They agree with SciPy to about 1e-11.
- The `fisk`, `burr`, `burr12`, `betaprime` and `genpareto` moments have
  a closed form. They agree with SciPy to about 1e-13.
- The `invweibull` moments of the shape 1.5 do not exist. The variance
  of that law is negative, so the package gives `ERR_UNDEFINED`. SciPy
  gives the same negative variance and a NaN for the skewness and for
  the kurtosis.
- The `gengamma` moments have a closed form. They agree with SciPy to
  about 1e-14. A shape `c` below 0 gives a moment that exists only up
  to a limit. The package gives `ERR_UNDEFINED` above that limit.
  SciPy integrates the law and gives a number.
- The `genexpon` moments come from an integration of the density of the
  law. SciPy integrates the quantile function instead. The two routes
  differ by about 1e-11 on the mean and by about 1e-10 on the kurtosis.
- The `halfgennorm` moments have a closed form. They agree with SciPy to
  about 1e-15.
- The `johnsonsb` moments come from an integration of the density of the
  law. SciPy integrates the density too and loses digits. The two routes
  differ by about 4e-7 on the variance. They differ by about 1e-4 on the
  kurtosis of the shape 1 and the shape 2. The package gives the value
  of the law.
- The `johnsonsu` moments have a closed form. They agree with SciPy to
  about 1e-15.
- The `mielke` moments have a closed form. The third moment of the
  shape 2 and the shape 3 does not exist. The package gives
  `ERR_UNDEFINED` for the skewness and for the kurtosis. SciPy gives an
  infinity for the skewness and a NaN for the kurtosis.
- The `truncpareto` law has a bounded support, so all of its moments
  exist for every shape. SciPy gives `ERR_UNDEFINED` for a shape that
  is not above the order of the moment. The package gives the value.
- The `tukeylambda` cumulative function searches the root of the
  quantile function. SciPy solves the same root with a different
  method. The two agree to about 1e-12.
- The `gennorm` and the `halfgennorm` moments have a closed form in the
  gamma function. They agree with SciPy to about 1e-15.
- The `chi`, `dgamma`, `invgamma` and `loggamma` moments have a closed
  form. They agree with SciPy to about 1e-13.
- The `invgauss` quantile searches the root of the logarithm of the
  cumulative function or of the survival function. A tiny target keeps
  its digits. SciPy solves the same root with the code of Boost. The
  two agree to about 1e-14 in the body of the law.
- SciPy loses the `invgauss` root in the deep tail. At a target of
  1e-300 and a mean of 0.05, SciPy gives 1.1187e+248 and reports that
  it found no solution. The package gives 0.000707432903319254979.
  The value of mpmath 1.4.1 confirms the answer of the package.
- The `invgauss` cumulative function and survival function use the
  logarithm of the normal law. The far tail keeps its digits. SciPy uses
  the same route.
- The `kappa3` moments have a closed form. SciPy integrates the quantile
  function and loses a few digits. The two agree to about 1e-12.
- The moment of order n of the `kappa3` law exists only when the shape
  is above n. The package gives `ERR_UNDEFINED` above that limit. SciPy
  gives a NaN.
- The `pearson3` law is the law of the normal form when the skewness
  falls below 1.6e-5. The shape of the law grows without bound as the
  skewness falls to zero. The package uses the law of the normal form
  there. The moments of the law have no fault there. SciPy gives a NaN
  for the moments and a value below 1e-300 for the density.
- The density of the `pearson3` law keeps its digits at a small
  skewness. The direct form takes the difference of two large numbers.
  The package uses the ratio of the point to the shape. SciPy loses
  about 6 digits at a skewness of 1e-4.
- The `genextreme` law with the shape 0 is the law of Gumbel. The
  package uses the law of Gumbel for the density, the cumulative
  function, the survival function and the quantile. SciPy uses the
  general form and loses the shape of 0.
- A small shape gives the skewness and the kurtosis of the Gumbel law
  for `genextreme`. The package uses the same cutoffs as SciPy. The two
  agree.
- The moments of the `genextreme` law need a shape above a bound. The
  package gives `ERR_UNDEFINED` below that bound. SciPy gives a NaN.
- The `skewnorm` cumulative function uses the T function of Owen. SciPy
  uses the code of Boost. The two agree to about 1e-14 in the body of
  the law. In the far left tail both forms cancel the same two terms.
  The relative accuracy of both falls there. It falls below a cumulative
  value of about 1e-7.
- The `skewnorm` far left tail uses a different form in the package.
  It keeps the digits of the tail. SciPy loses a few of them.
- The `skewnorm` quantile searches the root of the cumulative function
  or of the survival function. The two agree to about 1e-14 in the body
  of the law.
- SciPy loses the `skewnorm` root in the deep tail. At a target of
  1e-300 and a shape of 0.5, SciPy gives -36.231821510534026. The
  package gives -33.059565133730004. The value of mpmath 1.4.1 confirms
  the answer of the package.

- The `kstwo_pdf` call takes a central difference of the cumulative
  function. SciPy uses an exact derivative. The two agree to about 1e-8.
- The `ksone` and `kstwo` calls accept a sample size up to 1000. SciPy
  has no bound. The `ksone` sum cancels in its far tail for a large
  size.
- The `landau` calls give a density of 0 and a cumulative value of 0
  below -2.5. The integral cancels its terms in the far left tail, and
  its absolute error would pass 1e-8 there. The true value at -2.5 is
  about 1e-5. The answer has an absolute error up to 1e-5 below that
  point. SciPy holds a rational approximation there.
- The `nchypergeom_wallenius` call accepts a total up to 500. Each mass
  needs one integral. The Fisher form accepts a total up to 100,000.
- The `ncx2` and `ncf` calls sum a Poisson mixture over twelve standard
  deviations. They refuse a non-centrality above 100,000.
- The `nct` call integrates over the scale of a chi law. Its accuracy
  falls for a degrees of freedom above 1e6.

- The `levy_stable` calls follow the method of Nolan. The integral of
  Nolan needs a split at the peak of the density, where the helper
  equals one. The package makes that split. A rule without it loses
  about 8 digits.
- The rows of the stable law come from the Fourier form of the law, at 30
  digits of mpmath. SciPy 1.18.1 checks the rows. A point joins the rows
  only where the two agree to 1e-9. The package and SciPy agree to about
  1e-14 at every other point of the grid.
- The far tail of the stable law falls below the resolution of double
  precision. At a shape of 0.8 and a skewness of 1, the density of SciPy
  1.18.1 is not monotone. It gives 5.09e-18 at the point 0.0, 0.0 at 0.3,
  3.57e-139 at 0.4 and 2.26e-65 at 0.5.
- The package gives 1.84e-65 for the density at a shape of 0.8, a
  skewness of 1 and the point 0.5. It gives 1.07e-16 for the cumulative
  function there, where SciPy gives 2.26e-65 and 2.0e-204. The Fourier
  form of the law gives a negative density at that point. No method of
  double precision gives the value. The rows of the law leave that region
  out.
- The `random_correlation` call needs a draw that reaches a diagonal of
  1. A shape that no draw can reach gives `ERR_FIT`.
- The `random_table` call draws each cell from the hypergeometric law
  of what is left. The draw matches the law of Patefield in its
  distribution, but not its values, because SciPy uses another
  generator.
- The `vonmises_fisher` call needs a unit vector for its point and for
  its mean. A vector that differs from a unit vector by more than 1e-8
  gives `ERR_DOMAIN`.
- The `unitary_group` call gives the real parts of a matrix first, then
  its imaginary parts, each in row order. SciPy gives a complex matrix.
- The `histogram` calls need the content of n bins and the n+1 edges.
  The edges must rise. A bin with no mass gives `ERR_DOMAIN`.
- The `mixture` calls take continuous components only, as SciPy does.
  A discrete component gives `ERR_DOMAIN`.
- The `dist` calls of an object take a plain struct. A kind that the
  package does not know gives `ERR_DOMAIN`.
- The binomial object follows the new object classes of SciPy 1.18.1.
  Its cumulative function is the continuous extension of the law: the
  regularized incomplete beta function `I(1 - p)` over `(n - x, x + 1)`.
  The legacy `binom` family of SciPy floors the point instead. The two
  disagree at a point that is not a whole count.
- The draws of the multivariate laws use the inverse of the cumulative
  function of the chi law for the scales. SciPy draws the sum of the
  squares of the normals instead. One seed does not give the same
  values as SciPy, but the same law.
- The `covariance` calls take a kind and a flat buffer. The kind 1
  holds the diagonal. The kind 2 holds the precision matrix. The kind 3
  holds the lower Cholesky factor. The kind 4 holds the eigenvalues and
  then the eigenvectors. A buffer of another length gives `ERR_LENGTH`.

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
- The Hommel form of `multipletests` costs time in the square of the
  count of the list. A list above 10,000 values gives the fault
  `ERR_TOO_LONG` before any work. R gives no bound. A caller with a
  longer list uses another form, or splits the list.
- `combine_pvalues` gives the fault `ERR_OVERFLOW` when the Fisher or
  the Stouffer statistic would be infinite. That happens at a p-value
  of 0, and at a p-value of 1 for Stouffer. SciPy returns an infinite
  statistic there. A caller who ignores codes sees no infinity here.
- `combine_pvalues` rejects a negative weight and a weight that is not
  finite. SciPy accepts a negative weight. A negative weight has no
  meaning for a combination of evidence.
- `multipletests` and `false_discovery_control` give `ERR_EMPTY` for an
  empty list. SciPy returns an empty array. A caller who has no p-values
  checks the length first.
- `combine_pvalues` gives `ERR_DOMAIN` for weights with the Fisher or the
  Tippett method. SciPy ignores them there. A caller who passes weights
  by mistake sees the fault and not a silent change.
- The `invgauss` and the `skewnorm` quantile search a root of a
  continuous function. The search has a bound on its steps and on the
  range that it covers. A hostile shape cannot make it run forever.
- The `kappa3` and the `genextreme` moments give the fault
  `ERR_UNDEFINED` when the moment does not exist. A caller sees the
  fault and not a silent NaN.
- The `pearson3` moments have no fault at a small skewness. The law is
  the law of the normal form there. The moments stay finite and the
  caller sees a value.
- The resampling calls bound their work before any draw. A count of
  resamples above 1,000,000 gives `ERR_TOO_LONG`. A product of the
  resamples and the sample size above 500,000,000 gives the same code.
  SciPy has no such bound. A caller with a larger job splits it.
- `goodness_of_fit` refits and sorts every resample. The product of the
  resamples and the sample size stops at 2,000,000, or the call gives
  `ERR_TOO_LONG`. The other resampling calls stop at 500,000,000.
- The resampling calls draw from a generator that the caller passes.
  The package never makes its own seed. One seed gives one answer on
  every run and on both Goose releases.
- A block that returns NaN or an infinity gives `ERR_NAN` or
  `ERR_INFINITY`. The call never turns such a value into a p-value.
- A confidence level outside 0 to 1 gives `ERR_DOMAIN`. SciPy does not
  check it in the resampling calls.
- `sobol_indices` refuses more than 32 dimensions and more than
  1,048,576 points with `ERR_TOO_LONG`.
- Each exact test has a bound. The table total is at most 400. The
  nuisance grid holds at most 1,024 points. The Page test takes at most
  12 columns. The Ansari table holds at most 400,000 cells.
  A larger input gives `ERR_TOO_LONG` or `ERR_OVERFLOW`.
- SciPy enumerates the exact tests without a bound. A caller with a
  larger table uses the asymptotic form.

- `tukey_hsd` and `dunnett` refuse more than 64 groups with
  `ERR_TOO_LONG`. A group of one value gives `ERR_TOO_SHORT`.
- A search that does not end gives `ERR_FIT`. The fits stop after 200
  steps, and the interval search of `boxcox_interval` stops after 500.
- A result of any new call that would be NaN or infinite gives a fault
  code. A caller who ignores the code never sees an invalid value with
  the code `OK`. SciPy returns NaN or an infinity there.
- `bws_test` does not copy the unseeded Monte Carlo of SciPy. The
  caller passes the generator. The call draws the split with its own
  draw, which rejects the low values of the word. The shuffle of the
  standard library draws the low values more often on the release
  v2026.3. The two releases would give different answers for one seed.
- `combine_pvalues` rejects a list of weights that holds only zeros.
  SciPy returns a value that is not a number there.
- `combine_pvalues` scales the weights by the largest one before it
  squares them. A weight near 1e-160 keeps its share of the evidence.
  SciPy loses the digits of such a weight.
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
- The moments of a family of a power check the order of the moment
  against the shape of the law. The law of Burr has no mean below a
  shape of 1. The law of the generalized Pareto has no mean below a
  shape of 1. The package gives `ERR_UNDEFINED` there.
  SciPy returns NaN or an infinity.
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
- The Rice law gives its mean and its variance to full digits. Its
  skewness and kurtosis lose digits for a shape of 10 or more, because
  the central moments cancel. At the shape 30 the kurtosis has a relative
  error of about 2e-3. SciPy 1.18.1 loses digits there too. The shape must
  be at most 3000. A larger shape gives the domain fault.
- The families of a truncation check their bounds before any work. The
  bounds of the truncated normal law must lie in order. The end of the
  support of the truncated Pareto law has to be above 1. The package
  gives `ERR_DOMAIN` there.
- The families that come from the gamma law check their shapes before
  any work. A shape at or below 0 gives `ERR_DOMAIN`. The shape `c` of
  the generalized gamma law must not be 0. The package gives
  `ERR_DOMAIN` there. SciPy returns NaN.
- The families of the Johnson kind check the shape `b` before any work.
  A shape at or below 0 gives `ERR_DOMAIN`. The shape `a` may be any
  finite value.
- The moments of the new families check the order of the moment against
  the shape of the law. The check runs before any division. The law of
  Mielke has no third moment when its shape `s` is 3 or below. The law
  of the inverse gamma law has no mean when its shape is 1 or below.
- The skew t family of the Johnson kind has no moment when a shape is
  not above half the order. The package gives `ERR_UNDEFINED` there.
  SciPy returns an infinity or a NaN.

These risks stay.

- The package does not check the meaning of the values. A value that
  holds a price, a weight or a count is a number. The package leaves the
  choice to the caller.
- The package does not check the units. A sample of lengths and a sample
  of weights give the same mean.
- The tests compare the values with the reference libraries within a
  tolerance in units in the last place. A tolerance of 8 units in the
  last place is not an exact match.

- The binned calls bound their work before any table is built. A bin
  count above `MAX_BINS` of 1,048,576 gives `ERR_TOO_LONG`. A count of
  cells above `MAX_CELLS` of 4,194,304 gives the same code. A count of
  dimensions above `MAX_DIMS` of 16 gives the same code. The check
  comes before the allocation of the cells.
- The `ecdf_censored` and the `ecdf_censored_interval` calls sort the
  times of the events and of the censored values. The sort is the sort
  of the standard library, so the work grows with n log n. A sample
  above `MAX_INPUT` gives `ERR_TOO_LONG` before the sort.
- The `expectile` call searches on the range of the sample. The search
  takes at most `MAX_EXPECTILE_STEPS` of 200 steps. A search that does
  not end gives `ERR_FIT`. A caller with a hard bound checks the level
  before the call.
- The `entropy` call rejects a negative probability with
  `ERR_DOMAIN`. SciPy gives -inf there. A probability of zero is
  allowed, as in SciPy.
- The `relative_entropy` call gives `ERR_UNDEFINED` when a positive
  probability has no match in the second set. SciPy gives inf there. The call also gives `ERR_LENGTH` for vectors
  with another length. SciPy uses the shorter length and gives a silent
  value.
- The `differential_entropy` call gives `ERR_UNDEFINED` for a sample
  with ties. SciPy gives NaN or an infinity there.
- The `tmin` and the `tmax` calls give `ERR_UNDEFINED` when no value
  of the sample lies inside the limits. SciPy raises there.
- The `trimboth` and the `trim1` calls give `ERR_DOMAIN` for a
  proportion that SciPy refuses. The cut count is the count of the
  values at each end. A caller cannot lose more than the share that the
  proportion names.
- The `circmean` call gives `ERR_UNDEFINED` when the resultant of the
  sample is zero. SciPy gives NaN there.
- The `directional_stats` call gives `ERR_UNDEFINED` for a row of
  length zero and for a mean of length zero. SciPy gives NaN there.
- The `describe`, `gzscore`, `zmap` and `lmoment` calls give a fault
  code where SciPy gives NaN. A sample that is too short, a constant
  sample and a zero scale each give `ERR_TOO_SHORT` or `ERR_UNDEFINED`.
  The result never holds a silent NaN.
- The association and density calls cap their work against a hostile
  input. A table holds at most `MAX_TABLE_CELLS` of 1,048,576
  cells and `MAX_TABLE_SIDE` of 1,024 rows or columns. A density holds at
  most `MAX_KDE_POINTS` of 65,536 points and `MAX_KDE_DIMS` of 16
  dimensions. Its sample points times its evaluation points stay under
  `MAX_KDE_WORK` of 67,108,864. A Siegel estimate holds at most
  `MAX_SIEGEL_POINTS` of 2,048 points. A larger input gives
  `ERR_TOO_LONG` before any work starts.
- SciPy has no such caps.
- A multiscale graph correlation holds its pairs squared times one plus
  its repetitions under `MAX_MGC_WORK`. SciPy sets no such bound.
- A kernel density of a sample with a singular covariance gives
  `ERR_DOMAIN`. SciPy raises an error from the linear algebra there.
- A result that would overflow gives `ERR_OVERFLOW`. It never gives a
  wrong number with `OK`. This holds for weighted distances, for
  Dirichlet parameters, for a kernel density and for a contingency
  table. SciPy returns a wrong value or NaN there.
- A table count above 2^53 gives `ERR_DOMAIN`, because a larger count is
  not exact. A conditional odds ratio with more than 1,048,576 counts in
  its support gives `ERR_TOO_LONG`.
- `weightedtau_ranked` checks its ranks for NaN and infinity.
- A draw from a kernel estimate carries the noise of the band width. A
  draw without the noise would repeat the sample points.
- A multiscale graph correlation divides by the two roots of the
  variances. A tiny scale and a large scale keep the statistic.
- A `goodness_of_fit` call holds at most `MAX_GOF_REPS` of 100,000
  Monte Carlo draws. The draws times the sample size stay under
  `MAX_FIT_WORK` of 2,000,000. A `kde_resample` call holds at most
  `MAX_KDE_POINTS` draws.
  The seeded stream is the one of the caller, so a run repeats.
- The fit of the chi-square, the student and the F families searches
  with the method of Nelder and Mead. The search stops after
  `FIT_NM_STEPS` of 200 steps. It stops earlier when the simplex is
  small against the scale of the parameters and of the objective.
- The degrees of freedom of a fit stay below `FIT_DF_CAP` of 1,000,000.
  The laws are flat above it. SciPy returns its own bound there. The
  fits need as many distinct values as they have parameters.
- The `ppcc_max` call searches for the shape with the best correlation.
  Its objective is flat near its minimum. The shape that the search
  stops at therefore differs between the platforms by about 1e-7 of the
  value. The call agrees with SciPy to that bound.

- The `random_correlation` call checks every eigenvalue for NaN, for an
  infinity and for a negative value before any work. A sum that differs
  from the count of the eigenvalues gives `ERR_DOMAIN`.
- The `random_table` call checks that the row sums and the column sums
  are equal. It checks that each sum is a whole number. It checks that
  each sum is at least 0. A total above 1,000,000 gives `ERR_DOMAIN`.
- The `covariance` calls check the buffer for NaN and for an infinity.
  A factor that is not lower triangular gives `ERR_DOMAIN`.
- The `histogram` calls check every bin for NaN, for an infinity and
  for a negative value. A histogram with no mass gives `ERR_DOMAIN`.
- The `dist` calls of an object check every point for NaN and for an
  infinity. A probability outside its range gives `ERR_DOMAIN`.
- The `mixture` calls check every weight for NaN, for an infinity and
  for a negative value. They also check that the weights sum to 1.
- The searches of a quantile bound their growth and their bisection. A
  count whose search cannot end gives `ERR_FIT`. SciPy searches without
  a bound. The draw of the multivariate laws checks the degrees of
  freedom of the inverse and faults when the draw fails. SciPy gives a
  degenerate matrix there.
- A mean direction of a von Mises-Fisher draw holds only finite values.
  A NaN direction gives `ERR_DOMAIN`. The value scales against the
  infinity silently elsewhere in SciPy.
- The densities of the Weibull family work in the logarithm. A point
  far out gives 0 and not a NaN. SciPy returns a NaN at the same point.
- The moments of the families of a power give `ERR_OVERFLOW` when the
  moment overflows. SciPy returns a NaN or an infinity there.

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
goose run examples/variance_tests.goose
goose run examples/transforms.goose
goose run examples/exact_tests.goose
goose run examples/resampling.goose
goose run examples/comparisons.goose
goose run examples/survival.goose
goose run examples/families.goose
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

The script `test/make_family_vectors.py` writes the vectors of every
group of families. It needs Python 3, NumPy and SciPy 1.18.1. The
moments of the families with no closed form come from an integration of
the density of the law. The script uses that route for the same reason
as the package.

The environment variable `STATISTICS_SUITES` names one group of the test
suites. The groups follow the files of the package. The older calls have
the groups `vectors`, `assoc`, `accept`, `hostile`, `special` and
`distributions`. The newer calls have the groups `hypothesis`,
`normality`, `multipletests`, `variance`, `transform`, `exact`,
`resample`, `multicomp`, `fitting` and `time`. The name `tests` runs the
whole set of the newer ones. Without the variable every suite runs.

The families have one group for each file. The groups are
`distributions_elementary`, `distributions_bounded`, `distributions_weibull` and
`distributions_exponent`. The groups `distributions_power`, `distributions_growth`,
`distributions_truncated` and `distributions_gamma` follow. The groups
`distributions_generalized`, `distributions_symmetric` and `distributions_johnson`
follow them. The last groups are `distributions_special`, `distributions_extreme`
and `distributions_skew`. The groups `distributions_discrete_a`,
`distributions_discrete_b` and `distributions_rice` follow them.

The script `test/mutate.sh` names the group that covers the file of each
probe. A probe then costs only that group. A probe that no group of its
list catches runs the full suite as a fallback. The script reports such a
probe as a failure, so a wrong group never hides a missing test.

A change can make a call loop forever, so the script stops each probe
after ten minutes. A probe that stops counts as caught.

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

The multiple-test corrections take 3 times their input.
`combine_pvalues` takes 1 time.

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
