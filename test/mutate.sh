#!/usr/bin/env bash
# Mutation probes. Each probe flips one rule in a copy of the package.
# The test program must then fail. A probe that passes means a test is
# missing. Run it from the package directory with the compiler on PATH.
# The script exits 1 when a probe survives.
set -u
here="$(cd "$(dirname "$0")/.." && pwd)"
work="${TMPDIR:-/tmp}/statistics-mutate"
# The probes check the rules, so they run a smaller hostile count.
export STATISTICS_HOSTILE_N="${STATISTICS_HOSTILE_N:-2000}"
# The round trip of the quantiles costs an inverse per point. The probes
# run it on fewer points than the test program does.
export STATISTICS_ROUND_N="${STATISTICS_ROUND_N:-200}"
# A change can make a call loop forever. The probe gives up after this
# many seconds, and a probe that stops still counts as caught.
PROBE_TIMEOUT="${PROBE_TIMEOUT:-600}"
rm -rf "$work"
mkdir -p "$work"

# file | text in the code | replacement | name of the rule
probes=$(cat <<'PROBES'
exactsum.goose|if lo != 0.0 {|if false {|exact sum drops the rounding errors
exactsum.goose|if y == x - hi { hi = x; }|if false { hi = x; }|exact sum rounds half to even across partials
exactsum.goose|if m >= SCALE_LIMIT { return SCALE_DOWN; }|if m >= 1e308 { return SCALE_DOWN; }|scale for large values
exactsum.goose|if same { return xs[0], 0; }|if false { return xs[0], 0; }|mean of equal values is exact
spread.goose|var s = acc_total(s2) - a * a / n;|var s = acc_total(s2);|variance correction term
spread.goose|return s / (n - ddof), t, 0;|return s / n, t, 0;|sample variance divisor
spread.goose|let r = sqrt(v) * t;|let r = sqrt(v * t * t);|deviation from the scaled root
spread.goose|acc_add(a, log(x / pow(2.0, e)));|acc_add(a, log(x));|geometric mean splits the exponent
spread.goose|if x == 0.0 { return 0.0, 0; }|if x == 0.0 { return nan(), 5; }|harmonic mean of a zero
spread.goose|if x < 0.0 { return nan(), 5; }|if x < -1.0 { return nan(), 5; }|harmonic mean of a negative value
order.goose|sort(out) { a, b => a < b };|sort(out) { a, b => a > b };|sort order
order.goose|return mid(s[n / 2 - 1], s[n / 2]), 0;|return s[n / 2], 0;|median of an even sample
order.goose|if count > best {|if count >= best {|mode tie goes to the first value
order.goose|let r = (v - width / 2.0) + width * (n / 2.0 - below) / equal;|let r = (v - width / 2.0) + width * (n / 2.0 - below) / (equal + 1);|grouped median
order.goose|if method == 6 {|if method == 99 {|linear method
order.goose|if method == 3 { a = 0.0; b = 1.0; }|if method == 3 { a = 0.0; b = 0.5; }|interpolated inverted cdf constants
order.goose|if method == 4 { a = 0.5; b = 0.5; }|if method == 4 { a = 0.4; b = 0.5; }|hazen constants
order.goose|if method == 5 { a = 0.0; b = 0.0; }|if method == 5 { a = 0.0; b = 0.1; }|weibull constants
order.goose|if method == 7 { a = 1.0 / 3.0; b = 1.0 / 3.0; }|if method == 7 { a = 0.3; b = 1.0 / 3.0; }|median unbiased constants
order.goose|if method == 8 { a = 0.375; b = 0.375; }|if method == 8 { a = 0.375; b = 0.4; }|normal unbiased constants
order.goose|if idx - p == 0.0 { return s[clamp_index(p, n)]; }|if idx - p == 1.0 { return s[clamp_index(p, n)]; }|inverted cdf on a whole index
order.goose|if v - floor(v) == 0.0 { g = 0.5; } else { g = 1.0; }|if v - floor(v) == 0.0 { g = 0.4; } else { g = 1.0; }|averaged inverted cdf weight
order.goose|if idx - p == 0.0 \&\& (p as i64) % 2 == 1|if idx - p == 0.0 \&\& (p as i64) % 2 == 0|closest observation parity
order.goose|if t >= 0.5 { return b - d * (1.0 - t); }|if t >= 2.0 { return b - d * (1.0 - t); }|interpolation exact at the upper end
order.goose|out.push(quantile_sorted(s[..], (i as f64) / groups, method));|out.push(quantile_sorted(s[..], (i as f64) / (groups + 1), method));|cut point fractions
shape.goose|return finish(g1 * sqrt((n * (n - 1)) as f64) / (n - 2));|return finish(g1);|unbiased skewness factor
shape.goose|let g2 = pearson - 3.0;|let g2 = pearson - 2.0;|Fisher kurtosis constant
shape.goose|total += binom * p * r.m[j];|total += binom * p * r.m[j] * 1.0000001;|central moment shift to the true mean
shape.goose|if order == 1 { return 0.0, 0; }|if order == 1 { return 1.0, 0; }|first central moment
weighted.goose|let cut = floor(xs.len * proportion) as i64;|let cut = ceil(xs.len * proportion) as i64;|trim cut size
weighted.goose|if x >= lower \&\& x <= upper|if x > lower \&\& x <= upper|limits are inclusive
weighted.goose|let denom = total - 1.0;|let denom = total;|frequency weight divisor
weighted.goose|let denom = wsum - w2sum / wsum;|let denom = wsum;|reliability weight divisor
weighted.goose|if p.w > 0.0 \&\& cum >= half|if p.w > 0.0 \&\& cum > half|weighted median at the half way weight
check.goose|if xs.len > MAX_INPUT { return 9; }|if xs.len > MAX_INPUT + 1 { return 9; }|input length bound
check.goose|if xs.len < need { return 2; }|if xs.len < 0 { return 2; }|sample too short
check.goose|if w < 0.0 { return 5; }|if w < -1.0 { return 5; }|negative weight
check.goose|if ws.len != n { return 6; }|if false { return 6; }|weight length
cross.goose|let e = ((ah * bh - p) + ah * bl + al * bh) + al * bl;|let e = 0.0;|two product drops the error term
cross.goose|if n_xx == n_yy { den = n_xx; } else { den = sqrt(n_xx) * sqrt(n_yy); }|var den: f64 = sqrt(n_xx) * sqrt(n_yy);|correlation shares the denominator of equal variances
cross.goose|if v > 1.0 { return 1.0, 0; }|if v > 2.0 { return 1.0, 0; }|correlation clamps to one
cross.goose|if e > 400.0 { return 400 - e as i64; }|if e > 1000.0 { return 400 - e as i64; }|column scale bounds the column magnitude
cross.goose|if form == WEIGHT_FREQUENCY { return r.w * (r.w - 1.0 / r.top); }|if form == WEIGHT_FREQUENCY { return r.w * r.w - 1.0 / r.top; }|reliability weight divisor
cross.goose|if form == WEIGHT_POPULATION { return r.w * r.w; }|if form == WEIGHT_POPULATION { return r.w * (r.w - 1.0 / r.top); }|frequency weight divisor
ranks.goose|stable_sort(idx) { a, b => a.v < b.v };|sort(idx) { a, b => a.v < b.v };|rank sort is stable
ranks.goose|if method == TIE_DENSE { rank = dense + 1; }|if method == TIE_DENSE { rank = first; }|dense rank numbers groups
ranks.goose|if method == TIE_ORDINAL { rank = first + (m - k); }|if method == TIE_ORDINAL { rank = first; }|ordinal rank follows the input order
ranks.goose|if method == TIE_MAX { rank = first + (count - 1); }|if method == TIE_MAX { rank = first; }|max rank ends at the group end
ranks.goose|if a[i] <= a[j] {|if a[i] < a[j] {|merge counts no swap for a tie in y
ranks.goose|return n0 - t.n1 - t.n2 + t.n3 - 2 * swaps, t;|return n0 - t.n1 - t.n2 + t.n3, t;|kendall counts the discordant pairs
ranks.goose|let den = nf * nf * (mf - 1.0) / mf;|let den = nf * nf * (mf - 1.0) / (mf - 1.0);|tau c denominator scale
line.goose|let slope = scale_pow2(r.n_xy / r.n_xx, r.exp_x - r.exp_y);|let slope = scale_pow2(r.s_xy / r.s_xx, r.exp_x - r.exp_y);|line slope uses the centered sums
line.goose|let slope = scale_pow2(r.s_xy / r.s_xx, r.exp_x - r.exp_y);|let slope = scale_pow2(r.n_xy / r.n_xx, r.exp_x - r.exp_y);|proportional slope goes through the origin
line.goose|let intercept_err = scale_pow2(se * sqrt(r.s_xx / n), -r.exp_y);|let intercept_err = scale_pow2(se, -r.exp_y);|intercept error carries the sum of squares
line.goose|if xs.len > MAX_THEIL { return nan(), nan(), 9; }|if xs.len > MAX_THEIL * 2 { return nan(), nan(), 9; }|theil sen pair bound
line.goose|let intercept = my - slope * mx;|let intercept = my - slope * mx + 1.0;|theil sen intercept uses the median of x
matrix.goose|if rows < 2 { return 2; }|if rows < 1 { return 2; }|matrix needs two rows
special_erf.goose|return 1.0 - erf_series(x), OK;|return 1.0 - erf_series(abs(x)), OK;|erfc keeps the sign near zero
special_erf.goose|if x > 27.25 { return 0.0, OK; }|if x > 26.0 { return 0.0, OK; }|erfc underflow bound
special_erf.goose|if p <= 0.5 { return erfcinv_core(p), OK; }|if p <= 0.0 { return erfcinv_core(p), OK; }|erfcinv split of the range
special_erf.goose|if a == 1.0 { return nan(), ERR_OVERFLOW; }|if a == 1.0 { return 1.0, OK; }|erfinv at one
special_gamma.goose|if x > 171.6243769563027 { return nan(), 8; }|if x > 170.0 { return nan(), 8; }|gamma overflow bound
special_gamma.goose|let sg = gamma_sign(n);|let sg = 1.0;|sign of the gamma and the log gamma on the negative side
special_erf.goose|if x == -1.0 { return nan(), ERR_OVERFLOW; }|if x == -1.0 { return nan(), ERR_DOMAIN; }|log1p at minus one
special_erf.goose|if x < -40.0 { return -1.0, OK; }|if x < -20.0 { return -1.0, OK; }|expm1 saturates at minus one
special_incomplete.goose|if a > MAX_SHAPE|if a > 2.0 * MAX_SHAPE|shape bound
special_incomplete.goose|if delta > 2.5 {|if delta > -0.9 {|fraction start of the large shapes
special_incomplete.goose|if delta < -0.5 {|if delta < 5.0 {|series start of the large shapes
special_incomplete.goose|if a >= 2500.0 {|if a >= 250000.0 {|saddle point branch start
special_gamma.goose|if hi.hi < 1e8 {|if hi.hi < 1e300 {|log beta of the large shapes
special_gamma.goose|if lo.hi >= 20.0 {|if lo.hi >= -1.0 {|log beta of two large shapes
special_inverse.goose|if target > 1e-290 {|if false {|polish of the gamma inverse
special_inverse.goose|var xn = xr - (v - target) / dens;|var xn = xr + (v - target) / dens;|direction of the polish step
special_incomplete.goose|if a >= 1e100 { return dd_from(0.0); }|if a >= 1e3 { return dd_from(0.0); }|Stirling sum of a huge argument
special_gamma.goose|let n = (floor(STIRLING_SHIFT - z.hi) + 1.0) as i64;|let n = floor(STIRLING_SHIFT - z.hi) as i64;|shift of the log gamma reaches the Stirling range
special_inverse.goose|if p == 1.0 { return 1.0, 0; }|if p == 1.0 { return 0.5, 0; }|betaincinv at one
special_incomplete.goose|        if upper { return 1.0, 0; }|        if upper { return 0.5, 0; }|betainc at zero
special_incomplete.goose|        if upper { return 0.0, 0; }|        if upper { return 0.5, 0; }|betainc at one
special_gamma.goose|if is_nan(x) { return nan(), 3; }|if is_nan(x) { return 0.0, 3; }|betainc with a NaN
special.goose|if is_inf(p) { return ERR_DOMAIN; }|if is_inf(p) { return ERR_INFINITY; }|probability that is infinite
special_incomplete.goose|while k <= MAX_ITER_SERIES {|while k <= 3 {|bound of the power series
special_incomplete.goose|while i <= MAX_ITER_CF {|while i <= 3 {|bound of the gamma fraction
special_incomplete.goose|while m <= MAX_ITER_CF {|while m <= 3 {|bound of the beta fraction
special_incomplete.goose|let w = sqrt(2.0 * bigd);|let w = sqrt(2.0 * bigd) * 1.0000001;|saddle point variable
special_incomplete.goose|ssum += bn * y[0] * e;|ssum += bn * y[0];|saddle point leading term
special_incomplete.goose|if !upper \&\& n % 2 == 0 { bn = -bn; }|if upper \&\& n % 2 == 0 { bn = -bn; }|saddle point series sign of the lower tail
special_incomplete.goose|if upper \&\& a < 1e-3 {|if upper \&\& a < 1e-300 {|branch of the tiny shape
special_incomplete.goose|if abs(u.hi) >= 1e-3 { return 0.0, false; }|if abs(u.hi) >= 1e-300 { return 0.0, false; }|guard of the tiny shape
special_incomplete.goose|return dd_mul_d(s, a);|return dd_mul_d(s, 1.0);|series of the log gamma of a tiny shape
special_incomplete.goose|while k <= 100 {|while k <= 1 {|series of the tiny shape tail
special_incomplete.goose|if direct > 1.0 - 1e-15 {|if direct > 2.0 {|betainc upper side of a tiny shape
special_incomplete.goose|if abs(part.hi) <= 1e-33 * abs(ssum.hi) { done = true; break; }|if abs(part.hi) <= 1e300 * abs(ssum.hi) { done = true; break; }|stop of the beta series
special_incomplete.goose|while n <= MAX_ITER_SERIES {|while n <= 3 {|bound of the beta series
distributions.goose|if z < 0.7 {|if z < 0.0 {|branch of the normal tail
distributions.goose|return 0.5 * q * expf_dd(dd_mul_d(dd_mul(zd, zd), -0.5));|return 0.5 * q * exp(-0.5 * z * z);|square of the normal argument
distributions.goose|return norm_upper(dd_neg(norm_zdd(x, mu, sigma))), OK;|return 1.0 - norm_upper(norm_zdd(x, mu, sigma)), OK;|the normal cumulative function is the tail itself
distributions.goose|return norm_upper(norm_zdd(x, mu, sigma)), OK;|return 1.0 - norm_upper(dd_neg(norm_zdd(x, mu, sigma))), OK;|the normal survival function is the tail itself
distributions.goose|let ph, pl = two_prod(p, wh);|let ph, pl = two_prod(p, 1.0);|the uniform quantile uses the width
distributions.goose|return s + (e + pl + p * wl), OK;|return s + e, OK;|the uniform quantile keeps the error term
distributions.goose|return expf_dd(dd_neg(expon_arg(rate, x))), OK;|return exp(-rate * x), OK;|the exponential survival uses the pair
distributions.goose distributions_discrete.goose distributions_shape.goose distributions_t.goose|return exp_or_overflow(l);|return exp(l.hi), OK;|the density of a family keeps its digits
distributions.goose|let ph, pl = two_prod(sigma, zh);|let ph, pl = two_prod(sigma, 1.0);|the lognormal quantile uses the scale
distributions_shape.goose|var m3 = dd_sub(g3, dd_mul_d(dd_mul(g2, g1), 3.0));|var m3 = dd_sub(g3, dd_mul_d(dd_mul(g2, g1), 2.0));|third central moment of the Weibull
distributions_shape.goose|m4 = dd_add(m4, dd_mul_d(dd_mul(g2, g11), 6.0));|m4 = dd_add(m4, dd_mul_d(dd_mul(g2, g11), 3.0));|fourth central moment of the Weibull
distributions_shape.goose|return g1, g2, dd_sub(g2, dd_mul(g1, g1)), OK;|return g1, g2, dd_sub(g2, dd_mul(g1, dd_from(1.0))), OK;|variance of the Weibull
distributions_shape.goose|return betainc_upper(a, b, x);|return betainc(b, a, 1.0 - x);|the beta survival uses the upper tail
distributions_t.goose|betainc_upper(0.5, 0.5 * df, y.hi)|betainc(0.5 * df, 0.5, y.hi)|the Student t tail uses the upper beta
distributions_t.goose|if x >= 0.0 { return 1.0 - tail, tail, OK; }|if x >= 0.0 { return tail, 1.0 - tail, OK; }|the Student t tail at a positive value
distributions_shape.goose|let r = dd_to_f64(dd_div(m3, dd_mul(v, dd_sqrt(v))));|let r = dd_to_f64(dd_div(m3, dd_mul(v, dd_from(1.0))));|skewness of the Weibull
distributions_shape.goose|let r = dd_to_f64(dd_div(dd_sub(m4, dd_mul_d(v2, 3.0)), v2));|let r = dd_to_f64(dd_div(m4, v2));|excess kurtosis of the Weibull
distributions_discrete.goose|return betainc_upper(kf + 1.0, n - kf, p);|return betainc(n - kf, kf + 1.0, 1.0 - p);|the binomial cumulative function uses the upper beta
distributions_discrete.goose|return betainc(kf + 1.0, n - kf, p);|return betainc_upper(n - kf, kf + 1.0, 1.0 - p);|the binomial survival uses the upper beta
distributions_discrete.goose|return gammaincc(floor(k) + 1.0, rate);|return gammainc(floor(k) + 1.0, rate);|the Poisson cumulative function uses the upper gamma
distributions_discrete.goose|return gammainc(floor(k) + 1.0, rate);|return gammaincc(floor(k) + 1.0, rate);|the Poisson survival uses the lower gamma
distributions_discrete.goose|if k < 1.0 { return 0.0, OK; }|if k < 0.0 { return 0.0, OK; }|the geometric support starts at 1
distributions_discrete.goose|return betainc(r, floor(k) + 1.0, p);|return betainc_upper(r, floor(k) + 1.0, p);|the negative binomial cumulative function uses the lower beta
distributions_discrete.goose|var r = dd_mul_d(dd_mul_d(t, n - k), d - k);|var r = dd_mul_d(dd_mul_d(t, n - k), d - k + 1.0);|ratio of the hypergeometric masses
distributions_discrete.goose|while j >= kmin {|while j > kmin {|the hypergeometric sum starts at the mode
distributions_discrete.goose|if v >= p { r = mid; } else { l = mid + 1.0; }|if v >= p { r = mid + 1.0; } else { l = mid + 1.0; }|the search of a discrete quantile
distributions_discrete.goose|if is_inf(hi) {|if false {|the search doubles its bound
distributions_discrete.goose|floor(k) != k { return ERR_DOMAIN; }|false { return ERR_DOMAIN; }|a count with a fraction is out of range
distributions_t.goose|if d2 <= 4.0 { return nan(), ERR_UNDEFINED; }|if d2 <= 2.0 { return nan(), ERR_UNDEFINED; }|the variance of the F family needs df2 above 4
distributions_t.goose|if df <= 4.0 { return nan(), ERR_UNDEFINED; }|if df <= 2.0 { return nan(), ERR_UNDEFINED; }|the kurtosis of the Student t needs df above 4
distributions_discrete.goose|return (2.0 - p) / sqrt(1.0 - p), OK;|return (2.0 - p) / sqrt(p), OK;|skewness of the geometric family
distributions_discrete.goose|return 6.0 / r + p * p / (r * (1.0 - p)), OK;|return 6.0 / r + p * p / (1.0 - p), OK;|kurtosis of the negative binomial
distributions_discrete.goose|let den = dd_mul_d(dd_mul_d(dd_from(m), m), m - 1.0);|let den = dd_mul_d(dd_mul_d(dd_from(m), m), m);|variance of the hypergeometric
distributions_discrete.goose|if kmax - kmin + 1.0 > MAX_DIST_SUPPORT {|if kmax - kmin + 1.0 > 1e18 {|bound of the hypergeometric support
distributions_shape.goose|if z.hi > 170.0 { return dd_from(nan()), ERR_OVERFLOW; }|if z.hi > 17.0 { return dd_from(nan()), ERR_OVERFLOW; }|bound of the Weibull moment
distributions_t.goose|betainc_upper(0.5 * d2, 0.5 * d1, w.hi)|betainc(0.5 * d1, 0.5 * d2, u.hi)|the F survival uses the upper beta
distributions.goose|if z > 38.6 { return 0.0; }|if z > 37.0 { return 0.0; }|bound of the normal tail
distributions_shape.goose|return base + sign * pdf * u.lo;|return base + sign * pdf * u.hi;|the rounding step of a beta ratio uses the low part
distributions_discrete.goose|if p == 0.0 { return lo, OK; }|if p == 0.0 { return lo + 1.0, OK; }|the probability 0 gives the smallest count
distributions_discrete.goose|return n * d / m, OK;|return n * n / m, OK;|mean of the hypergeometric
distributions_t.goose|if p <= 0.5 {|if p <= 0.0 {|the F quantile uses the upper tail above the middle
@tests_multiplicity.goose|let v = ps[order[k]] * (nf - k);|let v = ps[order[k]] * nf;|the Holm multiplier
@tests_multiplicity.goose|let frac = ps[order[k]] * nf / (k + 1) * factor;|let frac = ps[order[k]] * nf / (k + 2) * factor;|the Benjamini-Hochberg divisor
@tests_multiplicity.goose|for k in 0..n { out[k] = ps[k] * nf; }|for k in 0..n { out[k] = ps[k]; }|the Bonferroni multiplier
@statistics.goose|const MAX_EXACT_MW: i64 = 8;|const MAX_EXACT_MW: i64 = 7;|the exact Mann-Whitney size bound
@statistics.goose|const MAX_MW_CELLS: i64 = 100000;|const MAX_MW_CELLS: i64 = 1000;|the cell bound of the Mann-Whitney table
@statistics.goose|const MAX_EXACT_WILCOXON: i64 = 50;|const MAX_EXACT_WILCOXON: i64 = 49;|the exact Wilcoxon pair bound
@statistics.goose|const MAX_ENUM_WILCOXON: i64 = 13;|const MAX_ENUM_WILCOXON: i64 = 12;|the sign pattern bound of the Wilcoxon test
@statistics.goose|const MAX_EXACT_KENDALL: i64 = 33;|const MAX_EXACT_KENDALL: i64 = 32;|the exact Kendall pair bound
@statistics.goose|const MAX_EXACT_KS: i64 = 10000;|const MAX_EXACT_KS: i64 = 9999;|the exact Kolmogorov-Smirnov product bound
@tests.goose|alt > 2 { return ERR_DOMAIN; }|alt > 3 { return ERR_DOMAIN; }|an unknown alternative gives the domain fault
@tests.goose|return 2.0 * p, c;|return p, c;|the two-sided t p-value doubles the tail
@tests.goose|df = 1.0 / (ra * ra / (n1 - 1.0) + rb * rb / (n2 - 1.0));|df = 1.0 / (ra * ra / (n1 - 1.0) + rb * rb / n2);|Welch degrees of freedom
@tests.goose|let mag = min(0.5, abs(diff));|let mag = min(0.4, abs(diff));|Yates correction size
@tests.goose|let target = d * (1.0 + 1e-7);|let target = d * (1.0 + 1e-3);|binomial two-sided tolerance
@tests.goose|let gamma = 1.0 + 1e-14;|let gamma = 1.0 + 1e-2;|Fisher two-sided tolerance
@tests_wilcoxon.goose|if zero_method == ZERO_ZSPLIT {|if zero_method == ZERO_PRATT {|the zsplit zero method
@tests_wilcoxon.goose|if zero_method == ZERO_PRATT {|if zero_method == ZERO_ZSPLIT {|the pratt zero method
@tests_ks.goose|let h = round(d * lcm) as i64;|let h = floor(d * lcm) as i64;|the lattice height rounds
@tests_ks.goose|return 2.0 * p;|return p;|the square method doubles the tail
@tests_ks.goose|if product <= MAX_EXACT_KS \&\& !tied {|if product <= MAX_EXACT_KS {|a tie in a Kolmogorov-Smirnov sample gives the normal form
@tests_ks.goose|let up = (k + 1) as f64 / n - f;|let up = k as f64 / n - f;|the upper step of the one-sample statistic
@tests_ks.goose|let down = f - k as f64 / n;|let down = f - (k + 1) as f64 / n;|the lower step of the one-sample statistic
@tests_ks.goose|if alt == ALT_GREATER { d = dplus; }|if alt == ALT_GREATER { d = dminus; }|the one-sided two-sample statistic
@tests_correlation.goose|a.push(sqrt(0.5));|a.push(sqrt(0.4));|the Shapiro-Wilk coefficient of three values
@tests_correlation.goose|norm_std_ppf((i - 0.375) / an25)|norm_std_ppf((i - 0.4) / an25)|the Shapiro-Wilk normal scores
@tests_correlation.goose|if alt == ALT_TWO_SIDED { x = 0.5 * (1.0 - abs(r)); }|if alt == ALT_TWO_SIDED { x = 0.4 * (1.0 - abs(r)); }|the Pearson p-value argument
@tests_correlation.goose|if n == 2 { return good_result(r, 1.0, df, nan()); }|if n == 2 { return good_result(r, 0.5, df, nan()); }|the Pearson test of two points
@tests_correlation.goose|if n <= MAX_EXACT_KENDALL|if n <= 20|the exact Kendall bound in the call
@normality.goose|return exp(1.2937 - 5.709 * adj + 0.0186 * adj * adj);|return exp(1.2937 - 5.709 * adj + 0.0186 * adj * adj * 1.01);|the top branch of the Anderson-Darling p-value
@normality.goose|return exp(0.9177 - 4.279 * adj - 1.38 * adj * adj);|return exp(0.9177 - 4.379 * adj - 1.38 * adj * adj);|the middle branch of the Anderson-Darling p-value
@normality.goose|return 1.0 - exp(-8.318 + 42.796 * adj - 59.938 * adj * adj);|return 1.0 - exp(-8.418 + 42.796 * adj - 59.938 * adj * adj);|the lower middle branch of the Anderson-Darling p-value
@normality.goose|return 1.0 - exp(-13.436 + 101.14 * adj - 223.73 * adj * adj);|return 1.0 - exp(-13.436 + 101.5 * adj - 223.73 * adj * adj);|the bottom branch of the Anderson-Darling p-value
@normality.goose|let adj = a2 * (1.0 + 0.75 / n + 2.25 / (n * n));|let adj = a2 * (1.0 + 0.75 / n);|the sample size adjustment of the Anderson-Darling p-value
@normality.goose|let w2 = -1.0 + sqrt(2.0 * (beta2 - 1.0));|let w2 = -1.0 + sqrt(beta2 - 1.0);|the skewness z transform
@normality.goose|var y = b1 * sqrt(((n + 1.0) * (n + 3.0)) / (6.0 * (n - 2.0)));|var y = b1 * sqrt(((n + 1.0) * (n + 3.0)) / (5.0 * (n - 2.0)));|the skewness z scaling
@normality.goose|let sqrtbeta1 = 6.0 * (n * n - 5.0 * n + 2.0) / ((n + 7.0) * (n + 9.0))|let sqrtbeta1 = 6.0 * (n * n - 5.0 * n + 3.0) / ((n + 7.0) * (n + 9.0))|the kurtosis z skewness link
@normality.goose|let q = (1.0 - 2.0 / a) / denom;|let q = (1.0 - 3.0 / a) / denom;|the kurtosis z cube root term
@normality.goose|return 2.0 * upper, 0;|return 4.0 * upper, 0;|the two-sided p-value doubles one tail
@normality.goose|let k2 = zs * zs + zk * zk;|let k2 = zs * zs - zk * zk;|the K-squared statistic adds the squares
@normality.goose|return chi_result(k2, 2.0);|return chi_result(k2, 3.0);|the K-squared chi-square tail
@normality.goose|let jb = n / 6.0 * (b1 * b1 + b2 * b2 / 4.0);|let jb = n / 6.0 * (b1 * b1 + b2 * b2 / 3.0);|the Jarque-Bera statistic
@normality.goose|return chi_result(jb, 2.0);|return chi_result(jb, 3.0);|the Jarque-Bera chi-square tail
@normality.goose|let scale2, cs = variance_core(xs, 1);|let scale2, cs = variance_core(xs, 0);|the Anderson-Darling scale divides by n minus one
@normality.goose|let w: f64 = 2 * (i + 1) - 1;|let w: f64 = 2 * (i + 1);|the Anderson-Darling sum weight
@normality.goose|let c = check_sample(xs, 8);|let c = check_sample(xs, 5);|the skewness test and the omnibus test need 8 values
@normality.goose|let c = check_sample(xs, 5);|let c = check_sample(xs, 4);|the kurtosis test needs 5 values and the Anderson-Darling test too
@normality.goose|let c = check_sample(xs, 2);|let c = check_sample(xs, 1);|the Jarque-Bera test needs 2 values
@normality.goose|if scale2 <= 0.0 { return fault_result(ERR_UNDEFINED); }|if scale2 <= 0.0 { return fault_result(ERR_DOMAIN); }|the Anderson-Darling test faults on a constant sample
@normality.goose|if y == 0.0 { y = 1.0; }|if y == 0.0 { y = 0.0; }|the skewness z replaces a zero input by one
@normality.goose|if q < 0.0 { term2 = -term2; }|if q < 0.0 { term2 = term2; }|the kurtosis z restores the sign of the cube root
@normality.goose|sum += w * (log_tail(z1, false) + log_tail(z2, true));|sum += w * (log_tail(z1, true) + log_tail(z2, true));|the Anderson-Darling lower tail is the lower tail
@normality.goose|if adj >= 153.4 { return 0.0; }|if adj >= 402.0 { return 0.0; }|the p-value rule clamps above the turn of its quadratic
@normality.goose|return log(0.5 * scaled_q(a * INV_SQRT2)) - 0.5 * a * a;|return log(0.5 * scaled_q(a * INV_SQRT2)) - 0.4 * a * a;|the logarithm of the scaled far tail
@normality.goose|if upper != (z > 0.0) { return 0.0; }|if upper == (z > 0.0) { return 0.0; }|the far tail on the other side is one

@tests_multiplicity.goose|method > ADJUST_NONE |method > ADJUST_HOMMEL |the bound of the adjustment methods
@tests_multiplicity.goose|method < 0 { return 0, ERR_DOMAIN; }|method < -1 { return 0, ERR_DOMAIN; }|the lower bound of the adjustment methods
@tests_multiplicity.goose|if method == ADJUST_HOMMEL && ps.len > HOMMEL_MAX { return 0, ERR_TOO_LONG; }|if method == ADJUST_HOMMEL && ps.len > 100000 { return 0, ERR_TOO_LONG; }|the bound of the Hommel form
@tests_multiplicity.goose|let lim = ps[order[k]] * (nf - k);|let lim = ps[order[k]] * nf;|the Hochberg multiplier
@tests_multiplicity.goose|h += 1.0 / i;|h += 1.0 / (i + 1);|the harmonic sum of the Benjamini-Yekutieli form
@tests_multiplicity.goose|let frac = ps[order[k]] * nf / (k + 1) * factor;|let frac = ps[order[k]] * nf / (k + 1);|the factor of the Benjamini-Yekutieli form
@tests_multiplicity.goose|for k in 0..n { out[k] = ps[k]; }|for k in 0..n { out[k] = ps[k] * nf; }|the identity method keeps the values
@tests_multiplicity.goose|if out[k] > 1.0 { out[k] = 1.0; }|if false { out[k] = 1.0; }|the adjusted p-value clamps to 1
@tests_multiplicity.goose|if m != ADJUST_BH && m != ADJUST_BY { m = -1; }|if m != ADJUST_BH { m = -1; }|the methods of the false discovery control
@tests_multiplicity.goose|if weights.len != ps.len { return fault_result(ERR_LENGTH); }|if weights.len == 0 { return fault_result(ERR_LENGTH); }|the length of the weights
@tests_multiplicity.goose|if w < 0.0 { return fault_result(ERR_DOMAIN); }|if w < -1.0 { return fault_result(ERR_DOMAIN); }|a negative weight
@tests_multiplicity.goose|let u = w / wmax;|let u = w;|the weights are scaled by the largest one
@tests_multiplicity.goose|sq += u * u;|sq += u;|the sum of squares of the weights
@tests_multiplicity.goose|if wmax <= 0.0 { return fault_result(ERR_DOMAIN); }|if false { return fault_result(ERR_DOMAIN); }|the weights of zero
@tests_multiplicity.goose|if method != COMBINE_STOUFFER { return fault_result(ERR_DOMAIN); }|if false { return fault_result(ERR_DOMAIN); }|the weights apply to the Stouffer method only
@tests_multiplicity.goose|if p == 0.0 { return fault_result(ERR_OVERFLOW); }|if p < 0.0 { return fault_result(ERR_OVERFLOW); }|the Fisher overflow fault
@tests_multiplicity.goose|if ps[i] == 1.0 { return fault_result(ERR_OVERFLOW); }|if ps[i] < 0.0 { return fault_result(ERR_OVERFLOW); }|the Stouffer overflow fault
@tests_multiplicity.goose|acc_add(a, log(p));|acc_add(a, p);|the logarithm of the Fisher statistic
@tests_multiplicity.goose|return chi_result(-2.0 * sum, 2.0 * k);|return chi_result(-2.0 * sum, k);|the degrees of freedom of the Fisher statistic
@tests_multiplicity.goose|let e, ce = expm1(k * l);|let e, ce = expm1(l);|the Tippett power of the count
@tests_multiplicity.goose|let e, ce = expm1(k * l);|let e, ce = expm1(k * l * 1.01);|the Tippett factor of the count
@tests_multiplicity.goose|top -= u * z;|top += u * z;|the sign of the Stouffer sum
@tests_multiplicity.goose|let zs = top / sqrt(sq);|let zs = top / sq;|the divisor of the Stouffer statistic
@tests_multiplicity.goose|if lo >= 1.0 { return good_result(lo, 1.0, nan(), nan()); }|if false { return good_result(lo, 1.0, nan(), nan()); }|the Tippett edge of one
@tests_multiplicity.goose|let v = sp[n - m + 1 + t] * m / (t + 2);|let v = sp[n - m + 1 + t] * m / (t + 1);|the Hommel divisor
@tests_multiplicity.goose|var v = sp[k] * m;|var v = sp[k];|the Hommel multiplier of the head
@tests_multiplicity.goose|if q1 < v { v = q1; }|if false { v = q1; }|the Hommel minimum of the head
@tests_multiplicity.goose|if v > pa[k] { pa[k] = v; }|if v < pa[k] { pa[k] = v; }|the Hommel running maximum
@tests_multiplicity.goose|if sp[k] > v { v = sp[k]; }|if false { v = sp[k]; }|the Hommel floor of the p-value
@check.goose|p > 1.0 { return 5; }|p > 2.0 { return 5; }|the range of the p-values
@variance.goose|let dfn = (k - 1) as f64;|let dfn = k as f64;|the degrees of freedom of the Levene test
@variance.goose|return vr_finish_chi(stat, k - 1);|return vr_finish_chi(stat, k);|the degrees of freedom of the Bartlett test
@variance.goose|acc_add(pooled, (n - 1) * v);|acc_add(pooled, n * v);|the pooled variance uses the sample divisor
@variance.goose|acc_add(recips, 1.0 / (n - 1));|acc_add(recips, 1.0 / n);|the reciprocal sum of the Bartlett test
@variance.goose|for k in 1..tot { if idx[k].v == idx[k - 1].v { tied = true; } }|for k in 1..tot { if idx[k].v == idx[k - 1].v { tied = false; } }|the tie correction of the Mood test
@variance.goose|let z = (acc_total(s) - expected) / sqrt(n1 * n2 * (n1 + n2 + 1.0) / 12.0);|let z = (acc_total(s) + expected) / sqrt(n1 * n2 * (n1 + n2 + 1.0) / 12.0);|the rank sum difference of the rank-sum test
@variance.goose|let expected = n1 * (n1 + n2 + 1.0) / 2.0;|let expected = n1 * (n1 + n2) / 2.0;|the expected rank sum of the rank-sum test
@variance_multi.goose|if method == MW_EXACT {|if false {|the exact form of the Ansari test
@transform.goose|if x <= 0.0 { return nan(), ERR_DOMAIN; }|if x < 0.0 { return nan(), ERR_DOMAIN; }|the Box-Cox form needs a positive value
@transform.goose|if abs(l) < 1e-19 { return g; }|if abs(l) < 1e-9 { return g; }|the Box-Cox form at a near zero lambda
@transform.goose|if method != NORMMAX_PEARSONR \&\& method != NORMMAX_MLE { return nan(), ERR_DOMAIN; }|if method != NORMMAX_PEARSONR \&\& method != NORMMAX_MLE { return nan(), ERR_UNDEFINED; }|an unknown method of the best fit
@transform_yeo.goose|la >= lb|la > lb|the bounds of the probability plot
@numeric.goose|if fa * fb > 0.0 { return nan(), ERR_DOMAIN; }|if fa * fb > 1.0 { return nan(), ERR_DOMAIN; }|the root finder needs a sign change
@resample.goose|INTERVAL_BCA { return rs_boot_fault(ERR_DOMAIN); }|INTERVAL_BCA { return rs_boot_fault(ERR_UNDEFINED); }|an unknown interval gives the domain fault
@resample.goose|if method.kind != RESAMPLE_BOOTSTRAP { return rs_boot_fault(ERR_DOMAIN); }|if method.kind != RESAMPLE_MONTE_CARLO { return rs_boot_fault(ERR_DOMAIN); }|the kind of the resampling method
@resample.goose|if method.interval == INTERVAL_BCA && sizes.len > 1 { return rs_boot_fault(ERR_DOMAIN); }|if false { return rs_boot_fault(ERR_DOMAIN); }|the BCa form refuses several samples
@fitting.goose|if !(m > 0.0) { return nan(), nan(), ERR_FIT; }|if !(m < 0.0) { return nan(), nan(), ERR_FIT; }|a fit needs a positive scale
@fitting.goose|if c != OK { return fault_result(c); }|if c != OK { return fault_result(ERR_FIT); }|a censored data fault keeps its own code
exact.goose|if grid < 1 { return ERR_DOMAIN; }|if grid < 1 { return ERR_UNDEFINED; }|the least count of the grid
exact_cvm.goose|if params.len != family_arity(family) { return fault_result(ERR_LENGTH); }|if params.len != family_arity(family) { return fault_result(ERR_DOMAIN); }|the arity of the family
exact_goodness.goose|if total < 4 { return fault_result(ERR_TOO_SHORT); }|if total < 4 { return fault_result(ERR_DOMAIN); }|the least count of the k-sample test
exact_page.goose|if columns < 3 { return fault_result(ERR_TOO_SHORT); }|if columns < 3 { return fault_result(ERR_DOMAIN); }|the least count of the columns of the Page test
exact_page.goose|if data.len % columns != 0 { return fault_result(ERR_LENGTH); }|if data.len % columns != 0 { return fault_result(ERR_DOMAIN); }|the length of the data of the Page test
multicomp.goose|if s < 2 { return ERR_TOO_SHORT; }|if s < 2 { return ERR_DOMAIN; }|the least count of a group
multicomp.goose|control >= k { return 0, ERR_DOMAIN; }|control >= k { return 0, ERR_UNDEFINED; }|the index of the control group
descriptive.goose|if per < 0.0 || per > 100.0 { return nan(), ERR_DOMAIN; }|if per < 0.0 || per > 101.0 { return nan(), ERR_DOMAIN; }|the bounds of the percentile
descriptive.goose|if kind == PCT_MEAN { return ((left + right)) * 50.0 / n, OK; }|if kind == PCT_MEAN { return ((left + right)) * 50.0 / (n + 1), OK; }|the mean share of a score
descriptive.goose|if kind == PCT_STRICT { return (left) * 100.0 / n, OK; }|if kind == PCT_STRICT { return (left) * 100.0 / (n + 1), OK; }|the strict share of a score
descriptive.goose|if scale <= 0.0 { return nan(), ERR_DOMAIN; }|if scale < 0.0 { return nan(), ERR_DOMAIN; }|the scale of the deviation
descriptive.goose|if p == 0.0 {|if p == 1.0 {|the exponent zero of the power mean
descriptive.goose|return m0 * (1.0 + 0.5 * p * (ln2_avg - ln_avg * ln_avg)), OK;|return m0 * (1.0 + 0.6 * p * (ln2_avg - ln_avg * ln_avg)), OK;|the linearized form of the power mean
descriptive.goose|if x <= 0.0 { return nan(), ERR_DOMAIN; }|if x < 0.0 { return nan(), ERR_DOMAIN; }|the logarithm of the geometric deviation
entropy.goose|let v = s2 * t * t / (nf - 1.0);|let v = s2 * t * t / nf;|the divisor of the second k-statistic
entropy.goose|let v = nf * s3 * t * t * t / ((nf - 1.0) * (nf - 2.0));|let v = nf * s3 * t * t * t / ((nf - 1.0) * (nf - 3.0));|the divisor of the third k-statistic
descriptive.goose|if standardize && max_order >= 3 {|if standardize && max_order >= 4 {|the standardization of the L-moments
descriptive.goose|if alpha < 0.0 || alpha > 1.0 { return nan(), ERR_DOMAIN; }|if alpha < 0.0 || alpha > 1.5 { return nan(), ERR_DOMAIN; }|the bounds of the expectile level
descriptive.goose|let denom = ex_choose(n - 1, r) * n;|let denom = ex_choose(n - 1, r) * (n + 1);|the denominator of the b coefficient
descriptive.goose|if (r - k) % 2 != 0 { return -v; }|if (r - k) % 2 != 1 { return -v; }|the sign of the L-moment coefficient
descriptive.goose|if sd == 0.0 { return 0, ERR_UNDEFINED; }|if sd < 0.0 { return 0, ERR_UNDEFINED; }|a constant reference for the z-scores
frequency.goose|let d = acc_total(acc);|let d = acc_total(acc) + 1.0;|the tie correction of the rank test
frequency.goose|probs.push(j / (n as f64));|probs.push(j / (n + 1) as f64);|the step probability of the empirical distribution
frequency.goose|return count / (xs.len as f64), OK;|return count / ((xs.len + 1) as f64), OK;|the share of the values at or below a point
frequency.goose|return ((jp1) + d) / (n), OK;|return ((jp1) + d) / (n + 1), OK;|the discontinuous estimated distribution
frequency.goose|let p = ((jp1) + delta - a) / ((n) + 1.0 - a - b);|let p = ((jp1) + delta - a) / ((n) + 2.0 - a - b);|the continuous estimated distribution
frequency.goose|let s = (b - a) / (2.0 * (nbins - 1));|let s = (b - a) / (2.0 * nbins);|the growth of the range of a histogram
binned.goose|if x >= edges[mid] { lo_i = mid + 1; } else { hi_i = mid; }|if x > edges[mid] { lo_i = mid + 1; } else { hi_i = mid; }|the bin of a value on an edge
check.goose|if level <= 0.0 || level >= 1.0 { return 5; }|if level <= 0.0 || level > 1.0 { return 5; }|the bounds of the level of the band
frequency.goose|let zse = z * sf * sqrt(cum);|let zse = z * sqrt(cum);|the width of the Greenwood band
entropy.goose|if p < 0.0 { return ERR_DOMAIN; }|if p <= 0.0 { return ERR_DOMAIN; }|a negative probability of the entropy
entropy.goose|if base == 1.0 { return ERR_DOMAIN; }|if base == 1.5 { return ERR_DOMAIN; }|the base of the entropy
entropy.goose|if order < 1 || order > 4 { return nan(), ERR_DOMAIN; }|if order < 1 || order > 5 { return nan(), ERR_DOMAIN; }|the orders of the k-statistic
entropy.goose|if order < 1 || order > 2 { return nan(), ERR_DOMAIN; }|if order < 1 || order > 3 { return nan(), ERR_DOMAIN; }|the orders of the variance of the k-statistic
entropy.goose|if m >= n { return nan(), ERR_DOMAIN; }|if false { return nan(), ERR_DOMAIN; }|a huge window never wraps the doubling
entropy.goose|let term1 = acc_total(acc) / ((n - m));|let term1 = acc_total(acc) / ((n - m + 1));|the divisor of the van Es estimator
entropy.goose|let v = num / (nf * (nf - 1.0) * (nf - 2.0) * (nf - 3.0));|let v = num / (nf * (nf - 1.0) * (nf - 2.0));|the divisor of the fourth k-statistic
trimmed.goose|if proportion < 0.0 || proportion >= 0.5 { return 0, ERR_DOMAIN; }|if proportion < 0.0 || proportion >= 0.6 { return 0, ERR_DOMAIN; }|the bound of the proportion of trimboth
trimmed.goose|if proportion < 0.0 || proportion >= 1.0 { return 0, ERR_DOMAIN; }|if proportion < 0.0 || proportion >= 1.1 { return 0, ERR_DOMAIN; }|the bound of the proportion of trim1
trimmed.goose|let cut = (proportion * n) as i64;|let cut = (proportion * n) as i64 + 1;|the cut count of a trim
trimmed.goose|if !((x > low) || (x == low && inclusive)) { continue; }|if !((x >= low) || (x == low && inclusive)) { continue; }|the lower limit of tmin
trimmed.goose|if keep.len <= ddof { return nan(), ERR_TOO_SHORT; }|if keep.len < ddof { return nan(), ERR_TOO_SHORT; }|the least count of the trimmed variance
circular.goose|if high <= low { return ERR_DOMAIN; }|if high < low { return ERR_DOMAIN; }|the range of the angles
circular.goose|let scale = TWO_PI / period;|let scale = TWO_PI / (period + 1.0);|the scale of the angles
circular.goose|return 1.0 - r, OK;|return 1.0 + r, OK;|the circular variance
circular.goose|var v = sqrt(-2.0 * log(rr));|var v = sqrt(2.0 * log(rr));|the circular deviation
circular.goose|if norm == 0.0 && normalize { return nan(), ERR_UNDEFINED; }|if norm < 0.0 && normalize { return nan(), ERR_UNDEFINED; }|a zero row of the directions
check.goose|if nbins < 1 { return ERR_DOMAIN; }|if nbins < 0 { return ERR_DOMAIN; }|the least count of the bins
binned.goose|edges.push(lo + (hi - lo) * (i) / (nbins));|edges.push(lo + (hi - lo) * (i) / (nbins + 1));|the size of a bin
binned.goose|index = index * (shape[d] + 1) + slots[d];|index = index * shape[d] + slots[d];|the row-major number of a bin
binned.goose|index = index * (shape[d] - 1) + (s - 1);|index = index * shape[d] + (s - 1);|the row-major number of a cell
binned.goose|if counts[d] == 0 {|if counts[d] == 1 {|the value of an empty cell
binned.goose|if stat == STAT_MEAN {|if stat == STAT_STD {|the statistic of the values of a bin
binned.goose|if dims > MAX_DIMS { return t, ERR_TOO_LONG; }|if dims > MAX_DIMS + 1 { return t, ERR_TOO_LONG; }|the bound of the dimensions
binned.goose|if total > MAX_CELLS { return t, ERR_TOO_LONG; }|if total > MAX_CELLS { return t, ERR_DOMAIN; }|the bound of the cells
PROBES
)

# The groups of the suites that cover one file. A probe runs the groups in
# order and stops at the first group that catches the change. The fallback
# group runs when no group of the list catches it. A probe then never
# survives on a wrong guess about the group.
groups_for() {
  case "$1" in
    distributions.goose) echo "distributions" ;;
    special.goose) echo "special distributions" ;;
    exactsum.goose) echo "accept distributions" ;;
    spread.goose) echo "accept vectors" ;;
    order.goose) echo "accept vectors" ;;
    shape.goose) echo "accept" ;;
    weighted.goose) echo "assoc accept vectors" ;;
    cross.goose) echo "assoc accept" ;;
    ranks.goose) echo "assoc accept" ;;
    line.goose) echo "assoc accept" ;;
    matrix.goose) echo "assoc accept" ;;
    tests.goose) echo "hypothesis multipletests" ;;
    statistics.goose) echo "hypothesis descriptive frequency entropy trimmed circular binned" ;;
    check.goose) echo "hypothesis vectors assoc hostile descriptive frequency entropy trimmed circular binned" ;;
    variance.goose) echo "variance" ;;
    transform.goose) echo "transform" ;;
    exact.goose) echo "exact" ;;
    resample.goose) echo "resample" ;;
    multicomp.goose) echo "multicomp" ;;
    fitting.goose) echo "fitting" ;;
    numeric.goose) echo "transform fitting" ;;
    descriptive.goose) echo "descriptive" ;;
    frequency.goose) echo "frequency" ;;
    entropy.goose) echo "entropy" ;;
    trimmed.goose) echo "trimmed" ;;
    circular.goose) echo "circular" ;;
    binned.goose) echo "binned" ;;
    distributions_discrete.goose) echo "distributions" ;;
    distributions_shape.goose) echo "distributions" ;;
    distributions_t.goose) echo "distributions" ;;
    exact_cvm.goose) echo "exact" ;;
    exact_goodness.goose) echo "exact" ;;
    exact_page.goose) echo "exact" ;;
    resample_permutation.goose) echo "resample" ;;
    resample_sobol.goose) echo "resample" ;;
    special_erf.goose) echo "special distributions" ;;
    special_gamma.goose) echo "special distributions" ;;
    special_incomplete.goose) echo "special distributions" ;;
    special_inverse.goose) echo "special distributions" ;;
    tests_correlation.goose) echo "hypothesis" ;;
    tests_ks.goose) echo "hypothesis" ;;
    tests_multiplicity.goose) echo "multipletests" ;;
    tests_wilcoxon.goose) echo "hypothesis" ;;
    transform_yeo.goose) echo "transform" ;;
    variance_multi.goose) echo "variance" ;;
    *) echo "" ;;
  esac
}

# One probe of the list. The helper writes its verdict into the work
# directory, because a parallel run cannot share a counter.
run_probe() {
  local n="$1" file="$2" from="$3" to="$4" name="$5"
  dir="$work/p$n"
  # A file name with a leading @ holds a call of the hypothesis tests.
  case "$file" in @*) file="${file#@}" ;; esac
  # A probe may name several files that hold the same rule.
  local first="${file%% *}"
  local only
  only="$(groups_for "$first")"
  if [ -z "$only" ]; then only="full"; fi
  mkdir -p "$dir"
  (cd "$here" && tar cf - --exclude=.git --exclude=mutate-work --exclude=tmp --exclude=test/oracle-rust/target --exclude=test/oracle-node/node_modules .) | (cd "$dir" && tar xf -)
  local f
  for f in $file; do
    python3 "$here/test/mutate_one.py" "$dir/$f" "$from" "$to" || {
      echo "BROKEN PROBE $n: $name"
      echo broken > "$work/verdict$n"
      rm -rf "$dir"
      return
    }
  done
  local g
  for g in $only; do
    # The word full names every suite of the package.
    local pick="$g"
    if [ "$g" = "full" ]; then pick=""; fi
    if ! (cd "$dir" && STATISTICS_SUITES="$pick" timeout "$PROBE_TIMEOUT" goose run test/run.goose >"$dir/out.txt" 2>&1); then
      # The group caught the change.
      rm -rf "$dir"
      return
    fi
  done
  if [ "$only" != "full" ]; then
    # The groups did not catch the change. The full suite may still catch it.
    if ! (cd "$dir" && STATISTICS_SUITES="" timeout "$PROBE_TIMEOUT" goose run test/run.goose >"$dir/full.txt" 2>&1); then
      echo "NEEDS THE FULL SUITE $n: $name"
      echo full > "$work/verdict$n"
      rm -rf "$dir"
      return
    fi
  fi
  echo "SURVIVED $n: $name"
  echo survived > "$work/verdict$n"
  rm -rf "$dir"
}

# The probes run at most four at a time. A copy of the package holds at
# most a few hundred megabytes, and the machine needs the rest.
jobs="${MUTATE_JOBS:-4}"
total=0

while IFS='|' read -r file from to name; do
  [ -z "$file" ] && continue
  # The variable MUTATE_ONLY keeps the probes whose line matches the pattern.
  if [ -n "${MUTATE_ONLY:-}" ]; then
    echo "$file|$from|$to|$name" | grep -Eq -- "$MUTATE_ONLY" || continue
  fi
  echo "$file|$from|$to|$name"

done <<< "$probes" > "$work/list.txt"

while IFS='|' read -r file from to name; do
  while [ "$(jobs -rp | wc -l)" -ge "$jobs" ]; do sleep 1; done
  total=$((total + 1))
  run_probe "$total" "$file" "$from" "$to" "$name" &
done < "$work/list.txt"
wait

failures=$(cat "$work"/verdict* 2>/dev/null | wc -l)
echo "$total probes, $failures failures"
[ "$failures" -eq 0 ]
