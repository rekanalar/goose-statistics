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
@statistics.goose|if alt == ALT_TWO_SIDED { x = 0.5 * (1.0 - abs(r)); }|if alt == ALT_TWO_SIDED { x = 0.4 * (1.0 - abs(r)); }|the Pearson p-value argument
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
@numeric.goose|if (fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0) {|if (fa > 0.0 && fb > 0.0) {|the root finder needs a sign change
@resample.goose|INTERVAL_BCA { return rs_boot_fault(ERR_DOMAIN); }|INTERVAL_BCA { return rs_boot_fault(ERR_UNDEFINED); }|an unknown interval gives the domain fault
@resample.goose|if method.kind != RESAMPLE_BOOTSTRAP { return rs_boot_fault(ERR_DOMAIN); }|if method.kind != RESAMPLE_MONTE_CARLO { return rs_boot_fault(ERR_DOMAIN); }|the kind of the resampling method
@resample.goose|if method.interval == INTERVAL_BCA && sizes.len > 1 { return rs_boot_fault(ERR_DOMAIN); }|if false { return rs_boot_fault(ERR_DOMAIN); }|the BCa form refuses several samples
@fitting.goose|if !(m > 0.0) { return nan(), nan(), nan(), ERR_FIT; }|if !(m < 0.0) { return nan(), nan(), nan(), ERR_FIT; }|a fit needs a positive mean
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
distributions_location.goose|if z < -1.0 { return atan(-1.0 / z) / PI; }|if false { return atan(-1.0 / z) / PI; }|the Cauchy lower tail keeps its digits
distributions_location.goose|q = -1.0 / tan(PI * p);|q = -2.0 / tan(PI * p);|the Cauchy quantile of a low probability
distributions_location.goose|if z < 0.0 { return 0.5 * exp(z); }|if z < 0.0 { return 1.0 - 0.5 * exp(-z); }|the Laplace tail below the location
distributions_location.goose|if p < 0.5 { q = log(2.0 * p); }|if p < 0.5 { q = log(p); }|the Laplace quantile below the median
distributions_location.goose|let l, c2 = log1p(-p);\n    return loc_result(loc + scale * sqrt(-2.0 * l));|let l = log(1.0 - p); let c2 = 0;\n    return loc_result(loc + scale * sqrt(-2.0 * l));|the Rayleigh quantile of a tiny probability keeps its digits
distributions_location.goose|return 2.0 * SQRT_PI * (PI - 3.0) / (d * sqrt(d)), OK;|return 2.0 * SQRT_PI * (PI - 2.0) / (d * sqrt(d)), OK;|the Rayleigh skewness
distributions_location.goose|return 12.0 * sqrt(6.0) * ZETA3 / (PI * PI * PI), OK;|return 12.0 * sqrt(5.0) * ZETA3 / (PI * PI * PI), OK;|the Gumbel skewness
distributions_location.goose|return 2.4, OK;|return 2.0, OK;|the Gumbel excess kurtosis
distributions_location.goose|return -v, c;|return v, c;|the left Gumbel skewness changes sign
distributions_location.goose|if scale <= 0.0 { return ERR_DOMAIN; }|if scale < 0.0 { return ERR_DOMAIN; }|a scale of 0 is a domain fault
distributions_location.goose|fn loc_result(v: f64) -> f64, i64 {\n    if is_nan(v) { return nan(), ERR_UNDEFINED; }\n    if is_inf(v) { return nan(), ERR_OVERFLOW; }|fn loc_result(v: f64) -> f64, i64 {\n    if is_nan(v) { return nan(), ERR_UNDEFINED; }|an overflowing result is a fault
distributions_location.goose|if z <= 0.0 { return 0.0, OK; }\n    if is_inf(z) { return 0.0, OK; }|if z <= 0.0 { return 0.0, OK; }|the Rayleigh density at an infinite score
distributions_location.goose|if z < 0.0 { return nan(), ERR_UNDEFINED; }|if z < 0.0 { return 0.0, OK; }|the right Gumbel density has no value at an overflowed score
distributions_location.goose|if z > 0.0 { return nan(), ERR_UNDEFINED; }|if z > 0.0 { return 0.0, OK; }|the left Gumbel density has no value at an overflowed score
contingency.goose|while d < shape.len { s *= shape[d]; d += 1; }|while d < shape.len - 1 { s *= shape[d]; d += 1; }|the stride of an axis
contingency.goose|e *= shares[offsets[d] + idx];|e *= shares[offsets[d]];|the expected count uses the index of each axis
contingency.goose|if total == 0.0 { return ERR_UNDEFINED; }|if total == 0.0 { return OK; }|an empty table has no expected counts
contingency.goose|let yates = correction \&\& (rows - 1) * (cols - 1) == 1;|let yates = correction;|the continuity correction of a 2 by 2 table only
contingency.goose|let mag = min(0.5, abs(diff));|let mag = abs(diff);|the correction moves a count by at most one half
contingency.goose|value = phi2 / (min(rows, cols) - 1);|value = phi2 / (max(rows, cols) - 1);|the Cramer divisor
contingency.goose|value = phi2 / sqrt(((rows - 1) * (cols - 1)) as f64);|value = phi2 / ((rows - 1) * (cols - 1));|the Tschuprow divisor
contingency.goose|else { value = phi2 / (1.0 + phi2); }|else { value = phi2; }|the Pearson coefficient
contingency.goose|if v != floor(v) { return nan(), ERR_DOMAIN; }|if false { return nan(), ERR_DOMAIN; }|association needs whole counts
contingency.goose|if shape.len > MAX_DIMS { return ERR_TOO_LONG; }|if shape.len > MAX_DIMS + 5 { return ERR_TOO_LONG; }|the dimension bound
divergence.goose|if lam == 0.0 { return 2.0 * xlogy_term(o, o / e); }|if lam == 0.0 { return 2.0 * xlogy_term(e, o / e); }|the log-likelihood term
divergence.goose|if lam == -1.0 { return 2.0 * xlogy_term(e, e / o); }|if lam == -1.0 { return 2.0 * xlogy_term(o, e / o); }|the modified log-likelihood term
divergence.goose|/ (lam * (lam + 1.0));|/ (lam * lam + 1.0);|the general term divisor
divergence.goose|if x == 0.0 \&\& !is_nan(y) { return 0.0; }|if false { return 0.0; }|a zero count times a log is 0
divergence.goose|return chi_result(s, (k - 1 - ddof));|return chi_result(s, (k - 1));|the degrees of freedom
divergence.goose|if abs(sum_obs - sum_exp) > 1e-8 * max(sum_obs, sum_exp) { return fault_result(ERR_DOMAIN); }|if false { return fault_result(ERR_DOMAIN); }|the sums of the two lists agree
divergence.goose|if is_nan(lam) { return fault_result(ERR_NAN); }|if false { return fault_result(ERR_NAN); }|a NaN power
divergence.goose|if is_inf(lam) { return fault_result(ERR_INFINITY); }|if false { return fault_result(ERR_INFINITY); }|an infinite power
clip.goose|lo = m - sd * low;|lo = m - sd * high;|the low bound uses the low factor
clip.goose|hi = m + sd * high;|hi = m + sd * low;|the high bound uses the high factor
clip.goose|if x >= lo \&\& x <= hi { cur[kept] = x; kept += 1; }|if x > lo \&\& x < hi { cur[kept] = x; kept += 1; }|the bounds are inclusive
clip.goose|if dropped == 0 { break; }|break;|the passes repeat until nothing drops
clip.goose|let sd, cs = stdev_core(cur[..], 0);|let sd, cs = stdev_core(cur[..], 1);|the deviation divides by n
contingency.goose|if s > MAX_TABLE_SIDE|if s > MAX_TABLE_SIDE + 1000|the table side bound
contingency.goose|if rows < 2 |if rows < 1 |association needs two rows
contingency.goose|cols < 2 { return nan(), ERR_TOO_SHORT; }|cols < 1 { return nan(), ERR_TOO_SHORT; }|association needs two columns
contingency.goose|v > MAX_COUNT { return ERR_DOMAIN; }|v > MAX_COUNT * 1e10 { return ERR_DOMAIN; }|a count above the exact limit
divergence.goose|if is_nan(t)|if false|a NaN term gives the undefined fault
divergence.goose|is_inf(t) { return nan(), ERR_UNDEFINED; }|is_inf(t) { return nan(), ERR_OVERFLOW; }|an infinite term gives the undefined fault
divergence.goose|ddof < 0 |ddof < -9 |the lower bound of ddof
clip.goose|if low < 0.0 |if low < -9.0 |a negative low factor
clip.goose|high < 0.0 { return fault, ERR_DOMAIN; }|high < -9.0 { return fault, ERR_DOMAIN; }|a negative high factor
assoc_rank.goose|return 1.0 / (1.0 + rank);|return 1.0 / (2.0 + rank);|the weigher is the hyperbolic one
assoc_rank.goose|out[perm[k]] = n - 1 - k;|out[perm[k]] = k;|the rank counts down from the top
assoc_rank.goose|if additive { return s * (count - 1.0); }|if additive { return s * count; }|the additive group weight
assoc_rank.goose|out[k] = out[k] - 1.0;|out[k] = out[k] - 0.5;|the dense ranks start at zero
assoc_rank.goose|if w.ys[perm[offset + j]] <= w.ys[perm[middle + k]] {|if w.ys[perm[offset + j]] < w.ys[perm[middle + k]] {|a tie in the merge is not discordant
assoc_rank.goose|w.exchanges += wr * left;|w.exchanges += left;|the product form multiplies the weights
assoc_rank.goose|if w.additive { w.exchanges += wr * (length0 - j) + left; }|if w.additive { w.exchanges += wr * (length0 - j); }|the additive form adds the left weight
assoc_rank.goose|return (a + b) / 2.0, OK;|return a, OK;|the two directions are averaged
assoc_rank.goose|let r, c3 = correlation(rx[..], ry[..]);|let r, c3 = correlation(xs, ys);|the Spearman rho correlates the ranks
assoc_rank.goose|if yb <= ya { up += 1.0; }|if yb < ya { up += 1.0; }|the upper rank counts ties
assoc_rank.goose|if yb >= ya { down += 1.0; }|if yb > ya { down += 1.0; }|the lower rank counts ties
assoc_rank.goose|stat = 1.0 - 3.0 * num / (n * n - 1);|stat = 1.0 - 3.0 * num / (n * n);|the continuous xi normaliser
assoc_rank.goose|stat = 1.0 - n * num / den;|stat = 1.0 - num / den;|the tied xi uses n
assoc_rank.goose|den *= 2.0;|den *= 1.0;|the tied xi denominator doubles
assoc_rank.goose|return sqrt(0.4) / sqrt(nf);|return sqrt(0.5) / sqrt(nf);|the continuous xi deviation
assoc_rank.goose|let c = 2.0 * nf - 2.0 * i + 1.0;|let c = 2.0 * nf - 2.0 * i;|the xi weight of rank i
assoc_rank.goose|an /= nf * nf * nf * nf;|an /= nf * nf * nf;|the an scale of the xi deviation
distance.goose|if p == 1 { sum += abs(d) * gap; }|if p == 1 { sum += d * gap; }|the Wasserstein gap uses the absolute value
distance.goose|else { sum += d * d * gap; }|else { sum += d * gap; }|the energy gap uses the square
distance.goose|let d = cu / tu - cv / tv;|let d = cu / tv - cv / tu;|the gap is the difference of the shares
distance.goose|let gap = x2 - x;|let gap = x2 + x;|the length of the interval
distance.goose|while i < u.len && u[i].x == x { cu += u[i].w; i += 1; }|while i < u.len && u[i].x == x { cu += 1.0; i += 1; }|the sweep adds the weight of the value
distance.goose|if cw != OK { return nan(), cw; }|if false { return nan(), cw; }|the weight check of the first list
distance.goose|let r = ((nf - 1.5) * nf * d * d - 0.5 * s2) / den;|let r = ((nf - 1.5) * nf * d * d + 0.5 * s2) / den;|the transform subtracts the spread
distance.goose|let den = (nf - 1.0) * (nf - 2.0);|let den = (nf - 1.0) * nf;|the denominator of the transform
distance_nd.goose|s += d * d;|s += d;|the cost is the Euclidean distance
distance_nd.goose|let r = sqrt(s);|let r = s;|the root of the squared distance
distance_nd.goose|dist[n + j] < bt { bt|dist[n + j] > bt { bt|the closest sink with demand
distance_nd.goose|if f < amt { amt = f; }|if false { amt = f; }|a reverse arc caps the move
distance_nd.goose|for w in uw { sup.push(w / su); }|for w in uw { sup.push(w); }|the supplies sum to 1
distance_nd.goose|v.len % dims != 0 { return ERR_LENGTH; }|false { return ERR_LENGTH; }|the length check of the second set
distance_nd.goose|if dims > MAX_NDW_DIMS { return ERR_TOO_LONG; }|if dims > MAX_NDW_DIMS + 1 { return ERR_TOO_LONG; }|the bound of the dimensions
distance_nd.goose|if is_inf(r) { return ERR_OVERFLOW; }|if false { return ERR_OVERFLOW; }|a distance that overflows
distance_nd.goose|total += cost[k] * flow[k];|total += cost[k] * flow[k] * 2.0;|the total cost of the flow
assoc_mgc.goose|let nm = (n - 1) as f64;|let nm = n as f64;|the divisor of the column centering
assoc_mgc.goose|col[t].v != col[t - 1].v { rank += 1; }|col[t].v == col[t - 1].v { rank += 1; }|equal distances share a rank
assoc_mgc.goose|if i != j { a[i * n + j] = d[i * n + j] - mean; }|a[i * n + j] = d[i * n + j] - mean;|the diagonal of a centered table is 0
assoc_mgc.goose|let at = ord[j] * n + ord[i];|let at = ord[i] * n + ord[j];|the second table is read transposed
assoc_mgc.goose|ex[k] * ey[l] / nn|ex[k] * ey[l] / n as f64|the mean term divides by n squared
assoc_mgc.goose|if vx[k] <= 0.0 |if vx[k] < 0.0 |a variance of 0 gives a correlation of 0
assoc_mgc.goose|let t = s * (s - 3.0) / 4.0 - 0.5;|let t = s * (s - 3.0) / 4.0;|the shape of the beta threshold
assoc_mgc.goose|1.0 - 0.02 / s|1.0 - 0.05 / s|the level of the beta threshold
assoc_mgc.goose|if size as f64 < big { return global, nx, ny; }|if false { return global, nx, ny; }|a small region keeps the global scale
assoc_mgc.goose|if l < ny - 1 && lab[c + 1] == 0 && map[c + 1] > thr { lab[c + 1] = label; stack.push(c + 1); }|if false { stack.push(c + 1); }|the region joins its right neighbor
assoc_mgc.goose|if ps >= stat { count += 1; }|if ps > stat { count += 1; }|an equal statistic counts as extreme
assoc_mgc.goose|let p = (1 + count) as f64 / (1 + reps) as f64;|let p = count as f64 / (1 + reps) as f64;|the p-value adds one to the count
assoc_mgc.goose|if nu < 5 |if nu < 4 |each sample of the pooled test needs 5 points
assoc_mgc.goose|if ys.len / ydims != n { return mgc_fault(), ERR_LENGTH; }|if false { return mgc_fault(), ERR_LENGTH; }|the samples need the same number of points
assoc_mgc.goose|if reps > MAX_MGC_REPS |if reps > MAX_MGC_REPS + 1 |the bound of the repetitions
assoc_mgc.goose|if n * n * (1 + reps) > MAX_MGC_WORK { return ERR_TOO_LONG; }|if false { return ERR_TOO_LONG; }|the bound of the work
assoc_mgc.goose|if reps < 0 { return ERR_DOMAIN; }|if false { return ERR_DOMAIN; }|a negative count of repetitions
line_siegel.goose|if dx == 0.0 { continue; }|if dx == 0.0 { }|a pair with the same x gives no slope
line_siegel.goose|let s = (ys[i] - ys[j]) / dx;|let s = (ys[j] - ys[i]) / dx;|the slope of a pair
line_siegel.goose|inter.push(ys[i] - s * xs[i]);|inter.push(ys[i] + s * xs[i]);|the intercept of a pair
line_siegel.goose|res.push(ys[i] - slope * xs[i]);|res.push(ys[i] + slope * xs[i]);|the hierarchical residual
line_siegel.goose|if slopes.len == 0 { return nan(), nan(), ERR_UNDEFINED; }|if false { return nan(), nan(), ERR_UNDEFINED; }|no slope at a point is undefined
line_siegel.goose|if xs.len > MAX_SIEGEL_POINTS { return nan(), nan(), ERR_TOO_LONG; }|if false { return nan(), nan(), ERR_TOO_LONG; }|the point bound
laws.goose|if fa > 0.5 { mass = sfa - sfb; } else { mass = fb - fa; }|if fa > 0.5 { mass = sfa - sfb; } else { mass = fb + fa; }|the mass of the range
laws.goose|x >= ub { return 0.0, OK; }|x > ub { return 0.0, OK; }|the density at the high bound is zero
laws.goose|if lb >= ub { return ERR_DOMAIN; }|if lb > ub { return ERR_DOMAIN; }|equal bounds give the domain fault
laws.goose|if !(mass > 0.0) { return ERR_UNDEFINED; }|if !(mass >= 0.0) { return ERR_UNDEFINED; }|a range with no mass is undefined
laws.goose|if x >= ub { return 1.0, OK; }|if x >= ub { return 0.0, OK; }|the cdf at the high bound
laws.goose|v = (sfa - s) / mass;|v = (sfa + s) / mass;|the cdf of an upper range
laws.goose|fam == LAW_CAUCHY && (is_inf(lb)|fam == LAW_CAUCHY && (is_nan(lb)|a Cauchy bound that is infinite
laws.goose|t += (p1 - 1.0) * log(lo);|t += p1 * log(lo);|the lower power of the order density
laws.goose|t += (p2 - p1) * log(hi);|t += (p2 - p1 + 1.0) * log(hi);|the upper power of the order density
laws.goose|return betainc(r, n - r + 1, f);|return betainc(r, n - r, f);|the cdf of an order statistic
laws.goose|let q, cq = betaincinv(r, n - r + 1, p);|let q, cq = betaincinv(r, n - r, p);|the quantile of an order statistic
laws.goose|r > n { return ERR_DOMAIN; }|r >= n { return ERR_DOMAIN; }|the rank may equal n
laws.goose|if n > MAX_ORDER_N { return ERR_TOO_LONG; }|if false { return ERR_TOO_LONG; }|the size bound
laws.goose|if need_counts && x.len != alpha.len { return ERR_LENGTH; }|if false { return ERR_LENGTH; }|the length check
laws.goose|if total != n { return ERR_DOMAIN; }|if false { return ERR_DOMAIN; }|the counts sum to n
laws.goose|v != floor(v) { return ERR_DOMAIN; }|false { return ERR_DOMAIN; }|the counts are whole numbers
laws.goose|for v in alpha { if v <= 0.0 { return ERR_DOMAIN; } }|for v in alpha { if v < 0.0 { return ERR_DOMAIN; } }|alpha is positive
laws.goose|s -= log(x[i]) + t;|s -= t;|the log of a count
laws.goose|var s = log(n as f64) + head;|var s = head;|the log of n
laws.goose|if n == 0 { return 0.0, OK; }|if n == 0 { return 1.0, OK; }|no draws give log 1
laws.goose|else { out.push(-n * alpha[i] * alpha[j] * k); }|else { out.push(n * alpha[i] * alpha[j] * k); }|the sign of the covariance
laws.goose|for v in alpha { out.push(n * v * (a0 - v) * k); }|for v in alpha { out.push(n * v * (a0 + v) * k); }|the variance of a category
laws.goose|        lo = cut;|        lo = lb;|an infinite low bound is cut at the tail
fitting_mle.goose|let sd = log(mean) - mlog;|let sd = log(mean) + mlog;|the gamma equation constant
fitting_mle.goose|return log(a) - d - sd;|return log(a) + d - sd;|the gamma equation sign
fitting_mle.goose|let scale = mean / a;|let scale = mean * a;|the gamma scale
fitting_mle.goose|return d1 - d2 - mlog1;|return d1 - d2 - mlog;|the inner beta equation
fitting_mle.goose|let scale = exp(lmax + log(sp / s.len) / m);|let scale = exp(lmax + log(sp) / m);|the Weibull scale
fitting_shapes.goose|let a = (k + 1) as f64;|let a = (k + 2) as f64;|the Filliben median shape
fitting_shapes.goose|sum = sum + head - a * log1p(z * z / df);|sum = sum + head - a * log1p(z * z) / df;|the shape of the student law
fitting_mle.goose|if family == FAMILY_BETA { return betadist_cdf(x, a, b); }|if family == FAMILY_BETA { return betadist_sf(x, a, b); }|the beta cdf
fitting_mle.goose|if family == FAMILY_GAMMA { return gammadist_sf(x, a, b); }|if family == FAMILY_GAMMA { return gammadist_cdf(x, a, b); }|the gamma survival
fitting_mle.goose|else { v, c = weibull_ppf(u, a, b); }|else { v, c = weibull_ppf(u, b, a); }|the Weibull draw
fitting.goose|if family >= FAMILY_WEIBULL { return mle_cdf(family, x, a, b); }|if family > FAMILY_WEIBULL { return mle_cdf(family, x, a, b); }|the Weibull cdf is routed
fitting.goose|else if family == FAMILY_BETA { out.push(a); out.push(b); out.push(0.0); out.push(1.0); }|else if family == FAMILY_BETA { out.push(a); out.push(b); out.push(0.0); out.push(0.0); }|the beta parameters list
fitting.goose|for x in xs { if !(x < 1.0) { return fault_result(ERR_DOMAIN); } }|for x in xs { if !(x < 2.0) { return fault_result(ERR_DOMAIN); } }|the beta values are below one
fitting.goose|sum += (2.0 * (i + 1) - 1.0) / nf * (lo[i] + hi[n - 1 - i]);|sum += (2.0 * (i + 1) - 1.0) / nf * (lo[i] + hi[i]);|the statistic pairs the lower and the upper tail
anderson_families.goose|return -nf - ad_sum(lo[..], hi[..]);|return -nf + ad_sum(lo[..], hi[..]);|the sign of the statistic
anderson_families.goose|hi = log_one_minus_exp_neg(e);|hi = -e;|the Gumbel right survival
anderson_families.goose|lo = log_one_minus_exp_neg(e);|lo = -e;|the Gumbel left cdf
anderson_families.goose|hi = logistic_log_cdf(-z);|hi = logistic_log_cdf(z);|the logistic survival
anderson_families.goose|if z >= 0.0 { return -log(1.0 + exp(-z)); }|if z >= 0.0 { return -log(1.0 + exp(z)); }|the logistic cdf
anderson_families.goose|if a2 <= crit[0] { return level[0]; }|if a2 <= crit[0] { return level[1]; }|the p-value below the table
anderson_families.goose|let t = (a2 - crit[i - 1]) / (crit[i] - crit[i - 1]);|let t = (a2 - crit[i]) / (crit[i] - crit[i - 1]);|the interpolation weight
anderson_families.goose|return level[crit.len - 1];|return level[0];|the p-value above the table
anderson_families.goose|loc = mn - s * log(w / n);|loc = mn + s * log(w / n);|the Gumbel location
anderson_families.goose|return s - mean + num / den;|return s + mean + num / den;|the Gumbel scale equation
anderson_families.goose|if f == ANDERSON_GUMBEL_L { p1 = -p1; }|if f == ANDERSON_GUMBEL_R { p1 = -p1; }|the Gumbel left location
anderson_families.goose|round3(v / (1.0 + 0.6 / n))|round3(v / (1.0 + 0.5 / n))|the exponential critical scale
anderson_families.goose|round3(v / (1.0 + 0.25 / n))|round3(v / (1.0 + 0.2 / n))|the logistic critical scale
anderson_families.goose|round3(v / (1.0 + 0.2 / sqrt(n)))|round3(v / (1.0 + 0.2 / n))|the Gumbel critical scale
anderson_families.goose|round3(v + 0.0005)|round3(v)|the Weibull critical value rounds up
anderson_families.goose|else if d == 0.5 && (f - 2.0 * floor(f / 2.0)) != 0.0 { r = f + 1.0; }|else if d == 0.5 { r = f + 1.0; }|a tie rounds to even
anderson_families.goose|if pos > 10.0 { pos = 10.0; }|if pos > 9.0 { pos = 9.0; }|the last Weibull row
anderson_families.goose|var d = 1e-7 * range;|var d = 1e-3 * range;|the scan starts near the minimum
anderson_families.goose|if c != OK { return fault_result(c); }|if false { return fault_result(c); }|the sample check
anderson_families.goose|if f != ANDERSON_EXPON &&|if f == 99 &&|the family check
probplot.goose|out[0] = 1.0 - last;|out[0] = last;|the first median
probplot.goose|let last = pow(0.5, 1.0 / nf);|let last = pow(0.5, 1.0 / (nf + 1.0));|the last median
probplot.goose|out[k] = ((k + 1) - 0.3175) / (nf + 0.365);|out[k] = ((k + 1) - 0.3) / (nf + 0.4);|the inner medians
probplot.goose|if lam == 0.0 { return log(p) - log(1.0 - p); }|if lam == 0.0 { return 0.0; }|the Tukey lambda of zero
probplot.goose|return (pow(p, lam) - pow(1.0 - p, lam)) / lam;|return (pow(p, lam) + pow(1.0 - p, lam)) / lam;|the Tukey lambda quantile
probplot.goose|if family == PLOT_RAYLEIGH { return rayleigh_ppf(p); }|if family == PLOT_RAYLEIGH { return expon_ppf(p); }|the Rayleigh quantile
probplot.goose|if family == PLOT_GUMBEL_L { return gumbel_l_ppf(p); }|if family == PLOT_GUMBEL_L { return gumbel_r_ppf(p); }|the Gumbel left quantile
probplot.goose|if family == PLOT_WEIBULL { return weibull_ppf(p, shape); }|if family == PLOT_WEIBULL { return weibull_ppf(p, 1.0); }|the Weibull shape
probplot.goose|if n < 1 || !(b > a) { return ERR_DOMAIN; }|if n < 1 || !(b >= a) { return ERR_DOMAIN; }|equal ends of the grid
probplot.goose|if n > 1 { s = a + (b - a) * k / (n - 1); }|if n > 1 { s = a + (b - a) * k / n; }|the grid step
probplot.goose|return 1.0 - r;|return 1.0 + r;|the loss is one minus r
probplot.goose|if xs.len > MAX_PPCC_POINTS { return nan(), ERR_TOO_LONG; }|if false { return nan(), ERR_TOO_LONG; }|the sample bound of the search
probplot.goose|if cl != OK { osm.resize(0); osr.resize(0); return none, cl; }|if cl != OK { return none, cl; }|a fault leaves the plot empty
probplot.goose|n > MAX_PPCC_GRID { return ERR_TOO_LONG; }|false { return ERR_TOO_LONG; }|the grid bound
laws.goose|if alpha.len > MAX_DIRMULT_K { return ERR_TOO_LONG; }|if false { return ERR_TOO_LONG; }|the category bound
laws.goose|if fam == LAW_LOGISTIC { return logi_pdf(x, a, b); }|if fam == LAW_LOGISTIC { return logi_cdf(x, a, b); }|the logistic density
laws.goose|let e = exp(-abs((x - a) / b));|let e = exp(-((x - a) / b));|the logistic density is even
laws.goose|if a > lo && a < hi { cuts.push(a); }|if false { cuts.push(a); }|the location splits the range

assoc_somers.goose|let d = (p - q) / den;|let d = (p + q) / den;|the sign of the statistic
assoc_somers.goose|let den = n * n - sri2;|let den = n * n;|the denominator of the statistic
assoc_somers.goose|p += a * conc;|p += a * disc;|the concordant count
assoc_somers.goose|q += a * disc;|q += a * conc;|the discordant count
assoc_somers.goose|s += a * (conc - disc) * (conc - disc);|s += a * (conc - disc);|the sum of the squares
assoc_somers.goose|let var_s = s - (p - q) * (p - q) / n;|let var_s = s + (p - q) * (p - q) / n;|the null variance
assoc_somers.goose|let z = (p - q) / sqrt(4.0 * var_s);|let z = (p - q) / sqrt(2.0 * var_s);|the test statistic
assoc_somers.goose|if alt == ALT_LESS { pv = lower; }|if alt == ALT_LESS { pv = upper; }|the lower alternative
assoc_somers.goose|else if alt == ALT_GREATER { pv = upper; }|else if alt == ALT_GREATER { pv = lower; }|the upper alternative
assoc_somers.goose|else { pv = 2.0 * min(lower, upper); }|else { pv = min(lower, upper); }|the two-sided p-value
assoc_somers.goose|if nonzero < 2 { return corr_fault(ERR_DOMAIN), ERR_DOMAIN; }|if nonzero < 1 { return corr_fault(ERR_DOMAIN), ERR_DOMAIN; }|one nonzero cell
assoc_somers.goose|if rows < 2 || cols < 2 { return corr_fault(ERR_UNDEFINED), ERR_UNDEFINED; }|if rows < 2 || cols < 2 { return corr_fault(ERR_DOMAIN), ERR_DOMAIN; }|a single row has no order
assoc_somers.goose|if p == q { return corr_fault(ERR_UNDEFINED), ERR_UNDEFINED; }|if false { return corr_fault(ERR_UNDEFINED), ERR_UNDEFINED; }|a zero variance
assoc_somers.goose|if lx.len * ly.len > MAX_TABLE_CELLS { return corr_fault(ERR_TOO_LONG), ERR_TOO_LONG; }|if lx.len * ly.len > MAX_TABLE_CELLS * 2 { return corr_fault(ERR_TOO_LONG), ERR_TOO_LONG; }|the table bound
assoc_somers.goose|if out.len == 0 || out[out.len - 1] != v { out.push(v); }|if out.len == 0 { out.push(v); }|the distinct levels
assoc_somers.goose|for v in table { if v != floor(v) { return corr_fault(ERR_DOMAIN), ERR_DOMAIN; } }|for v in table { if false { return corr_fault(ERR_DOMAIN), ERR_DOMAIN; } }|the counts are whole numbers
contingency.goose|if cells != table.len { return ERR_LENGTH; }|if false { return ERR_LENGTH; }|the length of the table
contingency_risk.goose|let or = a * d / (b * cc);|let or = a * cc / (b * d);|the ratio of the diagonals
contingency_risk.goose|let se = sqrt(1.0 / a + 1.0 / b + 1.0 / cc + 1.0 / d);|let se = sqrt(1.0 / a + 1.0 / b + 1.0 / cc);|the standard error of the log ratio
contingency_risk.goose|let z, c = norm_ppf(0.5 * level + 0.5);|let z, c = norm_ppf(level);|the z value of the two-sided interval
contingency_risk.goose|let l = log(or) - z * se;|let l = log(or) + z * se;|the lower end of the interval
contingency_risk.goose|let h = log(or) + z * se;|let h = log(or) - z * se;|the upper end of the interval
contingency_risk.goose|if alt == ALT_LESS { return 0.0; }|if alt == ALT_GREATER { return 0.0; }|the lower end of the less interval
contingency_risk.goose|if alt == ALT_GREATER { return 1.0 / 0.0; }|if alt == ALT_LESS { return 1.0 / 0.0; }|the upper end of the greater interval
contingency_risk.goose|if x == lo { return odds_ratio_conditional_at(x, 0.0, draws, type1, pop, alt, level, logs[..]); }|if x == lo { return odds_ratio_conditional_at(x, 1.0, draws, type1, pop, alt, level, logs[..]); }|the ratio at the low end of the support
contingency_risk.goose|if x == hi { return interval_fault(ERR_OVERFLOW), ERR_OVERFLOW; }|if false { return interval_fault(ERR_OVERFLOW), ERR_OVERFLOW; }|the ratio at the high end of the support
contingency_risk.goose|let t, c = nchyp_tilt_mean(x, logs[..], lo);|let t, c = nchyp_tilt_mean(x + 1.0, logs[..], lo);|the mean equation of the tilt
contingency_risk.goose|num += (lo + j) * e;|num += e;|the mean counts the support
contingency_risk.goose|if lo + j <= k { num += e; }|if lo + j < k { num += e; }|the share of the count
contingency_risk.goose|let se = sqrt(1.0 / exposed_cases - 1.0 / exposed_total + 1.0 / control_cases - 1.0 / control_total);|let se = sqrt(1.0 / exposed_cases + 1.0 / exposed_total + 1.0 / control_cases + 1.0 / control_total);|the standard error of the log risk
contingency_risk.goose|if exposed_cases == 0.0 && control_cases == 0.0 {|if false {|the risk of two empty groups
contingency_risk.goose|if table.len != 4 { return ERR_LENGTH; }|if false { return ERR_LENGTH; }|the length of the table
kde.goose|if bw == KDE_SCOTT { f = pow(neff, -1.0 / (dd + 4.0)); }|if bw == KDE_SCOTT { f = pow(neff, -1.0 / (dd + 3.0)); }|the factor of Scott
kde.goose|else if bw == KDE_SILVERMAN { f = pow(neff * (dd + 2.0) / 4.0, -1.0 / (dd + 4.0)); }|else if bw == KDE_SILVERMAN { f = pow(neff * (dd + 3.0) / 4.0, -1.0 / (dd + 4.0)); }|the factor of Silverman
kde.goose|cov[a * dims + b] = cov[a * dims + b] + w[i] * da * db;|cov[a * dims + b] = cov[a * dims + b] + w[i] * da;|the covariance of the points
kde.goose|let denom = 1.0 - ss;|let denom = 1.0;|the correction of the weights
kde.goose|for s in sums { out.push(s / norm); }|for s in sums { out.push(s * norm); }|the normalization of the density
kde.goose|else { out.push(log(s) - norm); }|else { out.push(log(s) + norm); }|the logarithm of the density
kde.goose|value = value + kde_weight(k, i) * (hi - lo);|value = value + (hi - lo);|the weights of the box
kde.goose|let stdev = kde_cho(k, 0, 0);|let stdev = kde_cho(k, 0, 0) * 2.0;|the band width of the kernel
kde.goose|for i in 0..j + 1 { v = v + kde_cho(k, j, i) * z[i]; }|for i in 0..j { v = v + kde_cho(k, j, i) * z[i]; }|the diagonal of the noise
kde.goose|acc = acc + kde_weight(k, i);|acc = acc + 1.0 / k.n as f64;|the draw of the index
kde.goose|sum = sum + exp(-0.5 * q) * kde_weight(k, i);|sum = sum + exp(-0.5 * q);|the weights of the density
kde.goose|if dims > n { return kde_fault(ERR_DOMAIN), ERR_DOMAIN; }|if false { return kde_fault(ERR_DOMAIN), ERR_DOMAIN; }|more dimensions than points
distributions_elementary.goose|return cos(2.0 * x), OK;|return cos(x), OK;|the density of the anglit law
distributions_elementary.goose|let s = sin(x + PI / 4.0);\n    return s * s, OK;|let s = sin(x + PI / 4.0);\n    return s, OK;|the cdf of the anglit law
distributions_elementary.goose|return asin(sqrt(p)) - PI / 4.0, OK;|return sqrt(p) - PI / 4.0, OK;|the quantile of the anglit law
distributions_elementary.goose|return 1.0 / (PI * sqrt(x * (1.0 - x))), OK;|return 1.0 / (PI * sqrt(x * (1.0 + x))), OK;|the density of the arcsine law
distributions_elementary.goose|return 2.0 / PI * asin(sqrt(x)), OK;|return 2.0 / PI * asin(x), OK;|the cdf of the arcsine law
distributions_elementary.goose|let s = sin(PI / 2.0 * p);|let s = sin(PI * p);|the quantile of the arcsine law
distributions_elementary.goose|return (1.0 + cos(x)) / (2.0 * PI), OK;|return (1.0 - cos(x)) / (2.0 * PI), OK;|the density of the cosine law
distributions_elementary.goose|return (PI + x + sin(x)) / (2.0 * PI), OK;|return (PI + x - sin(x)) / (2.0 * PI), OK;|the cdf of the cosine law
distributions_elementary.goose|return (PI - x - sin(x)) / (2.0 * PI), OK;|return (PI - x + sin(x)) / (2.0 * PI), OK;|the survival function of the cosine law
distributions_elementary.goose|if p == 0.5 { return 0.0, OK; }|if p == 0.25 { return 0.0, OK; }|the median of the cosine law
distributions_elementary.goose|return exp(-0.5 * l * l) / (x * SQRT_2PI), OK;|return exp(-0.5 * l) / (x * SQRT_2PI), OK;|the density of the gibrat law
distributions_elementary.goose|return norm_cdf(log(x));|return norm_cdf(x);|the cdf of the gibrat law
distributions_elementary.goose|let r = exp(norm_std_ppf(p));|let r = norm_std_ppf(p);|the quantile of the gibrat law
distributions_elementary.goose|return 2.0 / (PI * (1.0 + x * x)), OK;|return 2.0 / (PI * (1.0 + x)), OK;|the density of the half Cauchy law
distributions_elementary.goose|return 2.0 / PI * atan(x), OK;|return 1.0 / PI * atan(x), OK;|the cdf of the half Cauchy law
distributions_elementary.goose|return tan(PI / 2.0 * p), OK;|return tan(PI / 4.0 * p), OK;|the quantile of the half Cauchy law
distributions_elementary.goose|return 2.0 * e / ((1.0 + e) * (1.0 + e)), OK;|return e / ((1.0 + e) * (1.0 + e)), OK;|the density of the half logistic law
distributions_elementary.goose|return -m / (1.0 + exp(-x)), cm;|return m / (1.0 + exp(-x)), cm;|the cdf of the half logistic law
distributions_elementary.goose|return 2.0 * exp(-0.5 * x * x) / SQRT_2PI, OK;|return exp(-0.5 * x * x) / SQRT_2PI, OK;|the density of the half normal law
distributions_elementary.goose|return erf(x * INV_SQRT2);|return erf(x);|the cdf of the half normal law
distributions_elementary.goose|return norm_std_ppf(0.5 + 0.5 * p), OK;|return norm_std_ppf(0.5 + 0.25 * p), OK;|the quantile of the half normal law
distributions_elementary.goose|return 2.0 / (PI * (exp(x) + exp(-x))), OK;|return 1.0 / (PI * (exp(x) + exp(-x))), OK;|the density of the hypsecant law
distributions_elementary.goose|return 2.0 / PI * atan(exp(x)), OK;|return 1.0 / PI * atan(exp(x)), OK;|the cdf of the hypsecant law
distributions_elementary.goose|return log(tan(PI * p / 2.0)), OK;|return log(tan(PI * p / 4.0)), OK;|the quantile of the hypsecant law
distributions_elementary.goose|fn halfcauchy_mean() -> f64, i64 { return nan(), ERR_UNDEFINED; }|fn halfcauchy_mean() -> f64, i64 { return 0.0, OK; }|the mean of the half Cauchy law
distributions_elementary.goose|return PI * PI / 16.0 - 0.5, OK;|return PI * PI / 8.0 - 0.5, OK;|the variance of the anglit law
distributions_bounded.goose|return 2.0 / PI * sqrt(1.0 - x * x), OK;|return 2.0 / PI * sqrt(1.0 + x * x), OK;|the density of the semicircular law
distributions_bounded.goose|return 0.5 + (x * sqrt(1.0 - x * x) + asin(x)) / PI, OK;|return 0.5 + (x * sqrt(1.0 - x * x) - asin(x)) / PI, OK;|the cdf of the semicircular law
distributions_bounded.goose|return 2.0 * b - 1.0, OK;|return 2.0 * b + 1.0, OK;|the quantile of the semicircular law
distributions_bounded.goose|if x < c { return 2.0 * x / c, OK; }|if x < c { return 2.0 * x / (1.0 - c), OK; }|the density of the triangular law
distributions_bounded.goose|if x < c { return x * x / c, OK; }|if x < c { return x * x / (1.0 - c), OK; }|the cdf of the triangular law
distributions_bounded.goose|if p < c { return sqrt(c * p), OK; }|if p < c { return sqrt((1.0 - c) * p), OK; }|the quantile of the triangular law
distributions_bounded.goose|return (c + 1.0) / 3.0, OK;|return (c + 2.0) / 3.0, OK;|the mean of the triangular law
distributions_bounded.goose|return (1.0 - c + c * c) / 18.0, OK;|return (1.0 - c + c * c) / 9.0, OK;|the variance of the triangular law
distributions_bounded.goose|if x < c { return u * x / c, OK; }|if x < c { return u * x / d, OK; }|the density of the trapezoid law
distributions_bounded.goose|if x < c { return x * x / c / (d - c + 1.0), OK; }|if x < c { return x * x / d / (d - c + 1.0), OK; }|the cdf of the trapezoid law
distributions_bounded.goose|if p < qc { return sqrt(p * c * (1.0 + d - c)), OK; }|if p < qc { return sqrt(p * d * (1.0 + d - c)), OK; }|the quantile of the trapezoid law
distributions_bounded.goose|if c > 0.0 { t = t + u * pow(c, n + 1) / (n + 2.0); }|if c > 0.0 { t = t + u * pow(c, n + 2) / (n + 2.0); }|the first moment of the trapezoid law
distributions_bounded.goose|t = t + u * (pow(d, n + 1) - pow(c, n + 1)) / (n + 1.0);|t = t + u * (pow(d, n + 1) - pow(c, n + 1)) / (n + 2.0);|the moment of the plateau of the trapezoid law
distributions_bounded.goose|return c / (c * x + 1.0) / l, OK;|return c / (c * x + 1.0) * l, OK;|the density of the Bradford law
distributions_bounded.goose|return m / l, OK;|return m * l, OK;|the cdf of the Bradford law
distributions_bounded.goose|return m / c, OK;|return m * c, OK;|the quantile of the Bradford law
distributions_bounded.goose|return (c - l) / (c * l), OK;|return (c + l) / (c * l), OK;|the mean of the Bradford law
distributions_bounded.goose|return exp(-0.5 * d * d / x) / (SQRT_2PI * x * sqrt(x)), OK;|return exp(-0.5 * d * d / x) / (SQRT_2PI * sqrt(x)), OK;|the density of the wald law
distributions_bounded.goose|return a + e * b, OK;|return a + b, OK;|the cdf of the wald law
distributions_bounded.goose|return sqrt(2.0 / PI) * x * x * exp(-0.5 * x * x), OK;|return sqrt(2.0 / PI) * x * exp(-0.5 * x * x), OK;|the density of the maxwell law
distributions_bounded.goose|return sqrt(2.0 * v), OK;|return sqrt(v), OK;|the quantile of the maxwell law
distributions_bounded.goose|return exp(-0.5 * (x + e)) / SQRT_2PI, OK;|return exp(-0.5 * (x + e)) * SQRT_2PI, OK;|the density of the moyal law
distributions_bounded.goose|return -log(2.0 * v * v), OK;|return -log(2.0 * v), OK;|the quantile of the moyal law
distributions_weibull.goose|return c / 2.0 * pow(a, c - 1.0) * exp(-pow(a, c)), OK;|return c * pow(a, c - 1.0) * exp(-pow(a, c)), OK;|the density of the double Weibull law
distributions_weibull.goose|return exp(-pow(x, -c)), OK;|return exp(-pow(x, c)), OK;|the cdf of the inverse Weibull law
distributions_weibull.goose|return -pow(-log(p), 1.0 / c), OK;|return pow(-log(p), 1.0 / c), OK;|the quantile of the Weibull law of the negative line
distributions_weibull.goose|let s = n as f64 / c + 1.0;|let s = n as f64 / c + 2.0;|the moment of the truncated Weibull law
distributions_weibull.goose|return (ea - e) / d, OK;|return (e - ea) / d, OK;|the cdf of the truncated Weibull law
distributions_exponent.goose|return c * exp(-0.5 * x * x) / SQRT_2PI * pow(p, c - 1.0);|return c * exp(-0.5 * x * x) / SQRT_2PI * pow(p, c);|the density of the power of the normal law
distributions_exponent.goose|let r = -norm_std_ppf(pow(1.0 - p, 1.0 / c));|let r = norm_std_ppf(pow(1.0 - p, 1.0 / c));|the quantile of the power of the normal law
distributions_exponent.goose|return exp(-norm_std_ppf(pow(1.0 - p, 1.0 / c)) * s), OK;|return exp(norm_std_ppf(pow(1.0 - p, 1.0 / c)) * s), OK;|the quantile of the power of the log-normal law
distributions_exponent.goose|return b * pow(x, b - 1.0) * exp(1.0 + xb - exp(xb));|return b * pow(x, b - 1.0) * exp(xb - exp(xb));|the density of the generalized exponential law
distributions_exponent.goose|return a * c * pow(m, a - 1.0) * exp(-xc) * pow(x, c - 1.0);|return a * c * pow(m, a - 1.0) * exp(-xc) * pow(x, c);|the density of the Weibull law of the exponential kind
distributions_numeric.goose|gauss_legendre(gx, gw, FAM_NODES);|gauss_legendre(gx, gw, FAM_NODES - 1);|the nodes of the rule of the moments
distributions_numeric.goose|total += exp(l - lstar) * gw[j] * width / 2.0;|total += exp(l - lstar) * gw[j] * width;|the scale of the rule of the moments
distributions_numeric.goose|total += d * pow(x, n as f64) * gw[j] * width / 2.0;|total += d * pow(x, n as f64) * gw[j] * width;|the scale of the rule of the moments of a law on the whole line
distributions_numeric.goose|uhi = uhi * 2.0;|uhi = uhi * 1.01;|the widening of the scan for the peak
distributions_exponent.goose|fam_raws(1, s) { x => powerlognorm_dens(x, c, s) }|fam_raws(1, 1.0) { x => powerlognorm_dens(x, c, s) }|the scale of the rule of the moments of the power of the log-normal law
distributions_discrete_a.goose|if k == 0.0 { return 1.0 - p, OK; }|if k == 0.0 { return p, OK; }|the mass of the count 0 of the Bernoulli law
distributions_discrete_a.goose|if k < 1.0 { return 1.0 - p, OK; }|if k < 1.0 { return p, OK; }|the cumulative function of the Bernoulli law
distributions_discrete_a.goose|if q <= 1.0 - p { return 0.0, OK; }|if q <= p { return 0.0, OK; }|the quantile of the Bernoulli law
distributions_discrete_a.goose|let lb1, lcb1 = lbeta(k + a, n - k + b);|let lb1, lcb1 = lbeta(k + a, n - k + a);|the density of the beta binomial law
distributions_discrete_a.goose|if k < 0.0 || k > n { return 0.0, OK; }|if k < 0.0 || k > n + 1.0 { return 0.0, OK; }|the support of the beta binomial law
distributions_discrete_a.goose|top = top * 2.0;|top = top * 1.5;|the growth of the bound of the search of the quantile
distributions_discrete_a.goose|if F(mid) >= q { r = mid; } else { l = mid + 1.0; }|if F(mid) > q { r = mid; } else { l = mid + 1.0; }|the rule of the search of the quantile
distributions_discrete_a.goose|return (1.0 - exp(-a)) / (1.0 + exp(-a)) * exp(-a * abs(k)), OK;|return (1.0 - exp(-a)) / (1.0 + exp(-a)) * exp(-a * k), OK;|the mass of the discrete Laplace law
distributions_discrete_a.goose|if x != floor(x) { return ERR_DOMAIN; }|if x != floor(x) { return OK; }|the check of a count that must be a whole number
distributions_discrete_b.goose|return -m * exp(-lambda * k), OK;|return -m * exp(-lambda * (k + 1.0)), OK;|the mass of the Planck law
distributions_discrete_b.goose|let m, cm = expm1(-lambda * (floor(k) + 1.0));|let m, cm = expm1(-lambda * floor(k));|the cumulative function of the Planck law
distributions_discrete_b.goose|next.push(old * (1.0 - p));|next.push(old * p);|the masses of the Poisson binomial law
distributions_discrete_b.goose|next[k + 1] = next[k + 1] + old * p;|next[k + 1] = next[k + 1] + old * (1.0 - p);|the masses of the Poisson binomial law
distributions_discrete_b.goose|return alpha * b, OK;|return b, OK;|the mass of the Yule-Simon law
distributions_discrete_b.goose|let v = 1.0 - k * b;|let v = 1.0 - (k + 1.0) * b;|the cumulative function of the Yule-Simon law
distributions_discrete_b.goose|total += d[0] * weight;|total += d[0] * weight * 2.0;|the Euler transform of the zeta function
distributions_discrete_b.goose|let den = 1.0 - pow(2.0, 1.0 - a);|let den = 1.0 - pow(2.0, a - 1.0);|the factor of the series of the zeta function
distributions_discrete_b.goose|total += pow(k, -a);|total += pow(k, a);|the harmonic number of the Zipf law of the finite kind
distributions_discrete_b.goose|return 1.0 / (h * pow(k, a)), OK;|return 1.0 / (h * pow(k, -a)), OK;|the mass of the Zipf law of the finite kind
distributions_discrete_b.goose|return zipfian_harmonic(floor(k), a) / h, OK;|return zipfian_harmonic(floor(k) + 1.0, a) / h, OK;|the cumulative function of the Zipf law of the finite kind
distributions_discrete_b.goose|while skellam_cdf(r, mu1, mu2) < q { r = r * 2.0; }|while skellam_cdf(r, mu1, mu2) < q { r = r * 3.0; }|the bound of the search of the quantile of the Skellam law
distributions_rice.goose|return x * exp(-0.5 * d * d) * rice_i0e(x * b);|return x * exp(-0.5 * d * d) * rice_i0e(x);|the density of the Rice law
distributions_rice.goose|return total / sqrt(2.0 * PI * z);|return total / sqrt(PI * z);|the asymptotic series of the Bessel function
distributions_rice.goose|let q = 0.25 * z * z;|let q = 0.5 * z * z;|the power series of the Bessel function
distributions_rice.goose|let lw = -lam + j * log(lam) - lg;|let lw = -lam + j * log(lam);|the Poisson weight of the Marcum sum
distributions_rice.goose|let half = 0.5 * x * x;|let half = x * x;|the argument of the incomplete gamma of the Marcum sum
distributions_rice.goose|g, cc = gammaincc(j + 1.0, half);|g, cc = gammainc(j + 1.0, half);|the upper tail of the Marcum sum
distributions_rice.goose|return -m, OK;|return m, OK;|the cumulative function of the Rice law with the shape 0
distributions_rice.goose|const RICE_MAX_B: f64 = 3000.0;|const RICE_MAX_B: f64 = 30000.0;|the bound of the shape of the Rice law
distributions_rice.goose|let q = 1.0 - p;|let q = p;|the upper half of the quantile of the Rice law
distributions_search_a.goose|if x <= 0.0 || x >= 1.0 { return 0.0, OK; }|if x <= 0.0 || x > 1.0 { return 0.0, OK; }|the support of the Argus law
distributions_search_a.goose|let num = argus_phi(chi * sqrt(1.0 - x * x));|let num = argus_phi(chi * sqrt(1.0 + x * x));|the survival function of the Argus law
distributions_search_a.goose|return 0.5 * g;|return g;|the helper of the Argus law
distributions_search_a.goose|l = dpl_log_sub(t1, t2);|l = dpl_log_add(t1, t2);|the cumulative function of the double Pareto log-normal law
distributions_search_a.goose|let sum = dpl_log_add(r1, r2);|let sum = dpl_log_add(r1, r2) + 1.0;|the sum of the terms of the double Pareto log-normal law
distributions_search_a.goose|let l = (a - 1.0) * log(x) - x - fam_lngamma(a);|let l = (a - 1.0) * log(x) - x + fam_lngamma(a);|the density of the Erlang law
distributions_search_a.goose|let l = -log(2.0) - 0.5 * x * x - log(k) + log(exponnorm_erfcx(u));|let l = -log(2.0) - 0.5 * x * x - log(k) - log(exponnorm_erfcx(u));|the density of the sum of a normal and an exponential law
distributions_search_a.goose|return (1.0 / (1.0 + d1 * d1) + 1.0 / (1.0 + d2 * d2)) / PI, OK;|return (1.0 / (1.0 + d1 * d1) - 1.0 / (1.0 + d2 * d2)) / PI, OK;|the density of the folded Cauchy law
distributions_search_a.goose|return (atan(d1) + atan(d2)) / PI, OK;|return (atan(d1) - atan(d2)) / PI, OK;|the cumulative function of the folded Cauchy law
distributions_search_a.goose|return 0.5 * (e1 + e2), OK;|return 0.5 * (e1 - e2), OK;|the cumulative function of the folded normal law
distributions_search_b.goose|return 1.0 / (x * log(b / a)), OK;|return 1.0 / (x * log(a / b)), OK;|the density of the log-uniform law
distributions_search_b.goose|if b <= a { return ERR_DOMAIN; }|if b < a { return ERR_DOMAIN; }|the bounds of the log-uniform law
distributions_search_b.goose|next.push((a + b) / (m - 1) as f64);|next.push((a + b) / m as f64);|the recursion of the cardinal B-spline
distributions_search_b.goose|for j in 0..min(fl + 1, v.len) { cdf += v[j]; }|for j in 0..min(fl, v.len) { cdf += v[j]; }|the sum of the steps of the Irwin-Hall law
distributions_search_b.goose|return exp((a - 1.0) * log(x) + (b - 1.0) * log1p_core(-x) - c * log1p_core(z * x));|return exp((a - 1.0) * log(x) + (b - 1.0) * log1p_core(-x) + c * log1p_core(z * x));|the raw density of the Gauss hypergeometric law
distributions_search_b.goose|return (p - 1.0) * log(x) - 0.5 * b * (x + 1.0 / x) + b - log(2.0) - lkve;|return (p - 1.0) * log(x) - 0.5 * b * (x + 1.0 / x) + b + log(2.0) - lkve;|the log density of the generalized inverse Gaussian law
distributions_search_b.goose|let d2 = (a + b) * (a - b);|let d2 = (a + b) * (a + b);|the width of the Bessel function of the generalized hyperbolic law
distributions_search_c.goose|return log(a / PI) + lk - a * s + b * x + g - log(s);|return log(a / PI) + lk - a * s + b * x - g - log(s);|the log density of the normal inverse Gaussian law
distributions_search_c.goose|let r2 = 1.0 + 1.0 / (rho * rho);|let r2 = 1.0 + rho * rho;|the constant of the relative Breit-Wigner law
distributions_search_c.goose|return sqrt(2.0 * r2 / (1.0 + sqrt(r2))) * 2.0 / PI;|return sqrt(2.0 * r2 / (1.0 + sqrt(r2))) * PI;|the constant of the relative Breit-Wigner density
distributions_search_c.goose|inner += zw[j] * ph * pow(d, k - 1.0);|inner += zw[j] * ph * pow(d, k);|the power of the inner integral of the studentized range law
distributions_search_c.goose|total += sw[i] * exp(lc + (df - 1.0) * log(s) - 0.5 * df * s * s) * inner;|total += sw[i] * exp(lc + (df - 1.0) * log(s) - 0.5 * df * s * s) * 2.0 * inner;|the weight of the outer integral of the studentized range law
distributions_search_c.goose|let v = fam_quad(0.0, y) { t => exp(kappa * (cos(t) - 1.0)) / (2.0 * PI * vonmises_i0e(kappa)) };|let v = fam_quad(0.0, y) { t => exp(kappa * (cos(t) - 1.0)) / (PI * vonmises_i0e(kappa)) };|the integral of the von Mises law
distributions_closed.goose|return tl_out(tl_ncdf(a - 1.0 / x) / tl_ncdf(a));|return tl_out(tl_ncdf(a - 1.0 / x));|the normalising divisor of the alpha law
distributions_closed.goose|return tl_out(1.0 / (a - tl_nppf(q * tl_ncdf(a))));|return tl_out(1.0 / (a - tl_nppf(q)));|the quantile of the alpha law
distributions_closed.goose|return tl_out(tl_ncdf((s - 1.0 / s) / c));|return tl_out(tl_ncdf((s - 1.0 / s)));|the scale of the fatigue life law
distributions_closed.goose|return tl_out(2.0 * tl_ncdf(-1.0 / sqrt(x)));|return tl_out(tl_ncdf(-1.0 / sqrt(x)));|the factor of the Levy cumulative function
distributions_closed.goose|return tl_out(1.0 - 2.0 * tl_ncdf(-1.0 / sqrt(-x)));|return tl_out(2.0 * tl_ncdf(-1.0 / sqrt(-x)));|the mirror of the Levy cumulative function
distributions_closed.goose|if x < 1.0 { return tl_out(0.5 * pow(x, c)); }|if x < 1.0 { return tl_out(pow(x, c)); }|the left half of the log Laplace law
distributions_closed.goose|return tl_out(-expm1_raw(-c * log1p_raw(x)));|return tl_out(-expm1_raw(-log1p_raw(x)));|the shape of the Lomax law
distributions_closed.goose|return tl_out(-expm1_raw(-b * log(x)));|return tl_out(-expm1_raw(-log(x)));|the shape of the Pareto law
distributions_closed.goose|return tl_out(pow(q, 1.0 / a));|return tl_out(pow(q, a));|the quantile of the power law
distributions_closed.goose|return betainc(0.5 * c, 0.5 * c, 0.5 * (x + 1.0));|return betainc(0.5 * c, 0.5 * c, x);|the beta argument of the rdist law
distributions_closed.goose|if x <= 0.0 { return tl_out(0.5 * (1.0 - a) + (1.0 - a) / PI * atan(x / (1.0 - a))); }|if x <= 0.0 { return tl_out(0.5 * (1.0 - a) + (1.0 + a) / PI * atan(x / (1.0 - a))); }|the left branch of the skew Cauchy law
distributions_closed.goose|if x < PI { return tl_out(atan(r * tan(0.5 * x)) / PI); }|if x < PI { return tl_out(atan(r * tan(x)) / PI); }|the left branch of the wrapped Cauchy law
distributions_closed.goose|let u = 1.0 - h * p;|let u = 1.0 + h * p;|the inner term of the kappa4 law
distributions_noncentral.goose|return m / beta / (m - 1.0) * exp(-0.5 * beta * beta) + SQRT_2PI * tl_ncdf(beta);|return m / beta / m * exp(-0.5 * beta * beta) + SQRT_2PI * tl_ncdf(beta);|the constant of the crystal ball law
distributions_noncentral.goose|cdf += t;|cdf += 2.0 * t;|the theta series of the Kolmogorov limit
distributions_noncentral.goose|let v = fam_quad(0.0, x) { t => exp(kappa * (cos(t) - 1.0)) * c };|let v = fam_quad(0.0, x) { t => exp(kappa * (cos(t) - 1.0)) * 2.0 * c };|the integral of the von Mises line law
distributions_noncentral.goose|return exp(-lam + j * log(lam) - fam_lngamma(j + 1.0));|return exp(-lam + j * log(lam) - fam_lngamma(j));|the Poisson weight of the mixtures
distributions_noncentral.goose|if cp == OK { s += nc_weight(j, lam) * r * p; }|if cp == OK { s += nc_weight(j, lam) * p; }|the scale factor of the non-central F density
distributions_noncentral.goose|total += fam_quad(edges[i], edges[i + 1]) { s => tl_ncdf(t * s - nc) * exp(nct_lscale(s, df)) };|total += fam_quad(edges[i], edges[i + 1]) { s => tl_ncdf(t * s + nc) * exp(nct_lscale(s, df)) };|the sign of the non-centrality of the t law
distributions_noncentral.goose|sf += pow(1.0 - x, n);|sf += pow(1.0 - x, n - 1.0);|the first term of the one sided sum
distributions_noncentral.goose|if s > 7.24 || (s > 3.76 && n > 99.0) {|if s > 7.24 {|the asymptotic branch of the two sided law
distributions_nch.goose|let l = w + k * log(odds);|let l = w;|the odds weight of the Fisher law
distributions_nch.goose|let big = odds * (n - k) + (m - n - (d - k));|let big = odds * (n - k) + (m - n - d);|the divisor of the Wallenius integral
distributions_landau.goose|return min(max(v / PI, 0.0), 1.0);|return min(max(v, 0.0), 1.0);|the factor of the Landau cdf
distributions_landau.goose|return 1.0 - sf, OK;|return sf, OK;|the survival complement of the Landau law
distributions_search_d.goose|n.zeta = -beta * tan(0.5 * PI * alpha);|n.zeta = beta * tan(0.5 * PI * alpha);|the shift of the stable law
distributions_search_d.goose|if alpha != 1.0 { x0 = x + zeta; }|if alpha != 1.0 { x0 = x - zeta; }|the change to the form of Nolan
distributions_search_d.goose|n.c3 = 1.0 / PI;|n.c3 = 2.0 / PI;|the factor of the cumulative function of the stable law
distributions_search_d.goose|if n.alpha < 1.0 { return 0.0; }|if n.alpha < 1.0 { return inf(); }|the lower end of the helper of the stable law
distributions_search_d.goose|let bound = LEVY_X_TOL * pow(alpha, 1.0 / alpha);|let bound = LEVY_X_TOL * alpha;|the bound of the point near zeta of the stable law
distributions_search_d.goose|n.c2 = alpha / (PI * abs(alpha - 1.0) * (x0 - n.zeta));|n.c2 = alpha / (PI * abs(alpha - 1.0) * (x0 + n.zeta));|the factor of the density of the stable law
distributions_search_numeric.goose|hi = hi * 2.0 + 1.0;|hi = hi * 1.5 + 1.0;|the growth of the bracket of the search of a quantile
distributions_numeric.goose|return m2 - m1 * m1;|return m2 - m1;|the variance of a family
distributions_numeric.goose|return (m3 - 3.0 * m1 * m2 + 2.0 * m1 * m1 * m1) / (v * sqrt(v));|return (m3 - 3.0 * m1 * m2) / (v * sqrt(v));|the skewness of a family
distributions_power.goose|return 1.0 / (1.0 + pow(x, -c)), OK;|return 1.0 / (1.0 + pow(x, c)), OK;|the cdf of the log-logistic law
distributions_power.goose|return pow(1.0 + pow(x, -c), -d), OK;|return pow(1.0 + pow(x, c), -d), OK;|the cdf of the law of Burr
distributions_power.goose|return 1.0 - pow(1.0 + pow(x, c), -d), OK;|return pow(1.0 + pow(x, c), -d), OK;|the cdf of the law of Burr of the twelfth kind
distributions_power.goose|let l = (a - 1.0) * lx - (a + b) * log1p(x) - log(bt);|let l = (a) * lx - (a + b) * log1p(x) - log(bt);|the density of the beta prime law
distributions_power.goose|let b, cb = beta(1.0 + nc, 1.0 - nc);|let b, cb = beta(1.0 + nc, 1.0 + nc);|the raw moment of the log-logistic law
distributions_power.goose|let top, ct = beta(a + n as f64, b - n as f64);|let top, ct = beta(a + n as f64, b + n as f64);|the raw moment of the beta prime law
distributions_power.goose|return 1.0 / (1.0 - c) / (1.0 - c) / (1.0 - 2.0 * c), OK;|return 1.0 / (1.0 - c) / (1.0 - 2.0 * c), OK;|the variance of the generalized Pareto law
distributions_growth.goose|return exp(log(c) - mult * a - (c + 1.0) * m), OK;|return exp(log(c) - mult * a - (c + 2.0) * m), OK;|the density of the generalized logistic law
distributions_growth.goose|return exp(-c * m), OK;|return exp(-c * m) * 2.0, OK;|the cdf of the generalized logistic law
distributions_growth.goose|return 2.0 * t / ((1.0 + t2) * (1.0 + t2)), OK;|return t / ((1.0 + t2) * (1.0 + t2)), OK;|the density of the generalized half logistic law
distributions_growth.goose|return (1.0 - t2) / (1.0 + t2), OK;|return (1.0 + t2) / (1.0 - t2), OK;|the cdf of the generalized half logistic law
distributions_growth.goose|return (1.0 - r) / c, OK;|return (1.0 + r) / c, OK;|the quantile of the generalized half logistic law
distributions_growth.goose|total += c * pow(l, n as f64) * exp(-c * y) * gw[j] * width / 2.0;|total += c * pow(l, n as f64) * exp(-c * y) * gw[j] * width;|the raw moment of the Gompertz law
distributions_truncated.goose|return exp(-x) / (-d), OK;|return exp(-x) / d, OK;|the density of the truncated exponential law
distributions_truncated.goose|return e / d, OK;|return e / -d, OK;|the cdf of the truncated exponential law
distributions_truncated.goose|return (1.0 - pow(x, -b)) / d, OK;|return (1.0 - pow(x, b)) / d, OK;|the cdf of the truncated Pareto law
distributions_truncated.goose|return p / m, OK;|return p * m, OK;|the density of the truncated normal law
distributions_gamma.goose|return exp(l), OK;|return exp(-l), OK;|the density of the chi law
distributions_gamma.goose|return 0.5 + 0.5 * v, OK;|return 0.5 + v, OK;|the cdf of the double gamma law
distributions_gamma.goose|return exp((a - 1.0) * log(ax) - ax - lg) / 2.0, OK;|return exp((a - 1.0) * log(ax) - ax - lg), OK;|the density of the double gamma law
distributions_gamma.goose|return 1.0 / (a - 1.0), OK;|return 1.0 / (a + 1.0), OK;|the mean of the inverse gamma law
distributions_gamma.goose|return exp(c * x - exp(x) - lg), OK;|return exp(c * x + exp(x) - lg), OK;|the density of the log gamma law
distributions_generalized.goose|let p = exp((c * a - 1.0) * lx - pow(x, c)) * abs(c) / g;|let p = exp((c * a - 1.0) * lx - pow(x, c)) / (abs(c) * g);|the density of the generalized gamma law
distributions_generalized.goose|return pow(x, k) / pow(1.0 + pow(x, s), k / s), OK;|return pow(x, k) / pow(1.0 + pow(x, s), 2.0 * k / s), OK;|the cdf of the law of Mielke
distributions_generalized.goose|let u = -e;|let u = e;|the density of the generalized exponential law
distributions_generalized.goose|let t = a + n as f64 / c;|let t = a + n as f64 * c;|the raw moment of the generalized gamma law
distributions_symmetric.goose|return p, OK;|return 2.0 * p, OK;|the density of the generalized normal law
distributions_symmetric.goose|return v / 2.0, OK;|return v, OK;|the cdf of the generalized normal law
distributions_symmetric.goose|return 1.0 - v, OK;|return 1.0 + v, OK;|the cdf of the law of Tukey and of lambda
distributions_johnson.goose|return b * p / (x * (1.0 - x)), OK;|return b * p / (x + (1.0 - x)), OK;|the density of the Johnson law of the B kind
distributions_johnson.goose|return b * p / sqrt(x * x + 1.0), OK;|return b * p / sqrt(x * x), OK;|the density of the Johnson law of the U kind
distributions_johnson.goose|return 1.0 / (1.0 + e), OK;|return 1.0 / (1.0 - e), OK;|the quantile of the Johnson law of the B kind
distributions_johnson.goose|return johnson_sinh((z - a) / b), OK;|return johnson_sinh((z + a) / b), OK;|the quantile of the Johnson law of the U kind
distributions_special.goose|let lp = -0.5 * log(2.0 * PI) - 1.5 * log(x) - d * d / (2.0 * x);|let lp = -0.5 * log(2.0 * PI) - 1.5 * log(x) - d * d / (x * x);|the density of the inverse Gaussian law
distributions_special.goose|return max(la, lb) + l;|return max(la, lb) - l;|the cdf of the inverse Gaussian law
distributions_special.goose|return sqrt(v / nu), OK;|return sqrt(v * nu), OK;|the quantile of the Nakagami law
distributions_special.goose|return m * (1.0 - 4.0 * nu * v) / 2.0 / nu / (v * sqrt(v)), OK;|return m * (1.0 - 4.0 * nu * v) / 4.0 / nu / (v * sqrt(v)), OK;|the skewness of the Nakagami law
distributions_special.goose|let lp = log(a) - (1.0 / a + 1.0) * kappa3_log_sum(x, a);|let lp = log(a) - (1.0 / a - 1.0) * kappa3_log_sum(x, a);|the density of the Kappa law of the third kind
distributions_special.goose|return v / beta + zeta, OK;|return v * beta + zeta, OK;|the quantile of the Pearson law of the third kind
distributions_extreme.goose|return exp(-exp(l / c)), OK;|return exp(-exp(l * c)), OK;|the cdf of the generalized extreme value law
distributions_extreme.goose|let lp = -exp(t) + (1.0 / c - 1.0) * l;|let lp = -exp(t) + (1.0 / c + 1.0) * l;|the density of the generalized extreme value law
distributions_extreme.goose|return exp(-exp(l / c)), OK;|return exp(exp(l / c)), OK;|the cdf of the generalized extreme value law
distributions_skew.goose|return 2.0 * p * fam_ndtr(a * x), OK;|return p * fam_ndtr(a * x), OK;|the density of the skew normal law
distributions_skew.goose|let v = fam_ndtr(x) - 2.0 * fam_owens_t(x, a);|let v = fam_ndtr(x) + 2.0 * fam_owens_t(x, a);|the cdf of the skew normal law
distributions_skew.goose|return -log((1.0 - q) * tot * kappa) * ki, OK;|return -log((1.0 - q) * tot * kappa) * kappa, OK;|the quantile of the asymmetric Laplace law
distributions_skew.goose|return sign * (extra - fam_owens_integral(aa * h, 1.0 / aa));|return sign * (extra + fam_owens_integral(aa * h, 1.0 / aa));|the T function of Owen
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
    fitting.goose) echo "fitting fitting_mle fitting_shapes" ;;
    numeric.goose) echo "transform fitting" ;;
    descriptive.goose) echo "descriptive" ;;
    frequency.goose) echo "frequency" ;;
    entropy.goose) echo "entropy" ;;
    trimmed.goose) echo "trimmed" ;;
    circular.goose) echo "circular" ;;
    binned.goose) echo "binned" ;;
    distributions_location.goose) echo "distributions_location" ;;
    clip.goose) echo "clip" ;;
    assoc_rank.goose) echo "assoc_rank" ;;
    distance.goose) echo "distance" ;;
    distance_nd.goose) echo "distance_nd" ;;
    assoc_mgc.goose) echo "assoc_mgc" ;;
    line_siegel.goose) echo "line_siegel" ;;
    laws.goose) echo "laws" ;;
    probplot.goose) echo "probplot" ;;
    anderson_families.goose) echo "anderson_families" ;;
    assoc_somers.goose) echo "assoc_somers" ;;
    contingency_risk.goose) echo "contingency_risk" ;;
    kde.goose) echo "kde" ;;
    fitting_mle.goose) echo "fitting_mle" ;;
    fitting_shapes.goose) echo "fitting_shapes" ;;
    distributions_elementary.goose) echo "distributions_elementary" ;;
    distributions_bounded.goose) echo "distributions_bounded" ;;
    distributions_weibull.goose) echo "distributions_weibull" ;;
    distributions_numeric.goose) echo "distributions_exponent" ;;
    distributions_exponent.goose) echo "distributions_exponent" ;;
    distributions_discrete_a.goose) echo "distributions_discrete_a" ;;
    distributions_discrete_b.goose) echo "distributions_discrete_b" ;;
    distributions_rice.goose) echo "distributions_rice" ;;
    distributions_search_numeric.goose) echo "distributions_search_a" ;;
    distributions_search_a.goose) echo "distributions_search_a" ;;
    distributions_search_b.goose) echo "distributions_search_b" ;;
    distributions_search_c.goose) echo "distributions_search_c" ;;
    distributions_closed.goose) echo "distributions_closed" ;;
    distributions_noncentral.goose) echo "distributions_noncentral" ;;
    distributions_nch.goose) echo "distributions_nch" ;;
    distributions_landau.goose) echo "distributions_landau" ;;
    distributions_search_d.goose) echo "distributions_search_d" ;;
    distributions_power.goose) echo "distributions_power" ;;
    distributions_growth.goose) echo "distributions_growth" ;;
    distributions_truncated.goose) echo "distributions_truncated" ;;
    distributions_gamma.goose) echo "distributions_gamma" ;;
    distributions_generalized.goose) echo "distributions_generalized" ;;
    distributions_symmetric.goose) echo "distributions_symmetric" ;;
    distributions_johnson.goose) echo "distributions_johnson" ;;
    distributions_special.goose) echo "distributions_special" ;;
    distributions_extreme.goose) echo "distributions_extreme" ;;
    distributions_skew.goose) echo "distributions_skew" ;;
    divergence.goose) echo "divergence" ;;
    contingency.goose) echo "contingency" ;;
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
    normality.goose) echo "normality" ;;
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
