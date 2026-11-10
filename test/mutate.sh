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
# The round trip of the quantiles costs an inverse per point, so the
# probes run it on fewer points than the test program does.
export STATISTICS_ROUND_N="${STATISTICS_ROUND_N:-200}"
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
check.goose|if xs.len > 4194304 { return 9; }|if xs.len > 4194305 { return 9; }|input length bound
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
ranks.goose|let cd = n0 - t.n1 - t.n2 + t.n3 - 2 * swaps;|let cd = n0 - t.n1 - t.n2 + t.n3;|kendall counts the discordant pairs
ranks.goose|let den = nf * nf * (mf - 1.0) / mf;|let den = nf * nf * (mf - 1.0) / (mf - 1.0);|tau c denominator scale
line.goose|let slope = scale_pow2(r.n_xy / r.n_xx, r.exp_x - r.exp_y);|let slope = scale_pow2(r.s_xy / r.s_xx, r.exp_x - r.exp_y);|line slope uses the centered sums
line.goose|let slope = scale_pow2(r.s_xy / r.s_xx, r.exp_x - r.exp_y);|let slope = scale_pow2(r.n_xy / r.n_xx, r.exp_x - r.exp_y);|proportional slope goes through the origin
line.goose|let intercept_err = scale_pow2(se * sqrt(r.s_xx / n), -r.exp_y);|let intercept_err = scale_pow2(se, -r.exp_y);|intercept error carries the sum of squares
line.goose|if xs.len > MAX_THEIL { return nan(), nan(), 9; }|if xs.len > MAX_THEIL * 2 { return nan(), nan(), 9; }|theil sen pair bound
line.goose|let intercept = my - slope * mx;|let intercept = my - slope * mx + 1.0;|theil sen intercept uses the median of x
matrix.goose|if rows < 2 { return 2; }|if rows < 1 { return 2; }|matrix needs two rows
special.goose|return 1.0 - erf_series(x), OK;|return 1.0 - erf_series(abs(x)), OK;|erfc keeps the sign near zero
special.goose|if x > 27.25 { return 0.0, OK; }|if x > 26.0 { return 0.0, OK; }|erfc underflow bound
special.goose|if p <= 0.5 { return erfcinv_core(p), OK; }|if p <= 0.0 { return erfcinv_core(p), OK; }|erfcinv split of the range
special.goose|if a == 1.0 { return nan(), ERR_OVERFLOW; }|if a == 1.0 { return 1.0, OK; }|erfinv at one
special.goose|if x > 171.6243769563027 { return nan(), 8; }|if x > 170.0 { return nan(), 8; }|gamma overflow bound
special.goose|let sg = gamma_sign(n);|let sg = 1.0;|sign of the gamma and the log gamma on the negative side
special.goose|if x == -1.0 { return nan(), ERR_OVERFLOW; }|if x == -1.0 { return nan(), ERR_DOMAIN; }|log1p at minus one
special.goose|if x < -40.0 { return -1.0, OK; }|if x < -20.0 { return -1.0, OK; }|expm1 saturates at minus one
special.goose|if a > MAX_SHAPE|if a > 2.0 * MAX_SHAPE|shape bound
special.goose|if delta > 2.5 {|if delta > -0.9 {|fraction start of the large shapes
special.goose|if delta < -0.5 {|if delta < 5.0 {|series start of the large shapes
special.goose|if a >= 2500.0 {|if a >= 250000.0 {|saddle point branch start
special.goose|if hi.hi < 1e8 {|if hi.hi < 1e300 {|log beta of the large shapes
special.goose|if lo.hi >= 20.0 {|if lo.hi >= -1.0 {|log beta of two large shapes
special.goose|if target > 1e-290 {|if false {|polish of the gamma inverse
special.goose|var xn = xr - (v - target) / dens;|var xn = xr + (v - target) / dens;|direction of the polish step
special.goose|if a >= 1e100 { return dd_from(0.0); }|if a >= 1e3 { return dd_from(0.0); }|Stirling sum of a huge argument
special.goose|let n = (floor(STIRLING_SHIFT - z.hi) + 1.0) as i64;|let n = floor(STIRLING_SHIFT - z.hi) as i64;|shift of the log gamma reaches the Stirling range
special.goose|if p == 1.0 { return 1.0, 0; }|if p == 1.0 { return 0.5, 0; }|betaincinv at one
special.goose|        if upper { return 1.0, 0; }|        if upper { return 0.5, 0; }|betainc at zero
special.goose|        if upper { return 0.0, 0; }|        if upper { return 0.5, 0; }|betainc at one
special.goose|if is_nan(x) { return nan(), 3; }|if is_nan(x) { return 0.0, 3; }|betainc with a NaN
special.goose|if is_inf(p) { return ERR_DOMAIN; }|if is_inf(p) { return ERR_INFINITY; }|probability that is infinite
special.goose|while k <= MAX_ITER_SERIES {|while k <= 3 {|bound of the power series
special.goose|while i <= MAX_ITER_CF {|while i <= 3 {|bound of the gamma fraction
special.goose|while m <= MAX_ITER_CF {|while m <= 3 {|bound of the beta fraction
special.goose|let w = sqrt(2.0 * bigd);|let w = sqrt(2.0 * bigd) * 1.0000001;|saddle point variable
special.goose|ssum += bn * y[0] * e;|ssum += bn * y[0];|saddle point leading term
special.goose|if !upper \&\& n % 2 == 0 { bn = -bn; }|if upper \&\& n % 2 == 0 { bn = -bn; }|saddle point series sign of the lower tail
special.goose|if upper \&\& a < 1e-3 {|if upper \&\& a < 1e-300 {|branch of the tiny shape
special.goose|if abs(u.hi) >= 1e-3 { return 0.0, false; }|if abs(u.hi) >= 1e-300 { return 0.0, false; }|guard of the tiny shape
special.goose|return dd_mul_d(s, a);|return dd_mul_d(s, 1.0);|series of the log gamma of a tiny shape
special.goose|while k <= 100 {|while k <= 1 {|series of the tiny shape tail
special.goose|if direct > 1.0 - 1e-15 {|if direct > 2.0 {|betainc upper side of a tiny shape
special.goose|if abs(part.hi) <= 1e-33 * abs(ssum.hi) { done = true; break; }|if abs(part.hi) <= 1e300 * abs(ssum.hi) { done = true; break; }|stop of the beta series
special.goose|while n <= MAX_ITER_SERIES {|while n <= 3 {|bound of the beta series
distributions.goose|if z < 0.7 {|if z < 0.0 {|branch of the normal tail
distributions.goose|return 0.5 * q * expf_dd(dd_mul_d(dd_mul(zd, zd), -0.5));|return 0.5 * q * exp(-0.5 * z * z);|square of the normal argument
distributions.goose|return norm_upper(dd_neg(norm_zdd(x, mu, sigma))), OK;|return 1.0 - norm_upper(norm_zdd(x, mu, sigma)), OK;|the normal cumulative function is the tail itself
distributions.goose|return norm_upper(norm_zdd(x, mu, sigma)), OK;|return 1.0 - norm_upper(dd_neg(norm_zdd(x, mu, sigma))), OK;|the normal survival function is the tail itself
distributions.goose|let ph, pl = two_prod(p, wh);|let ph, pl = two_prod(p, 1.0);|the uniform quantile uses the width
distributions.goose|return s + (e + pl + p * wl), OK;|return s + e, OK;|the uniform quantile keeps the error term
distributions.goose|return expf_dd(dd_neg(expon_arg(rate, x))), OK;|return exp(-rate * x), OK;|the exponential survival uses the pair
distributions.goose|return exp_or_overflow(l);|return exp(l.hi), OK;|the lognormal density keeps its digits
distributions.goose|let ph, pl = two_prod(sigma, zh);|let ph, pl = two_prod(sigma, 1.0);|the lognormal quantile uses the scale
distributions.goose|var m3 = dd_sub(g3, dd_mul_d(dd_mul(g2, g1), 3.0));|var m3 = dd_sub(g3, dd_mul_d(dd_mul(g2, g1), 2.0));|third central moment of the Weibull
distributions.goose|m4 = dd_add(m4, dd_mul_d(dd_mul(g2, g11), 6.0));|m4 = dd_add(m4, dd_mul_d(dd_mul(g2, g11), 3.0));|fourth central moment of the Weibull
distributions.goose|return g1, g2, dd_sub(g2, dd_mul(g1, g1)), OK;|return g1, g2, dd_sub(g2, dd_mul(g1, dd_from(1.0))), OK;|variance of the Weibull
distributions.goose|return betainc_upper(a, b, x);|return betainc(b, a, 1.0 - x);|the beta survival uses the upper tail
distributions.goose|betainc_upper(0.5, 0.5 * df, y.hi)|betainc(0.5 * df, 0.5, y.hi)|the Student t tail uses the upper beta
distributions.goose|if x >= 0.0 { return 1.0 - tail, tail, OK; }|if x >= 0.0 { return tail, 1.0 - tail, OK; }|the Student t tail at a positive value
distributions.goose|let r = dd_to_f64(dd_div(m3, dd_mul(v, dd_sqrt(v))));|let r = dd_to_f64(dd_div(m3, dd_mul(v, dd_from(1.0))));|skewness of the Weibull
distributions.goose|let r = dd_to_f64(dd_div(dd_sub(m4, dd_mul_d(v2, 3.0)), v2));|let r = dd_to_f64(dd_div(m4, v2));|excess kurtosis of the Weibull
distributions.goose|return betainc_upper(kf + 1.0, n - kf, p);|return betainc(n - kf, kf + 1.0, 1.0 - p);|the binomial cumulative function uses the upper beta
distributions.goose|return betainc(kf + 1.0, n - kf, p);|return betainc_upper(n - kf, kf + 1.0, 1.0 - p);|the binomial survival uses the upper beta
distributions.goose|return gammaincc(floor(k) + 1.0, rate);|return gammainc(floor(k) + 1.0, rate);|the Poisson cumulative function uses the upper gamma
distributions.goose|return gammainc(floor(k) + 1.0, rate);|return gammaincc(floor(k) + 1.0, rate);|the Poisson survival uses the lower gamma
distributions.goose|if k < 1.0 { return 0.0, OK; }|if k < 0.0 { return 0.0, OK; }|the geometric support starts at 1
distributions.goose|return betainc(r, floor(k) + 1.0, p);|return betainc_upper(r, floor(k) + 1.0, p);|the negative binomial cumulative function uses the lower beta
distributions.goose|var r = dd_mul_d(dd_mul_d(t, n - k), d - k);|var r = dd_mul_d(dd_mul_d(t, n - k), d - k + 1.0);|ratio of the hypergeometric masses
distributions.goose|while j >= kmin {|while j > kmin {|the hypergeometric sum starts at the mode
distributions.goose|if v >= p { r = mid; } else { l = mid + 1.0; }|if v >= p { r = mid + 1.0; } else { l = mid + 1.0; }|the search of a discrete quantile
distributions.goose|if is_inf(hi) {|if false {|the search doubles its bound
distributions.goose|floor(k) != k { return ERR_DOMAIN; }|false { return ERR_DOMAIN; }|a count with a fraction is out of range
distributions.goose|if d2 <= 4.0 { return nan(), ERR_UNDEFINED; }|if d2 <= 2.0 { return nan(), ERR_UNDEFINED; }|the variance of the F family needs df2 above 4
distributions.goose|if df <= 4.0 { return nan(), ERR_UNDEFINED; }|if df <= 2.0 { return nan(), ERR_UNDEFINED; }|the kurtosis of the Student t needs df above 4
distributions.goose|return (2.0 - p) / sqrt(1.0 - p), OK;|return (2.0 - p) / sqrt(p), OK;|skewness of the geometric family
distributions.goose|return 6.0 / r + p * p / (r * (1.0 - p)), OK;|return 6.0 / r + p * p / (1.0 - p), OK;|kurtosis of the negative binomial
distributions.goose|let den = dd_mul_d(dd_mul_d(dd_from(m), m), m - 1.0);|let den = dd_mul_d(dd_mul_d(dd_from(m), m), m);|variance of the hypergeometric
distributions.goose|if kmax - kmin + 1.0 > MAX_DIST_SUPPORT {|if kmax - kmin + 1.0 > 1e18 {|bound of the hypergeometric support
distributions.goose|if z.hi > 170.0 { return dd_from(nan()), ERR_OVERFLOW; }|if z.hi > 17.0 { return dd_from(nan()), ERR_OVERFLOW; }|bound of the Weibull moment
distributions.goose|betainc_upper(0.5 * d2, 0.5 * d1, w.hi)|betainc(0.5 * d1, 0.5 * d2, u.hi)|the F survival uses the upper beta
distributions.goose|if z > 38.6 { return 0.0; }|if z > 37.0 { return 0.0; }|bound of the normal tail
distributions.goose|return base + sign * pdf * u.lo;|return base + sign * pdf * u.hi;|the rounding step of a beta ratio uses the low part
distributions.goose|if p == 0.0 { return lo, OK; }|if p == 0.0 { return lo + 1.0, OK; }|the probability 0 gives the smallest count
distributions.goose|return n * d / m, OK;|return n * n / m, OK;|mean of the hypergeometric
distributions.goose|if p <= 0.5 {|if p <= 0.0 {|the F quantile uses the upper tail above the middle
PROBES
)

# One probe of the list. The helper writes its verdict into the work
# directory, because a parallel run cannot share a counter.
run_probe() {
  local n="$1" file="$2" from="$3" to="$4" name="$5"
  dir="$work/p$n"
  mkdir -p "$dir"
  (cd "$here" && tar cf - --exclude=.git --exclude=mutate-work .) | (cd "$dir" && tar xf -)
  python3 "$here/test/mutate_one.py" "$dir/$file" "$from" "$to" || {
    echo "BROKEN PROBE $n: $name"
    echo broken > "$work/verdict$n"
    rm -rf "$dir"
    return
  }
  if (cd "$dir" && goose run test/run.goose >"$dir/out.txt" 2>&1); then
    echo "SURVIVED $n: $name"
    echo survived > "$work/verdict$n"
  fi
  rm -rf "$dir"
}

# The probes run at most four at a time. A copy of the package holds at
# most a few hundred megabytes, and the machine needs the rest.
jobs="${MUTATE_JOBS:-4}"
total=0

while IFS='|' read -r file from to name; do
  [ -z "$file" ] && continue
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
