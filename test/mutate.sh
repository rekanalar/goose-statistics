#!/usr/bin/env bash
# Mutation probes. Each probe flips one rule in a copy of the package.
# The test program must then fail. A probe that passes means a test is
# missing. Run it from the package directory with the compiler on PATH.
# The script exits 1 when a probe survives.
set -u
here="$(cd "$(dirname "$0")/.." && pwd)"
work="${TMPDIR:-/tmp}/statistics-mutate"
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
