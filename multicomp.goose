// multicomp: the comparisons of many group means and the Bayesian intervals.
//
// tukey_hsd compares every pair of groups with the studentized range.
// dunnett compares each group with one control group. bayes_mvs and
// mvsdist give the posterior of the mean, the variance and the deviation
// of one sample. The package computes each probability by deterministic
// quadrature with fixed panels, so a call never draws a random number.
namespace statistics;

import std;
import math;

// The mean, the variance and the deviation of one sample, each with the
// mean of its posterior and an interval. The interval holds the equal
// tails of the posterior.
struct MvsResult {
    mean: f64,
    mean_interval: Interval,
    variance: f64,
    variance_interval: Interval,
    deviation: f64,
    deviation_interval: Interval,
}

// The parameters of the three posteriors of mvsdist. SciPy returns three
// distribution objects. The package has no distribution object, so this
// struct holds their parameters instead.
//
// Up to 1000 values the mean follows a Student t law with mean_df degrees
// of freedom, the location mean and the scale mean_scale. The variance
// follows an inverse gamma law with the shape variance_shape and the scale
// variance_scale. The deviation follows a generalized gamma law of the
// shape deviation_shape, the power -2 and the scale deviation_scale.
// Above 1000 values each posterior is a normal law. Then gaussian is true,
// the three shape fields and mean_df hold 0, and each scale field is the
// deviation of that law. The mean, variance and deviation fields hold the
// mean of each posterior.
struct MvsDist {
    gaussian: bool,
    mean: f64,
    mean_scale: f64,
    mean_df: f64,
    variance: f64,
    variance_scale: f64,
    variance_shape: f64,
    deviation: f64,
    deviation_scale: f64,
    deviation_shape: f64,
}

const MC_NODES: i64 = 24;          // the points of one quadrature panel
const MC_MAX_GROUPS: i64 = 64;     // the most groups of one call
const MC_MAX_PANELS: i64 = 256;    // the most panels of one integral
const MC_MAX_MVS: i64 = 1000;      // the most values of the exact posterior

// ---- the quadrature kernels -------------------------------------------

// The upper tail of the standard normal law.
fn mc_nsf(x: f64) -> f64 {
    let e, c = erfc(x * 0.7071067811865476);
    return 0.5 * e;
}

// One minus the survival factor of k terms that each miss with
// probability m: 1 - (1 - m)^k, kept exact for a small m.
fn mc_miss_all(logsum: f64) -> f64 {
    if logsum < -40.0 { return 1.0; }
    return -expm1_core(logsum);
}

// The upper tail of the range of k standard normal values at w. The
// tail is the integral of k phi(z) sf(z)^(k-1) (1 - (1 - r)^(k-1)) with
// r = sf(z + w) / sf(z). Every factor stays positive, so a tail of 1e-30
// keeps its digits. SciPy forms one minus the cdf and loses them.
fn mc_range_tail(w: f64, k: i64, gx: const f64[:], gw: const f64[:]) -> f64 {
    if w <= 0.0 { return 1.0; }
    let km1 = (k - 1) as f64;
    let kf = k as f64;
    let lo = -(0.5 * w + 12.0);
    let hi = 12.0;
    let np = ceil((hi - lo) / 6.0) as i64;
    let h = (hi - lo) / np;
    var total: f64 = 0.0;
    for p in 0..np {
        let mid = lo + h * (p + 0.5);
        for i in 0..gx.len {
            let z = mid + 0.5 * h * gx[i];
            let sz = mc_nsf(z);
            let sw = mc_nsf(z + w);
            let r = sw / sz;
            var f: f64 = 1.0;
            if r < 1.0 { f = mc_miss_all(km1 * log1p_core(-r)); }
            let dens = 0.3989422804014327 * exp(-0.5 * z * z);
            total += gw[i] * 0.5 * h * kf * dens * pow(sz, km1) * f;
        }
    }
    return total;
}

// The nodes and the weights for the law of u = chi / sqrt(df), where chi
// is a chi variable with df degrees of freedom. The weight holds the
// density and the width of the panel. The grid covers 14 deviations around
// 1. When that range reaches 0, the grid runs over log(u) from the point
// where the density is below 1e-20 of its peak, so a small df keeps its
// accuracy.
fn mc_ugrid(df: f64, gx: const f64[:], gw: const f64[:], us: f64[>..]&, wu: f64[>..]&) {
    us.resize(0);
    wu.resize(0);
    let lg, sg, cl = lgamma(0.5 * df);
    let logc = 0.5 * df * log(df) + (1.0 - 0.5 * df) * log(2.0) - lg;
    let sd = 1.0 / sqrt(2.0 * df);
    let lo = 1.0 - 14.0 * sd;
    let hi = 1.0 + 14.0 * sd;
    let logs = lo <= 0.0;
    var a = lo;
    var b = hi;
    if logs {
        a = -50.0 / df;
        b = log(hi);
    }
    var np = ceil((b - a) / (5.0 * sd)) as i64;
    if np > MC_MAX_PANELS { np = MC_MAX_PANELS; }
    if np < 1 { np = 1; }
    let h = (b - a) / np;
    for p in 0..np {
        let mid = a + h * (p + 0.5);
        for i in 0..gx.len {
            let x = mid + 0.5 * h * gx[i];
            var u = x;
            var jac: f64 = 1.0;
            if logs { u = exp(x); jac = u; }
            let wt = gw[i] * 0.5 * h * jac * exp(logc + (df - 1.0) * log(u) - 0.5 * df * u * u);
            if wt > 0.0 { us.push(u); wu.push(wt); }
        }
    }
}

// The upper tail of the studentized range at q for k groups.
fn mc_range_sf(q: f64, k: i64, us: const f64[:], wu: const f64[:], gx: const f64[:], gw: const f64[:]) -> f64 {
    if q <= 0.0 { return 1.0; }
    var total: f64 = 0.0;
    for i in 0..us.len {
        total += wu[i] * mc_range_tail(q * us[i], k, gx, gw);
    }
    return min(1.0, max(0.0, total));
}

// The point q where the upper tail of the studentized range equals 1 -
// level. The call brackets q by doubling and then runs Brent, as SciPy does.
fn mc_range_ppf(level: f64, k: i64, us: const f64[:], wu: const f64[:], gx: const f64[:], gw: const f64[:]) -> f64, i64 {
    let alpha = 1.0 - level;
    return mc_critical(4.0, alpha, 1e-14, false) { x => mc_range_sf(x, k, us, wu, gx, gw) };
}

// The critical value of a test: the point where the chance of the
// exceedance equals alpha. The call doubles an upper bound until the tail
// lies below alpha, and then runs Brent on the tail. The block gives the
// chance of the exceedance at a point.
fn mc_critical<F>(hi0: f64, alpha: f64, xtol: f64, one_sided: bool) -> f64, i64 {
    var hi = hi0;
    var steps: i64 = 0;
    while F(hi) > alpha {
        hi *= 2.0;
        steps += 1;
        if steps > 40 { return nan(), ERR_FIT; }
    }
    var lo: f64 = 0.0;
    if one_sided { lo = -hi; }
    return brent_root(lo, hi, xtol, BRENT_RTOL) { x => F(x) / alpha - 1.0 };
}

// The checks that both group calls share: the groups, the size of each
// group and the count of groups.
fn mc_check_groups(groups: const f64[:], sizes: const i64[:]) -> i64 {
    let c = check_groups(groups, sizes);
    if c != OK { return c; }
    for s in sizes {
        if s < 2 { return ERR_TOO_SHORT; }
    }
    if sizes.len > MC_MAX_GROUPS { return ERR_TOO_LONG; }
    return OK;
}

// The mean of each group and the variance of each group with n - 1.
fn mc_moments(groups: const f64[:], sizes: const i64[:], means: f64[>..]&, vars: f64[>..]&) -> i64 {
    means.resize(0);
    vars.resize(0);
    var at: i64 = 0;
    for n in sizes {
        let g = groups[at..at + n];
        let m, cm = mean(g);
        if cm != OK { return cm; }
        let v, cv = variance_core(g, 1);
        if cv != OK { return cv; }
        means.push(m);
        vars.push(v);
        at += n;
    }
    return OK;
}

// The Tukey honest significant difference test of all pairs of groups.
// The groups sit one after the other in groups, as in f_oneway. Each pair
// i < j in the order (0,1), (0,2), ..., (k-2,k-1) gives one row of four
// numbers in out: the difference of the means i minus j, the p-value, and
// the low and the high end of the interval of the confidence level. The
// call returns the count of rows. With equal_var false the call runs the
// Games-Howell test: it uses the variance of each pair and the Welch
// degrees of freedom. Each group needs two values. A sample without any
// spread gives the undefined fault.
//
// SciPy: tukey_hsd(*samples, equal_var=True) and the method
// confidence_interval(confidence_level). SciPy computes the p-value as one
// minus the cdf, so a p-value below 1e-12 loses its digits. This call
// integrates the upper tail directly and keeps them.
fn tukey_hsd(out: f64[>..]&, groups: const f64[:], sizes: const i64[:], level: f64 = 0.95, equal_var: bool = true) -> i64, i64 {
    out.resize(0);
    var c = mc_check_groups(groups, sizes);
    if c != OK { return 0, c; }
    c = check_level(level);
    if c != OK { return 0, c; }
    var means: f64[>..] = [];
    var vars: f64[>..] = [];
    c = mc_moments(groups, sizes, means, vars);
    if c != OK { return 0, c; }
    let k = sizes.len;
    var gx: f64[>..] = [];
    var gw: f64[>..] = [];
    c = gauss_legendre(gx, gw, MC_NODES);
    if c != OK { return 0, c; }
    var us: f64[>..] = [];
    var wu: f64[>..] = [];
    var mse: f64 = 0.0;
    var crit: f64 = 0.0;
    if equal_var {
        var total: i64 = 0;
        for i in 0..k {
            mse += (sizes[i] - 1) * vars[i];
            total += sizes[i];
        }
        let df = (total - k) as f64;
        mse /= df;
        if is_inf(mse) { return 0, ERR_OVERFLOW; }
        if mse <= 0.0 { return 0, ERR_UNDEFINED; }
        mc_ugrid(df, gx[..], gw[..], us, wu);
        let q, cq = mc_range_ppf(level, k, us[..], wu[..], gx[..], gw[..]);
        if cq != OK { return 0, cq; }
        crit = q;
    }
    var rows: i64 = 0;
    for i in 0..k - 1 {
        for j in i + 1..k {
            var se: f64 = 0.0;
            var q = crit;
            let ni = sizes[i];
            let nj = sizes[j];
            if equal_var {
                se = sqrt((1.0 / ni + 1.0 / nj) * mse / 2.0);
            } else {
                let a = vars[i] / ni;
                let b = vars[j] / nj;
                let v = a + b;
                if v <= 0.0 { out.resize(0); return 0, ERR_UNDEFINED; }
                let den = a * a / (ni - 1) + b * b / (nj - 1);
                let df = v * v / den;
                if is_bad(df) || df <= 0.0 { out.resize(0); return 0, ERR_OVERFLOW; }
                se = sqrt(v / 2.0);
                mc_ugrid(df, gx[..], gw[..], us, wu);
                let qq, cq = mc_range_ppf(level, k, us[..], wu[..], gx[..], gw[..]);
                if cq != OK { out.resize(0); return 0, cq; }
                q = qq;
            }
            let stat = means[i] - means[j];
            let p = mc_range_sf(abs(stat) / se, k, us[..], wu[..], gx[..], gw[..]);
            let radius = q * se;
            out.push(stat);
            out.push(p);
            out.push(stat - radius);
            out.push(stat + radius);
            rows += 1;
        }
    }
    for v in out {
        if is_bad(v) { out.resize(0); return 0, ERR_OVERFLOW; }
    }
    return rows, OK;
}

// ---- Dunnett ----------------------------------------------------------

// The chance that at least one comparison passes the bound b. The
// statistics share one chi factor u and one normal factor y, because the
// correlation of two statistics has the form l_i l_j. Given u and y the
// statistics are independent, so the chance is a two-dimensional
// integral. The lists ls, ss and ms hold the distinct l, sqrt(1 - l^2)
// and the count of each l.
fn mc_dunnett_p(b: f64, alt: i64, ls: const f64[:], ss: const f64[:], ms: const f64[:], us: const f64[:], wu: const f64[:], gx: const f64[:], gw: const f64[:]) -> f64 {
    var bound = b;
    if alt == ALT_LESS { bound = -b; }
    let two = alt == ALT_TWO_SIDED;
    if two { bound = abs(b); }
    var total: f64 = 0.0;
    for iu in 0..us.len {
        let v = bound * us[iu];
        if v > 38.0 { continue; }
        if v < -38.0 { total += wu[iu]; continue; }
        let half = 12.0 + abs(v);
        let np = ceil(half / 3.0) as i64;
        let h = 2.0 * half / np;
        var inner: f64 = 0.0;
        for p in 0..np {
            let mid = -half + h * (p + 0.5);
            for i in 0..gx.len {
                let y = mid + 0.5 * h * gx[i];
                var acc: f64 = 0.0;
                var full = false;
                for j in 0..ls.len {
                    var miss = mc_nsf((v - ls[j] * y) / ss[j]);
                    if two { miss += mc_nsf((v + ls[j] * y) / ss[j]); }
                    if miss >= 1.0 { full = true; break; }
                    acc += ms[j] * log1p_core(-miss);
                }
                var f: f64 = 1.0;
                if !full { f = mc_miss_all(acc); }
                inner += gw[i] * 0.5 * h * 0.3989422804014327 * exp(-0.5 * y * y) * f;
            }
        }
        total += wu[iu] * inner;
    }
    return min(1.0, max(0.0, total));
}

// Dunnett's test of every group against the control group. The groups sit
// one after the other in groups, as in f_oneway, and control is the index
// of the control group in sizes. Each other group, in order, gives one row
// of four numbers in out: the t statistic, the p-value, and the low and
// the high end of the interval of the difference of the means, group
// minus control, at the level of 0.95. The call returns the count of
// rows. alt is an ALT_ constant. The interval of the alternative greater
// has +inf as its high end and the interval of the alternative less has
// -inf as its low end, as in SciPy. Each group needs two values. A sample
// without any spread gives the undefined fault.
//
// SciPy: dunnett(*samples, control=..., alternative=..., rng=...). SciPy
// draws a randomized quasi Monte Carlo with no fixed seed. Its p-value
// carries an error of about 1.5e-4 and changes from call to call. This
// call integrates the same probability deterministically, so the same
// input always gives the same output. SciPy leaves the sign of the
// critical value of the one-sided interval to its search. This call
// takes the critical value that gives the interval its meaning.
fn dunnett(out: f64[>..]&, groups: const f64[:], sizes: const i64[:], control: i64, alt: i64 = ALT_TWO_SIDED, level: f64 = 0.95) -> i64, i64 {
    out.resize(0);
    var c = mc_check_groups(groups, sizes);
    if c != OK { return 0, c; }
    let k = sizes.len;
    if control < 0 || control >= k { return 0, ERR_DOMAIN; }
    if check_alt(alt) != OK { return 0, ERR_DOMAIN; }
    c = check_level(level);
    if c != OK { return 0, c; }
    var means: f64[>..] = [];
    var vars: f64[>..] = [];
    c = mc_moments(groups, sizes, means, vars);
    if c != OK { return 0, c; }
    var total: i64 = 0;
    var ss: f64 = 0.0;
    for i in 0..k {
        total += sizes[i];
        ss += (sizes[i] - 1) * vars[i];
    }
    let df = (total - k) as f64;
    let s2 = ss / df;
    if is_inf(s2) { return 0, ERR_OVERFLOW; }
    if s2 <= 0.0 { return 0, ERR_UNDEFINED; }
    let sd = sqrt(s2);
    let nc = sizes[control];
    var ts: f64[>..] = [];
    var scale: f64[>..] = [];
    var ls: f64[>..] = [];
    var rs: f64[>..] = [];
    var ms: f64[>..] = [];
    for i in 0..k {
        if i == control { continue; }
        let f = sqrt(1.0 / sizes[i] + 1.0 / nc);
        ts.push((means[i] - means[control]) / f / sd);
        scale.push(f * sd);
        let l = 1.0 / sqrt(1.0 * nc / sizes[i] + 1.0);
        var seen = false;
        for j in 0..ls.len {
            if ls[j] == l { ms[j] += 1.0; seen = true; break; }
        }
        if !seen { ls.push(l); rs.push(sqrt(1.0 - l * l)); ms.push(1.0); }
    }
    var gx: f64[>..] = [];
    var gw: f64[>..] = [];
    c = gauss_legendre(gx, gw, MC_NODES);
    if c != OK { return 0, c; }
    var us: f64[>..] = [];
    var wu: f64[>..] = [];
    mc_ugrid(df, gx[..], gw[..], us, wu);
    // The critical value of the interval: the bound where the chance of
    // the exceedance equals 1 - level. The less form uses the greater
    // form, because the two are mirror images.
    let alpha = 1.0 - level;
    var calt = alt;
    if alt == ALT_LESS { calt = ALT_GREATER; }
    let one_sided = calt != ALT_TWO_SIDED;
    let crit, cc = mc_critical(2.0, alpha, 1e-12, one_sided) { x => mc_dunnett_p(x, calt, ls[..], rs[..], ms[..], us[..], wu[..], gx[..], gw[..]) };
    if cc != OK { return 0, cc; }
    var rows: i64 = 0;
    var at: i64 = 0;
    for i in 0..k {
        if i == control { continue; }
        let t = ts[at];
        let p = mc_dunnett_p(t, alt, ls[..], rs[..], ms[..], us[..], wu[..], gx[..], gw[..]);
        let diff = means[i] - means[control];
        let half = crit * scale[at];
        out.push(t);
        out.push(p);
        if alt == ALT_LESS { out.push(-inf()); } else { out.push(diff - half); }
        if alt == ALT_GREATER { out.push(inf()); } else { out.push(diff + half); }
        rows += 1;
        at += 1;
    }
    for r in 0..rows {
        for j in 0..2 {
            let v = out[4 * r + j];
            if is_bad(v) { out.resize(0); return 0, ERR_OVERFLOW; }
        }
    }
    return rows, OK;
}

// ---- Bayesian mean, variance and deviation ----------------------------


fn mc_bad_dist() -> MvsDist {
    return MvsDist { gaussian: false, mean: nan(), mean_scale: nan(), mean_df: nan(), variance: nan(), variance_scale: nan(), variance_shape: nan(), deviation: nan(), deviation_scale: nan(), deviation_shape: nan() };
}

// The posterior parameters of a sample, after the call checked it.
fn mc_dist(xs: const f64[:]) -> MvsDist, i64 {
    let n = xs.len;
    let xbar, cm = mean(xs);
    if cm != OK { return mc_bad_dist(), cm; }
    let cv, cc = variance_core(xs, 0);
    if cc != OK { return mc_bad_dist(), cc; }
    if cv <= 0.0 { return mc_bad_dist(), ERR_UNDEFINED; }
    let nf = n as f64;
    var d = mc_bad_dist();
    d.mean = xbar;
    if n > MC_MAX_MVS {
        d.gaussian = true;
        d.mean_scale = sqrt(cv / nf);
        d.mean_df = 0.0;
        d.variance = cv;
        d.variance_scale = sqrt(2.0 / nf) * cv;
        d.variance_shape = 0.0;
        d.deviation = sqrt(cv);
        d.deviation_scale = sqrt(cv / (2.0 * nf));
        d.deviation_shape = 0.0;
    } else {
        // The mean of the variance needs 3 < n and the mean of the
        // deviation needs 2 < n. SciPy gives infinity or NaN below that.
        if n < 4 { return mc_bad_dist(), ERR_UNDEFINED; }
        let a = 0.5 * (nf - 1.0);
        let vscale = 0.5 * nf * cv;
        let l1, s1, c1 = lgamma(a - 0.5);
        let l2, s2, c2 = lgamma(a);
        d.mean_scale = sqrt(cv / (nf - 1.0));
        d.mean_df = nf - 1.0;
        d.variance_shape = a;
        d.variance_scale = vscale;
        d.variance = vscale / (a - 1.0);
        d.deviation_shape = a;
        d.deviation_scale = sqrt(vscale);
        d.deviation = d.deviation_scale * exp(l1 - l2);
    }
    if is_bad(d.mean) || is_bad(d.variance) || is_bad(d.deviation) { return mc_bad_dist(), ERR_OVERFLOW; }
    if is_bad(d.mean_scale) || is_bad(d.variance_scale) || is_bad(d.deviation_scale) { return mc_bad_dist(), ERR_OVERFLOW; }
    if d.mean_scale <= 0.0 || d.variance_scale <= 0.0 || d.deviation_scale <= 0.0 {
        return mc_bad_dist(), ERR_OVERFLOW;
    }
    return d, OK;
}

// The parameters of the posterior of the mean, the variance and the
// deviation of a sample, with the mean of each posterior. It needs four
// values: the mean of the variance does not exist below that, and the
// call gives the undefined fault there where SciPy gives infinity or NaN.
// A sample without any spread gives the undefined fault too.
//
// SciPy: mvsdist(data).
fn mvsdist(xs: const f64[:]) -> MvsDist, i64 {
    let c = check_sample(xs, 2);
    if c != OK { return mc_bad_dist(), c; }
    return mc_dist(xs);
}

// The interval of the equal tails of the level for a posterior of the
// kind that mc_dist describes. part is 0 for the mean, 1 for the
// variance and 2 for the deviation.
fn mc_interval(d: MvsDist, part: i64, level: f64) -> Interval, i64 {
    var qs: f64[2] = [0.5 * (1.0 - level), 0.5 * (1.0 + level)];
    var vals: f64[2] = [0.0, 0.0];
    for i in 0..2 {
        let q = qs[i];
        if d.gaussian {
            let z = norm_std_ppf(q);
            if part == 0 { vals[i] = d.mean + d.mean_scale * z; }
            else if part == 1 { vals[i] = d.variance + d.variance_scale * z; }
            else { vals[i] = d.deviation + d.deviation_scale * z; }
        } else if part == 0 {
            let t, ct = tdist_ppf(q, d.mean_df);
            if ct != OK { return Interval { low: nan(), high: nan() }, ct; }
            vals[i] = d.mean + d.mean_scale * t;
        } else {
            let shape = d.variance_shape;
            let g, cg = gammainccinv(shape, q);
            if cg != OK { return Interval { low: nan(), high: nan() }, cg; }
            if g <= 0.0 { return Interval { low: nan(), high: nan() }, ERR_OVERFLOW; }
            if part == 1 { vals[i] = d.variance_scale / g; }
            else { vals[i] = d.deviation_scale / sqrt(g); }
        }
        if is_bad(vals[i]) { return Interval { low: nan(), high: nan() }, ERR_OVERFLOW; }
    }
    return Interval { low: vals[0], high: vals[1] }, OK;
}

// The Bayesian estimate of the mean, the variance and the deviation of a
// sample: the mean of each posterior and the equal-tail interval of the
// level. The level lies strictly between 0 and 1. SciPy uses 0.90 as the
// default. The rules of mvsdist hold for the sample.
//
// SciPy: bayes_mvs(data, alpha=0.90) returns three results of the form
// (statistic, minmax).
fn bayes_mvs(xs: const f64[:], level: f64 = 0.9) -> MvsResult, i64 {
    let bad = MvsResult { mean: nan(), mean_interval: Interval { low: nan(), high: nan() }, variance: nan(), variance_interval: Interval { low: nan(), high: nan() }, deviation: nan(), deviation_interval: Interval { low: nan(), high: nan() } };
    var c = check_sample(xs, 2);
    if c != OK { return bad, c; }
    c = check_level(level);
    if c != OK { return bad, c; }
    let d, cd = mc_dist(xs);
    if cd != OK { return bad, cd; }
    let im, c1 = mc_interval(d, 0, level);
    if c1 != OK { return bad, c1; }
    let iv, c2 = mc_interval(d, 1, level);
    if c2 != OK { return bad, c2; }
    let ds, c3 = mc_interval(d, 2, level);
    if c3 != OK { return bad, c3; }
    return MvsResult { mean: d.mean, mean_interval: im, variance: d.variance, variance_interval: iv, deviation: d.deviation, deviation_interval: ds }, OK;
}
