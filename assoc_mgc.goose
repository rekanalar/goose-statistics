// The multiscale graph correlation of two samples.
//
// The call finds a dependence of any shape between two samples that hold
// the same number of points. For each pair of neighbor counts (k, l), the
// local correlation of two distance tables says how well the k nearest
// points of the first sample match the l nearest points of the second one.
// The statistic is the largest local correlation inside the largest
// connected region of the large correlations. The p-value comes from the
// permutations of the second sample.
//
// Each distance table is centered by the sum of its column over n - 1,
// and its diagonal is 0. The ranks of a column are dense: equal distances
// share a rank. The local covariance of the scale (k, l) sums the product
// of A[i][j] and B[j][i] over the pairs with a rank of at most k in the
// first table and a rank of at most l in the second table. It subtracts
// the product of the two sums over n squared. The local variances use the
// same sum of a table with itself.
//
// A permutation of the second sample permutes the rows and the columns of
// its table. The ranks and the centering move with the permutation, so
// one pass over the pairs gives each permuted map. The local variances do
// not change.
//
// The result holds the map of the local correlations. The rows follow the
// ranks of the first table, and the columns follow the ranks of the second
// table. Both have n rows and columns when no distance repeats.
namespace statistics;

import std;
import math;

// A value of one column with its row. The sort uses the value.
struct MgcCell {
    v: f64,
    i: i64,
}

// The result of the test. The scale is the pair of neighbor counts at the
// statistic. The map has `rows` rows and `cols` columns.
struct MgcResult {
    statistic: f64,
    pvalue: f64,
    kx: i64,
    ky: i64,
    rows: i64,
    cols: i64,
}

// A result with a fault code. Every field is empty.
fn mgc_fault() -> MgcResult {
    return MgcResult { statistic: nan(), pvalue: nan(), kx: 0, ky: 0, rows: 0, cols: 0 };
}

// The Euclidean distances between the points of a sample. The points go in
// row by row. The result is the fault code.
fn mgc_distances(out: f64[>..]&, xs: const f64[:], n: i64, dims: i64) -> i64 {
    out.resize(0);
    for i in 0..n {
        for j in 0..n {
            var s: f64 = 0.0;
            for k in 0..dims {
                let d = xs[i * dims + k] - xs[j * dims + k];
                s += d * d;
            }
            let r = sqrt(s);
            if is_inf(r) { return ERR_OVERFLOW; }
            out.push(r);
        }
    }
    return OK;
}

// The centered table with a zero diagonal, and the dense ranks of each
// column. The result is the largest rank.
fn mgc_center(a: f64[>..]&, ranks: i64[>..]&, d: const f64[:], n: i64) -> i64 {
    a.resize(0);
    ranks.resize(0);
    for k in 0..n * n {
        a.push(0.0);
        ranks.push(0);
    }
    let nm = (n - 1) as f64;
    var col: MgcCell[>..] = [];
    var top: i64 = 0;
    for j in 0..n {
        var s: f64 = 0.0;
        col.resize(0);
        for i in 0..n {
            s += d[i * n + j];
            col.push(MgcCell { v: d[i * n + j], i: i });
        }
        let mean = s / nm;
        for i in 0..n {
            if i != j { a[i * n + j] = d[i * n + j] - mean; }
        }
        sort(col) { p, q => p.v < q.v };
        var rank: i64 = 0;
        for t in 0..n {
            if t == 0 || col[t].v != col[t - 1].v { rank += 1; }
            ranks[col[t].i * n + j] = rank;
        }
        if rank > top { top = rank; }
    }
    return top;
}

// The local covariances of two centered tables. The pair (i, j) holds A[i][j]
// and B[j][i], with the ranks of the first and of the transposed second
// table. The order permutes the second sample. The map has nx rows and ny
// columns.
fn mgc_cov_map(out: f64[>..]&, a: const f64[:], b: const f64[:], ra: const i64[:], rb: const i64[:], ord: const i64[:], n: i64, nx: i64, ny: i64) {
    out.resize(0);
    var ex: f64[>..] = [];
    var ey: f64[>..] = [];
    for c in 0..nx * ny { out.push(0.0); }
    for k in 0..nx { ex.push(0.0); }
    for l in 0..ny { ey.push(0.0); }
    for i in 0..n {
        for j in 0..n {
            let at = ord[j] * n + ord[i];
            let k = ra[i * n + j] - 1;
            let l = rb[at] - 1;
            let x = a[i * n + j];
            let y = b[at];
            out[k * ny + l] += x * y;
            ex[k] += x;
            ey[l] += y;
        }
    }
    for k in 0..nx {
        for l in 0..ny {
            var v = out[k * ny + l];
            if k > 0 { v += out[(k - 1) * ny + l]; }
            if l > 0 { v += out[k * ny + l - 1]; }
            if k > 0 && l > 0 { v -= out[(k - 1) * ny + l - 1]; }
            out[k * ny + l] = v;
        }
    }
    for k in 1..nx { ex[k] += ex[k - 1]; }
    for l in 1..ny { ey[l] += ey[l - 1]; }
    let nn = (n * n) as f64;
    for k in 0..nx {
        for l in 0..ny { out[k * ny + l] -= ex[k] * ey[l] / nn; }
    }
}

// The diagonal of the local covariances of a table with itself.
fn mgc_variances(out: f64[>..]&, a: const f64[:], ra: const i64[:], ident: const i64[:], n: i64, nx: i64) {
    var m: f64[>..] = [];
    mgc_cov_map(m, a, a, ra, ra, ident, n, nx, nx);
    out.resize(0);
    for k in 0..nx { out.push(m[k * nx + k]); }
}

// The local correlations. A scale with a variance of 0 or less gives 0,
// and a value above 1 gives 1.
fn mgc_corr_map(map: f64[>..]&, cov: const f64[:], vx: const f64[:], vy: const f64[:], nx: i64, ny: i64) {
    map.resize(0);
    for k in 0..nx {
        for l in 0..ny {
            var c = cov[k * ny + l] / (sqrt(vx[k]) * sqrt(vy[l]));
            if c > 1.0 { c = 1.0; }
            if vx[k] <= 0.0 || vy[l] <= 0.0 { c = 0.0; }
            map.push(c);
        }
    }
}

// The marks of the largest connected region of the cells above the
// threshold. A cell joins the four cells next to it. A tie goes to the
// region that starts first in row order. The result is the size.
fn mgc_component(keep: i64[>..]&, map: const f64[:], nx: i64, ny: i64, thr: f64) -> i64 {
    keep.resize(0);
    var lab: i64[>..] = [];
    for c in 0..nx * ny {
        lab.push(0);
        keep.push(0);
    }
    var stack: i64[>..] = [];
    var label: i64 = 0;
    var best: i64 = 0;
    var best_label: i64 = 0;
    for start in 0..nx * ny {
        if map[start] > thr && lab[start] == 0 {
            label += 1;
            var size: i64 = 0;
            stack.resize(0);
            stack.push(start);
            lab[start] = label;
            while stack.len > 0 {
                let c = stack[stack.len - 1];
                stack.resize(stack.len - 1);
                size += 1;
                let k = c / ny;
                let l = c % ny;
                if k > 0 && lab[c - ny] == 0 && map[c - ny] > thr { lab[c - ny] = label; stack.push(c - ny); }
                if k < nx - 1 && lab[c + ny] == 0 && map[c + ny] > thr { lab[c + ny] = label; stack.push(c + ny); }
                if l > 0 && lab[c - 1] == 0 && map[c - 1] > thr { lab[c - 1] = label; stack.push(c - 1); }
                if l < ny - 1 && lab[c + 1] == 0 && map[c + 1] > thr { lab[c + 1] = label; stack.push(c + 1); }
            }
            if size > best { best = size; best_label = label; }
        }
    }
    for c in 0..nx * ny {
        if lab[c] == best_label && best > 0 { keep[c] = 1; }
    }
    return best;
}

// The statistic of a map and the scale of the statistic. The result is the
// statistic, the neighbor count of the first sample and the neighbor count
// of the second sample. The threshold is the beta rule of the sample size.
fn mgc_statistic(keep: i64[>..]&, map: const f64[:], nx: i64, ny: i64, base: f64) -> f64, i64, i64 {
    let global = map[(nx - 1) * ny + ny - 1];
    if nx == 1 || ny == 1 { return global, nx, ny; }
    var thr = base;
    if global > thr { thr = global; }
    let size = mgc_component(keep, map, nx, ny, thr);
    if size == 0 { return global, nx, ny; }
    var big = ceil(0.02 * max(nx, ny) as f64);
    big = big * min(nx, ny) as f64;
    if size as f64 < big { return global, nx, ny; }
    var top = -inf();
    for c in 0..nx * ny {
        if keep[c] == 1 && map[c] > top { top = map[c]; }
    }
    if top < global { return global, nx, ny; }
    var at: i64 = 0;
    for c in 0..nx * ny {
        if keep[c] == 1 && map[c] >= top { at = c; }
    }
    return top, at / ny + 1, at % ny + 1;
}

// The test on two distance tables of n points. The tables hold n by n
// values in row order. The map goes into `map`.
fn mgc_core(map: f64[>..]&, dx: const f64[:], dy: const f64[:], n: i64, reps: i64, r: rng&) -> MgcResult, i64 {
    var a: f64[>..] = [];
    var b: f64[>..] = [];
    var ra: i64[>..] = [];
    var rb: i64[>..] = [];
    let nx = mgc_center(a, ra, dx, n);
    let ny = mgc_center(b, rb, dy, n);
    var ident: i64[>..] = [];
    for k in 0..n { ident.push(k); }
    var vx: f64[>..] = [];
    var vy: f64[>..] = [];
    mgc_variances(vx, a[..], ra[..], ident[..], n, nx);
    mgc_variances(vy, b[..], rb[..], ident[..], n, ny);
    var cov: f64[>..] = [];
    mgc_cov_map(cov, a[..], b[..], ra[..], rb[..], ident[..], n, nx, ny);
    mgc_corr_map(map, cov[..], vx[..], vy[..], nx, ny);
    // the threshold of the sample size: the beta law with shape t twice
    let s = (n - 1) as f64;
    let t = s * (s - 3.0) / 4.0 - 0.5;
    let q, cq = betaincinv(t, t, 1.0 - 0.02 / s);
    if cq != OK { return mgc_fault(), cq; }
    let base = q * 2.0 - 1.0;
    var keep: i64[>..] = [];
    let stat, kx, ky = mgc_statistic(keep, map[..], nx, ny, base);
    if is_nan(stat) { return mgc_fault(), ERR_UNDEFINED; }
    var count: i64 = 0;
    var pmap: f64[>..] = [];
    for rep in 0..reps {
        rs_shuffle_index(ident[..], 0, n, r);
        mgc_cov_map(cov, a[..], b[..], ra[..], rb[..], ident[..], n, nx, ny);
        mgc_corr_map(pmap, cov[..], vx[..], vy[..], nx, ny);
        let ps, pkx, pky = mgc_statistic(keep, pmap[..], nx, ny, base);
        if ps >= stat { count += 1; }
    }
    let p = (1 + count) as f64 / (1 + reps) as f64;
    return MgcResult { statistic: stat, pvalue: p, kx: kx, ky: ky, rows: nx, cols: ny }, OK;
}

// The checks of the options and of the work. The order is the repetitions,
// the size and then the work of the permutations.
fn mgc_check_work(n: i64, reps: i64) -> i64 {
    if reps < 0 { return ERR_DOMAIN; }
    if reps > MAX_MGC_REPS || n > MAX_MGC_PAIRS { return ERR_TOO_LONG; }
    if n * n * (1 + reps) > MAX_MGC_WORK { return ERR_TOO_LONG; }
    return OK;
}

// The checks of a sample of points. The order is the bound, the empty
// sample, the length, a NaN and then an infinity.
fn mgc_check_points(xs: const f64[:], dims: i64) -> i64 {
    if dims < 1 { return ERR_DOMAIN; }
    if xs.len > MAX_INPUT { return ERR_TOO_LONG; }
    if xs.len == 0 { return ERR_EMPTY; }
    if xs.len % dims != 0 { return ERR_LENGTH; }
    return OK;
}

// The test on two samples of points. The points go in row by row, and the
// first sample has `xdims` values to a point and the second has `ydims`.
// The map of the local correlations goes into `map`, row by row.
fn multiscale_graphcorr(map: f64[>..]&, xs: const f64[:], ys: const f64[:], xdims: i64, ydims: i64, reps: i64, r: rng&) -> MgcResult, i64 {
    map.resize(0);
    let c1 = mgc_check_points(xs, xdims);
    if c1 != OK { return mgc_fault(), c1; }
    let c2 = mgc_check_points(ys, ydims);
    if c2 != OK { return mgc_fault(), c2; }
    let n = xs.len / xdims;
    if ys.len / ydims != n { return mgc_fault(), ERR_LENGTH; }
    let cw = mgc_check_work(n, reps);
    if cw != OK { return mgc_fault(), cw; }
    for v in xs { if is_nan(v) { return mgc_fault(), ERR_NAN; } }
    for v in ys { if is_nan(v) { return mgc_fault(), ERR_NAN; } }
    for v in xs { if is_inf(v) { return mgc_fault(), ERR_INFINITY; } }
    for v in ys { if is_inf(v) { return mgc_fault(), ERR_INFINITY; } }
    if n < 5 { return mgc_fault(), ERR_TOO_SHORT; }
    var dx: f64[>..] = [];
    var dy: f64[>..] = [];
    let cx = mgc_distances(dx, xs, n, xdims);
    if cx != OK { return mgc_fault(), cx; }
    let cy = mgc_distances(dy, ys, n, ydims);
    if cy != OK { return mgc_fault(), cy; }
    return mgc_core(map, dx[..], dy[..], n, reps, r);
}

// The test on two distance tables. Each table holds n by n values in row
// order, so its length is a square. The tables need not be symmetric.
fn multiscale_graphcorr_dist(map: f64[>..]&, dx: const f64[:], dy: const f64[:], reps: i64, r: rng&) -> MgcResult, i64 {
    map.resize(0);
    if dx.len > MAX_MGC_PAIRS * MAX_MGC_PAIRS || dy.len > MAX_MGC_PAIRS * MAX_MGC_PAIRS { return mgc_fault(), ERR_TOO_LONG; }
    if dx.len == 0 || dy.len == 0 { return mgc_fault(), ERR_EMPTY; }
    if dx.len != dy.len { return mgc_fault(), ERR_LENGTH; }
    var n: i64 = 0;
    while n * n < dx.len { n += 1; }
    if n * n != dx.len { return mgc_fault(), ERR_LENGTH; }
    let cw = mgc_check_work(n, reps);
    if cw != OK { return mgc_fault(), cw; }
    for v in dx { if is_nan(v) { return mgc_fault(), ERR_NAN; } }
    for v in dy { if is_nan(v) { return mgc_fault(), ERR_NAN; } }
    for v in dx { if is_inf(v) { return mgc_fault(), ERR_INFINITY; } }
    for v in dy { if is_inf(v) { return mgc_fault(), ERR_INFINITY; } }
    if n < 5 { return mgc_fault(), ERR_TOO_SHORT; }
    return mgc_core(map, dx, dy, n, reps, r);
}

// The test of two samples of points of the same width. The test pools the
// samples and labels each point with its sample. Each sample needs at
// least 5 points.
fn multiscale_graphcorr_twosample(map: f64[>..]&, u: const f64[:], v: const f64[:], dims: i64, reps: i64, r: rng&) -> MgcResult, i64 {
    map.resize(0);
    let c1 = mgc_check_points(u, dims);
    if c1 != OK { return mgc_fault(), c1; }
    let c2 = mgc_check_points(v, dims);
    if c2 != OK { return mgc_fault(), c2; }
    let nu = u.len / dims;
    let nv = v.len / dims;
    if nu < 5 || nv < 5 {
        for x in u { if is_nan(x) { return mgc_fault(), ERR_NAN; } }
        for x in v { if is_nan(x) { return mgc_fault(), ERR_NAN; } }
        for x in u { if is_inf(x) { return mgc_fault(), ERR_INFINITY; } }
        for x in v { if is_inf(x) { return mgc_fault(), ERR_INFINITY; } }
        return mgc_fault(), ERR_TOO_SHORT;
    }
    var pooled: f64[>..] = [];
    for x in u { pooled.push(x); }
    for x in v { pooled.push(x); }
    var label: f64[>..] = [];
    for k in 0..nu { label.push(0.0); }
    for k in 0..nv { label.push(1.0); }
    return multiscale_graphcorr(map, pooled[..], label[..], dims, 1, reps, r);
}
