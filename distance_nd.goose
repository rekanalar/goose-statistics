// The first Wasserstein distance of two sets of points in d dimensions.
//
// The cost of a move is the Euclidean distance. The call solves the
// transport problem exactly. The weights of each set are scaled to a sum
// of 1. The solver is the method of the shortest augmenting path: each
// step finds the cheapest path from a source with supply left to a sink
// with demand left, over the forward arcs and the reverse arcs of the
// flow so far. The step moves the largest amount that the path allows.
// The reduced costs stay at or above 0, so a dense Dijkstra search finds
// the path.
//
// Each set holds at most `MAX_NDW_POINTS` points of at most
// `MAX_NDW_DIMS` dimensions. Both sets hold points in the same layout,
// row by row. A supply or a demand below 1e-14 of the total counts as
// used up, so the distance can miss at most 1e-14 times the diameter of
// the points.
namespace statistics;

import std;
import math;

const NDW_EPS: f64 = 1e-14;

// The checks of two sets of points. The order is the bounds, the empty
// set, the length of a set, a NaN and then an infinity.
fn ndw_check(u: const f64[:], v: const f64[:], dims: i64) -> i64 {
    if dims < 1 { return ERR_DOMAIN; }
    if dims > MAX_NDW_DIMS { return ERR_TOO_LONG; }
    if u.len > MAX_NDW_POINTS * dims || v.len > MAX_NDW_POINTS * dims { return ERR_TOO_LONG; }
    if u.len == 0 || v.len == 0 { return ERR_EMPTY; }
    if u.len % dims != 0 || v.len % dims != 0 { return ERR_LENGTH; }
    for x in u { if is_nan(x) { return ERR_NAN; } }
    for x in v { if is_nan(x) { return ERR_NAN; } }
    for x in u { if is_inf(x) { return ERR_INFINITY; } }
    for x in v { if is_inf(x) { return ERR_INFINITY; } }
    return OK;
}

// The matrix of distances between the points of the first set and the
// points of the second set. The result is the fault code.
fn ndw_costs(out: f64[>..]&, u: const f64[:], v: const f64[:], dims: i64) -> i64 {
    out.resize(0);
    let n = u.len / dims;
    let m = v.len / dims;
    for i in 0..n {
        for j in 0..m {
            var s: f64 = 0.0;
            for k in 0..dims {
                let d = u[i * dims + k] - v[j * dims + k];
                s += d * d;
            }
            let r = sqrt(s);
            if is_inf(r) { return ERR_OVERFLOW; }
            out.push(r);
        }
    }
    return OK;
}

// The cheapest total cost of moving the supplies onto the demands. The
// supplies and the demands change while the solver runs.
fn ndw_solve(cost: const f64[:], n: i64, m: i64, sup: f64[>..]&, dem: f64[>..]&) -> f64, i64 {
    let nv = n + m;
    var flow: f64[>..] = [];
    for k in 0..n * m { flow.push(0.0); }
    var pot: f64[>..] = [];
    var dist: f64[>..] = [];
    var par: i64[>..] = [];
    var done: i64[>..] = [];
    for k in 0..nv { pot.push(0.0); dist.push(0.0); par.push(-1); done.push(0); }
    let cap = 20 * nv + 100;
    var steps: i64 = 0;
    while true {
        var any = false;
        for k in 0..nv { dist[k] = inf(); par[k] = -1; done[k] = 0; }
        for i in 0..n {
            if sup[i] > NDW_EPS { dist[i] = 0.0; any = true; }
        }
        if !any { break; }
        steps += 1;
        if steps > cap { return nan(), ERR_FIT; }
        // The search of the shortest paths, over the reduced costs.
        while true {
            var u: i64 = -1;
            var best = inf();
            for k in 0..nv {
                if done[k] == 0 && dist[k] < best { best = dist[k]; u = k; }
            }
            if u < 0 { break; }
            done[u] = 1;
            if u < n {
                for j in 0..m {
                    let w = n + j;
                    if done[w] == 0 {
                        var rc = cost[u * m + j] + pot[u] - pot[w];
                        if rc < 0.0 { rc = 0.0; }
                        if best + rc < dist[w] { dist[w] = best + rc; par[w] = u; }
                    }
                }
            } else {
                let j = u - n;
                for i in 0..n {
                    if done[i] == 0 && flow[i * m + j] > NDW_EPS {
                        var rc = pot[u] - pot[i] - cost[i * m + j];
                        if rc < 0.0 { rc = 0.0; }
                        if best + rc < dist[i] { dist[i] = best + rc; par[i] = u; }
                    }
                }
            }
        }
        // The sink with demand left that is the closest.
        var t: i64 = -1;
        var bt = inf();
        for j in 0..m {
            if dem[j] > NDW_EPS && dist[n + j] < bt { bt = dist[n + j]; t = n + j; }
        }
        if t < 0 { break; }
        for k in 0..nv {
            var d = dist[k];
            if d > bt { d = bt; }
            pot[k] += d;
        }
        // The amount that the path allows.
        var amt = dem[t - n];
        var node = t;
        while par[node] >= 0 {
            let p = par[node];
            if node < n {
                let f = flow[node * m + (p - n)];
                if f < amt { amt = f; }
            }
            node = p;
        }
        if sup[node] < amt { amt = sup[node]; }
        // The move along the path.
        node = t;
        while par[node] >= 0 {
            let p = par[node];
            if node >= n {
                flow[p * m + (node - n)] += amt;
            } else {
                let at = node * m + (p - n);
                flow[at] -= amt;
                if flow[at] < NDW_EPS { flow[at] = 0.0; }
            }
            node = p;
        }
        sup[node] -= amt;
        if sup[node] < NDW_EPS { sup[node] = 0.0; }
        dem[t - n] -= amt;
        if dem[t - n] < NDW_EPS { dem[t - n] = 0.0; }
    }
    var total: f64 = 0.0;
    for k in 0..n * m { total += cost[k] * flow[k]; }
    return total, OK;
}

// The shared body of the two calls. The weights are scaled to a sum of 1.
fn ndw_run(u: const f64[:], v: const f64[:], dims: i64, uw: const f64[:], vw: const f64[:], weighted: bool) -> f64, i64 {
    let c = ndw_check(u, v, dims);
    if c != OK { return nan(), c; }
    let n = u.len / dims;
    let m = v.len / dims;
    var sup: f64[>..] = [];
    var dem: f64[>..] = [];
    if weighted {
        let cu = check_weights(uw, n, 1);
        if cu != OK { return nan(), cu; }
        let cv = check_weights(vw, m, 1);
        if cv != OK { return nan(), cv; }
        var su: f64 = 0.0;
        for w in uw { su += w; }
        var sv: f64 = 0.0;
        for w in vw { sv += w; }
        if is_inf(su) || is_inf(sv) { return nan(), ERR_OVERFLOW; }
        for w in uw { sup.push(w / su); }
        for w in vw { dem.push(w / sv); }
    } else {
        for k in 0..n { sup.push(1.0 / n as f64); }
        for k in 0..m { dem.push(1.0 / m as f64); }
    }
    var cost: f64[>..] = [];
    let cc = ndw_costs(cost, u, v, dims);
    if cc != OK { return nan(), cc; }
    return ndw_solve(cost[..], n, m, sup, dem);
}

// The distance of two sets of points with equal weights. The points go
// in row by row, with `dims` values to a point.
fn wasserstein_distance_nd(u: const f64[:], v: const f64[:], dims: i64) -> f64, i64 {
    return ndw_run(u, v, dims, [], [], false);
}

// The distance of two sets of points with a weight for each point.
fn wasserstein_distance_nd_weighted(u: const f64[:], v: const f64[:], dims: i64, uw: const f64[:], vw: const f64[:]) -> f64, i64 {
    return ndw_run(u, v, dims, uw, vw, true);
}
