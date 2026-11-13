// The rank association calls of the package: the Spearman rho, the
// point biserial correlation, the weighted Kendall tau and the xi
// correlation of Chatterjee.
//
// Each call returns a CorrelationResult and then a fault code. A fault
// gives NaN in every field of the result.
namespace statistics;

import std;
import math;

// The weight of a rank in the weighted tau. The rank 0 has the weight 1.
fn wt_weigher(rank: f64) -> f64 {
    return 1.0 / (1.0 + rank);
}

// The order of the pairs by the first list, then by the second list. The
// sort is stable, so a tie keeps the order of the input. The result
// holds the index of every pair.
fn wt_perm(out: i64[>..]&, xs: const f64[:], ys: const f64[:]) {
    out.resize(0);
    var idx: Vi[>..] = [];
    for k in 0..xs.len { idx.push(Vi { v: xs[k], i: k }); }
    stable_sort(idx) { a, b => a.v < b.v };
    // The sort above orders by the first list alone. The pairs of one
    // value of the first list are then ordered by the second list.
    var out2: i64[>..] = [];
    var k: i64 = 0;
    while k < idx.len {
        var j = k + 1;
        while j < idx.len && idx[j].v == idx[k].v { j += 1; }
        var sub: Vi[>..] = [];
        for m in k..j { sub.push(Vi { v: ys[idx[m].i], i: idx[m].i }); }
        stable_sort(sub) { a, b => a.v < b.v };
        for m in 0..sub.len { out2.push(sub[m].i); }
        k = j;
    }
    for v in out2 { out.push(v); }
}

// The ranks of the pairs in the order of the permutation. A rank counts
// down from n - 1, so the last pair of the order has the rank 0.
fn wt_ranks(out: f64[>..]&, perm: const i64[:]) {
    out.resize(0);
    let n = perm.len;
    for k in 0..n { out.push(0.0); }
    for k in 0..n { out[perm[k]] = n - 1 - k; }
}

// The dense ranks of a list of ranks. The equal values share a rank,
// and the ranks count the distinct values in increasing order.
fn wt_dense_ranks(out: f64[>..]&, ranks: const f64[:]) -> i64 {
    out.resize(0);
    let n, c = rankdata(out, ranks, TIE_DENSE);
    if c != OK { return c; }
    for k in 0..n { out[k] = out[k] - 1.0; }
    return OK;
}

// The weight of the group of the pairs that share a value. The additive
// form gives the sum of the weights times the count of the other members,
// and the other form gives the sum of the products of the weights.
fn wt_group_weight(perm: const i64[:], ranks: const f64[:], first: i64, last: i64, additive: bool) -> f64 {
    var s: f64 = 0.0;
    var sq: f64 = 0.0;
    var k = first;
    while k < last {
        let w = wt_weigher(ranks[perm[k]]);
        s += w;
        sq += w * w;
        k += 1;
    }
    let count = (last - first);
    if additive { return s * (count - 1.0); }
    return (s * s - sq) / 2.0;
}

// The key of a group of pairs. The mode is 0 for the pairs that share a
// value of both lists, 1 for the pairs that share a value of the first
// list and 2 for the pairs that share a value of the second list.
fn wt_same_key(xs: const f64[:], ys: const f64[:], perm: const i64[:], a: i64, b: i64, mode: i64) -> bool {
    if mode == 1 { return xs[perm[a]] == xs[perm[b]]; }
    if mode == 2 { return ys[perm[a]] == ys[perm[b]]; }
    return xs[perm[a]] == xs[perm[b]] && ys[perm[a]] == ys[perm[b]];
}

// The sum of the weights of the groups of the pairs that share a value.
// The mode gives the key of a group.
fn wt_group_sum(xs: const f64[:], ys: const f64[:], perm: const i64[:], ranks: const f64[:], mode: i64, additive: bool) -> f64 {
    let n = perm.len;
    var total: f64 = 0.0;
    var first: i64 = 0;
    var k: i64 = 1;
    while k < n {
        if !wt_same_key(xs, ys, perm, first, k, mode) {
            total += wt_group_weight(perm, ranks, first, k, additive);
            first = k;
        }
        k += 1;
    }
    total += wt_group_weight(perm, ranks, first, n, additive);
    return total;
}

// The weight of the groups of the pairs that share a value of the first
// list and of the second list, of the pairs that share a value of the
// first list only, and of the pairs that share a value of the second
// list only. The three sums go into the three outputs.
fn wt_ties(xs: const f64[:], ys: const f64[:], perm: const i64[:], ranks: const f64[:], joint: f64&, same_x: f64&, same_y: f64&, additive: bool) {
    joint = wt_group_sum(xs, ys, perm, ranks, 0, additive);
    same_x = wt_group_sum(xs, ys, perm, ranks, 1, additive);
    same_y = wt_group_sum(xs, ys, perm, ranks, 2, additive);
}

// The work of the merge sort of the weighted tau. The buffers of the
// merge sort are passed by reference, because a recursive call cannot
// hold a growable value.
struct WtWork {
    xs: const f64[:],
    ys: const f64[:],
    ranks: const f64[:],
    exchanges: f64,
    additive: bool,
}

// The weight of the range of the permutation and the weight of the
// discordant pairs of the range. The merge sort orders the range by the
// second list. The weight of a discordant pair is the sum of the weights
// of its members in the additive form, and the product of them in the
// other form.
recursive fn wt_weigh(w: WtWork&, perm: i64[>..]&, temp: i64[>..]&, offset: i64, length: i64) -> f64 {
    if length == 1 { return wt_weigher(w.ranks[perm[offset]]); }
    let length0 = length / 2;
    let length1 = length - length0;
    let middle = offset + length0;
    let residual = wt_weigh(w, perm, temp, offset, length0);
    let weight = wt_weigh(w, perm, temp, middle, length1) + residual;
    if w.ys[perm[middle - 1]] < w.ys[perm[middle]] { return weight; }
    var i: i64 = 0;
    var j: i64 = 0;
    var k: i64 = 0;
    var left = residual;
    while j < length0 && k < length1 {
        if w.ys[perm[offset + j]] <= w.ys[perm[middle + k]] {
            temp[i] = perm[offset + j];
            left -= wt_weigher(w.ranks[temp[i]]);
            j += 1;
        } else {
            temp[i] = perm[middle + k];
            let wr = wt_weigher(w.ranks[temp[i]]);
            if w.additive { w.exchanges += wr * (length0 - j) + left; }
            else { w.exchanges += wr * left; }
            k += 1;
        }
        i += 1;
    }
    while j < length0 {
        temp[i] = perm[offset + j];
        left -= wt_weigher(w.ranks[temp[i]]);
        j += 1;
        i += 1;
    }
    while k < length1 {
        temp[i] = perm[middle + k];
        k += 1;
        i += 1;
    }
    for m in 0..i { perm[offset + m] = temp[m]; }
    return weight;
}

// The weighted Kendall tau of two samples. The ranks come from the order
// of the pairs, or from the list of ranks of the caller. The additive
// form adds the weights of the two members of a pair, and the other form
// multiplies them. The call gives no p-value, because the null
// distribution of the statistic is unknown.
fn wt_one(xs: const f64[:], ys: const f64[:], ranks: const f64[:], use_ranks: bool, additive: bool) -> f64, i64 {
    let c = check_pair(xs, ys, 2);
    if c != OK { return nan(), c; }
    let n = xs.len;
    var perm: i64[>..] = [];
    wt_perm(perm, xs, ys);
    var own: f64[>..] = [];
    if use_ranks {
        if ranks.len != n { return nan(), ERR_LENGTH; }
        let cr = wt_dense_ranks(own, ranks);
        if cr != OK { return nan(), cr; }
    } else {
        wt_ranks(own, perm);
    }
    var temp: i64[>..] = [];
    for k in 0..n { temp.push(0); }
    var w = WtWork { xs: xs, ys: ys, ranks: own, exchanges: 0.0, additive: additive };
    var joint: f64 = 0.0;
    var same_x: f64 = 0.0;
    var same_y: f64 = 0.0;
    // The joint and first-list groups use the order before the merge sort.
    // The second-list groups use the order after it.
    wt_ties(xs, ys, perm, own, joint, same_x, same_y, additive);
    wt_weigh(w, perm, temp, 0, n);
    var ignore_joint: f64 = 0.0;
    var ignore_x: f64 = 0.0;
    wt_ties(xs, ys, perm, own, ignore_joint, ignore_x, same_y, additive);
    var s: f64 = 0.0;
    var sq: f64 = 0.0;
    for k in 0..n {
        let ww = wt_weigher(own[perm[k]]);
        s += ww;
        sq += ww * ww;
    }
    var tot: f64 = 0.0;
    if additive { tot = s * (n - 1); }
    else { tot = (s * s - sq) / 2.0; }
    let num = (tot - (same_y + same_x - joint)) - 2.0 * w.exchanges;
    let den = sqrt(tot - same_x) * sqrt(tot - same_y);
    if den == 0.0 { return nan(), ERR_UNDEFINED; }
    var tau = num / den;
    if is_nan(tau) { return nan(), ERR_UNDEFINED; }
    if tau > 1.0 { tau = 1.0; }
    if tau < -1.0 { tau = -1.0; }
    return tau, OK;
}

// The weighted Kendall tau with the ranks of the caller is one direction.
// Without ranks the call averages the two directions of the lexical order.
fn weightedtau_core(xs: const f64[:], ys: const f64[:], ranks: const f64[:], use_ranks: bool, additive: bool) -> f64, i64 {
    if use_ranks { return wt_one(xs, ys, ranks, true, additive); }
    let a, ca = wt_one(xs, ys, ranks, false, additive);
    if ca != OK { return nan(), ca; }
    let b, cb = wt_one(ys, xs, ranks, false, additive);
    if cb != OK { return nan(), cb; }
    return (a + b) / 2.0, OK;
}

// The weighted Kendall tau of two samples. The ranks come from the order
// of the pairs.
fn weightedtau(xs: const f64[:], ys: const f64[:], additive: bool = true) -> f64, i64 {
    let none: f64[0] = [];
    return weightedtau_core(xs, ys, none[..], false, additive);
}

// The weighted Kendall tau with the ranks of the caller. The equal ranks
// share a rank, and the ranks count the distinct values in increasing
// order.
fn weightedtau_ranked(xs: const f64[:], ys: const f64[:], ranks: const f64[:], additive: bool = true) -> f64, i64 {
    return weightedtau_core(xs, ys, ranks, true, additive);
}

// The Spearman rho of two samples. The p-value comes from the beta law
// of the correlation of the ranks.
fn spearmanrho(xs: const f64[:], ys: const f64[:], alt: i64 = ALT_TWO_SIDED) -> CorrelationResult, i64 {
    let ca = check_alt(alt);
    let c = check_pair(xs, ys, 2);
    if c != OK { return corr_fault(c), c; }
    if ca != OK { return corr_fault(ca), ca; }
    var rx: f64[>..] = [];
    var ry: f64[>..] = [];
    let n1, c1 = rankdata(rx, xs, TIE_AVERAGE);
    if c1 != OK { return corr_fault(c1), c1; }
    let n2, c2 = rankdata(ry, ys, TIE_AVERAGE);
    if c2 != OK { return corr_fault(c2), c2; }
    let r, c3 = correlation(rx[..], ry[..]);
    if c3 != OK { return corr_fault(c3), c3; }
    let n = xs.len;
    if n == 2 { return corr_good(r, 1.0, n - 2), OK; }
    let p, cb = corr_beta_p(r, n, alt);
    if cb != OK { return corr_fault(cb), cb; }
    return corr_good(r, p, n - 2), OK;
}

// The point biserial correlation of a list with two values and a list of
// measurements. The p-value comes from the beta law of the correlation.
fn pointbiserialr(xs: const f64[:], ys: const f64[:], alt: i64 = ALT_TWO_SIDED) -> CorrelationResult, i64 {
    let ca = check_alt(alt);
    let c = check_pair(xs, ys, 2);
    if c != OK { return corr_fault(c), c; }
    if ca != OK { return corr_fault(ca), ca; }
    // The first list holds two values, so one of them is 0 and the other
    // one is not 0.
    var low: f64 = 0.0;
    var high: f64 = 0.0;
    var seen: i64 = 0;
    for x in xs {
        if seen == 0 { low = x; high = x; seen = 1; }
        else if x < low { low = x; }
        else if x > high { high = x; }
    }
    if high == low { return corr_fault(ERR_UNDEFINED), ERR_UNDEFINED; }
    let r, c2 = correlation(xs, ys);
    if c2 != OK { return corr_fault(c2), c2; }
    let n = xs.len;
    if n == 2 { return corr_good(r, 1.0, n - 2), OK; }
    let p, cb = corr_beta_p(r, n, alt);
    if cb != OK { return corr_fault(cb), cb; }
    return corr_good(r, p, n - 2), OK;
}

// The xi correlation of Chatterjee and the p-value of the test of
// independence. The statistic uses the ranks of the second list in the
// order of the first list.
fn chatterjeexi(xs: const f64[:], ys: const f64[:], y_continuous: bool = false) -> CorrelationResult, i64 {
    let c = check_pair(xs, ys, 2);
    if c != OK { return corr_fault(c), c; }
    let n = xs.len;
    // The pairs in the order of the first list.
    var idx: Vi[>..] = [];
    for k in 0..n { idx.push(Vi { v: xs[k], i: k }); }
    stable_sort(idx) { a, b => a.v < b.v };
    // The rank of the second value of a pair is the number of the pairs
    // whose second value is not larger. The rank from below is the
    // number of the pairs whose second value is not smaller.
    var r: f64[>..] = [];
    var l: f64[>..] = [];
    for k in 0..n { r.push(0.0); l.push(0.0); }
    for a in 0..n {
        let ya = ys[idx[a].i];
        var up: f64 = 0.0;
        var down: f64 = 0.0;
        for b in 0..n {
            let yb = ys[idx[b].i];
            if yb <= ya { up += 1.0; }
            if yb >= ya { down += 1.0; }
        }
        r[a] = up;
        l[a] = down;
    }
    var num: f64 = 0.0;
    for k in 1..n { num += abs(r[k] - r[k - 1]); }
    var stat: f64 = 0.0;
    if y_continuous {
        stat = 1.0 - 3.0 * num / (n * n - 1);
    } else {
        var den: f64 = 0.0;
        for k in 0..n { den += (n - l[k]) * l[k]; }
        den *= 2.0;
        if den == 0.0 { return corr_fault(ERR_UNDEFINED), ERR_UNDEFINED; }
        stat = 1.0 - n * num / den;
    }
    let std = xi_std(r, l, y_continuous, n);
    if std == 0.0 { return corr_fault(ERR_UNDEFINED), ERR_UNDEFINED; }
    let p = norm_sf(stat / std);
    if is_nan(p) { return corr_fault(ERR_UNDEFINED), ERR_UNDEFINED; }
    return corr_good(stat, p, nan()), OK;
}

// The asymptotic deviation of the xi correlation under the law of
// independence.
fn xi_std(r: const f64[:], l: const f64[:], y_continuous: bool, n: i64) -> f64 {
    let nf: f64 = n;
    if y_continuous { return sqrt(0.4) / sqrt(nf); }
    // The ranks in increasing order, and their partial sums.
    var u: f64[>..] = [];
    for k in 0..n { u.push(r[k]); }
    stable_sort(u) { a, b => a < b };
    var v: f64[>..] = [];
    var sum: f64 = 0.0;
    for k in 0..n { sum += u[k]; v.push(sum); }
    var an: f64 = 0.0;
    var bn: f64 = 0.0;
    var cn: f64 = 0.0;
    var dn: f64 = 0.0;
    for k in 0..n {
        let i = (k + 1);
        let c = 2.0 * nf - 2.0 * i + 1.0;
        an += c * u[k] * u[k];
        bn += (v[k] + (nf - i) * u[k]) * (v[k] + (nf - i) * u[k]);
        cn += c * u[k];
        dn += l[k] * (nf - l[k]);
    }
    an /= nf * nf * nf * nf;
    bn /= nf * nf * nf * nf * nf;
    cn /= nf * nf * nf;
    dn /= nf * nf * nf;
    if dn == 0.0 { return 0.0; }
    let tau2 = (an - 2.0 * bn + cn * cn) / (dn * dn);
    if tau2 < 0.0 { return 0.0; }
    return sqrt(tau2) / sqrt(nf);
}
