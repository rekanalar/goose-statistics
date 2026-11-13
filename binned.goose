// The binned statistics: the statistic of the values of each bin of one
// or many dimensions.
//
// The bin numbers follow the rule of NumPy. In each dimension the slot
// 0 holds the values below the first edge and the slot nbin holds the
// values above the last edge, where nbin is the count of the edges plus
// one. The values that round to the last edge stay in the last bin.
// The bin number of a point is the row-major index over these slots.
namespace statistics;

import std;
import math;

// The digits of the rounding of the right edge rule. The rule rounds
// the values and the edge to this many digits.
fn bd_digits(gap: f64) -> i64 {
    if gap <= 0.0 { return -1; }
    var d = (-log10(gap)) as i64 + 6;
    if d < 0 { return 0; }
    // 10^300 is the largest power that stays finite for any value.
    if d > 300 { return 300; }
    return d;
}

// The rounding of a value to a number of digits.
fn bd_round(x: f64, digits: i64) -> f64 {
    var scale: f64 = 1.0;
    var i: i64 = 0;
    while i < digits { scale *= 10.0; i += 1; }
    return round(x * scale) / scale;
}

// The edges of one dimension. The list holds the edges. The edges come
// from the range, or from the sample when the range is not given. The
// edges are appended to the list, one dimension after the other. A
// range with equal ends widens by 0.5 on each side. The edges must
// ascend.
fn bd_edges(edges: f64[>..]&, xs: const f64[:], nbins: i64, low: f64, high: f64, use_low: bool, use_high: bool) -> i64, i64 {
    let c = check_sample(xs, 1);
    if c != 0 { return 0, c; }
    let bc = check_bins(nbins);
    if bc != 0 { return 0, bc; }
    var lo = low;
    var hi = high;
    if !use_low && !use_high {
        let a, c2 = minimum(xs);
        if c2 != 0 { return 0, c2; }
        let b, c3 = maximum(xs);
        if c3 != 0 { return 0, c3; }
        lo = a;
        hi = b;
    } else {
        if use_low {
            let lc = check_scalar(low);
            if lc != OK { return 0, lc; }
        }
        if use_high {
            let hc = check_scalar(high);
            if hc != OK { return 0, hc; }
        }
        if use_low && use_high {
            let lc = check_limits(low, high);
            if lc != 0 { return 0, lc; }
        }
    }
    if lo == hi {
        lo = lo - 0.5;
        hi = hi + 0.5;
    }
    let n = nbins + 1;
    let start = edges.len;
    var i: i64 = 0;
    while i < n {
        edges.push(lo + (hi - lo) * (i) / (nbins));
        i += 1;
    }
    i = start + 1;
    while i < edges.len {
        if edges[i] <= edges[i - 1] { return 0, ERR_DOMAIN; }
        i += 1;
    }
    return nbins, OK;
}

// The edges of one dimension from the caller. The list holds the edges
// and they must ascend.
fn bd_edges_given(edges: f64[>..]&, given: const f64[:]) -> i64, i64 {
    edges.resize(0);
    if given.len < 2 { return 0, ERR_DOMAIN; }
    var i: i64 = 1;
    while i < given.len {
        if given[i] <= given[i - 1] { return 0, ERR_DOMAIN; }
        i += 1;
    }
    for v in given { edges.push(v); }
    return given.len - 1, OK;
}

// The slot of a value in one dimension. The slot 0 is below the first
// edge and the slot nbin is above the last edge, where nbin is the
// count of the edges. The values that round to the last edge stay in
// the last bin, as in SciPy.
fn bd_slot_of(x: f64, edges: const f64[:]) -> i64 {
    let n = edges.len;
    // The count of the edges at or below x, by binary search.
    var lo_i: i64 = 0;
    var hi_i: i64 = n;
    while lo_i < hi_i {
        let mid = (lo_i + hi_i) / 2;
        if x >= edges[mid] { lo_i = mid + 1; } else { hi_i = mid; }
    }
    let b = lo_i;
    if b == 0 { return 0; }
    if b == n {
        let digits = bd_digits(bd_min_gap(edges));
        if digits >= 0 {
            if bd_round(x, digits) == bd_round(edges[n - 1], digits) { return n - 1; }
        }
        return n;
    }
    return b;
}

// The smallest gap between two edges.
fn bd_min_gap(edges: const f64[:]) -> f64 {
    var gap: f64 = 0.0;
    var i: i64 = 1;
    while i < edges.len {
        let d = edges[i] - edges[i - 1];
        if d > 0.0 && (gap == 0.0 || d < gap) { gap = d; }
        i += 1;
    }
    return gap;
}

// The statistic of the values of one bin.
fn bd_statistic(values: const f64[:], stat: i64) -> f64, i64 {
    if stat == STAT_MEAN {
        let m, c = exact_mean(values);
        return m, c;
    }
    if stat == STAT_STD {
        // The deviation with the divisor n, as in SciPy.
        if values.len == 0 { return nan(), ERR_UNDEFINED; }
        let v, c = variance_core(values, 0);
        if c != 0 { return nan(), c; }
        return sqrt(v), OK;
    }
    if stat == STAT_MEDIAN {
        let m, c = median(values);
        return m, c;
    }
    return nan(), ERR_DOMAIN;
}

// The slots of a point in each dimension. The list shape holds the
// count of the edges of each dimension and the list all_edges the edges
// of each dimension, one dimension after the other.
fn bd_slots_of(point: const f64[:], all_edges: const f64[:], shape: const i64[:]) -> i64[>..] {
    var slots: i64[>..] = [];
    var start: i64 = 0;
    var d: i64 = 0;
    while d < shape.len {
        let e = all_edges[start..(start + shape[d])];
        slots.push(bd_slot_of(point[d], e));
        start += shape[d];
        d += 1;
    }
    return slots;
}

// The row-major index of the slots over the count of the slots of each
// dimension. The count of the slots of a dimension is the count of its
// edges plus one.
fn bd_index_of(slots: const i64[:], shape: const i64[:]) -> i64 {
    var index: i64 = 0;
    var d: i64 = 0;
    while d < shape.len {
        index = index * (shape[d] + 1) + slots[d];
        d += 1;
    }
    return index;
}

// The cell of a point in the table. The cell is the row-major index
// over the bins of the table, where the bins exclude the outer slots.
// A point outside the edges of any dimension has no cell, which the
// result marks with -1.
fn bd_cell_of(slots: const i64[:], shape: const i64[:]) -> i64 {
    var index: i64 = 0;
    var d: i64 = 0;
    while d < shape.len {
        let s = slots[d];
        if s == 0 || s == shape[d] { return -1; }
        index = index * (shape[d] - 1) + (s - 1);
        d += 1;
    }
    return index;
}

// The number of cells of a table.
fn bd_cell_count(shape: const i64[:]) -> i64 {
    var total: i64 = 1;
    var d: i64 = 0;
    while d < shape.len {
        let side = shape[d] - 1;
        // A product above MAX_CELLS is reported as MAX_CELLS + 1, so the
        // product never wraps.
        if side > 0 && total > MAX_CELLS / side { return MAX_CELLS + 1; }
        total *= side;
        d += 1;
    }
    return total;
}

// The statistic of each cell of a table. The list cells holds the cell
// of each point, and the list values the value of each point. The list
// stats holds the statistic of each cell and the list empty the flag of
// each cell.
fn bd_fill_cells(stats: f64[>..]&, empty: f64[>..]&, cells: const f64[:], values: const f64[:], stat: i64, total: i64) -> i64 {
    var counts: i64[>..] = [];
    var sums: f64[>..] = [];
    var mins: f64[>..] = [];
    var maxs: f64[>..] = [];
    var d: i64 = 0;
    while d < total {
        counts.push(0);
        sums.push(0.0);
        mins.push(inf());
        maxs.push(-inf());
        d += 1;
    }
    for i in 0..values.len {
        let cell = cells[i] as i64;
        if cell < 0 || cell >= total { continue; }
        counts[cell] += 1;
        sums[cell] += values[i];
        if values[i] < mins[cell] { mins[cell] = values[i]; }
        if values[i] > maxs[cell] { maxs[cell] = values[i]; }
    }
    d = 0;
    while d < total {
        var v: f64 = 0.0;
        if counts[d] == 0 {
            v = 0.0;
            empty.push(1.0);
        } else {
            if stat == STAT_COUNT { v = counts[d]; }
            else if stat == STAT_SUM { v = sums[d]; }
            else if stat == STAT_MIN { v = mins[d]; }
            else if stat == STAT_MAX { v = maxs[d]; }
            else { v = nan(); }
            empty.push(0.0);
        }
        stats.push(v);
        d += 1;
    }
    return OK;
}

// The values of the points grouped by cell. The list sorted holds the
// values of the points that have a cell, cell after cell. The list
// starts holds one more entry than the cells: the values of cell d lie
// in sorted[starts[d]..starts[d + 1]]. One counting pass does the work,
// so the cost grows with the cells plus the points.
fn bd_group(sorted: f64[>..]&, starts: i64[>..]&, cells: const f64[:], values: const f64[:], total: i64) {
    sorted.resize(0);
    starts.resize(0);
    var d: i64 = 0;
    while d <= total { starts.push(0); d += 1; }
    var kept: i64 = 0;
    for i in 0..values.len {
        let cell = cells[i] as i64;
        if cell < 0 || cell >= total { continue; }
        starts[cell + 1] += 1;
        kept += 1;
    }
    d = 0;
    while d < total { starts[d + 1] += starts[d]; d += 1; }
    var fill: i64[>..] = [];
    d = 0;
    while d < total { fill.push(starts[d]); d += 1; }
    var k: i64 = 0;
    while k < kept { sorted.push(0.0); k += 1; }
    for i in 0..values.len {
        let cell = cells[i] as i64;
        if cell < 0 || cell >= total { continue; }
        sorted[fill[cell]] = values[i];
        fill[cell] += 1;
    }
}

// The statistic of each cell of a table, for the statistics that need
// the values of the cell. The empty cells keep the value 0.
fn bd_refine_cells(stats: f64[>..]&, sorted: const f64[:], starts: const i64[:], stat: i64) -> i64 {
    if stat == STAT_COUNT || stat == STAT_SUM || stat == STAT_MIN || stat == STAT_MAX { return OK; }
    var d: i64 = 0;
    while d < stats.len {
        if starts[d + 1] > starts[d] {
            let v, c = bd_statistic(sorted[starts[d]..starts[d + 1]], stat);
            if c != 0 { return c; }
            stats[d] = v;
        }
        d += 1;
    }
    return OK;
}

// The statistic of each cell of a table for a block of the caller. The
// block takes the values of a cell and gives one number. A block result
// of NaN or infinity gives ERR_UNDEFINED.
fn bd_apply_block<F>(stats: f64[>..]&, sorted: const f64[:], starts: const i64[:]) -> i64 {
    var d: i64 = 0;
    while d < stats.len {
        if starts[d + 1] > starts[d] {
            let v = F(sorted[starts[d]..starts[d + 1]]);
            if is_bad(v) { return ERR_UNDEFINED; }
            stats[d] = v;
        }
        d += 1;
    }
    return OK;
}

// The slots and the cell of each point. The list shape holds the count
// of the edges of each dimension.
fn bd_points(cells: f64[>..]&, bins: f64[>..]&, rows: const f64[:], dims: i64, all_edges: const f64[:], shape: const i64[:]) -> i64, i64 {
    cells.resize(0);
    bins.resize(0);
    let nrows = rows.len / dims;
    var r: i64 = 0;
    while r < nrows {
        let point = rows[r * dims..(r + 1) * dims];
        let slots = bd_slots_of(point, all_edges, shape);
        bins.push(bd_index_of(slots, shape));
        cells.push(bd_cell_of(slots, shape));
        r += 1;
    }
    return cells.len, OK;
}

// The table of one dimension. The list shape holds the count of the
// edges of the dimension.
fn bd_table(stats: f64[>..]&, edges_out: f64[>..]&, shape_out: f64[>..]&, binnumber: f64[>..]&, empty: f64[>..]&, xs: const f64[:], values: const f64[:], stat: i64, shape: const i64[:], all_edges: const f64[:], dims: i64) -> BinnedTable, i64 {
    var t = BinnedTable { .. };
    let total = bd_cell_count(shape);
    if total > MAX_CELLS { return t, ERR_TOO_LONG; }
    var cells: f64[>..] = [];
    var bins: f64[>..] = [];
    let _np, r = bd_points(cells, bins, xs, dims, all_edges, shape);
    if r != 0 { return t, r; }
    stats.resize(0);
    empty.resize(0);
    edges_out.resize(0);
    shape_out.resize(0);
    binnumber.resize(0);
    let c2 = bd_fill_cells(stats, empty, cells, values, stat, total);
    if c2 != 0 { return t, c2; }
    var sorted: f64[>..] = [];
    var starts: i64[>..] = [];
    bd_group(sorted, starts, cells, values, total);
    let c3 = bd_refine_cells(stats, sorted[..], starts[..], stat);
    if c3 != 0 { return t, c3; }
    for v in all_edges { edges_out.push(v); }
    for v in shape { shape_out.push(v); }
    for v in bins { binnumber.push(v); }
    t.dims = dims;
    t.cells = total;
    t.points = values.len;
    return t, OK;
}

// Signature: binned_statistic(out, xs, values, stat, nbins) -> BinnedTable, code
// The statistic of the values of each bin of one dimension. The list
// xs holds the position of each point and the list values the value of
// each point. The statistic is one of the STAT_ constants. The output
// lists of the struct hold the statistic of each bin, the edges, the
// count of the edges, the bin of each point and the flag of the empty
// bins.
fn binned_statistic(stats: f64[>..]&, edges_out: f64[>..]&, shape_out: f64[>..]&, binnumber: f64[>..]&, empty: f64[>..]&, xs: const f64[:], values: const f64[:], stat: i64, nbins: i64) -> BinnedTable, i64 {
    var t = BinnedTable { .. };
    if stat < STAT_MEAN || stat > STAT_MAX { return t, ERR_DOMAIN; }
    let c = check_pair(xs, values, 0);
    if c != 0 { return t, c; }
    var all_edges: f64[>..] = [];
    let nb, c2 = bd_edges(all_edges, xs, nbins, 0.0, 0.0, false, false);
    if c2 != 0 { return t, c2; }
    var shape: i64[>..] = [];
    shape.push(nb + 1);
    return bd_table(stats, edges_out, shape_out, binnumber, empty, xs, values, stat, shape, all_edges, 1);
}

// Signature: binned_statistic_edges(out, xs, values, stat, edges) -> BinnedTable, code
// The statistic of the values of each bin of one dimension, with the
// edges given by the caller. The edges must ascend.
fn binned_statistic_edges(stats: f64[>..]&, edges_out: f64[>..]&, shape_out: f64[>..]&, binnumber: f64[>..]&, empty: f64[>..]&, xs: const f64[:], values: const f64[:], stat: i64, given: const f64[:]) -> BinnedTable, i64 {
    var t = BinnedTable { .. };
    if stat < STAT_MEAN || stat > STAT_MAX { return t, ERR_DOMAIN; }
    let c = check_pair(xs, values, 0);
    if c != 0 { return t, c; }
    var all_edges: f64[>..] = [];
    let nb, c2 = bd_edges_given(all_edges, given);
    if c2 != 0 { return t, c2; }
    var shape: i64[>..] = [];
    shape.push(nb + 1);
    return bd_table(stats, edges_out, shape_out, binnumber, empty, xs, values, stat, shape, all_edges, 1);
}

// Signature: binned_statistic_with<F>(stats, edges, shape, binnumber, empty, xs, values, nbins) -> BinnedTable, code
// The statistic of the values of each bin of one dimension, where the
// caller gives the statistic. The block takes the values of a bin and
// gives one number. A block result of NaN or infinity gives
// ERR_UNDEFINED. The edges come from the range of the sample.
fn binned_statistic_with<F>(stats: f64[>..]&, edges_out: f64[>..]&, shape_out: f64[>..]&, binnumber: f64[>..]&, empty: f64[>..]&, xs: const f64[:], values: const f64[:], nbins: i64) -> BinnedTable, i64 {
    var t = BinnedTable { .. };
    let c = check_pair(xs, values, 0);
    if c != 0 { return t, c; }
    var all_edges: f64[>..] = [];
    let nb, c2 = bd_edges(all_edges, xs, nbins, 0.0, 0.0, false, false);
    if c2 != 0 { return t, c2; }
    var shape: i64[>..] = [];
    shape.push(nb + 1);
    return bd_table_with(stats, edges_out, shape_out, binnumber, empty, xs, values, shape, all_edges, 1) { v => F(v) };
}

// The table of one dimension where the caller gives the statistic.
fn bd_table_with<F>(stats: f64[>..]&, edges_out: f64[>..]&, shape_out: f64[>..]&, binnumber: f64[>..]&, empty: f64[>..]&, xs: const f64[:], values: const f64[:], shape: const i64[:], all_edges: const f64[:], dims: i64) -> BinnedTable, i64 {
    var t = BinnedTable { .. };
    let total = bd_cell_count(shape);
    if total > MAX_CELLS { return t, ERR_TOO_LONG; }
    var cells: f64[>..] = [];
    var bins: f64[>..] = [];
    let _np, r = bd_points(cells, bins, xs, dims, all_edges, shape);
    if r != 0 { return t, r; }
    stats.resize(0);
    empty.resize(0);
    edges_out.resize(0);
    shape_out.resize(0);
    binnumber.resize(0);
    let c2 = bd_fill_cells(stats, empty, cells, values, STAT_COUNT, total);
    if c2 != 0 { return t, c2; }
    var sorted: f64[>..] = [];
    var starts: i64[>..] = [];
    bd_group(sorted, starts, cells, values, total);
    let c3 = bd_apply_block(stats, sorted[..], starts[..]) { v => F(v) };
    if c3 != 0 { return t, c3; }
    for v in all_edges { edges_out.push(v); }
    for v in shape { shape_out.push(v); }
    for v in bins { binnumber.push(v); }
    t.dims = dims;
    t.cells = total;
    t.points = values.len;
    return t, OK;
}

// Signature: binned_statistic_2d(out, xs, ys, values, stat, nx, ny) -> BinnedTable, code
// The statistic of the values of each bin of two dimensions. The list
// xs and the list ys hold the position of each point and the list
// values the value of each point.
fn binned_statistic_2d(stats: f64[>..]&, edges_out: f64[>..]&, shape_out: f64[>..]&, binnumber: f64[>..]&, empty: f64[>..]&, xs: const f64[:], ys: const f64[:], values: const f64[:], stat: i64, nx: i64, ny: i64) -> BinnedTable, i64 {
    var t = BinnedTable { .. };
    if stat < STAT_MEAN || stat > STAT_MAX { return t, ERR_DOMAIN; }
    let c = check_pair(xs, values, 0);
    if c != 0 { return t, c; }
    let c2 = check_pair(ys, values, 0);
    if c2 != 0 { return t, c2; }
    var all_edges: f64[>..] = [];
    let a, c3 = bd_edges(all_edges, xs, nx, 0.0, 0.0, false, false);
    if c3 != 0 { return t, c3; }
    let b, c4 = bd_edges(all_edges, ys, ny, 0.0, 0.0, false, false);
    if c4 != 0 { return t, c4; }
    var shape: i64[>..] = [];
    shape.push(a + 1);
    shape.push(b + 1);
    // The rows of the table hold the two coordinates of each point.
    if xs.len != ys.len { return t, ERR_LENGTH; }
    var rows: f64[>..] = [];
    var i: i64 = 0;
    while i < xs.len {
        rows.push(xs[i]);
        rows.push(ys[i]);
        i += 1;
    }
    return bd_table(stats, edges_out, shape_out, binnumber, empty, rows, values, stat, shape, all_edges, 2);
}

// Signature: binned_statistic_2d_with<F>(stats, edges, shape, binnumber, empty, xs, ys, values, nx, ny) -> BinnedTable, code
// The statistic of the values of each bin of two dimensions, where the
// caller gives the statistic. The block takes the values of a bin and
// gives one number.
fn binned_statistic_2d_with<F>(stats: f64[>..]&, edges_out: f64[>..]&, shape_out: f64[>..]&, binnumber: f64[>..]&, empty: f64[>..]&, xs: const f64[:], ys: const f64[:], values: const f64[:], nx: i64, ny: i64) -> BinnedTable, i64 {
    var t = BinnedTable { .. };
    let c = check_pair(xs, values, 0);
    if c != 0 { return t, c; }
    let c2 = check_pair(ys, values, 0);
    if c2 != 0 { return t, c2; }
    var all_edges: f64[>..] = [];
    let a, c3 = bd_edges(all_edges, xs, nx, 0.0, 0.0, false, false);
    if c3 != 0 { return t, c3; }
    let b, c4 = bd_edges(all_edges, ys, ny, 0.0, 0.0, false, false);
    if c4 != 0 { return t, c4; }
    var shape: i64[>..] = [];
    shape.push(a + 1);
    shape.push(b + 1);
    if xs.len != ys.len { return t, ERR_LENGTH; }
    var rows: f64[>..] = [];
    var i: i64 = 0;
    while i < xs.len {
        rows.push(xs[i]);
        rows.push(ys[i]);
        i += 1;
    }
    return bd_table_with(stats, edges_out, shape_out, binnumber, empty, rows, values, shape, all_edges, 2) { v => F(v) };
}

// Signature: binned_statistic_dd(out, sample, dims, values, stat, nbins) -> BinnedTable, code
// The statistic of the values of each bin of many dimensions. The list
// sample holds the rows of the points one after the other, and dims is
// the count of the dimensions. The bin count is the same in each
// dimension.
fn binned_statistic_dd(stats: f64[>..]&, edges_out: f64[>..]&, shape_out: f64[>..]&, binnumber: f64[>..]&, empty: f64[>..]&, sample: const f64[:], dims: i64, values: const f64[:], stat: i64, nbins: i64) -> BinnedTable, i64 {
    var t = BinnedTable { .. };
    if stat < STAT_MEAN || stat > STAT_MAX { return t, ERR_DOMAIN; }
    let c = check_sample(sample, 1);
    if c != 0 { return t, c; }
    let vc = check_sample(values, 1);
    if vc != 0 { return t, vc; }
    if dims < 1 { return t, ERR_DOMAIN; }
    if dims > MAX_DIMS { return t, ERR_TOO_LONG; }
    if sample.len % dims != 0 || sample.len / dims != values.len { return t, ERR_LENGTH; }
    let edges = bd_edges_dd(sample, dims, nbins);
    if edges.len == 0 { return t, ERR_DOMAIN; }
    var shape: i64[>..] = [];
    for d in 0..dims { shape.push(nbins + 1); }
    return bd_table(stats, edges_out, shape_out, binnumber, empty, sample, values, stat, shape, edges, dims);
}

// Signature: binned_statistic_dd_with<F>(stats, edges, shape, binnumber, empty, sample, dims, values, nbins) -> BinnedTable, code
// The statistic of the values of each bin of many dimensions, where the
// caller gives the statistic. The block takes the values of a bin and
// gives one number.
fn binned_statistic_dd_with<F>(stats: f64[>..]&, edges_out: f64[>..]&, shape_out: f64[>..]&, binnumber: f64[>..]&, empty: f64[>..]&, sample: const f64[:], dims: i64, values: const f64[:], nbins: i64) -> BinnedTable, i64 {
    var t = BinnedTable { .. };
    let c = check_sample(sample, 1);
    if c != 0 { return t, c; }
    let vc = check_sample(values, 1);
    if vc != 0 { return t, vc; }
    if dims < 1 { return t, ERR_DOMAIN; }
    if dims > MAX_DIMS { return t, ERR_TOO_LONG; }
    if sample.len % dims != 0 || sample.len / dims != values.len { return t, ERR_LENGTH; }
    let edges = bd_edges_dd(sample, dims, nbins);
    if edges.len == 0 { return t, ERR_DOMAIN; }
    var shape: i64[>..] = [];
    for d in 0..dims { shape.push(nbins + 1); }
    return bd_table_with(stats, edges_out, shape_out, binnumber, empty, sample, values, shape, edges, dims) { v => F(v) };
}

// The edges of each dimension of a table of many dimensions. The list
// holds the edges of each dimension, one dimension after the other. The
// edges of a dimension come from the range of its column.
fn bd_edges_dd(sample: const f64[:], dims: i64, nbins: i64) -> f64[>..] {
    var all_edges: f64[>..] = [];
    let rows = sample.len / dims;
    var d: i64 = 0;
    while d < dims {
        var column: f64[>..] = [];
        var r: i64 = 0;
        while r < rows {
            column.push(sample[r * dims + d]);
            r += 1;
        }
        let nb, c = bd_edges(all_edges, column[..], nbins, 0.0, 0.0, false, false);
        if c != OK { return []; }
        if nb < 1 { return []; }
        d += 1;
    }
    return all_edges;
}
