// Tests against the reference vectors of the association calls. The
// expected values come from rational arithmetic. The tolerance of a call
// is in units of the relative machine epsilon of the expected value.
namespace suite;

import std;
import statistics from "../..";

// Counts a failure of a call that must give a fault.
fn expect_fault(bad: i64&, code: i64, name: u8[:], k: i64) {
    if code != 0 { return; }
    bad += 1;
    if bad <= 12 {
        print("FAIL ", name, " sample ", k, " gave no fault\n");
    }
}

// The covariance and the correlation of every pair.
fn run_assoc_pairs() -> i64 {
    var bad: i64 = 0;
    let count = A_OFFS.len - 1;
    for k in 0..count {
        let n = A_NCOLS[k];
        let start = A_OFFS[k];
        let xs = A_FLAT[start..start + n];
        let ys = A_FLAT[start + n..A_OFFS[k + 1]];
        let fit = A_FIT[k] == 1;
        let constant = A_CONST[k] == 1;
        if fit {
            let cov, c1 = statistics::covariance(xs, ys);
            expect(bad, within_unit(cov, E_COV[k], 4.0), "covariance", k, cov, E_COV[k], c1);
            let pcov, c2 = statistics::pcovariance(xs, ys);
            expect(bad, within_unit(pcov, E_PCOV[k], 4.0), "pcovariance", k, pcov, E_PCOV[k], c2);
        } else {
            let _cov, c1 = statistics::covariance(xs, ys);
            expect_fault(bad, c1, "covariance of a large sample", k);
            let _pcov, c2 = statistics::pcovariance(xs, ys);
            expect_fault(bad, c2, "pcovariance of a large sample", k);
        }
        if constant {
            let _r, c3 = statistics::correlation(xs, ys);
            expect_fault(bad, c3, "correlation of a constant column", k);
        } else if fit {
            let r, c3 = statistics::correlation(xs, ys);
            expect(bad, within_unit(r, E_CORR[k], 4.0), "correlation", k, r, E_CORR[k], c3);
        }
    }
    return bad;
}

// The weighted covariance and the weighted correlation of every pair. The
// weights are whole numbers, so the repeated sample is the second source.
fn run_assoc_weighted() -> i64 {
    var bad: i64 = 0;
    let count = A_OFFS.len - 1;
    for k in 0..count {
        let n = A_NCOLS[k];
        let start = A_OFFS[k];
        let xs = A_FLAT[start..start + n];
        let ys = A_FLAT[start + n..A_OFFS[k + 1]];
        let ws = A_WFLAT[A_WOFFS[k]..A_WOFFS[k + 1]];
        if A_WFIT[k] != 1 { continue; }
        let pop, c1 = statistics::weighted_covariance(xs, ys, ws, statistics::WEIGHT_POPULATION);
        expect(bad, within_unit(pop, E_WCOV_POP[k], 4.0), "weighted covariance population", k, pop, E_WCOV_POP[k], c1);
        let freq, c2 = statistics::weighted_covariance(xs, ys, ws, statistics::WEIGHT_FREQUENCY);
        expect(bad, within_unit(freq, E_WCOV_FREQ[k], 4.0), "weighted covariance frequency", k, freq, E_WCOV_FREQ[k], c2);
        let rel, c3 = statistics::weighted_covariance(xs, ys, ws, statistics::WEIGHT_RELIABILITY);
        expect(bad, within_unit(rel, E_WCOV_REL[k], 4.0), "weighted covariance reliability", k, rel, E_WCOV_REL[k], c3);
        let corr, c4 = statistics::weighted_correlation(xs, ys, ws);
        expect(bad, within_unit(corr, E_WCOV_CORR[k], 4.0), "weighted correlation", k, corr, E_WCOV_CORR[k], c4);
    }
    return bad;
}

// The ranks of every sample and the two rank correlations.
fn run_assoc_ranks() -> i64 {
    var bad: i64 = 0;
    let count = R_OFFS.len - 1;
    for k in 0..count {
        let start = R_OFFS[k];
        let end = R_OFFS[k + 1];
        let xs = R_FLAT[start..end];
        let ys = R_FLAT_Y[start..end];
        var avg: f64[>..] = [];
        var mn: f64[>..] = [];
        var mx: f64[>..] = [];
        var dense: f64[>..] = [];
        var ord: f64[>..] = [];
        let n1, c1 = statistics::rankdata(avg, xs, statistics::TIE_AVERAGE);
        expect(bad, c1 == 0, "rankdata average", k, 0.0, 0.0, c1);
        let n2, c2 = statistics::rankdata(mn, xs, statistics::TIE_MIN);
        expect(bad, c2 == 0, "rankdata min", k, 0.0, 0.0, c2);
        let n3, c3 = statistics::rankdata(mx, xs, statistics::TIE_MAX);
        expect(bad, c3 == 0, "rankdata max", k, 0.0, 0.0, c3);
        let n4, c4 = statistics::rankdata(dense, xs, statistics::TIE_DENSE);
        expect(bad, c4 == 0, "rankdata dense", k, 0.0, 0.0, c4);
        let n5, c5 = statistics::rankdata(ord, xs, statistics::TIE_ORDINAL);
        expect(bad, c5 == 0, "rankdata ordinal", k, 0.0, 0.0, c5);
        if n1 == 0 || n2 == 0 || n3 == 0 || n4 == 0 || n5 == 0 { continue; }
        let same = avg.len == xs.len && mn.len == xs.len && mx.len == xs.len && dense.len == xs.len && ord.len == xs.len;
        let got_len: f64 = avg.len;
        let want_len: f64 = xs.len;
        expect(bad, same, "rankdata count", k, got_len, want_len, 0);
        if !same { continue; }
        var ok = true;
        for i in 0..xs.len {
            if avg[i] != E_RAVG[start + i] { ok = false; }
            if mn[i] != E_RMIN[start + i] { ok = false; }
            if mx[i] != E_RMAX[start + i] { ok = false; }
            if dense[i] != E_RDENSE[start + i] { ok = false; }
            if ord[i] != E_RORD[start + i] { ok = false; }
        }
        expect(bad, ok, "rankdata values", k, avg[0], E_RAVG[start], 0);
        // A flag of 0 marks a pair whose rank correlation has no value:
        // a constant column, or an exact correlation of zero. The fault
        // of the constant column has its own test in the accept suite.
        if R_FIT[k] != 1 { continue; }
        let rho, c6 = statistics::spearman(xs, ys);
        expect(bad, within_unit(rho, E_SPEARMAN[k], 4.0), "spearman", k, rho, E_SPEARMAN[k], c6);
        let tb, c7 = statistics::kendall(xs, ys, statistics::KENDALL_B);
        expect(bad, within_unit(tb, E_KENDALL_B[k], 4.0), "kendall b", k, tb, E_KENDALL_B[k], c7);
        let tc, c8 = statistics::kendall(xs, ys, statistics::KENDALL_C);
        expect(bad, within_unit(tc, E_KENDALL_C[k], 4.0), "kendall c", k, tc, E_KENDALL_C[k], c8);
    }
    return bad;
}

// The line of every pair and the Theil-Sen estimate.
fn run_assoc_line() -> i64 {
    var bad: i64 = 0;
    let count = L_OFFS.len - 1;
    for k in 0..count {
        let start = L_OFFS[k];
        let end = L_OFFS[k + 1];
        let n = (end - start) / 2;
        let xs = L_FLAT[start..start + n];
        let ys = L_FLAT[start + n..end];
        if L_FIT[k] != 1 {
            let _l, c1 = statistics::linear_regression(xs, ys, false);
            expect_fault(bad, c1, "line of a constant column", k);
        } else {
            let l, c1 = statistics::linear_regression(xs, ys, false);
            expect(bad, within_unit(l.slope, E_SLOPE[k], 4.0), "slope", k, l.slope, E_SLOPE[k], c1);
            expect(bad, within_unit(l.intercept, E_INTERCEPT[k], 4.0), "intercept", k, l.intercept, E_INTERCEPT[k], c1);
        }
        if L_ERR_FIT[k] == 1 {
            let se, ie, c2 = statistics::regression_errors(xs, ys);
            expect(bad, within_unit(se, E_SLOPE_ERR[k], 4.0), "slope error", k, se, E_SLOPE_ERR[k], c2);
            expect(bad, within_unit(ie, E_INTERCEPT_ERR[k], 4.0), "intercept error", k, ie, E_INTERCEPT_ERR[k], c2);
        }
        if L_THEIL[k] == 1 {
            // The intercept is the median of y minus the slope times the
            // median of x. The subtraction cancels, and the rounding of
            // the medians differs between the compilers, so the
            // tolerance is wide.
            let ts, ti, c3 = statistics::theil_sen(xs, ys);
            expect(bad, within_unit(ts, E_THEIL_SLOPE[k], 64.0), "theil slope", k, ts, E_THEIL_SLOPE[k], c3);
            expect(bad, within_unit(ti, E_THEIL_INTERCEPT[k], 64.0), "theil intercept", k, ti, E_THEIL_INTERCEPT[k], c3);
        }
    }
    return bad;
}

// The covariance matrix and the correlation matrix of every sample.
fn run_assoc_matrix() -> i64 {
    var bad: i64 = 0;
    let count = M_OFFS.len - 1;
    // The cells of every matrix follow the cells of the matrix before it.
    var cell_at: i64 = 0;
    var cell_start: i64 = 0;
    for k in 0..count {
        let nc = M_COLS[k];
        let rows = M_ROWS[k];
        cell_start = cell_at;
        cell_at += nc * nc;
        let data = M_FLAT[M_OFFS[k]..M_OFFS[k + 1]];
        var cells: f64[>..] = [];
        if M_CONST[k] == 1 {
            let n1, c1 = statistics::correlation_matrix(cells, data, nc);
            expect_fault(bad, c1, "correlation matrix of a constant column", k);
            if cells.len != 0 {
                bad += 1;
                if bad <= 12 { print("FAIL correlation matrix wrote a cell of a fault sample ", k, "\n"); }
            }
            continue;
        }
        let n1, c1 = statistics::covariance_matrix(cells, data, nc, 1);
        let got_cells: f64 = n1;
        let want_cells: f64 = nc * nc;
        expect(bad, n1 == nc * nc, "covariance matrix count", k, got_cells, want_cells, c1);
        if c1 == 0 {
            var ok = true;
            for i in 0..nc {
                for j in 0..nc {
                    if !within_unit(cells[i * nc + j], E_MAT_COV[cell_start + i * nc + j], 4.0) { ok = false; }
                }
            }
            expect(bad, ok, "covariance matrix", k, cells[0], E_MAT_COV[cell_start], 0);
        }
        var ccells: f64[>..] = [];
        let n2, c2 = statistics::correlation_matrix(ccells, data, nc);
        let got_ccells: f64 = n2;
        let want_ccells: f64 = nc * nc;
        expect(bad, n2 == nc * nc, "correlation matrix count", k, got_ccells, want_ccells, c2);
        if c2 == 0 {
            var ok = true;
            for i in 0..nc {
                for j in 0..nc {
                    if !within_unit(ccells[i * nc + j], E_MAT_CORR[cell_start + i * nc + j], 4.0) { ok = false; }
                }
            }
            expect(bad, ok, "correlation matrix", k, ccells[0], E_MAT_CORR[cell_start], 0);
        }
    }
    return bad;
}

fn run_assoc() -> i64 {
    var bad: i64 = 0;
    bad += run_assoc_pairs();
    bad += run_assoc_weighted();
    bad += run_assoc_ranks();
    bad += run_assoc_line();
    bad += run_assoc_matrix();
    return bad;
}
