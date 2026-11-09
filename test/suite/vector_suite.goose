// Tests against the reference vectors. The tolerance of a call is in
// units of the relative machine epsilon of the expected value.
namespace suite;

import std;
import statistics from "../..";

const EPS: f64 = 2.220446049250313e-16;
const TINY: f64 = 5e-324;

// True when a equals e within ulps units.
fn within(a: f64, e: f64, ulps: f64) -> bool {
    if a == e { return true; }
    let tol = max(ulps * EPS * abs(e), ulps * TINY);
    return abs(a - e) <= tol;
}

// The same, but the scale is at least 1. The skewness and the kurtosis
// cancel to a value near zero, so a relative bound would be too strict.
fn within_unit(a: f64, e: f64, ulps: f64) -> bool {
    if a == e { return true; }
    return abs(a - e) <= ulps * EPS * max(abs(e), 1.0);
}

// Counts a failure and prints the first few.
fn expect(bad: i64&, ok: bool, name: u8[:], k: i64, got: f64, want: f64, code: i64) {
    if ok && code == 0 { return; }
    bad += 1;
    if bad <= 12 {
        print("FAIL ", name, " sample ", k, " got ", got, " want ", want, " code ", code, "\n");
    }
}

// A sorted copy of a sample, made by the test.
fn sorted_copy(xs: const f64[:]) -> f64[>..] {
    var out: f64[>..] = [];
    for x in xs { out.push(x); }
    sort(out) { a, b => a < b };
    return out;
}

// The weighted vectors. The whole number weights have the repeated sample
// as a second source.
fn run_weighted() -> i64 {
    var bad: i64 = 0;
    let count = W_OFFS.len - 1;
    for k in 0..count {
        let xs = W_FLAT[W_OFFS[k]..W_OFFS[k + 1]];
        let ws = W_WEIGHTS[W_OFFS[k]..W_OFFS[k + 1]];
        let m, c1 = statistics::weighted_mean(xs, ws);
        expect(bad, within(m, E_WMEAN[k], 4.0), "weighted_mean", k, m, E_WMEAN[k], c1);
        let p, c2 = statistics::weighted_variance(xs, ws, statistics::WEIGHT_POPULATION);
        expect(bad, within(p, E_WPOP[k], 4.0), "weighted variance population", k, p, E_WPOP[k], c2);
        let f, c3 = statistics::weighted_variance(xs, ws, statistics::WEIGHT_FREQUENCY);
        expect(bad, within(f, E_WFREQ[k], 4.0), "weighted variance frequency", k, f, E_WFREQ[k], c3);
        let r, c4 = statistics::weighted_variance(xs, ws, statistics::WEIGHT_RELIABILITY);
        expect(bad, within(r, E_WREL[k], 4.0), "weighted variance reliability", k, r, E_WREL[k], c4);
        let md, c5 = statistics::weighted_median(xs, ws);
        expect(bad, md == E_WMED[k], "weighted_median", k, md, E_WMED[k], c5);
    }
    return bad;
}

fn run_vectors() -> i64 {
    var bad: i64 = 0;
    let count = V_OFFS.len - 1;
    for k in 0..count {
        let xs = V_FLAT[V_OFFS[k]..V_OFFS[k + 1]];
        let c = V_CONST[k] == 1;
        let pos = V_POS[k] == 1;

        let s, c1 = statistics::sum(xs);
        expect(bad, within(s, E_SUM[k], 1.0), "sum", k, s, E_SUM[k], c1);
        let m, c2 = statistics::mean(xs);
        expect(bad, within(m, E_MEAN[k], 1.0), "mean", k, m, E_MEAN[k], c2);
        let fm, c2b = statistics::fmean(xs);
        expect(bad, within(fm, E_MEAN[k], 1.0), "fmean", k, fm, E_MEAN[k], c2b);

        let md, c3 = statistics::median(xs);
        expect(bad, within(md, E_MEDIAN[k], 1.0), "median", k, md, E_MEDIAN[k], c3);
        let ml, c4 = statistics::median_low(xs);
        expect(bad, ml == E_MLOW[k], "median_low", k, ml, E_MLOW[k], c4);
        let mh, c5 = statistics::median_high(xs);
        expect(bad, mh == E_MHIGH[k], "median_high", k, mh, E_MHIGH[k], c5);
        let mo, c6 = statistics::mode(xs);
        expect(bad, mo == E_MODE[k], "mode", k, mo, E_MODE[k], c6);

        let pv, c7 = statistics::pvariance(xs);
        expect(bad, within(pv, E_PVAR[k], 4.0), "pvariance", k, pv, E_PVAR[k], c7);
        let sv, c8 = statistics::variance(xs);
        expect(bad, within(sv, E_VAR[k], 4.0), "variance", k, sv, E_VAR[k], c8);
        let ps, c9 = statistics::pstdev(xs);
        expect(bad, within(ps, E_PSTDEV[k], 4.0), "pstdev", k, ps, E_PSTDEV[k], c9);
        let ss, c10 = statistics::stdev(xs);
        expect(bad, within(ss, E_STDEV[k], 4.0), "stdev", k, ss, E_STDEV[k], c10);
        let se, c11 = statistics::sem(xs);
        expect(bad, within(se, E_SEM[k], 4.0), "sem", k, se, E_SEM[k], c11);
        // A mean in the subnormal range has few bits, so a ratio with it
        // cannot keep a relative bound. The check covers normal means.
        if abs(E_MEAN[k]) > 1e-290 {
            let cv, c12 = statistics::variation(xs);
            expect(bad, within(cv, E_VARIATION[k], 4.0), "variation", k, cv, E_VARIATION[k], c12);
        }

        if pos {
            let g, c13 = statistics::geometric_mean(xs);
            expect(bad, within(g, E_GEO[k], 8.0), "geometric_mean", k, g, E_GEO[k], c13);
            let h, c14 = statistics::harmonic_mean(xs);
            expect(bad, within(h, E_HARM[k], 8.0), "harmonic_mean", k, h, E_HARM[k], c14);
        }

        // The quantile methods match NumPy within 2 units.
        for mth in 0..9 {
            for qi in 0..8 {
                let q = QS[qi];
                let v, cq = statistics::quantile(xs, q, mth);
                let e = E_QUANT[k * 72 + mth * 8 + qi];
                expect(bad, within(v, e, 2.0), "quantile", k, v, e, cq);
            }
        }
        let iq, c15 = statistics::iqr(xs, statistics::METHOD_LINEAR);
        expect(bad, within(iq, E_IQR[k], 4.0), "iqr", k, iq, E_IQR[k], c15);

        // The trimmed forms and the grouped median.
        let sorted = sorted_copy(xs);
        let lo = sorted[sorted.len / 4];
        let hi = sorted[3 * sorted.len / 4];
        let tm, c24 = statistics::tmean(xs, lo, hi);
        expect(bad, within(tm, E_TMEAN[k], 2.0), "tmean", k, tm, E_TMEAN[k], c24);
        let tr, c25 = statistics::trim_mean(xs, 0.1);
        expect(bad, within(tr, E_TRIM[k], 2.0), "trim_mean", k, tr, E_TRIM[k], c25);
        let gm, c26 = statistics::median_grouped(xs, 2.0);
        expect(bad, within(gm, E_GROUPED[k], 4.0), "median_grouped", k, gm, E_GROUPED[k], c26);

        // Shape. A constant sample has no skewness.
        if c {
            let _v, cs = statistics::skew(xs, true);
            if cs != statistics::ERR_UNDEFINED {
                bad += 1;
                print("FAIL skew of a constant sample gave code ", cs, "\n");
            }
            continue;
        }
        let sk, c16 = statistics::skew(xs, true);
        expect(bad, within_unit(sk, E_SKEW_B[k], 8.0), "skew", k, sk, E_SKEW_B[k], c16);
        let su, c17 = statistics::skew(xs, false);
        expect(bad, within_unit(su, E_SKEW_U[k], 8.0), "skew unbiased", k, su, E_SKEW_U[k], c17);
        let kf, c18 = statistics::kurtosis(xs, true, true);
        expect(bad, within_unit(kf, E_KURT_F_B[k], 8.0), "kurtosis fisher", k, kf, E_KURT_F_B[k], c18);
        let kp, c19 = statistics::kurtosis(xs, false, true);
        expect(bad, within_unit(kp, E_KURT_P_B[k], 8.0), "kurtosis pearson", k, kp, E_KURT_P_B[k], c19);
        let kfu, c20 = statistics::kurtosis(xs, true, false);
        expect(bad, within_unit(kfu, E_KURT_F_U[k], 16.0), "kurtosis fisher unbiased", k, kfu, E_KURT_F_U[k], c20);
        let kpu, c21 = statistics::kurtosis(xs, false, false);
        expect(bad, within_unit(kpu, E_KURT_P_U[k], 16.0), "kurtosis pearson unbiased", k, kpu, E_KURT_P_U[k], c21);
        let c3m, c22 = statistics::central_moment(xs, 3);
        expect(bad, within(c3m, E_CM3[k], 8.0) || abs(c3m - E_CM3[k]) <= 8.0 * EPS * E_PVAR[k] * sqrt(E_PVAR[k]), "central moment 3", k, c3m, E_CM3[k], c22);
        let c4m, c23 = statistics::central_moment(xs, 4);
        expect(bad, within(c4m, E_CM4[k], 8.0), "central moment 4", k, c4m, E_CM4[k], c23);
    }
    bad += run_quantile_grid();
    return bad;
}

// The nine quantile methods on 99 fractions of the first 200 samples. The
// expected values come from NumPy.
fn run_quantile_grid() -> i64 {
    var bad: i64 = 0;
    let samples = 200;
    for k in 0..samples {
        let xs = V_FLAT[V_OFFS[k]..V_OFFS[k + 1]];
        for mth in 0..9 {
            for i in 1..100 {
                let q = (i as f64) / 100.0;
                let v, cq = statistics::quantile(xs, q, mth);
                let e = E_Q99[k * 891 + mth * 99 + (i - 1)];
                expect(bad, within(v, e, 2.0), "quantile grid", k, v, e, cq);
            }
        }
    }
    return bad;
}
