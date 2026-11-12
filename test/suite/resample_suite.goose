// The tests of the resampling calls. The measured rows live in
// resample_vectors.goose. The rows here check the faults, the order of
// the checks, the seed repeat, the no-write rule and the Sobol indices.
namespace suite;

import std;
import math;
import statistics from "../..";

// The Ishigami function of a point in the unit cube, mapped to the cube
// from minus pi to pi.
fn rt_ishigami(p: const f64[:]) -> f64 {
    let pi = 3.141592653589793;
    let x1 = (2.0 * p[0] - 1.0) * pi;
    let x2 = (2.0 * p[1] - 1.0) * pi;
    let x3 = (2.0 * p[2] - 1.0) * pi;
    return sin(x1) + 7.0 * sin(x2) * sin(x2) + 0.1 * pow(x3, 4.0) * sin(x1);
}

fn rt_same(x: const f64[:], y: const f64[:]) -> bool {
    if x.len != y.len { return false; }
    for i in 0..x.len {
        if x[i] != y[i] { return false; }
    }
    return true;
}

fn run_resample_tests() -> i64 {
    var bad: i64 = 0;
    let a: f64[8] = [2.1, 3.4, 1.9, 5.6, 4.4, 3.8, 2.7, 4.9];
    let b: f64[8] = [3.9, 5.2, 4.8, 6.1, 5.5, 4.1, 6.6, 5.0];
    var ab: f64[16];
    for i in 0..8 { ab[i] = a[i]; ab[8 + i] = b[i]; }
    let one: i64[1] = [8];
    let two: i64[2] = [8, 8];
    var nul: f64[>..] = [];
    var r = rng { s: 99 };

    let bm = statistics::ResampleMethod { resamples: 200, .. };
    let pm = statistics::ResampleMethod { kind: statistics::RESAMPLE_PERMUTATION, resamples: 200, .. };
    let mm = statistics::ResampleMethod { kind: statistics::RESAMPLE_MONTE_CARLO, resamples: 200, .. };
    let np: f64[2] = [3.0, 1.5];

    // The defaults of the struct.
    let dm = statistics::ResampleMethod { .. };
    check(bad, dm.kind == statistics::RESAMPLE_BOOTSTRAP && dm.resamples == 9999
        && dm.interval == statistics::INTERVAL_BCA && dm.permutation == statistics::PERMUTE_INDEPENDENT, "the defaults of the method");

    // The seed repeat. The same seed gives the same interval, the same
    // null list and the same p-value.
    var r1 = rng { s: 12 };
    var r2 = rng { s: 12 };
    var n1: f64[>..] = [];
    var n2: f64[>..] = [];
    let x1, k1 = statistics::bootstrap(n1, a[..], one[..], bm, 0.9, statistics::ALT_TWO_SIDED, r1) { s => mean_of(s) };
    let x2, k2 = statistics::bootstrap(n2, a[..], one[..], bm, 0.9, statistics::ALT_TWO_SIDED, r2) { s => mean_of(s) };
    check(bad, k1 == 0 && k2 == 0 && x1.interval.low == x2.interval.low && x1.interval.high == x2.interval.high
        && x1.standard_error == x2.standard_error && rt_same(n1[..], n2[..]), "bootstrap, the same seed gives the same interval");
    check(bad, x1.interval.low < x1.interval.high && x1.interval.low > 2.0 && x1.interval.high < 5.0, "bootstrap, the interval brackets the mean");
    var r3 = rng { s: 13 };
    let x3, k3 = statistics::bootstrap(n2, a[..], one[..], bm, 0.9, statistics::ALT_TWO_SIDED, r3) { s => mean_of(s) };
    check(bad, k3 == 0 && !rt_same(n1[..], n2[..]), "bootstrap, another seed gives another list");
    var r4 = rng { s: 12 };
    var r5 = rng { s: 12 };
    let t1, k4 = statistics::permutation_test(n1, ab[..], two[..], pm, statistics::ALT_TWO_SIDED, r4) { s => mean_of(s[0..8]) - mean_of(s[8..16]) };
    let t2, k5 = statistics::permutation_test(n2, ab[..], two[..], pm, statistics::ALT_TWO_SIDED, r5) { s => mean_of(s[0..8]) - mean_of(s[8..16]) };
    check(bad, k4 == 0 && k5 == 0 && t1.pvalue == t2.pvalue && rt_same(n1[..], n2[..]) && n1.len == 200, "permutation, the same seed gives the same list");
    // A drawn p-value adds 1 to the count and to the total.
    check(bad, t1.pvalue >= 1.0 / 201.0 && t1.pvalue <= 1.0, "permutation, the drawn p-value is above 1 over total plus 1");
    var r6 = rng { s: 12 };
    var r7 = rng { s: 12 };
    let m1, k6 = statistics::monte_carlo_test(n1, a[..], statistics::FAMILY_NORMAL, np[..], mm, statistics::ALT_TWO_SIDED, r6) { s => mean_of(s) };
    let m2, k7 = statistics::monte_carlo_test(n2, a[..], statistics::FAMILY_NORMAL, np[..], mm, statistics::ALT_TWO_SIDED, r7) { s => mean_of(s) };
    check(bad, k6 == 0 && k7 == 0 && m1.pvalue == m2.pvalue && rt_same(n1[..], n2[..]), "monte carlo, the same seed gives the same list");
    var r8 = rng { s: 12 };
    var r9 = rng { s: 12 };
    let p1, k8 = statistics::power(statistics::FAMILY_NORMAL, np[..], 8, 0.05, 300, r8) { s => mean_of(s) / 100.0 };
    let p2, k9 = statistics::power(statistics::FAMILY_NORMAL, np[..], 8, 0.05, 300, r9) { s => mean_of(s) / 100.0 };
    check(bad, k8 == 0 && k9 == 0 && p1 == p2 && p1 == 1.0, "power, the same seed gives the same value");
    // A seed of 0 is not special, and a huge seed runs.
    var rz = rng { s: 0 };
    let z1, kz = statistics::bootstrap(n1, a[..], one[..], bm, 0.9, statistics::ALT_TWO_SIDED, rz) { s => mean_of(s) };
    check(bad, kz == 0 && z1.interval.low < z1.interval.high, "bootstrap, a seed of 0");
    var rf = rng { s: 18446744073709551615 };
    let z2, kf = statistics::bootstrap(n1, a[..], one[..], bm, 0.9, statistics::ALT_TWO_SIDED, rf) { s => mean_of(s) };
    check(bad, kf == 0, "bootstrap, the largest seed");

    // A one-sided request leaves the open end at infinity.
    let o1, ko1 = statistics::bootstrap(n1, a[..], one[..], bm, 0.9, statistics::ALT_LESS, r) { s => mean_of(s) };
    check(bad, ko1 == 0 && o1.interval.low < -1e300 && o1.interval.high < 10.0, "bootstrap, the less side");
    let o2, ko2 = statistics::bootstrap(n1, a[..], one[..], bm, 0.9, statistics::ALT_GREATER, r) { s => mean_of(s) };
    check(bad, ko2 == 0 && o2.interval.high > 1e300 && o2.interval.low > 0.0, "bootstrap, the greater side");

    // A sample of two samples, percentile and basic: the block sees both.
    let pm2 = statistics::ResampleMethod { resamples: 200, interval: statistics::INTERVAL_PERCENTILE, .. };
    let tw, ktw = statistics::bootstrap(n1, ab[..], two[..], pm2, 0.9, statistics::ALT_TWO_SIDED, r) { s => mean_of(s[0..8]) - mean_of(s[8..16]) };
    check(bad, ktw == 0 && tw.interval.high < 0.5 && tw.interval.low < tw.interval.high, "bootstrap, two samples");
    let _bx, kbc2 = statistics::bootstrap(n1, ab[..], two[..], bm, 0.9, statistics::ALT_TWO_SIDED, r) { s => mean_of(s[0..8]) - mean_of(s[8..16]) };
    check(bad, kbc2 == 5, "bootstrap, the bca form refuses two samples");

    // The faults of the bootstrap, in the order of the checks.
    let nanl: f64[8] = [2.1, 3.4, rt_nan(), 5.6, 4.4, 3.8, 2.7, 4.9];
    let infl: f64[8] = [2.1, 3.4, 1.0 / 0.0, 5.6, 4.4, 3.8, 2.7, 4.9];
    let em: f64[0] = [];
    let ez: i64[0] = [];
    let _f1, c1 = statistics::bootstrap(nul, em[..], ez[..], bm, 0.9, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    check(bad, c1 == statistics::ERR_EMPTY, "bootstrap, an empty sample");
    let _f2, c2 = statistics::bootstrap(nul, nanl[..], one[..], bm, 0.9, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    check(bad, c2 == statistics::ERR_NAN, "bootstrap, a NaN");
    let _f3, c3 = statistics::bootstrap(nul, infl[..], one[..], bm, 0.9, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    check(bad, c3 == statistics::ERR_INFINITY, "bootstrap, an infinity");
    let s7: i64[1] = [7];
    let _f4, c4 = statistics::bootstrap(nul, a[..], s7[..], bm, 0.9, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    check(bad, c4 == statistics::ERR_LENGTH, "bootstrap, sizes that miss the length");
    let s1: i64[2] = [1, 7];
    let _f5, c5 = statistics::bootstrap(nul, a[..], s1[..], bm, 0.9, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    check(bad, c5 == statistics::ERR_TOO_SHORT, "bootstrap, a sample of one value");
    let _f6, c6 = statistics::bootstrap(nul, a[..], one[..], bm, 0.0, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    let _f7, c7 = statistics::bootstrap(nul, a[..], one[..], bm, 1.0, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    let _f8, c8 = statistics::bootstrap(nul, a[..], one[..], bm, rt_nan(), statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    check(bad, c6 == statistics::ERR_DOMAIN && c7 == statistics::ERR_DOMAIN && c8 == statistics::ERR_NAN, "bootstrap, a level outside 0 and 1");
    let _f9, c9 = statistics::bootstrap(nul, a[..], one[..], bm, 0.9, 3, r) { s => mean_of(s) };
    check(bad, c9 == statistics::ERR_DOMAIN, "bootstrap, an unknown alternative");
    let bi = statistics::ResampleMethod { resamples: 200, interval: 7, .. };
    let _g1, d1 = statistics::bootstrap(nul, a[..], one[..], bi, 0.9, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    check(bad, d1 == statistics::ERR_DOMAIN, "bootstrap, an unknown interval constant");
    let bk = statistics::ResampleMethod { kind: 9, resamples: 200, .. };
    let _g2, d2 = statistics::bootstrap(nul, a[..], one[..], bk, 0.9, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    check(bad, d2 == statistics::ERR_DOMAIN, "bootstrap, an unknown kind constant");
    let b0 = statistics::ResampleMethod { resamples: 0, .. };
    let _g3, d3 = statistics::bootstrap(nul, a[..], one[..], b0, 0.9, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    check(bad, d3 == statistics::ERR_DOMAIN, "bootstrap, no resamples");
    let b1 = statistics::ResampleMethod { resamples: 1, .. };
    let _g4, d4 = statistics::bootstrap(nul, a[..], one[..], b1, 0.9, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    check(bad, d4 == statistics::ERR_TOO_SHORT, "bootstrap, one resample");
    let big = statistics::ResampleMethod { resamples: 1000001, .. };
    let _g5, d5 = statistics::bootstrap(nul, a[..], one[..], big, 0.9, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    check(bad, d5 == statistics::ERR_TOO_LONG && nul.len == 0, "bootstrap, above the resample bound");
    var wide: f64[>..] = [];
    for i in 0..600 { wide.push(i * 0.37 + (i % 7)); }
    let sw: i64[1] = [600];
    let work = statistics::ResampleMethod { resamples: 1000000, interval: statistics::INTERVAL_PERCENTILE, .. };
    let _g6, d6 = statistics::bootstrap(nul, wide[..], sw[..], work, 0.9, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    check(bad, d6 == statistics::ERR_TOO_LONG && nul.len == 0, "bootstrap, above the work bound");
    let kc: f64[8] = [2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0];
    let _g7, d7 = statistics::bootstrap(nul, kc[..], one[..], bm, 0.9, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    check(bad, d7 == statistics::ERR_UNDEFINED, "bootstrap, a constant sample");
    let _g8, d8 = statistics::bootstrap(nul, a[..], one[..], bm, 0.9, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) * (0.0 / 0.0) };
    check(bad, d8 == statistics::ERR_NAN && nul.len == 0, "bootstrap, a block that returns NaN");
    let _g9, d9 = statistics::bootstrap(nul, a[..], one[..], bm, 0.9, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) / 0.0 };
    check(bad, d9 == statistics::ERR_INFINITY, "bootstrap, a block that returns infinity");
    // A block that fails only inside the jackknife.
    let _h1, e1 = statistics::bootstrap(nul, a[..], one[..], bm, 0.9, statistics::ALT_TWO_SIDED, r) { s => if s.len < 8 { 0.0 / 0.0 } else { mean_of(s) } };
    check(bad, e1 == statistics::ERR_NAN, "bootstrap, a block that fails in the jackknife");
    // The check order: a bad length bound beats a bad level.
    let _h2, e2 = statistics::bootstrap(nul, nanl[..], one[..], bm, 7.0, 9, r) { s => mean_of(s) };
    check(bad, e2 == statistics::ERR_NAN, "bootstrap, the data check comes before the argument checks");

    // The faults of the permutation test.
    let _p1, q1 = statistics::permutation_test(nul, em[..], ez[..], pm, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    check(bad, q1 == statistics::ERR_EMPTY, "permutation, an empty sample");
    let _p2, q2 = statistics::permutation_test(nul, nanl[..], one[..], pm, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    check(bad, q2 == statistics::ERR_NAN, "permutation, a NaN");
    let _p3, q3 = statistics::permutation_test(nul, a[..], one[..], pm, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    check(bad, q3 == statistics::ERR_TOO_SHORT, "permutation, one sample of the independent type");
    let _p4, q4 = statistics::permutation_test(nul, ab[..], two[..], bm, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    check(bad, q4 == statistics::ERR_DOMAIN, "permutation, the kind of another call");
    let pt = statistics::ResampleMethod { kind: statistics::RESAMPLE_PERMUTATION, resamples: 200, permutation: 5, .. };
    let _p5, q5 = statistics::permutation_test(nul, ab[..], two[..], pt, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    check(bad, q5 == statistics::ERR_DOMAIN, "permutation, an unknown permutation type");
    let _p6, q6 = statistics::permutation_test(nul, ab[..], two[..], pm, 4, r) { s => mean_of(s) };
    check(bad, q6 == statistics::ERR_DOMAIN, "permutation, an unknown alternative");
    let pbig = statistics::ResampleMethod { kind: statistics::RESAMPLE_PERMUTATION, resamples: 2000000, .. };
    let _p7, q7 = statistics::permutation_test(nul, ab[..], two[..], pbig, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    check(bad, q7 == statistics::ERR_TOO_LONG && nul.len == 0, "permutation, above the resample bound");
    let sw2: i64[2] = [300, 300];
    var wide2: f64[>..] = [];
    for i in 0..600 { wide2.push(i * 0.37 + (i % 7)); }
    let pwork = statistics::ResampleMethod { kind: statistics::RESAMPLE_PERMUTATION, resamples: 1000000, .. };
    let _p8, q8 = statistics::permutation_test(nul, wide2[..], sw2[..], pwork, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    check(bad, q8 == statistics::ERR_TOO_LONG && nul.len == 0, "permutation, above the work bound");
    let kc2: f64[16] = [2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0];
    let _pa, qa = statistics::permutation_test(nul, kc2[..], two[..], pm, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    check(bad, qa == statistics::ERR_UNDEFINED, "permutation, a constant sample");
    let _pb, qb = statistics::permutation_test(nul, ab[..], two[..], pm, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) * (0.0 / 0.0) };
    check(bad, qb == statistics::ERR_NAN, "permutation, a block that returns NaN");
    let _pc, qc = statistics::permutation_test(nul, ab[..], two[..], pm, statistics::ALT_TWO_SIDED, r) { s => 1.0 / 0.0 };
    check(bad, qc == statistics::ERR_INFINITY, "permutation, a block that returns infinity");
    let sm: i64[2] = [8, 4];
    let sp = statistics::ResampleMethod { kind: statistics::RESAMPLE_PERMUTATION, resamples: 200, permutation: statistics::PERMUTE_SAMPLES, .. };
    var ab12: f64[12];
    for i in 0..12 { ab12[i] = ab[i]; }
    let _pd, qd = statistics::permutation_test(nul, ab12[..], sm[..], sp, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    check(bad, qd == statistics::ERR_LENGTH, "permutation, paired samples of unequal length");
    // The drawn forms of the three types give a list of the resample count.
    let ds = statistics::ResampleMethod { kind: statistics::RESAMPLE_PERMUTATION, resamples: 100, permutation: statistics::PERMUTE_SAMPLES, .. };
    let _pe, qe = statistics::permutation_test(nul, ab[..], two[..], ds, statistics::ALT_TWO_SIDED, r) { s => mean_of(s[0..8]) - mean_of(s[8..16]) };
    check(bad, qe == 0 && nul.len == 100, "permutation, drawn swaps of pairs");
    let dg = statistics::ResampleMethod { kind: statistics::RESAMPLE_PERMUTATION, resamples: 100, permutation: statistics::PERMUTE_PAIRINGS, .. };
    let _pf, qf = statistics::permutation_test(nul, ab[..], two[..], dg, statistics::ALT_TWO_SIDED, r) { s => mean_of(s[0..8]) * mean_of(s[8..16]) };
    check(bad, qf == 0 && nul.len == 100, "permutation, drawn pairings");
    let dsf = statistics::ResampleMethod { kind: statistics::RESAMPLE_PERMUTATION, resamples: 100, permutation: statistics::PERMUTE_SAMPLES, .. };
    let _pg, qg = statistics::permutation_test(nul, ab[..], one[..], dsf, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    check(bad, qg == statistics::ERR_LENGTH, "permutation, the sizes miss the length");

    // The faults of the Monte Carlo test.
    let _m1, u1 = statistics::monte_carlo_test(nul, em[..], statistics::FAMILY_NORMAL, np[..], mm, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    check(bad, u1 == statistics::ERR_EMPTY, "monte carlo, an empty sample");
    let _m2, u2 = statistics::monte_carlo_test(nul, nanl[..], statistics::FAMILY_NORMAL, np[..], mm, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    check(bad, u2 == statistics::ERR_NAN, "monte carlo, a NaN");
    let _m3, u3 = statistics::monte_carlo_test(nul, infl[..], statistics::FAMILY_NORMAL, np[..], mm, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    check(bad, u3 == statistics::ERR_INFINITY, "monte carlo, an infinity");
    let _m4, u4 = statistics::monte_carlo_test(nul, a[..], 12, np[..], mm, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    check(bad, u4 == statistics::ERR_DOMAIN, "monte carlo, an unknown family");
    let one_p: f64[1] = [3.0];
    let _m5, u5 = statistics::monte_carlo_test(nul, a[..], statistics::FAMILY_NORMAL, one_p[..], mm, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    check(bad, u5 == statistics::ERR_LENGTH, "monte carlo, the wrong count of parameters");
    let badp: f64[2] = [3.0, rt_nan()];
    let _m6, u6 = statistics::monte_carlo_test(nul, a[..], statistics::FAMILY_NORMAL, badp[..], mm, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    check(bad, u6 == statistics::ERR_NAN, "monte carlo, a NaN parameter");
    let negp: f64[2] = [3.0, -1.0];
    let _m7, u7 = statistics::monte_carlo_test(nul, a[..], statistics::FAMILY_NORMAL, negp[..], mm, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    check(bad, u7 != 0 && nul.len == 0, "monte carlo, a negative scale");
    let _m8, uu8 = statistics::monte_carlo_test(nul, a[..], statistics::FAMILY_NORMAL, np[..], bm, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    check(bad, uu8 == statistics::ERR_DOMAIN, "monte carlo, the kind of another call");
    let mbig = statistics::ResampleMethod { kind: statistics::RESAMPLE_MONTE_CARLO, resamples: 1000001, .. };
    let _m9, u9 = statistics::monte_carlo_test(nul, a[..], statistics::FAMILY_NORMAL, np[..], mbig, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    check(bad, u9 == statistics::ERR_TOO_LONG, "monte carlo, above the resample bound");
    let _ma, ua = statistics::monte_carlo_test(nul, a[..], statistics::FAMILY_NORMAL, np[..], mm, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) * (0.0 / 0.0) };
    check(bad, ua == statistics::ERR_NAN, "monte carlo, a block that returns NaN");
    let _mb, ub = statistics::monte_carlo_test(nul, a[..], statistics::FAMILY_NORMAL, np[..], mm, statistics::ALT_TWO_SIDED, r) { s => 1.0 / 0.0 };
    check(bad, ub == statistics::ERR_INFINITY, "monte carlo, a block that returns infinity");
    // Every family draws without a fault.
    var every = true;
    let fam1: f64[1] = [2.0];
    let fam2: f64[2] = [2.0, 3.0];
    for f in 0..10 {
        var pr = fam2[..];
        if f == statistics::FAMILY_EXPONENTIAL || f == statistics::FAMILY_CHI2 || f == statistics::FAMILY_T { pr = fam1[..]; }
        let _mf, uf = statistics::monte_carlo_test(nul, a[..], f, pr, rt_mc(20), statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
        if uf != 0 || nul.len != 20 { every = false; }
    }
    check(bad, every, "monte carlo, every family draws");

    // The faults of the power.
    let _w1, v1 = statistics::power(statistics::FAMILY_NORMAL, np[..], 0, 0.05, 100, r) { s => 0.5 };
    check(bad, v1 == statistics::ERR_EMPTY, "power, no observations");
    let _w2, v2 = statistics::power(statistics::FAMILY_NORMAL, np[..], 8, 0.0, 100, r) { s => 0.5 };
    let _w3, v3 = statistics::power(statistics::FAMILY_NORMAL, np[..], 8, 1.0, 100, r) { s => 0.5 };
    check(bad, v2 == statistics::ERR_DOMAIN && v3 == statistics::ERR_DOMAIN, "power, a level outside 0 and 1");
    let _w4, v4 = statistics::power(statistics::FAMILY_NORMAL, np[..], 8, 0.05, 0, r) { s => 0.5 };
    check(bad, v4 == statistics::ERR_DOMAIN, "power, no resamples");
    let _w5, v5 = statistics::power(statistics::FAMILY_NORMAL, np[..], 8, 0.05, 1000001, r) { s => 0.5 };
    check(bad, v5 == statistics::ERR_TOO_LONG, "power, above the resample bound");
    let _w6, v6 = statistics::power(statistics::FAMILY_NORMAL, np[..], 600, 0.05, 1000000, r) { s => 0.5 };
    check(bad, v6 == statistics::ERR_TOO_LONG, "power, above the work bound");
    let _w7, v7 = statistics::power(statistics::FAMILY_NORMAL, np[..], 4194305, 0.05, 10, r) { s => 0.5 };
    check(bad, v7 == statistics::ERR_TOO_LONG, "power, above the input bound");
    let _w8, v8 = statistics::power(statistics::FAMILY_NORMAL, np[..], 8, 0.05, 10, r) { s => 1.5 };
    check(bad, v8 == statistics::ERR_DOMAIN, "power, a p-value above 1");
    let _w9, v9 = statistics::power(statistics::FAMILY_NORMAL, np[..], 8, 0.05, 10, r) { s => 0.0 / 0.0 };
    check(bad, v9 == statistics::ERR_NAN, "power, a block that returns NaN");
    let _wa, va = statistics::power(statistics::FAMILY_NORMAL, np[..], 8, 0.05, 10, r) { s => 1.0 / 0.0 };
    check(bad, va == statistics::ERR_INFINITY, "power, a block that returns infinity");
    let _wb, vb = statistics::power(99, np[..], 8, 0.05, 10, r) { s => 0.5 };
    check(bad, vb == statistics::ERR_DOMAIN, "power, an unknown family");
    // The tolerance of the level: a p-value equal to the level is no hit.
    let pe, ke = statistics::power(statistics::FAMILY_NORMAL, np[..], 8, 0.05, 50, r) { s => 0.05 };
    check(bad, ke == 0 && pe == 0.0, "power, a p-value at the level is no hit");

    // The faults of the Sobol call.
    var fi: f64[>..] = [];
    var ft: f64[>..] = [];
    let sa, sk = statistics::sobol_indices(fi, ft, 3, 4096, r) { p => rt_ishigami(p) };
    check(bad, sk == 0 && sa == 3 && fi.len == 3 && ft.len == 3, "sobol, three indices of each kind");
    // The first-order indices of the Ishigami function are 0.3139, 0.4424
    // and 0. The total indices are 0.5576, 0.4424 and 0.2437.
    check(bad, abs(fi[0] - 0.3139) < 0.05 && abs(fi[1] - 0.4424) < 0.05 && abs(fi[2]) < 0.05, "sobol, the first-order indices of the Ishigami function");
    check(bad, abs(ft[0] - 0.5576) < 0.05 && abs(ft[1] - 0.4424) < 0.05 && abs(ft[2] - 0.2437) < 0.05, "sobol, the total indices of the Ishigami function");
    var rs1 = rng { s: 31 };
    var rs2 = rng { s: 31 };
    var fj: f64[>..] = [];
    var tj: f64[>..] = [];
    let _s1, sk1 = statistics::sobol_indices(fi, ft, 3, 512, rs1) { p => rt_ishigami(p) };
    let _s2, sk2 = statistics::sobol_indices(fj, tj, 3, 512, rs2) { p => rt_ishigami(p) };
    check(bad, sk1 == 0 && sk2 == 0 && rt_same(fi[..], fj[..]) && rt_same(ft[..], tj[..]), "sobol, the same seed gives the same indices");
    // An additive model has equal first and total indices.
    var ra = rng { s: 41 };
    let ad, ka = statistics::sobol_indices(fi, ft, 2, 4096, ra) { p => p[0] + 2.0 * p[1] };
    check(bad, ka == 0 && abs(fi[0] - 0.2) < 0.05 && abs(fi[1] - 0.8) < 0.05 && abs(ft[1] - 0.8) < 0.05, "sobol, an additive model");
    let _s3, k3x = statistics::sobol_indices(fi, ft, 33, 100, r) { p => p[0] };
    check(bad, k3x == statistics::ERR_TOO_LONG && fi.len == 0, "sobol, above the dimension bound");
    let _s4, k4x = statistics::sobol_indices(fi, ft, 3, 1048577, r) { p => p[0] };
    check(bad, k4x == statistics::ERR_TOO_LONG, "sobol, above the point bound");
    let _s5, k5x = statistics::sobol_indices(fi, ft, 0, 100, r) { p => 1.0 };
    check(bad, k5x == statistics::ERR_DOMAIN, "sobol, no variable");
    let _s6, k6x = statistics::sobol_indices(fi, ft, 3, 1, r) { p => p[0] };
    check(bad, k6x == statistics::ERR_TOO_SHORT, "sobol, one point");
    let _s7, k7x = statistics::sobol_indices(fi, ft, 3, 100, r) { p => 4.0 };
    check(bad, k7x == statistics::ERR_UNDEFINED, "sobol, a constant model");
    let _s8, k8x = statistics::sobol_indices(fi, ft, 3, 100, r) { p => 0.0 / 0.0 };
    check(bad, k8x == statistics::ERR_NAN, "sobol, a model that returns NaN");
    let _s9, k9x = statistics::sobol_indices(fi, ft, 3, 100, r) { p => 1.0 / 0.0 };
    check(bad, k9x == statistics::ERR_INFINITY, "sobol, a model that returns infinity");
    // The faults of the estimators on the evaluation lists.
    let ea: f64[4] = [1.0, 2.0, 3.0, 4.5];
    let eb: f64[4] = [2.0, 1.0, 4.0, 3.0];
    let eab: f64[8] = [1.5, 2.5, 2.0, 3.0, 2.0, 1.0, 5.0, 3.5];
    let _e1, x1k = statistics::sobol_from_evaluations(fi, ft, ea[..], eb[..], eab[..], 2);
    check(bad, x1k == 0 && fi.len == 2, "sobol, the estimators accept a good list");
    let _e2, x2k = statistics::sobol_from_evaluations(fi, ft, ea[..], eb[..], eab[..], 3);
    check(bad, x2k == statistics::ERR_LENGTH, "sobol, the page list misses the length");
    let ebs: f64[3] = [2.0, 1.0, 4.0];
    let _e3, x3k = statistics::sobol_from_evaluations(fi, ft, ea[..], ebs[..], eab[..], 2);
    check(bad, x3k == statistics::ERR_LENGTH, "sobol, the lists A and B differ in length");
    let ean: f64[4] = [1.0, rt_nan(), 3.0, 4.5];
    let _e4, x4k = statistics::sobol_from_evaluations(fi, ft, ean[..], eb[..], eab[..], 2);
    check(bad, x4k == statistics::ERR_NAN, "sobol, a NaN in an evaluation");
    let eak: f64[4] = [2.0, 2.0, 2.0, 2.0];
    let _e5, x5k = statistics::sobol_from_evaluations(fi, ft, eak[..], eak[..], eab[..], 2);
    check(bad, x5k == statistics::ERR_UNDEFINED, "sobol, constant evaluations");
    let _e6, x6k = statistics::sobol_from_evaluations(fi, ft, ea[..], eb[..], eab[..], 33);
    check(bad, x6k == statistics::ERR_TOO_LONG, "sobol, the estimators refuse above the dimension bound");
    let _e7, x7k = statistics::sobol_from_evaluations(fi, ft, em[..], em[..], em[..], 2);
    check(bad, x7k == statistics::ERR_EMPTY, "sobol, empty evaluation lists");

    // No call writes into the list of the caller.
    var keep: f64[8] = [2.1, 3.4, 1.9, 5.6, 4.4, 3.8, 2.7, 4.9];
    var keep2: f64[16];
    for i in 0..16 { keep2[i] = ab[i]; }
    let _k1, kk1 = statistics::bootstrap(nul, keep[..], one[..], bm, 0.9, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    let _k2, kk2 = statistics::permutation_test(nul, keep2[..], two[..], pm, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    let _k3, kk3 = statistics::permutation_test(nul, keep2[..], two[..], ds, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    let _k4, kk4 = statistics::permutation_test(nul, keep2[..], two[..], dg, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    let _k5, kk5 = statistics::monte_carlo_test(nul, keep[..], statistics::FAMILY_NORMAL, np[..], mm, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    var same = kk1 == 0 && kk2 == 0 && kk3 == 0 && kk4 == 0 && kk5 == 0;
    for i in 0..8 { if keep[i] != a[i] { same = false; } }
    for i in 0..16 { if keep2[i] != ab[i] { same = false; } }
    check(bad, same, "no call writes into the list of the caller");

    // A fault leaves no null list behind.
    let _n1, nk1 = statistics::bootstrap(nul, a[..], one[..], bm, 0.9, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    let _n2, nk2 = statistics::bootstrap(nul, a[..], one[..], bm, 7.0, statistics::ALT_TWO_SIDED, r) { s => mean_of(s) };
    check(bad, nk1 == 0 && nk2 != 0 && nul.len == 0, "a fault clears the null list");
    return bad;
}

fn rt_nan() -> f64 { return 0.0 / 0.0; }

fn rt_mc(n: i64) -> statistics::ResampleMethod {
    return statistics::ResampleMethod { kind: statistics::RESAMPLE_MONTE_CARLO, resamples: n, .. };
}
