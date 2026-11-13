// The tests of the kernel density estimate. The density of a sample sums
// to one over a wide range. The logarithm of the density is the logarithm
// of the density. The factor of Silverman is above the factor of Scott.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_kde_tests() -> i64 {
    var bad: i64 = 0;
    let d1: f64[7] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
    let p1: f64[4] = [0.5, 2.5, 4.3, 7.5];
    var k, ck = statistics::kde_fit(d1[..], 1, [], statistics::KDE_SCOTT, 0.0);
    chk_code(bad, "kde code", ck, 0);
    // the density of a point of the sample is above the density far away
    var out: f64[>..] = [];
    let c1 = statistics::kde_pdf(k, out, p1[..]);
    chk_code(bad, "kde pdf code", c1, 0);
    check(bad, out.len == 4, "the density holds one value per point");
    check(bad, out[0] > 0.0 && out[3] > 0.0, "the density is above zero");
    check(bad, out[2] > out[0], "the density at the middle is above the density at the end");
    // the logarithm of the density is the logarithm of the density
    var logs: f64[>..] = [];
    let c2 = statistics::kde_logpdf(k, logs, p1[..]);
    chk_code(bad, "kde logpdf code", c2, 0);
    for i in 0..4 { chk_val(bad, "kde logpdf", logs[i], c2, log(out[i]), 1e-10); }
    // the density sums to one over a wide range
    let whole, cw = statistics::kde_integrate_box_1d(k, -40.0, 40.0);
    chk_code(bad, "kde box code", cw, 0);
    chk_val(bad, "kde whole", whole, cw, 1.0, 1e-6);
    // the factor of Silverman is above the factor of Scott
    var ks, cks = statistics::kde_fit(d1[..], 1, [], statistics::KDE_SILVERMAN, 0.0);
    chk_code(bad, "kde silverman code", cks, 0);
    check(bad, ks.factor > k.factor, "the factor of Silverman is above the factor of Scott");
    // the covariance of the estimate holds the band width
    var cv: f64[>..] = [];
    let ccv = statistics::kde_covariance(k, cv);
    chk_code(bad, "kde covariance code", ccv, 0);
    check(bad, cv.len == 1, "the covariance of a sample of one dimension holds one value");
    check(bad, cv[0] > 0.0, "the covariance is above zero");
    // the draws of the estimate keep the count and the order
    var draws: f64[>..] = [];
    var gen = rng { s: 7 };
    let cr = statistics::kde_resample(k, draws, 20, gen);
    chk_code(bad, "kde resample code", cr, 0);
    check(bad, draws.len == 20, "the draws hold one value per point of a sample of one dimension");
    var same = true;
    var gen2 = rng { s: 7 };
    var draws2: f64[>..] = [];
    let _cr2 = statistics::kde_resample(k, draws2, 20, gen2);
    for i in 0..20 { if draws[i] != draws2[i] { same = false; } }
    check(bad, same, "the same seed gives the same draws");
    // the draws carry the noise of the band width, so most of them lie
    // outside the values of the sample
    var off: i64 = 0;
    for i in 0..draws.len {
        var on: bool = false;
        for j in 0..d1.len { if draws[i] == d1[j] { on = true; } }
        if !on { off += 1; }
    }
    check(bad, off > draws.len / 2, "the draws carry the noise of the band width");
    // the weight of a point decides how often the draws come from it
    var heavy: f64[4] = [0.0, 10.0, 20.0, 30.0];
    var hw: f64[4] = [1.0, 1.0, 1.0, 100.0];
    let kh, ckh = statistics::kde_fit(heavy[..], 1, hw[..], statistics::KDE_SCALAR, 0.05);
    chk_code(bad, "the weighted fit code", ckh, 0);
    var hdraws: f64[>..] = [];
    var hgen = rng { s: 11 };
    let chr = statistics::kde_resample(kh, hdraws, 40, hgen);
    chk_code(bad, "the weighted draw code", chr, 0);
    var near: i64 = 0;
    for i in 0..hdraws.len { if abs(hdraws[i] - 30.0) < 2.0 { near += 1; } }
    check(bad, near > 30, "the draws come from the point with the large weight");
    // the faults
    let odd: f64[3] = [1.0, 2.0, 3.0];
    let _k1, ck1 = statistics::kde_fit(odd[..], 2, [], statistics::KDE_SCOTT, 0.0);
    chk_code(bad, "a sample that is not a multiple of the dimensions", ck1, 6);
    let one: f64[1] = [1.0];
    let _k2, ck2 = statistics::kde_fit(one[..], 1, [], statistics::KDE_SCOTT, 0.0);
    chk_code(bad, "a sample of one point", ck2, 2);
    let _k3, ck3 = statistics::kde_fit(d1[..], 1, [], 9, 0.0);
    chk_code(bad, "an unknown bandwidth", ck3, 5);
    let _k4, ck4 = statistics::kde_fit(d1[..], 1, [], statistics::KDE_SCALAR, 0.0);
    chk_code(bad, "a zero factor", ck4, 5);
    let flat: f64[4] = [2.0, 2.0, 2.0, 2.0];
    let _k5, ck5 = statistics::kde_fit(flat[..], 1, [], statistics::KDE_SCOTT, 0.0);
    chk_code(bad, "a sample without spread", ck5, 5);
    let zw: f64[7] = [1.0, 1.0, 2.0, 2.0, 3.0, 1.0, 0.0];
    let _k6, ck6 = statistics::kde_fit(d1[..], 1, zw[..], statistics::KDE_SCOTT, 0.0);
    chk_code(bad, "a weight of zero", ck6, 5);
    return bad;
}
