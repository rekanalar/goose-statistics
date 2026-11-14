// Compare two survival curves with the logrank test, then test whether
// a sample fits the normal family. Run with: goose run examples/survival.goose
// The program must start in the package directory.
import std;
import statistics from "..";

fn main() {
    // Times in weeks. An uncensored time is a death. A right time is a
    // subject who left the study with no death seen.
    let ux: f64[8] = [6.0, 13.0, 21.0, 30.0, 37.0, 38.0, 49.0, 50.0];
    let rx: f64[3] = [31.0, 47.0, 80.0];
    let uy: f64[8] = [10.0, 10.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0];
    let ry: f64[2] = [34.0, 40.0];
    var none: f64[>..] = [];
    let a = statistics::CensoredData { uncensored: ux[..], left: none[..], right: rx[..], interval: none[..] };
    let b = statistics::CensoredData { uncensored: uy[..], left: none[..], right: ry[..], interval: none[..] };
    let r, c = statistics::logrank(a, b);
    if c != 0 { print("logrank failed, code ", c, "\n"); return 1; }
    print("logrank z ", r.statistic, " p ", r.pvalue, "\n");

    // The test of fit draws from a generator that the caller seeds.
    let xs: f64[10] = [9.1, 10.4, 11.2, 8.7, 10.9, 9.6, 12.1, 10.2, 9.9, 11.0];
    var g = rng { s: 7 };
    var fit: f64[>..] = [];
    let t, ct = statistics::goodness_of_fit(fit, xs[..], statistics::FAMILY_NORMAL, statistics::GOF_ANDERSON, 999, g);
    if ct != 0 { print("fit failed, code ", ct, "\n"); return 1; }
    print("anderson ", t.statistic, " p ", t.pvalue, " mean ", fit[0], " sd ", fit[1], "\n");

    // A sample of the chi-square law fits its family. The Filliben
    // statistic is a correlation, so a value near one means a good fit.
    let cs: f64[10] = [3.1, 4.2, 5.8, 6.4, 7.1, 8.3, 9.0, 10.2, 11.4, 12.1];
    var g2 = rng { s: 7 };
    var cfit: f64[>..] = [];
    let t2, ct2 = statistics::goodness_of_fit(cfit, cs[..], statistics::FAMILY_CHI2, statistics::GOF_FILLIBEN, 99, g2);
    if ct2 != 0 { print("chi-square fit failed, code ", ct2, "\n"); return 1; }
    print("filliben ", t2.statistic, " p ", t2.pvalue, " df ", cfit[0], " scale ", cfit[2], "\n");
    return 0;
}
