// The exact and distance tests: two exact tests of a 2 by 2 table, two
// distance tests of two samples, the k-sample Anderson-Darling test, the
// Kolmogorov-Smirnov dispatchers, the Poisson means test and the Page
// trend test.
//
// Run it with: goose run examples/exact_tests.goose
import std;
import statistics from "..";

fn main() {
    let a: f64[8] = [2.1, 3.4, 1.9, 5.6, 4.4, 3.8, 2.7, 4.9];
    let b: f64[8] = [3.9, 5.2, 4.8, 6.1, 5.5, 4.1, 6.6, 5.0];
    let c: f64[8] = [1.2, 2.5, 2.2, 3.1, 1.8, 2.9, 2.4, 3.3];

    // The exact tests of a 2 by 2 table. The count of the first row of the
    // first column is 7, of the second column 3, and so on.
    let barn, c1 = statistics::barnard_exact(7.0, 3.0, 2.0, 8.0);
    print("barnard   statistic = ", barn.statistic, ", p = ", barn.pvalue, ", code ", c1, "\n");
    let bos, c2 = statistics::boschloo_exact(7.0, 3.0, 2.0, 8.0);
    print("boschloo  statistic = ", bos.statistic, ", p = ", bos.pvalue, ", code ", c2, "\n");

    // The distance tests of two samples.
    let cv, c3 = statistics::cramervonmises_2samp(a[..], b[..]);
    print("cvm 2samp statistic = ", cv.statistic, ", p = ", cv.pvalue, ", code ", c3, "\n");
    let es, c4 = statistics::epps_singleton_2samp(a[..], b[..]);
    print("epps      statistic = ", es.statistic, ", p = ", es.pvalue, ", df ", es.df, ", code ", c4, "\n");
    var g = rng { s: 1 };
    let bw, c5 = statistics::bws_test(a[..], b[..], g);
    print("bws       statistic = ", bw.statistic, ", p = ", bw.pvalue, ", code ", c5, "\n");

    // The Kolmogorov-Smirnov dispatchers. The list out holds the location
    // of the statistic and its sign.
    var out: f64[>..] = [];
    let ks, c6 = statistics::kstest_2samp(out, a[..], b[..]);
    print("ks 2samp  statistic = ", ks.statistic, ", p = ", ks.pvalue, ", location ", out[0], ", sign ", out[1], ", code ", c6, "\n");
    let np: f64[2] = [3.5, 1.5];
    let k1, c7 = statistics::kstest(out, a[..], statistics::FAMILY_NORMAL, np[..]);
    print("ks 1samp  statistic = ", k1.statistic, ", p = ", k1.pvalue, ", location ", out[0], ", sign ", out[1], ", code ", c7, "\n");

    // The k-sample Anderson-Darling test. Three samples sit in one buffer,
    // and a list of sizes gives the length of each sample.
    var all: f64[>..] = [];
    for v in a { all.push(v); }
    for v in b { all.push(v); }
    for v in c { all.push(v); }
    let sizes: i64[3] = [8, 8, 8];
    var crit: f64[>..] = [];
    let ad, c8 = statistics::anderson_ksamp(crit, all[..], sizes[..]);
    print("anderson  statistic = ", ad.statistic, ", p = ", ad.pvalue, ", code ", c8, "\n");
    print("anderson  critical values: ", crit[0], " ", crit[3], " ", crit[6], "\n");

    // The Poisson means test.
    let pm, c9 = statistics::poisson_means_test(10.0, 2.0, 18.0, 2.5);
    print("poisson   statistic = ", pm.statistic, ", p = ", pm.pvalue, ", code ", c9, "\n");

    // The Page trend test. The data hold one subject for each row, and
    // three values in a row.
    var data: f64[>..] = [];
    for i in 0..8 {
        data.push(a[i]);
        data.push(b[i]);
        data.push(c[i]);
    }
    let pg, c10 = statistics::page_trend_test(data[..], 3);
    print("page      statistic = ", pg.statistic, ", p = ", pg.pvalue, ", method ", pg.df, ", code ", c10, "\n");
    return 0;
}
