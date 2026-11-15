// Time rows for the distribution objects: a call on a large input must
// end in a bounded time.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_objects_time() -> i64 {
    var bad: i64 = 0;
    var edges: f64[>..] = [];
    var content: f64[>..] = [];
    for i in 0..20001 { edges.push(i); }
    for i in 0..20000 { content.push(1.0 + (i % 7)); }
    var s: f64 = 0.0;
    for k in 0..50 {
        let v, c = statistics::histogram_cdf(content[..], edges[..], true, 50.0 * k + 0.5);
        if c != 0 { bad += 1; }
        s = s + v;
    }
    let q, cq = statistics::histogram_ppf(content[..], edges[..], true, 0.37);
    check(bad, cq == 0 && q > 0.0, "the quantile of a large histogram");
    check(bad, s > 0.0, "the cdf of a large histogram");
    let b = statistics::Dist { kind: 4, a: 2000.0, b: 0.3 };
    let m, cm = statistics::dist_entropy(b);
    check(bad, cm == 0 && m > 0.0, "the entropy of a large binomial law");
    return bad;
}
