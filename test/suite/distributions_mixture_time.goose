// Time rows for the mixture of laws: a mixture with many components
// must end in a bounded time.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_distributions_mixture_time() -> i64 {
    var bad: i64 = 0;
    var comps: statistics::Dist[>..] = [];
    var w: f64[>..] = [];
    for i in 0..40 {
        comps.push(statistics::Dist { kind: 1, a: i as f64 - 20.0, b: 1.0 });
        w.push(0.025);
    }
    var s: f64 = 0.0;
    for k in 0..20 {
        let v, c = statistics::mixture_cdf(comps[..], w[..], k as f64 - 10.0);
        if c != 0 { bad += 1; }
        s = s + v;
    }
    check(bad, s > 0.0 && s <= 20.0, "the cdf of a large mixture");
    let q, cq = statistics::mixture_icdf(comps[..], w[..], 0.37);
    check(bad, cq == 0 && q > -20.0 && q < 20.0, "the quantile of a large mixture");
    return bad;
}
