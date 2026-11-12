// The power transforms: the Box-Cox and the Yeo-Johnson transform, the
// search for the best lambda and the normal plot table.
//
// Run it with: goose run examples/transforms.goose
import std;
import statistics from "..";

fn main() {
    let skewed: f64[8] = [0.8, 1.5, 2.2, 3.1, 4.7, 6.3, 9.8, 12.5];
    let mixed: f64[8] = [-2.1, -0.4, 0.3, 1.1, 1.9, 2.6, 4.8, 7.5];
    var out: f64[>..] = [];

    // The lambda that makes the data most normal, by two criteria.
    let by_mle, c1 = statistics::boxcox_normmax(skewed[..], statistics::NORMMAX_MLE);
    let by_r, c2 = statistics::boxcox_normmax(skewed[..], statistics::NORMMAX_PEARSONR);
    print("box-cox lambda: mle ", by_mle, " (code ", c1, "), pearsonr ", by_r, " (code ", c2, ")\n");

    // The transform for a given lambda, and with the fitted lambda.
    let l, c3 = statistics::boxcox(out, skewed[..], 0.3);
    print("box-cox, lambda ", l, ", code ", c3, ", first value ", out[0], "\n");
    let fit, c4 = statistics::boxcox_fit(out, skewed[..]);
    print("box-cox fit: lambda ", fit, ", code ", c4, ", last value ", out[7], "\n");

    // The confidence interval of the lambda.
    let lo, hi, mid, c5 = statistics::boxcox_interval(skewed[..], 0.05);
    print("95 percent interval: ", lo, " to ", hi, " around ", mid, ", code ", c5, "\n");

    // The log-likelihood of one lambda.
    let llf, c6 = statistics::boxcox_llf(0.3, skewed[..]);
    print("log-likelihood at 0.3: ", llf, ", code ", c6, "\n");

    // Yeo-Johnson takes values of any sign.
    let yfit, c7 = statistics::yeojohnson_fit(out, mixed[..]);
    print("yeo-johnson fit: lambda ", yfit, ", code ", c7, ", first value ", out[0], "\n");
    let yl, c8 = statistics::yeojohnson_llf(0.4, mixed[..]);
    print("yeo-johnson log-likelihood at 0.4: ", yl, ", code ", c8, "\n");

    // A normal plot table. Each row is a lambda and the correlation r.
    let rows, c9 = statistics::boxcox_normplot(out, skewed[..], -1.0, 1.0, 5);
    print("normal plot, ", rows, " rows, code ", c9, "\n");
    for i in 0..rows { print("  lambda ", out[2 * i], "  r ", out[2 * i + 1], "\n"); }

    // A value at or below 0 is a fault for Box-Cox.
    let bad: f64[3] = [1.0, 0.0, 2.0];
    let _v, c10 = statistics::boxcox(out, bad[..], 0.5);
    print("box-cox of a zero: code ", c10, "\n");
}
