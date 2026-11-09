// The simple-statistics oracle of the statistics package. It reads samples
// from a text file. Each line holds one sample of space separated values.
// It writes one line of space separated results for each sample.
import { readFileSync } from "node:fs";
import * as ss from "simple-statistics";

const path = process.argv[2];
const text = readFileSync(path, "utf8");
for (const line of text.split("\n")) {
    if (line.trim() === "") { continue; }
    const vals = line.trim().split(/\s+/).map(Number);
    const mean = ss.mean(vals);
    const variance = ss.sampleVariance(vals);
    const popVariance = ss.variance(vals);
    const stdDev = Math.sqrt(variance);
    // The package simple-statistics throws on a negative value in the
    // geometric mean. The crate statrs gives NaN. This oracle gives NaN.
    let geometric = NaN;
    let harmonic = NaN;
    // The package simple-statistics has two forms. The form geometricMean
    // multiplies the values, so it overflows on a large sample. The form
    // logAverage uses the logarithms, so it holds. This oracle uses the
    // second form.
    if (vals.every((x) => x >= 0)) { geometric = ss.logAverage(vals); }
    if (vals.every((x) => x > 0)) { harmonic = ss.harmonicMean(vals); }
    // The skewness of simple-statistics is the adjusted form, the one with
    // the factor of the sample size. The crate statrs gives the plain form.
    const skew = vals.length > 2 ? ss.sampleSkewness(vals) : NaN;
    const median = ss.median(vals);
    console.log([mean, variance, popVariance, stdDev, geometric, harmonic, skew, median].join(" "));
}
