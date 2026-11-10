// The distribution mode of the jStat oracle. Each input line holds the
// family number, the operation number and up to four arguments. Each
// output line holds one result, or the word nan when the library gives
// no number.
//
// The families are the same numbers as in the other oracles: 0 norm,
// 1 unif, 2 expon, 3 lognorm, 4 weibull, 5 gammadist, 6 betadist, 7 chi2,
// 8 tdist, 9 fdist, 10 binom, 11 poisson, 12 geom, 13 nbinom, 14 hypergeom.
// The operations are 0 pdf, 1 cdf, 2 sf, 3 ppf, 4 mean, 5 var, 6 skew,
// 7 kurt.
//
// The library holds no geometric family. It holds no quantile for the
// binomial family and for the Poisson family. It holds no skewness and
// no kurtosis for any family. Those answers give NaN.
import { readFileSync } from "node:fs";
import * as ss from "simple-statistics";
import jStatPackage from "jstat";
const { jStat } = jStatPackage;

const path = process.argv[2];
const text = readFileSync(path, "utf8");

// The distribution mode.
if (process.argv[3] === "distribution") {
    for (const line of text.split("\n")) {
        if (line.trim() === "") { continue; }
        const w = line.trim().split(/\s+/).map(Number);
        const fam = w[0];
        const op = w[1];
        const x = w.length > 2 ? w[2] : 0;
        const a = w.length > 3 ? w[3] : 0;
        const b = w.length > 4 ? w[4] : 0;
        const c = w.length > 5 ? w[5] : 0;
        let r = NaN;
        try {
            r = distribution(fam, op, x, a, b, c);
        } catch (e) { r = NaN; }
        console.log(String(r));
    }
    process.exit(0);
}

// The result of one operation of one family.
function distribution(fam, op, x, a, b, c) {
    switch (fam) {
        case 0: return cont(jStat.normal, op, x, a, b);
        case 1: return cont(jStat.uniform, op, x, a, b);
        case 2: return cont(jStat.exponential, op, x, a);
        case 3: return cont(jStat.lognormal, op, x, a, b);
        // the library takes the scale first and the shape last
        case 4: return cont(jStat.weibull, op, x, b, a);
        case 5: return cont(jStat.gamma, op, x, a, b);
        case 6: return cont(jStat.beta, op, x, a, b);
        case 7: return cont(jStat.chisquare, op, x, a);
        case 8: return cont(jStat.studentt, op, x, a);
        case 9: return cont(jStat.centralF, op, x, a, b);
        case 10: return disc(jStat.binomial, op, x, a, b);
        case 11: return disc(jStat.poisson, op, x, a);
        case 13: return disc(jStat.negbin, op, x, a, b);
        case 14: return disc(jStat.hypgeom, op, x, a, b, c);
    }
    return NaN;
}

// The result of one operation of a continuous family.
function cont(d, op, x, ...args) {
    if (op === 0) { return d.pdf(x, ...args); }
    if (op === 1) { return d.cdf(x, ...args); }
    if (op === 2) { return 1 - d.cdf(x, ...args); }
    if (op === 3) { return d.inv(x, ...args); }
    if (op === 4) { return d.mean(...args); }
    if (op === 5) { return d.variance(...args); }
    return NaN;
}

// The result of one operation of a discrete family.
function disc(d, op, x, ...args) {
    if (op === 0) { return d.pdf(x, ...args); }
    if (op === 1) { return d.cdf(x, ...args); }
    if (op === 2) { return 1 - d.cdf(x, ...args); }
    if (op === 4) { return d.mean ? d.mean(...args) : NaN; }
    if (op === 5) { return d.variance ? d.variance(...args) : NaN; }
    return NaN;
}

// The special function mode. Each line holds the function number and three
// arguments. The package holds the error function, the gamma function and
// its log. The other numbers give NaN.
if (process.argv[3] === "special") {
    for (const line of text.split("\n")) {
        if (line.trim() === "") { continue; }
        const [f, a, b, x] = line.trim().split(/\s+/).map(Number);
        let r = NaN;
        try {
            if (f === 0) { r = ss.errorFunction(x); }
            if (f === 4) { r = ss.gamma(x); }
            if (f === 5) { r = ss.gammaln(x); }
        } catch (e) { r = NaN; }
        console.log(String(r));
    }
    process.exit(0);
}
for (const line of text.split("\n")) {
    if (line.trim() === "") { continue; }
    // A line with a bar holds two samples. The oracle writes the
    // covariance and the correlation of the pair.
    if (line.includes("|")) {
        const [left, right] = line.split("|");
        const xs = left.trim().split(/\s+/).map(Number);
        const ys = right.trim().split(/\s+/).map(Number);
        const covariance = xs.length > 1 ? ss.sampleCovariance(xs, ys) : NaN;
        const correlation = xs.length > 1 ? ss.sampleCorrelation(xs, ys) : NaN;
        console.log([covariance, correlation].join(" "));
        continue;
    }
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
