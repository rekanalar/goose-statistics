// The statrs oracle of the statistics package. It reads samples from a
// text file. Each line holds one sample of space separated values. It
// writes one line of space separated results for each sample.
use statrs::distribution::{
    Beta, Binomial, ChiSquared, Continuous, ContinuousCDF, Discrete, DiscreteCDF, Exp,
    FisherSnedecor, Gamma, Geometric, Hypergeometric, LogNormal, NegativeBinomial, Normal,
    Poisson, StudentsT, Uniform, Weibull,
};
use statrs::function::{beta, erf, gamma};
use statrs::statistics::{Accumulate, Data, Distribution, OnlineMoments, OrderStatistics, Statistics};
use statrs::stats_tests::chisquare::chisquare;
use statrs::stats_tests::f_oneway::f_oneway;
use statrs::stats_tests::fisher::fishers_exact_with_odds_ratio;
use statrs::stats_tests::anderson_darling::anderson_darling;
use statrs::stats_tests::skewtest::skewtest;
use statrs::stats_tests::ks_test::{ks_onesample, ks_twosample, KSOneSampleAlternativeMethod, KSTwoSampleAlternativeMethod};
use statrs::stats_tests::mannwhitneyu::{mannwhitneyu, MannWhitneyUMethod};
use statrs::stats_tests::ttest_onesample::ttest_onesample;
use statrs::stats_tests::{Alternative, NaNPolicy};
use std::fs;

// The result of the special function number fn. The crate holds no
// function for the others, and those give NaN.
fn special(f: i64, a: f64, b: f64, x: f64) -> f64 {
    match f {
        0 => erf::erf(x),
        1 => erf::erfc(x),
        2 => erf::erf_inv(x),
        3 => erf::erfc_inv(x),
        4 => gamma::gamma(x),
        5 => gamma::ln_gamma(x),
        7 => gamma::digamma(x),
        8 => beta::beta(a, b),
        9 => beta::ln_beta(a, b),
        11 => gamma::gamma_lr(a, x),
        12 => gamma::gamma_ur(a, x),
        15 => beta::beta_reg(a, b, x),
        16 => beta::inv_beta_reg(a, b, x),
        _ => f64::NAN,
    }
}

// The special function mode. Each line holds the function number and
// three arguments. The crate panics on a bad domain, and a panic gives NaN.
fn run_special(path: &str) {
    std::panic::set_hook(Box::new(|_| {}));
    let text = fs::read_to_string(path).expect("a readable file");
    for line in text.lines() {
        let w: Vec<f64> = line.split_whitespace().map(|v| v.parse::<f64>().expect("a number")).collect();
        if w.len() != 4 {
            continue;
        }
        let f = w[0] as i64;
        let r = std::panic::catch_unwind(|| special(f, w[1], w[2], w[3])).unwrap_or(f64::NAN);
        println!("{}", r);
    }
}

// The result of one operation of one family. The families are the same
// numbers as in the other oracles: 0 norm, 1 unif, 2 expon, 3 lognorm,
// 4 weibull, 5 gammadist, 6 betadist, 7 chi2, 8 tdist, 9 fdist, 10 binom,
// 11 poisson, 12 geom, 13 nbinom, 14 hypergeom. The operations are 0 pdf,
// 1 cdf, 2 sf, 3 ppf, 4 mean, 5 var, 6 skew, 7 kurt.
//
// The crate panics on a bad domain, and a panic gives NaN. The crate
// holds no skewness and no kurtosis for its families, so those give NaN.
fn distribution(fam: i64, op: i64, x: f64, a: f64, b: f64, c: f64) -> f64 {
    macro_rules! cont {
        ($d:expr) => {{
            let d = $d.expect("a valid family");
            match op {
                0 => d.pdf(x),
                1 => d.cdf(x),
                2 => d.sf(x),
                3 => d.inverse_cdf(x),
                4 => d.mean().unwrap_or(f64::NAN),
                5 => d.variance().unwrap_or(f64::NAN),
                _ => f64::NAN,
            }
        }};
    }
    macro_rules! disc {
        ($d:expr) => {{
            let d = $d.expect("a valid family");
            let k = if x >= 0.0 { x as u64 } else { 0 };
            match op {
                0 => d.pmf(k),
                1 => d.cdf(k),
                2 => d.sf(k),
                3 => d.inverse_cdf(x) as f64,
                4 => d.mean().unwrap_or(f64::NAN),
                5 => d.variance().unwrap_or(f64::NAN),
                _ => f64::NAN,
            }
        }};
    }
    // The crate holds no mean and no variance for the negative binomial
    // family, so those answers give NaN.
    macro_rules! disc_nm {
        ($d:expr) => {{
            let d = $d.expect("a valid family");
            let k = if x >= 0.0 { x as u64 } else { 0 };
            match op {
                0 => d.pmf(k),
                1 => d.cdf(k),
                2 => d.sf(k),
                3 => d.inverse_cdf(x) as f64,
                _ => f64::NAN,
            }
        }};
    }
    match fam {
        0 => cont!(Normal::new(a, b)),
        1 => cont!(Uniform::new(a, b)),
        2 => cont!(Exp::new(a)),
        3 => cont!(LogNormal::new(a, b)),
        4 => cont!(Weibull::new(a, b)),
        // the crate holds the rate, and the package holds the scale
        5 => cont!(Gamma::new(a, 1.0 / b)),
        6 => cont!(Beta::new(a, b)),
        7 => cont!(ChiSquared::new(a)),
        8 => cont!(StudentsT::new(0.0, 1.0, a)),
        9 => cont!(FisherSnedecor::new(a, b)),
        // the crate holds the probability first
        10 => disc!(Binomial::new(b, a as u64)),
        11 => disc!(Poisson::new(a)),
        12 => disc!(Geometric::new(a)),
        13 => disc_nm!(NegativeBinomial::new(a, b)),
        14 => disc!(Hypergeometric::new(a as u64, b as u64, c as u64)),
        _ => f64::NAN,
    }
}

// The distribution mode. Each line holds the family number, the operation
// number and up to four arguments.
fn run_distribution(path: &str) {
    std::panic::set_hook(Box::new(|_| {}));
    let text = fs::read_to_string(path).expect("a readable file");
    for line in text.lines() {
        let w: Vec<f64> = line
            .split_whitespace()
            .map(|v| v.parse::<f64>().expect("a number"))
            .collect();
        if w.len() < 2 {
            continue;
        }
        let fam = w[0] as i64;
        let op = w[1] as i64;
        let v = [
            w.get(2).copied().unwrap_or(0.0),
            w.get(3).copied().unwrap_or(0.0),
            w.get(4).copied().unwrap_or(0.0),
            w.get(5).copied().unwrap_or(0.0),
        ];
        let r = std::panic::catch_unwind(|| distribution(fam, op, v[0], v[1], v[2], v[3]))
            .unwrap_or(f64::NAN);
        println!("{}", r);
    }
}

// The hypothesis test mode. Each line holds a test name and its numbers:
// ttest1 mu n x.., chi n counts.., anova k n1 x.. n2 y.., fisher a b c d,
// mwu n1 x.. n2 y.. (asymptotic with the continuity correction),
// ks1 n x.. (standard normal, asymptotic) and ks2 n1 x.. n2 y.. (asymptotic).
// Each result line holds the statistic and the p-value. The crate panics
// or fails on some inputs, and a failure gives NaN.
fn run_tests(path: &str) {
    std::panic::set_hook(Box::new(|_| {}));
    let text = fs::read_to_string(path).expect("a readable file");
    for line in text.lines() {
        let mut it = line.split_whitespace();
        let name = it.next().unwrap_or("").to_string();
        let w: Vec<f64> = it.map(|v| v.parse::<f64>().expect("a number")).collect();
        let r = std::panic::catch_unwind(|| one_test(&name, &w)).unwrap_or((f64::NAN, f64::NAN));
        println!("{} {}", r.0, r.1);
    }
}

fn one_test(name: &str, w: &[f64]) -> (f64, f64) {
    let nan = (f64::NAN, f64::NAN);
    let take = |at: usize| -> (Vec<f64>, usize) {
        let n = w[at] as usize;
        (w[at + 1..at + 1 + n].to_vec(), at + 1 + n)
    };
    match name {
        "ttest1" => {
            let (x, _) = take(1);
            ttest_onesample(x, w[0], Alternative::TwoSided, NaNPolicy::Error).unwrap_or(nan)
        }
        "chi" => {
            let obs: Vec<usize> = w[1..1 + w[0] as usize].iter().map(|v| *v as usize).collect();
            chisquare(&obs, None, None).unwrap_or(nan)
        }
        "anova" => {
            let mut groups = Vec::new();
            let mut at = 1;
            for _ in 0..w[0] as usize {
                let (g, next) = take(at);
                groups.push(g);
                at = next;
            }
            f_oneway(groups, NaNPolicy::Error).unwrap_or(nan)
        }
        "fisher" => {
            let t = [w[0] as u64, w[1] as u64, w[2] as u64, w[3] as u64];
            fishers_exact_with_odds_ratio(&t, Alternative::TwoSided).unwrap_or(nan)
        }
        "mwu" => {
            let (x, at) = take(0);
            let (y, _) = take(at);
            mannwhitneyu(&x, &y, MannWhitneyUMethod::AsymptoticInclContinuityCorrection, Alternative::TwoSided)
                .unwrap_or(nan)
        }
        "ks1" => {
            let (x, _) = take(0);
            let d = Normal::new(0.0, 1.0).unwrap();
            ks_onesample(x, &d, KSOneSampleAlternativeMethod::TwoSidedAsymptotic, NaNPolicy::Error).unwrap_or(nan)
        }
        "ks2" => {
            let (x, at) = take(0);
            let (y, _) = take(at);
            ks_twosample(x, y, KSTwoSampleAlternativeMethod::TwoSidedAsymptotic, NaNPolicy::Error).unwrap_or(nan)
        }
        // The Anderson-Darling test against a normal fitted with the
        // mean and the unbiased scale, the way the SciPy oracle fits.
        "ad" => {
            let (x, _) = take(0);
            let n = x.len() as f64;
            let mean = x.iter().sum::<f64>() / n;
            let sd = (x.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / (n - 1.0)).sqrt();
            let d = Normal::new(mean, sd).unwrap();
            anderson_darling(&x, &d).unwrap_or(nan)
        }
        "skew" => {
            let (x, _) = take(0);
            skewtest(x, Alternative::TwoSided, NaNPolicy::Error).unwrap_or(nan)
        }
        _ => nan,
    }
}

fn main() {
    let path = std::env::args().nth(1).expect("a path to a file of samples");
    if std::env::args().nth(2).as_deref() == Some("special") {
        run_special(&path);
        return;
    }
    if std::env::args().nth(2).as_deref() == Some("tests") {
        run_tests(&path);
        return;
    }
    if std::env::args().nth(2).as_deref() == Some("distribution") {
        run_distribution(&path);
        return;
    }
    let text = fs::read_to_string(&path).expect("a readable file");
    for line in text.lines() {
        // A line with a bar holds two samples. The oracle writes the
        // covariance of the pair.
        if let Some((left, right)) = line.split_once('|') {
            let xs: Vec<f64> = left
                .split_whitespace()
                .map(|w| w.parse::<f64>().expect("a number"))
                .collect();
            let ys: Vec<f64> = right
                .split_whitespace()
                .map(|w| w.parse::<f64>().expect("a number"))
                .collect();
            if xs.is_empty() || ys.is_empty() {
                continue;
            }
            let covariance = xs.iter().covariance(ys.iter());
            println!("{}", covariance);
            continue;
        }
        let vals: Vec<f64> = line
            .split_whitespace()
            .map(|w| w.parse::<f64>().expect("a number"))
            .collect();
        if vals.is_empty() {
            continue;
        }
        let mean = vals.iter().mean();
        let variance = vals.iter().variance();
        let pop_variance = vals.iter().population_variance();
        let std_dev = vals.iter().std_dev();
        let geometric = vals.iter().geometric_mean();
        let harmonic = vals.iter().harmonic_mean();
        let mut om = OnlineMoments::<3>::default();
        om = vals.iter().fold(om, |a, x| a.push(*x));
        let skew = om.skewness().unwrap_or(f64::NAN);
        let mut data = Data::new(vals.clone());
        let median = data.median();
        println!(
            "{} {} {} {} {} {} {} {}",
            mean, variance, pop_variance, std_dev, geometric, harmonic, skew, median
        );
    }
}
