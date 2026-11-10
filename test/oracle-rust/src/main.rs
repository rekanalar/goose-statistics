// The statrs oracle of the statistics package. It reads samples from a
// text file. Each line holds one sample of space separated values. It
// writes one line of space separated results for each sample.
use statrs::function::{beta, erf, gamma};
use statrs::statistics::{Accumulate, Data, OnlineMoments, OrderStatistics, Statistics};
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

fn main() {
    let path = std::env::args().nth(1).expect("a path to a file of samples");
    if std::env::args().nth(2).as_deref() == Some("special") {
        run_special(&path);
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
