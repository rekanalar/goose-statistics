// The statrs oracle of the statistics package. It reads samples from a
// text file. Each line holds one sample of space separated values. It
// writes one line of space separated results for each sample.
use statrs::statistics::{Accumulate, Data, OnlineMoments, OrderStatistics, Statistics};
use std::fs;

fn main() {
    let path = std::env::args().nth(1).expect("a path to a file of samples");
    let text = fs::read_to_string(&path).expect("a readable file");
    for line in text.lines() {
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
