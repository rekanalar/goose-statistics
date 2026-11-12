// The answers of SciPy 1.18.1 for the logrank test and the
// goodness-of-fit test. The generator test/make_fitting_vectors.py
// writes this file. It holds the numbers and the checks that use them.
namespace suite;

import std;
import math;
import statistics from "../..";

fn lrow(bad: i64&, ux: const f64[:], rx: const f64[:], uy: const f64[:], ry: const f64[:], alt: i64, stat: f64, p: f64, name: u8[:]) {
    var none: f64[>..] = [];
    let x = statistics::CensoredData { uncensored: ux, left: none[..], right: rx, interval: none[..] };
    let y = statistics::CensoredData { uncensored: uy, left: none[..], right: ry, interval: none[..] };
    let r, c = statistics::logrank(x, y, alt);
    if c == 0 && near_tol(r.statistic, stat, 1e-12) && near_tol(r.pvalue, p, 1e-12) { return; }
    bad += 1;
    print("FAIL ", name, " got ", r.statistic, " ", r.pvalue, " code ", c, " want ", stat, " ", p, "\n");
}

// The statistic matches at 1e-12. The p-value is a Monte Carlo value: it
// must lie within 5 combined standard errors of the SciPy value.
fn grow(bad: i64&, data: const f64[:], family: i64, stat: i64, want_stat: f64, want_p: f64, name: u8[:]) {
    var out: f64[>..] = [];
    var g = rng { s: 12345 };
    let r, c = statistics::goodness_of_fit(out, data, family, stat, 3000, g);
    let se = sqrt(want_p * (1.0 - want_p) * (1.0 / 3000.0 + 1.0 / 20000.0)) + 1e-4;
    if c == 0 && near_tol(r.statistic, want_stat, 1e-12) && abs(r.pvalue - want_p) <= 5.0 * se { return; }
    bad += 1;
    print("FAIL ", name, " got ", r.statistic, " ", r.pvalue, " code ", c, " want ", want_stat, " ", want_p, "\n");
}

fn run_fitting_vectors() -> i64 {
    var bad: i64 = 0;
    lrow(bad, [2.0, 3.0, 5.0, 8.0], [4.0, 6.0, 9.0], [1.0, 2.0, 2.0, 4.0], [3.0, 5.0, 7.0], statistics::ALT_TWO_SIDED, -0.8042527331835678, 0.42125102949556004, "logrank, the measured pair");
    lrow(bad, [2.0, 3.0, 5.0, 8.0], [4.0, 6.0, 9.0], [1.0, 2.0, 2.0, 4.0], [3.0, 5.0, 7.0], statistics::ALT_LESS, -0.8042527331835678, 0.21062551474778002, "logrank, less");
    lrow(bad, [2.0, 3.0, 5.0, 8.0], [4.0, 6.0, 9.0], [1.0, 2.0, 2.0, 4.0], [3.0, 5.0, 7.0], statistics::ALT_GREATER, -0.8042527331835678, 0.78937448525222, "logrank, greater");
    lrow(bad, [6.0, 13.0, 21.0, 30.0, 37.0, 38.0, 49.0, 50.0, 63.0, 79.0, 86.0, 98.0, 202.0, 219.0], [31.0, 47.0, 80.0, 82.0, 82.0, 149.0], [10.0, 10.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 20.0, 24.0, 24.0, 25.0, 28.0, 30.0, 33.0, 35.0, 37.0, 40.0, 40.0, 46.0, 48.0, 76.0, 81.0, 82.0, 91.0, 112.0, 181.0], [34.0, 40.0, 70.0], statistics::ALT_TWO_SIDED, -2.7379909000095908, 0.006181578637461755, "logrank, the glioma data");
    lrow(bad, [10.0, 10.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 20.0, 24.0, 24.0, 25.0, 28.0, 30.0, 33.0, 35.0, 37.0, 40.0, 40.0, 46.0, 48.0, 76.0, 81.0, 82.0, 91.0, 112.0, 181.0], [34.0, 40.0, 70.0], [6.0, 13.0, 21.0, 30.0, 37.0, 38.0, 49.0, 50.0, 63.0, 79.0, 86.0, 98.0, 202.0, 219.0], [31.0, 47.0, 80.0, 82.0, 82.0, 149.0], statistics::ALT_TWO_SIDED, 2.7379909000095908, 0.006181578637461755, "logrank, swapped groups");
    lrow(bad, [1.5, 2.5, 3.5, 6.0], [], [2.0, 4.5, 7.0, 8.0, 9.5], [], statistics::ALT_TWO_SIDED, 1.7293443254660494, 0.08374748793695098, "logrank, no censoring");
    lrow(bad, [3.0, 5.0, 5.0], [5.0, 8.0], [2.0, 5.0, 9.0], [5.0], statistics::ALT_GREATER, 0.1067128754403218, 0.4575083839063884, "logrank, a censored time equals an event");
    grow(bad, [8.414, 10.481, 6.207, 12.792, 11.277, 9.416, 9.376, 10.608, 9.465, 9.548, 11.44, 11.029, 9.872, 9.829, 10.322, 8.772, 9.192, 11.097, 9.739, 7.251, 9.045, 11.313, 9.535, 9.703, 11.284, 13.649, 8.574, 12.696, 7.54, 10.35], statistics::FAMILY_NORMAL, statistics::GOF_ANDERSON, 0.30738800475545247, 0.5716714164291785, "gof, normal, anderson");
    grow(bad, [8.414, 10.481, 6.207, 12.792, 11.277, 9.416, 9.376, 10.608, 9.465, 9.548, 11.44, 11.029, 9.872, 9.829, 10.322, 8.772, 9.192, 11.097, 9.739, 7.251, 9.045, 11.313, 9.535, 9.703, 11.284, 13.649, 8.574, 12.696, 7.54, 10.35], statistics::FAMILY_NORMAL, statistics::GOF_KOLMOGOROV, 0.09670991182846067, 0.6673666316684166, "gof, normal, kolmogorov");
    grow(bad, [8.414, 10.481, 6.207, 12.792, 11.277, 9.416, 9.376, 10.608, 9.465, 9.548, 11.44, 11.029, 9.872, 9.829, 10.322, 8.772, 9.192, 11.097, 9.739, 7.251, 9.045, 11.313, 9.535, 9.703, 11.284, 13.649, 8.574, 12.696, 7.54, 10.35], statistics::FAMILY_NORMAL, statistics::GOF_CRAMER, 0.05003223203166938, 0.5127743612819359, "gof, normal, cramer");
    grow(bad, [1.983, 1.867, 1.169, 2.281, 0.715, 3.746, 5.335, 0.322, 4.729, 0.41, 0.23, 0.27, 0.643, 5.499, 1.765, 3.009, 3.249, 1.132, 0.818, 2.792, 1.569, 2.656, 4.608, 3.103, 0.963, 6.161, 0.463, 0.61, 0.009, 6.002], statistics::FAMILY_NORMAL, statistics::GOF_ANDERSON, 1.0291666191046112, 0.0087995600219989, "gof, skewed data against normal");
    grow(bad, [1.983, 1.867, 1.169, 2.281, 0.715, 3.746, 5.335, 0.322, 4.729, 0.41, 0.23, 0.27, 0.643, 5.499, 1.765, 3.009, 3.249, 1.132, 0.818, 2.792, 1.569, 2.656, 4.608, 3.103, 0.963, 6.161, 0.463, 0.61, 0.009, 6.002], statistics::FAMILY_EXPONENTIAL, statistics::GOF_ANDERSON, 0.3040293480124241, 0.8188090595470227, "gof, exponential, anderson");
    grow(bad, [1.983, 1.867, 1.169, 2.281, 0.715, 3.746, 5.335, 0.322, 4.729, 0.41, 0.23, 0.27, 0.643, 5.499, 1.765, 3.009, 3.249, 1.132, 0.818, 2.792, 1.569, 2.656, 4.608, 3.103, 0.963, 6.161, 0.463, 0.61, 0.009, 6.002], statistics::FAMILY_EXPONENTIAL, statistics::GOF_CRAMER, 0.03815204492471706, 0.8653067346632668, "gof, exponential, cramer");
    grow(bad, [1.955, 3.571, 1.181, 0.913, 0.903, 0.922, 0.699, 0.953, 3.582, 1.155, 1.924, 0.794, 1.825, 0.58, 1.084, 6.377, 1.162, 3.228, 2.166, 1.504, 1.115, 3.569, 1.482, 4.122, 1.071, 1.706, 2.18, 2.062, 0.786, 1.107], statistics::FAMILY_LOGNORMAL, statistics::GOF_ANDERSON, 0.6257358593646138, 0.10149492525373731, "gof, lognormal, anderson");
    grow(bad, [1.955, 3.571, 1.181, 0.913, 0.903, 0.922, 0.699, 0.953, 3.582, 1.155, 1.924, 0.794, 1.825, 0.58, 1.084, 6.377, 1.162, 3.228, 2.166, 1.504, 1.115, 3.569, 1.482, 4.122, 1.071, 1.706, 2.18, 2.062, 0.786, 1.107], statistics::FAMILY_LOGNORMAL, statistics::GOF_KOLMOGOROV, 0.17367910086924, 0.021248937553122346, "gof, lognormal, kolmogorov");
    grow(bad, [2.917, 2.979, 2.638, 3.752, 1.7, 2.798, 3.443, 1.404, 2.847, 2.209, 3.298, 1.203, 2.82, 3.571, 2.885, 1.964, 2.967, 2.014, 3.006, 1.385, 1.84, 1.098, 1.252, 2.664, 1.695, 2.548, 2.979, 3.644, 2.065, 1.956], statistics::FAMILY_UNIFORM, statistics::GOF_KOLMOGOROV, 0.11441848781713149, 0.7950102494875256, "gof, uniform, kolmogorov");
    grow(bad, [2.917, 2.979, 2.638, 3.752, 1.7, 2.798, 3.443, 1.404, 2.847, 2.209, 3.298, 1.203, 2.82, 3.571, 2.885, 1.964, 2.967, 2.014, 3.006, 1.385, 1.84, 1.098, 1.252, 2.664, 1.695, 2.548, 2.979, 3.644, 2.065, 1.956], statistics::FAMILY_UNIFORM, statistics::GOF_CRAMER, 0.06677236844869948, 0.7915604219789011, "gof, uniform, cramer");
    return bad;
}
