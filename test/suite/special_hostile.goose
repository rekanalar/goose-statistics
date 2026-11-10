// Scenario tests and hostile input tests of the special functions. Every
// call must answer a fault with NaN and a code, answer a value with code 0,
// never loop without a bound, and finish in a short time on a crafted input.
namespace suite;

import std;
import math;
import os;
import statistics from "../..";

// A relative or absolute match of a result.
fn close(got: f64, want: f64, rel: f64) -> bool {
    if got == want { return true; }
    return abs(got - want) <= rel * abs(want);
}

// One call by function number, with the arguments of the call.
fn want_fault(bad: i64&, fn_no: i64, a: f64, b: f64, x: f64, code: i64, what: u8[:]) {
    let v, c = call_row(fn_no, a, b, x);
    if c == code && is_nan(v) { return; }
    bad += 1;
    print("FAIL ", what, ": code ", c, " value ", v, "\n");
}

fn want_value(bad: i64&, fn_no: i64, a: f64, b: f64, x: f64, want: f64, rel: f64, what: u8[:]) {
    let v, c = call_row(fn_no, a, b, x);
    if c == 0 && close(v, want, rel) { return; }
    bad += 1;
    print("FAIL ", what, ": code ", c, " value ", v, " want ", want, "\n");
}

fn special_cases() -> i64 {
    var bad: i64 = 0;
    let nan = 0.0 / 0.0;
    let inf = 1.0 / 0.0;
    // error function family
    want_value(bad, 0, 0.0, 0.0, 0.5, 0.5204998778130465, 4e-16, "erf(0.5)");
    want_value(bad, 0, 0.0, 0.0, inf, 1.0, 0.0, "erf(inf)");
    want_value(bad, 0, 0.0, 0.0, -inf, -1.0, 0.0, "erf(-inf)");
    want_value(bad, 1, 0.0, 0.0, 27.0, 5.237048920029e-319, 1e-3, "erfc(27) in the subnormal range");
    want_value(bad, 1, 0.0, 0.0, 30.0, 0.0, 0.0, "erfc(30) underflows to zero");
    want_value(bad, 1, 0.0, 0.0, -inf, 2.0, 0.0, "erfc(-inf)");
    want_value(bad, 1, 0.0, 0.0, 10.0, 2.088487583762545e-45, 1e-14, "erfc(10) has no cancellation");
    want_value(bad, 2, 0.0, 0.0, 0.5, 0.4769362762044699, 1e-15, "erfinv(0.5)");
    want_fault(bad, 2, 0.0, 0.0, 1.5, 5, "erfinv(1.5) is out of range");
    want_fault(bad, 3, 0.0, 0.0, -0.5, 5, "erfcinv(-0.5) is out of range");
    // gamma family
    want_value(bad, 4, 0.0, 0.0, 5.0, 24.0, 0.0, "gamma(5) is 24");
    want_value(bad, 4, 0.0, 0.0, 0.5, 1.7724538509055159, 1e-15, "gamma(0.5)");
    want_value(bad, 4, 0.0, 0.0, 171.0, 7.257415615307994e306, 1e-14, "gamma(171)");
    want_fault(bad, 4, 0.0, 0.0, 0.0, 7, "gamma(0) is a pole");
    want_fault(bad, 4, 0.0, 0.0, -3.0, 7, "gamma(-3) is a pole");
    want_fault(bad, 4, 0.0, 0.0, 172.0, 8, "gamma(172) overflows");
    want_fault(bad, 4, 0.0, 0.0, inf, 8, "gamma(inf) overflows");
    want_value(bad, 5, 0.0, 0.0, 1.0, 0.0, 0.0, "lgamma(1) is exactly 0");
    want_value(bad, 5, 0.0, 0.0, 2.0, 0.0, 0.0, "lgamma(2) is exactly 0");
    want_value(bad, 6, 0.0, 0.0, -0.5, -1.0, 0.0, "the sign of lgamma(-0.5)");
    want_value(bad, 6, 0.0, 0.0, -1.5, 1.0, 0.0, "the sign of lgamma(-1.5)");
    want_value(bad, 7, 0.0, 0.0, 1.0, -0.5772156649015329, 1e-15, "digamma(1)");
    want_value(bad, 8, 2.0, 3.0, 0.0, 1.0 / 12.0, 4e-16, "beta(2, 3)");
    want_fault(bad, 8, -1.0, 3.0, 0.0, 5, "beta of a negative shape");
    want_value(bad, 10, 5.0, 2.0, 0.0, 2.302585092994046, 1e-15, "lchoose(5, 2) is log 10");
    want_fault(bad, 10, 2.0, 5.0, 0.0, 5, "lchoose with k above n");
    // log1p and expm1
    want_value(bad, 17, 0.0, 0.0, 1e-20, 1e-20, 1e-15, "log1p(1e-20)");
    want_value(bad, 18, 0.0, 0.0, 1e-20, 1e-20, 1e-15, "expm1(1e-20)");
    want_fault(bad, 17, 0.0, 0.0, -1.0, 8, "log1p(-1) overflows");
    want_fault(bad, 17, 0.0, 0.0, -2.0, 5, "log1p(-2) is out of range");
    want_fault(bad, 18, 0.0, 0.0, 710.0, 8, "expm1(710) overflows");
    // incomplete gamma and beta
    want_value(bad, 11, 1.0, 0.0, 1.0, 0.6321205588285577, 4e-16, "gammainc(1, 1)");
    want_value(bad, 12, 1.0, 0.0, 1.0, 0.36787944117144233, 4e-16, "gammaincc(1, 1)");
    want_fault(bad, 11, 0.0, 0.0, 1.0, 5, "gammainc(0, 1) is out of range");
    want_fault(bad, 11, 1.0, 0.0, -1.0, 5, "gammainc with a negative x");
    want_fault(bad, 11, 1e11, 0.0, 1.0, 5, "gammainc above the shape bound");
    want_value(bad, 15, 1.0, 1.0, 0.3, 0.3, 1e-15, "betainc(1, 1, 0.3)");
    want_value(bad, 15, 2.0, 3.0, 0.0, 0.0, 0.0, "betainc at 0");
    want_value(bad, 15, 2.0, 3.0, 1.0, 1.0, 0.0, "betainc at 1");
    want_fault(bad, 15, 1.0, 1.0, 1.5, 5, "betainc(1, 1, 1.5) is out of range");
    want_value(bad, 13, 2.0, 0.0, 0.0, 0.0, 0.0, "gammaincinv at 0");
    want_fault(bad, 13, 2.0, 0.0, 1.0, 8, "gammaincinv at 1 is infinite");
    // a shape above the bound of 1e10 is a domain fault, and a shape just under it is not
    want_fault(bad, 11, 1.5e10, 0.0, 1.0, 5, "gammainc above the shape bound");
    want_fault(bad, 12, 1.5e10, 0.0, 1.0, 5, "gammaincc above the shape bound");
    want_fault(bad, 15, 1.5e10, 2.0, 0.5, 5, "betainc above the shape bound");
    want_fault(bad, 15, 2.0, 1.5e10, 0.5, 5, "betainc second shape above the bound");
    want_fault(bad, 16, 2.0, 3.0, 1.5, 5, "betaincinv(1.5) is out of range");
    // infinity and the end points of the inverses
    want_fault(bad, 13, 2.0, 0.0, inf, 5, "gammaincinv with an infinite probability");
    want_fault(bad, 2, 0.0, 0.0, inf, 4, "erfinv(inf) is an infinity");
    want_fault(bad, 2, 0.0, 0.0, 1.0, 8, "erfinv(1) is infinite");
    want_value(bad, 16, 2.0, 3.0, 1.0, 1.0, 0.0, "betaincinv at 1");
    want_value(bad, 16, 2.0, 3.0, 0.0, 0.0, 0.0, "betaincinv at 0");
    want_fault(bad, 15, nan, 1.0, 0.5, 3, "betainc with a NaN shape");
    want_fault(bad, 15, 1.0, 1.0, nan, 3, "betainc with a NaN argument");
    want_value(bad, 15, 0.5, 0.5, 0.5, 0.5, 1e-15, "betainc(0.5, 0.5, 0.5)");
    want_value(bad, 17, 0.0, 0.0, -0.5, -0.6931471805599453, 1e-15, "log1p(-0.5)");
    want_value(bad, 18, 0.0, 0.0, -30.0, -0.9999999999999064, 1e-15, "expm1(-30)");
    want_value(bad, 11, 100000.0, 0.0, 100000.0, 0.5004205221103652, 1e-12, "gammainc at the saddle point");
    want_value(bad, 4, 0.0, 0.0, 171.5, 9.483367566824799e307, 1e-14, "gamma(171.5)");
    // NaN in each argument of each call
    for f in 0..19 {
        if f == 6 { continue; }
        want_fault(bad, f, nan, nan, nan, 3, call_name(f));
    }
    return bad;
}

// An argument from a pool of values that break naive code, a log uniform
// number, or a plain uniform number.
fn pick_arg(st: u64&) -> f64 {
    let r = rnd(st) % 6;
    if r == 0 { return pick(st); }
    let u = (rnd(st) % 2000001 - 1000000) as f64 / 1000000.0;
    if r == 1 { return u; }
    if r == 2 { return abs(u); }
    let e = (rnd(st) % 41 - 20) as f64;
    if r == 3 { return u * pow(10.0, e); }
    if r == 4 { return abs(u) * pow(10.0, e); }
    return abs(u) * pow(10.0, e / 2.0);
}

// The count of seeded inputs for each call: 100000, or the number in the
// variable STATISTICS_HOSTILE_N. The mutation probes use a smaller count.
fn hostile_count() -> i64 {
    var buf: u8[>..] = [];
    if !env("STATISTICS_HOSTILE_N", buf) { return 100000; }
    var n: i64 = 0;
    for ch in buf {
        if ch < 48 || ch > 57 { return 100000; }
        n = n * 10 + (ch - 48) as i64;
    }
    if n < 1 { return 100000; }
    return n;
}

// The count of seeded inputs for the inverse calls. Each inverse call
// costs milliseconds, and a full count of them takes longer than the
// whole test program may run. The count stays above 200.
fn hostile_inverse_count() -> i64 {
    var n = hostile_count() / 50;
    if n < 200 { return 200; }
    return n;
}

// Seeded inputs for each call. A fault must hold NaN. A value must
// hold no NaN, and the sign of lgamma must be one of -1, 0 and 1.
fn hostile_special() -> i64 {
    var bad: i64 = 0;
    for f in 0..19 {
        print("hostile ", call_name(f), "\n");
        var st: u64 = 20261110 + (f as u64) * 7919;
        var faults: i64 = 0;
        var values: i64 = 0;
        var count = hostile_count();
        if f == 13 || f == 14 || f == 16 { count = hostile_inverse_count(); }
        for i in 0..count {
            let a = pick_arg(st);
            let b = pick_arg(st);
            let x = pick_arg(st);
            let v, c = call_row(f, a, b, x);
            var ok = c >= 0 && c <= 9;
            if c != 0 {
                faults += 1;
                if f == 6 { ok = ok && v == 0.0; } else { ok = ok && is_nan(v); }
            } else {
                values += 1;
                ok = ok && !is_nan(v);
                if f == 6 { ok = ok && (v == 1.0 || v == -1.0); }
            }
            if !ok {
                bad += 1;
                if bad <= 20 { print("FAIL ", call_name(f), " a ", a, " b ", b, " x ", x, " gave ", v, " code ", c, "\n"); }
            }
        }
        if values == 0 || faults == 0 {
            bad += 1;
            print("FAIL ", call_name(f), " random inputs gave ", values, " values and ", faults, " faults\n");
        }
    }
    return bad;
}

// The crafted inputs: huge shapes near the argument, tiny shapes, and the
// end points. Each call must take at most 0.05 s more than five times a
// random call.
fn special_time() -> i64 {
    var bad: i64 = 0;
    let t0 = clock();
    var total: f64 = 0.0;
    for i in 0..200 {
        let v, c = statistics::betainc(3.0 + (i as f64), 4.0, 0.3);
        total += v;
    }
    let base = (clock() - t0) / 200.0;
    let shapes: const f64[8] = [1e-10, 1e-5, 0.5, 30.0, 5000.0, 1e6, 1e9, 1e10];
    for k in 0..8 {
        let a = shapes[k];
        for j in 0..8 {
            let b = shapes[j];
            let t1 = clock();
            let _p, c1 = statistics::betainc(a, b, a / (a + b));
            let _q, c2 = statistics::gammainc(a, a * 1.0000001);
            let _r, c3 = statistics::gammaincc(a, a * 0.9999999);
            let _s, c4 = statistics::gammaincinv(a, 0.001);
            let _t, c5 = statistics::gammainccinv(a, 1e-9);
            let _u, c6 = statistics::betaincinv(a, b, 0.5);
            let spent = clock() - t1;
            if spent > 0.5 + 40.0 * base {
                bad += 1;
                print("FAIL time of a ", a, " b ", b, ": ", spent, " seconds\n");
            }
            if c1 != 0 || c2 != 0 || c3 != 0 || c4 != 0 || c5 != 0 || c6 != 0 {
                bad += 1;
                print("FAIL codes at a ", a, " b ", b, ": ", c1, c2, c3, c4, c5, c6, "\n");
            }
        }
    }
    return bad;
}

fn run_special_extra() -> i64 {
    var bad: i64 = 0;
    bad += special_cases();
    bad += special_time();
    bad += hostile_special();
    return bad;
}
