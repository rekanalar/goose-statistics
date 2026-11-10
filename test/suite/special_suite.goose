// Tests of the special functions against the reference vectors. Each row
// holds the function number, up to three arguments, the exact value and the
// tolerance in units in the last place of the value.
namespace suite;

import std;
import statistics from "../..";
import special_vectors;

// The distance of a result from the exact value, in units in the last place
// of the exact value. A result that is a fault gives a large number.
fn ulp_error(got: f64, want: f64) -> f64 {
    if got == want { return 0.0; }
    if is_nan(got) { return 1e300; }
    var sp = abs(want) * 2.220446049250313e-16 * 0.5;
    if abs(want) < 2.2250738585072014e-308 { sp = 5e-324; }
    // the spacing of the value is the power of two at or below it, times 2^-52
    var m = abs(want);
    var p: f64 = 1.0;
    if m >= 1.0 {
        while p * 2.0 <= m { p = p * 2.0; }
    } else if m >= 2.2250738585072014e-308 {
        while p > m { p = p * 0.5; }
    } else {
        p = 2.2250738585072014e-308;
    }
    sp = p * 2.220446049250313e-16;
    if abs(want) < 2.2250738585072014e-308 { sp = 5e-324; }
    return abs(got - want) / sp;
}

// The names of the calls, in the order of the function numbers.
fn call_name(fn_no: i64) -> u8[:] {
    if fn_no == 0 { return "erf"; }
    if fn_no == 1 { return "erfc"; }
    if fn_no == 2 { return "erfinv"; }
    if fn_no == 3 { return "erfcinv"; }
    if fn_no == 4 { return "gamma"; }
    if fn_no == 5 { return "lgamma"; }
    if fn_no == 6 { return "lgamma sign"; }
    if fn_no == 7 { return "digamma"; }
    if fn_no == 8 { return "beta"; }
    if fn_no == 9 { return "lbeta"; }
    if fn_no == 10 { return "lchoose"; }
    if fn_no == 11 { return "gammainc"; }
    if fn_no == 12 { return "gammaincc"; }
    if fn_no == 13 { return "gammaincinv"; }
    if fn_no == 14 { return "gammainccinv"; }
    if fn_no == 15 { return "betainc"; }
    if fn_no == 16 { return "betaincinv"; }
    if fn_no == 17 { return "log1p"; }
    return "expm1";
}

// The result of the call of one row and its fault code.
fn call_row(fn_no: i64, a: f64, b: f64, x: f64) -> f64, i64 {
    if fn_no == 0 { return statistics::erf(x); }
    if fn_no == 1 { return statistics::erfc(x); }
    if fn_no == 2 { return statistics::erfinv(x); }
    if fn_no == 3 { return statistics::erfcinv(x); }
    if fn_no == 4 { return statistics::gamma(x); }
    if fn_no == 5 {
        let v, s, c = statistics::lgamma(x);
        return v, c;
    }
    if fn_no == 6 {
        let v, s, c = statistics::lgamma(x);
        return s, c;
    }
    if fn_no == 7 { return statistics::digamma(x); }
    if fn_no == 8 { return statistics::beta(a, b); }
    if fn_no == 9 { return statistics::lbeta(a, b); }
    if fn_no == 10 { return statistics::lchoose(a, b); }
    if fn_no == 11 { return statistics::gammainc(a, x); }
    if fn_no == 12 { return statistics::gammaincc(a, x); }
    if fn_no == 13 { return statistics::gammaincinv(a, x); }
    if fn_no == 14 { return statistics::gammainccinv(a, x); }
    if fn_no == 15 { return statistics::betainc(a, b, x); }
    if fn_no == 16 { return statistics::betaincinv(a, b, x); }
    if fn_no == 17 { return statistics::log1p(x); }
    return statistics::expm1(x);
}

// True when the result is within the absolute floor of a value near a zero
// of the log gamma or the digamma, where the relative error has no meaning.
fn near_zero_ok(fn_no: i64, got: f64, want: f64) -> bool {
    if fn_no != 5 && fn_no != 7 && fn_no != 9 && fn_no != 10 { return false; }
    return abs(got - want) <= 2e-15;
}

// The round trip of an inverse row. The forward call must return the input
// probability within the tolerance, when the shape is small enough for the
// round trip to be defined within it.
fn round_trip_ulp(fn_no: i64, a: f64, b: f64, p: f64, root: f64) -> f64, i64 {
    var v: f64 = 0.0;
    var c: i64 = 0;
    if fn_no == 13 { v, c = statistics::gammainc(a, root); }
    if fn_no == 14 { v, c = statistics::gammaincc(a, root); }
    if fn_no == 16 { v, c = statistics::betainc(a, b, root); }
    return ulp_error(v, p), c;
}

// Run all the rows. With report set, print the largest error of each call
// and return 0.
fn run_special(report: bool) -> i64 {
    var bad: i64 = 0;
    var worst: f64[19];
    var worst_row: i64[19];
    var fails: i64[19];
    for k in 0..19 {
        worst[k] = 0.0;
        worst_row[k] = 0;
        fails[k] = 0;
    }
    for r in 0..S_ROWS {
        let fn_no = S_FN[r];
        let a = S_A[r];
        let b = S_B[r];
        let x = S_X[r];
        let e = S_E[r];
        let tol = S_TOL[r];
        let got, code = call_row(fn_no, a, b, x);
        var err: f64 = 0.0;
        var ok = true;
        if fn_no == 13 || fn_no == 14 || fn_no == 16 {
            // x is the probability, e is the exact root
            var small = a <= 10.0;
            if fn_no == 16 { small = a <= 100.0 && b <= 100.0; }
            if code != 0 {
                ok = false;
                err = 1e300;
            } else if small {
                let rt, c2 = round_trip_ulp(fn_no, a, b, x, got);
                err = rt;
                var lim = 20.0;
                if fn_no == 16 { lim = 16.0; }
                // a root within 2 ULP is the best answer; the round trip then holds only to the conditioning
                ok = c2 == 0 && (rt <= lim || ulp_error(got, e) <= 2.0);
            } else {
                err = ulp_error(got, e);
                ok = err <= 64.0;
            }
        } else if fn_no == 6 {
            err = abs(got - e);
            ok = code == 0 && got == e;
        } else {
            err = ulp_error(got, e);
            ok = code == 0 && (err <= tol || near_zero_ok(fn_no, got, e));
        }
        if err > worst[fn_no] && !near_zero_ok(fn_no, got, e) {
            worst[fn_no] = err;
            worst_row[fn_no] = r;
        }
        if !ok {
            bad += 1;
            fails[fn_no] += 1;
            if (bad <= 30 && !report) || (report && fails[fn_no] <= 12) {
                print("FAIL ", call_name(fn_no), " row ", r, " a ", a, " b ", b, " x ", x, " got ", got, " want ", e, " code ", code, "\n");
            }
        }
    }
    if report {
        for k in 0..19 {
            print(k, " ", call_name(k), " worst ", worst[k], " row ", worst_row[k], " fails ", fails[k], "\n");
        }
        print("bad ", bad, "\n");
        return 0;
    }
    return bad;
}
