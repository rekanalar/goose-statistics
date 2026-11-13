// The measured rows of the entropy calls: the answers of SciPy 1.18.1.
// The generator test/make_shape_vectors.py writes this file.
namespace suite;

import std;
import math;
import statistics from "../..";

fn run_entropy_vectors() -> i64 {
    var bad: i64 = 0;
    var out: f64[>..] = [];
    var out2: f64[>..] = [];
    let s2: f64[10] = [2.5, -1.0, 4.0, 3.25, 10.5, 0.5, 7.75, 3.0, 8.0, 6.0];
    let de: f64[12] = [0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3];
    let h2: f64[2] = [0.5, 0.5];
    let h3: f64[2] = [0.9, 0.1];
    let h4: f64[2] = [1.0, 3.0];
    let v1, c1 = statistics::entropy(h2[..], 0.0);
    chk_val(bad, "entropy", v1, c1, 0.6931471805599453, 1e-12);
    let v2, c2 = statistics::entropy(h2[..], 2.0);
    chk_val(bad, "entropy base 2", v2, c2, 1.0, 1e-12);
    let v3, c3 = statistics::entropy(h4[..], 0.0);
    chk_val(bad, "entropy unnormalized", v3, c3, 0.5623351446188083, 1e-12);
    let v4, c4 = statistics::relative_entropy(h2[..], h3[..], 0.0);
    chk_val(bad, "relative_entropy", v4, c4, 0.5108256237659906, 1e-12);
    let v5, c5 = statistics::differential_entropy(de[..], 1, 0, 0.0);
    chk_val(bad, "differential_entropy vasicek", v5, c5, 0.6619797631102224, 1e-12);
    let v6, c6 = statistics::differential_entropy(de[..], 1, 3, 0.0);
    chk_val(bad, "differential_entropy vasicek window 3", v6, c6, 0.6619797631102224, 1e-12);
    let v7, c7 = statistics::differential_entropy(de[..], 2, 0, 0.0);
    chk_val(bad, "differential_entropy van es", v7, c7, 1.0923850544446876, 1e-12);
    let v8, c8 = statistics::differential_entropy(de[..], 2, 3, 0.0);
    chk_val(bad, "differential_entropy van es window 3", v8, c8, 1.0923850544446876, 1e-12);
    let v9, c9 = statistics::differential_entropy(de[..], 3, 0, 0.0);
    chk_val(bad, "differential_entropy ebrahimi", v9, c9, 0.8754687373538999, 1e-12);
    let v10, c10 = statistics::differential_entropy(de[..], 3, 3, 0.0);
    chk_val(bad, "differential_entropy ebrahimi window 3", v10, c10, 0.8754687373538999, 1e-12);
    let v11, c11 = statistics::differential_entropy(de[..], 4, 0, 0.0);
    chk_val(bad, "differential_entropy correa", v11, c11, 0.7440341858318904, 1e-12);
    let v12, c12 = statistics::differential_entropy(de[..], 4, 3, 0.0);
    chk_val(bad, "differential_entropy correa window 3", v12, c12, 0.7440341858318904, 1e-12);
    let v13, c13 = statistics::differential_entropy(de[..], 0, 0, 0.0);
    chk_val(bad, "differential_entropy auto", v13, c13, 0.8754687373538999, 1e-12);
    let v14, c14 = statistics::differential_entropy(de[..], 1, 0, 2.0);
    chk_val(bad, "differential_entropy base 2", v14, c14, 0.9550349214079688, 1e-12);
    let v15, c15 = statistics::kstat(s2[..], 1);
    chk_val(bad, "kstat 1", v15, c15, 4.45, 1e-12);
    let v16, c16 = statistics::kstat(s2[..], 2);
    chk_val(bad, "kstat 2", v16, c16, 12.816666666666666, 1e-12);
    let v17, c17 = statistics::kstat(s2[..], 3);
    chk_val(bad, "kstat 3", v17, c17, 9.7296875, 1e-12);
    let v18, c18 = statistics::kstat(s2[..], 4);
    chk_val(bad, "kstat 4", v18, c18, -108.3248511904762, 1e-12);
    let v19, c19 = statistics::kstatvar(s2[..], 1);
    chk_val(bad, "kstatvar 1", v19, c19, 1.2816666666666667, 1e-12);
    let v20, c20 = statistics::kstatvar(s2[..], 2);
    chk_val(bad, "kstatvar 2", v20, c20, 21.0037748015873, 1e-12);
    return bad;
}
