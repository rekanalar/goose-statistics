package main

import (
	"encoding/json"
	"fmt"
	"os"
	"sort"

	"gonum.org/v1/gonum/stat"
)

// The test mode. The input is a list of cases. A case of kind "chi" holds the
// observed counts in a and the expected counts in b, and gives the chi-square
// statistic. A case of kind "ks" holds two samples and gives the two-sample
// Kolmogorov-Smirnov statistic. The library holds no p-value for either.
type testCase struct {
	Kind string    `json:"kind"`
	A    []float64 `json:"a"`
	B    []float64 `json:"b"`
}

func runTests(path string) {
	raw, err := os.ReadFile(path)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	var cases []testCase
	if err := json.Unmarshal(raw, &cases); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	out := make([]float64, 0, len(cases))
	for _, c := range cases {
		switch c.Kind {
		case "chi":
			out = append(out, stat.ChiSquare(c.A, c.B))
		case "ks":
			x := append([]float64(nil), c.A...)
			y := append([]float64(nil), c.B...)
			sort.Float64s(x)
			sort.Float64s(y)
			out = append(out, stat.KolmogorovSmirnov(x, nil, y, nil))
		default:
			out = append(out, 0)
		}
	}
	b, _ := json.Marshal(out)
	fmt.Println(string(b))
}
