// The gonum oracle of the statistics package. It reads a JSON list of
// samples and writes a JSON list of results from gonum.org/v1/gonum/stat.
// With the argument pairs it reads a JSON list of pairs of samples and
// writes the covariance and the correlation of each pair.
package main

import (
	"encoding/json"
	"fmt"
	"math"
	"os"

	"gonum.org/v1/gonum/stat"
)

type row struct {
	Mean     float64 `json:"mean"`
	Variance float64 `json:"variance"`
	StdDev   float64 `json:"stddev"`
	Skew     float64 `json:"skew"`
	ExKurt   float64 `json:"exkurt"`
}

type pairRow struct {
	Covariance  float64 `json:"covariance"`
	Correlation float64 `json:"correlation"`
}

type pairIn struct {
	X []float64 `json:"x"`
	Y []float64 `json:"y"`
}

// clean maps a NaN or an infinity to 0, which JSON can hold.
func clean(x float64) float64 {
	if math.IsNaN(x) || math.IsInf(x, 0) {
		return 0
	}
	return x
}

func main() {
	if len(os.Args) > 2 && os.Args[2] == "special" {
		runSpecial(os.Args[1])
		return
	}
	raw, err := os.ReadFile(os.Args[1])
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	if len(os.Args) > 2 && os.Args[2] == "pairs" {
		var pairs []pairIn
		if err := json.Unmarshal(raw, &pairs); err != nil {
			fmt.Fprintln(os.Stderr, err)
			os.Exit(1)
		}
		out := make([]pairRow, 0, len(pairs))
		for _, p := range pairs {
			out = append(out, pairRow{
				Covariance:  clean(stat.Covariance(p.X, p.Y, nil)),
				Correlation: clean(stat.Correlation(p.X, p.Y, nil)),
			})
		}
		b, _ := json.Marshal(out)
		fmt.Println(string(b))
		return
	}
	var samples [][]float64
	if err := json.Unmarshal(raw, &samples); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	out := make([]row, 0, len(samples))
	for _, xs := range samples {
		out = append(out, row{
			Mean:     clean(stat.Mean(xs, nil)),
			Variance: clean(stat.Variance(xs, nil)),
			StdDev:   clean(stat.StdDev(xs, nil)),
			Skew:     clean(stat.Skew(xs, nil)),
			ExKurt:   clean(stat.ExKurtosis(xs, nil)),
		})
	}
	b, _ := json.Marshal(out)
	fmt.Println(string(b))
}
