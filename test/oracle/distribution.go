// The distribution mode of the gonum oracle. Each input line holds the
// family number, the operation number and up to four arguments. Each
// output line holds one result, or the word nan when the library gives
// no number.
//
// The families: 0 norm, 1 unif, 2 expon, 3 lognorm, 4 weibull, 5 gammadist,
// 6 betadist, 7 chi2, 8 tdist, 9 fdist, 10 binom, 11 poisson, 12 geom,
// 13 nbinom, 14 hypergeom. The operations: 0 pdf, 1 cdf, 2 sf, 3 ppf,
// 4 mean, 5 var, 6 skew, 7 kurt.
//
// The library holds no geometric family, no negative binomial family and
// no hypergeometric family. It holds no quantile for the binomial family
// and for the Poisson family. Those answers give NaN.
package main

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"strconv"
	"strings"

	"gonum.org/v1/gonum/stat/distuv"
)

// The result of one operation of one family. The arguments are the value
// or the probability, and up to three parameters.
func distribution(fam int, op int, x float64, a float64, b float64, c float64) float64 {
	switch fam {
	case 0:
		d := distuv.Normal{Mu: a, Sigma: b}
		return continuous(d, op, x)
	case 1:
		d := distuv.Uniform{Min: a, Max: b}
		return continuous(d, op, x)
	case 2:
		d := distuv.Exponential{Rate: a}
		return continuous(d, op, x)
	case 3:
		d := distuv.LogNormal{Mu: a, Sigma: b}
		return continuous(d, op, x)
	case 4:
		d := distuv.Weibull{K: a, Lambda: b}
		return continuous(d, op, x)
	case 5:
		// The library holds the rate, and the package holds the scale.
		d := distuv.Gamma{Alpha: a, Beta: 1 / b}
		return continuous(d, op, x)
	case 6:
		d := distuv.Beta{Alpha: a, Beta: b}
		return continuous(d, op, x)
	case 7:
		d := distuv.ChiSquared{K: a}
		return continuous(d, op, x)
	case 8:
		d := distuv.StudentsT{Mu: 0, Sigma: 1, Nu: a}
		return continuous(d, op, x)
	case 9:
		d := distuv.F{D1: a, D2: b}
		return continuous(d, op, x)
	case 10:
		d := distuv.Binomial{N: a, P: b}
		return discrete(d, op, x)
	case 11:
		d := distuv.Poisson{Lambda: a}
		return discrete(d, op, x)
	}
	return math.NaN()
}

// A continuous family of the library.
type continuousDist interface {
	Prob(x float64) float64
	CDF(x float64) float64
	Survival(x float64) float64
	Quantile(p float64) float64
	Mean() float64
	Variance() float64
}

// A discrete family of the library.
type discreteDist interface {
	Prob(x float64) float64
	CDF(x float64) float64
	Survival(x float64) float64
	Mean() float64
	Variance() float64
}

// The result of one operation of a continuous family.
func continuous(d continuousDist, op int, x float64) float64 {
	switch op {
	case 0:
		return d.Prob(x)
	case 1:
		return d.CDF(x)
	case 2:
		return d.Survival(x)
	case 3:
		return d.Quantile(x)
	case 4:
		return d.Mean()
	case 5:
		return d.Variance()
	}
	return skewness(d, op)
}

// The result of one operation of a discrete family.
func discrete(d discreteDist, op int, x float64) float64 {
	switch op {
	case 0:
		return d.Prob(x)
	case 1:
		return d.CDF(x)
	case 2:
		return d.Survival(x)
	case 4:
		return d.Mean()
	case 5:
		return d.Variance()
	}
	return math.NaN()
}

// The skewness and the excess kurtosis of a family. The library gives
// them through the interfaces of the moments, and it holds them for some
// families only.
func skewness(d interface{}, op int) float64 {
	type moments interface {
		Skewness() float64
		ExKurtosis() float64
	}
	m, ok := d.(moments)
	if !ok {
		return math.NaN()
	}
	if op == 6 {
		return m.Skewness()
	}
	return m.ExKurtosis()
}

// runDistribution reads the lines of the file and writes one result for
// each. A panic of the library gives NaN, because it panics on a bad
// domain.
func runDistribution(path string) {
	f, err := os.Open(path)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	defer f.Close()
	sc := bufio.NewScanner(f)
	for sc.Scan() {
		w := strings.Fields(sc.Text())
		if len(w) < 2 {
			continue
		}
		fam, _ := strconv.Atoi(w[0])
		op, _ := strconv.Atoi(w[1])
		v := [4]float64{}
		for i := 0; i < 4 && i+2 < len(w); i++ {
			v[i], _ = strconv.ParseFloat(w[i+2], 64)
		}
		fmt.Println(guardDist(fam, op, v[0], v[1], v[2], v[3]))
	}
}

func guardDist(fam int, op int, x float64, a float64, b float64, c float64) (out string) {
	defer func() {
		if r := recover(); r != nil {
			out = "nan"
		}
	}()
	return strconv.FormatFloat(distribution(fam, op, x, a, b, c), 'g', -1, 64)
}
