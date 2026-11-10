// The special function mode of the gonum oracle. Each input line holds
// the function number and three arguments. Each output line holds one
// result, or the word nan when the library gives no number.
package main

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"strconv"
	"strings"

	"gonum.org/v1/gonum/mathext"
)

// special gives the result of the function number fn. A call that the
// libraries do not hold gives NaN.
func special(fn int, a, b, x float64) float64 {
	switch fn {
	case 0:
		return math.Erf(x)
	case 1:
		return math.Erfc(x)
	case 2:
		return math.Erfinv(x)
	case 3:
		return math.Erfcinv(x)
	case 4:
		return math.Gamma(x)
	case 5:
		v, _ := math.Lgamma(x)
		return v
	case 6:
		_, s := math.Lgamma(x)
		return float64(s)
	case 7:
		return mathext.Digamma(x)
	case 8:
		return mathext.Beta(a, b)
	case 9:
		return mathext.Lbeta(a, b)
	case 11:
		return mathext.GammaIncReg(a, x)
	case 12:
		return mathext.GammaIncRegComp(a, x)
	case 13:
		return mathext.GammaIncRegInv(a, x)
	case 14:
		return mathext.GammaIncRegCompInv(a, x)
	case 15:
		return mathext.RegIncBeta(a, b, x)
	case 16:
		return mathext.InvRegIncBeta(a, b, x)
	case 17:
		return math.Log1p(x)
	case 18:
		return math.Expm1(x)
	}
	return math.NaN()
}

// runSpecial reads the lines of the file and writes one result for each.
// A panic of the library gives NaN, because gonum panics on a bad domain.
func runSpecial(path string) {
	f, err := os.Open(path)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	defer f.Close()
	sc := bufio.NewScanner(f)
	for sc.Scan() {
		w := strings.Fields(sc.Text())
		if len(w) != 4 {
			continue
		}
		fn, _ := strconv.Atoi(w[0])
		v := [3]float64{}
		for i := 0; i < 3; i++ {
			v[i], _ = strconv.ParseFloat(w[i+1], 64)
		}
		fmt.Println(guard(fn, v[0], v[1], v[2]))
	}
}

func guard(fn int, a, b, x float64) (out string) {
	defer func() {
		if r := recover(); r != nil {
			out = "nan"
		}
	}()
	return strconv.FormatFloat(special(fn, a, b, x), 'g', -1, 64)
}
