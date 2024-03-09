package util

import (
	"math"
	"sort"

	"gonum.org/v1/gonum/mat"
)

func MaxIdx[T float64 | int](x []T) int {
    ix := make([] struct {
        ID int
        V T
    }, len(x))
    for i := 0; i < len(x); i++ {
        ix[i].ID = i
        ix[i].V = x[i]
    }

    sort.Slice(ix, func(i, j int) bool {
        return ix[i].V < ix[j].V
    })
    
    return ix[0].ID
}

func Integral(x, y []float64) float64 {
	var sum float64 = 0
	for i := 0; i < len(x) - 1; i++ {
		h := x[i + 1] - x[i]
		sum += (y[i + 1] + y[i]) * h / 2
	}

	return sum
}


func Linspace(a_start, a_end float64, a_n int) (ret []float64) {
	ret = make([]float64, a_n)
	if a_n == 1 {
		ret[0] = a_end
		return ret
	}
	delta := (a_end - a_start) / (float64(a_n) - 1)
	for i := 0; i < a_n; i++ {
		ret[i] = float64(a_start) + (delta * float64(i))
	}
	return
}

func Eye(n int) mat.Matrix {
	ones := make([]float64, n)
	for i := range ones {
		ones[i] = 1
	}
	eye := mat.NewDiagDense(n, ones)
	return eye
}


func DiagFlat(v float64, n, k int) mat.Matrix {
	zeros := make([]float64, n*n)
	diag := mat.NewDense(n, n, zeros)
	if k >= 0 {
		for i := 0; i < n-k; i++ {
			diag.Set(i, i+k, 1)
		}
	} else {
		for i := 0; i < n-k; i++ {
			diag.Set(i+k, i, 1)
		}
	}
	return diag
}

func D(n int, h float64) mat.Matrix {
	var d mat.Dense
	eye := Eye(n)
	d.Scale(-1, eye)
	diag := DiagFlat(1, n, 1)
	d.Add(&d, diag)
	d.Scale(1/h, &d)
	return &d
}


func D2(n int, h float64) mat.Matrix {
	var d2 mat.Dense
	d := D(n, h)
	var mdt mat.Dense
	mdt.Scale(-1, d.T())
	d2.Product(d, &mdt)
	return &d2
}


func Sin(x mat.Vector) mat.Vector {
	r, c := x.Dims()
	y := mat.NewDense(r, c, nil)
	for _r := 0; _r < r; _r++ {
		for _c := 0; _c < c; _c++ {
			y.Set(_r, _c, math.Sin(x.At(_r, _c)))
		}
	}
	return y.ColView(0)
}


func Square(x []float64) []float64 {
    y := make([]float64, len(x))
    for i := 0; i < len(x); i++ {
        y[i] = x[i] * x[i]
    }
    return y
}
