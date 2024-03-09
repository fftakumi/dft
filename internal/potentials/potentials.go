// https://www.jstage.jst.go.jp/article/xshinpo/47/0/47_5/_pdf/-char/ja
package potentials

import "math"

type Potential interface {
	V(x []float64) []float64
}

// 調和振動子
type HarmonicOscillator struct {}
func (h HarmonicOscillator) V(x []float64) []float64 {
	y := make([]float64, len(x))
	for i, v := range x {
		y[i] = v * v
	}
	return y
}

// 非調和振動子
type AnharmonicOscillator struct {
	Mu float64
	Lamda float64
}
func (a AnharmonicOscillator) V(x []float64) []float64 {
	y := make([]float64, len(x))
	for i, v := range x {
		v2 := v * v
		v4 := v2 * v2
		y[i] = a.Mu * v2 + a.Lamda * v4
	}
	return y
}

// Morseポテンシャル
type Morse struct {}

func (m Morse) V(x []float64) []float64 {
	y := make([]float64, len(x))
	for i, v := range x {
		y[i] = math.Exp(-2 * v) -2 * math.Exp(-v)
	}
	return y
}

// 変形Poschl Teller
type ModifiedPoschlTeller struct {
	V0 float64
}

func (m ModifiedPoschlTeller) V(x []float64) []float64 {
	y := make([]float64, len(x))
	for i, v := range x {
		cosh2 := math.Cosh(v) * math.Cosh(v)
		y[i] = -m.V0 / cosh2
	}
	return y
}
