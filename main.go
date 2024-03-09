package main

import (
	"errors"
	"fmt"
	"log"
	"math"

	"github.com/fftakumi/dft-go/internal/potentials"
	"github.com/fftakumi/dft-go/internal/util"
	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/plotutil"
	"gonum.org/v1/plot/vg"
)

func matPrint(X mat.Matrix) {
	fa := mat.Formatted(X, mat.Prefix(""), mat.Squeeze())
	fmt.Printf("%v\n", fa)
}

func ElectronDensity(nElec int, x, psi []float64) []float64 {
    numx := len(x)
    if (numx != len(psi)) {
        panic("長さか異なります。")
    }
    psi2 := make([]float64, numx)
    floats.MulTo(psi2, psi, psi)

    I := util.Integral(x, psi)
    normPsi := make([]float64, numx)
    floats.Add(psi, normPsi)
    floats.Scale(I, normPsi)
    
    isEven := (nElec % 2) == 0
    var fn []float64
    if (isEven) {
        fn = make([]float64, nElec / 2.)
    } else {
        fn = make([]float64, nElec / 2. + 1.)
    }
    floats.AddConst(2, fn)

    density := make([]float64, numx)
    for _, v := range fn {
        floats.AddScaledTo(density, nil, v, psi)
    }
    return density
}

func Hamiltonian(x []float64, p potentials.Potential) mat.Dense {
	numGrid := len(x)
	var term1 mat.Dense
	d2 := util.D2(numGrid, x[1]-x[0])
	term1.Scale(-1.0/2.0, d2)

	term2 := mat.NewDiagDense(len(x), p.V(x))

	term1.Add(&term1, term2)
	return term1
}

func SolveSchrodinger(x []float64, v potentials.Potential) (mat.Eigen, error) {
	h := Hamiltonian(x, v)
	var eig mat.Eigen
	ok := eig.Factorize(&h, mat.EigenRight)
	if ok {
		return eig, nil
	} else {
		log.Fatal("Eigendecomposition failed")
        return eig, errors.New("固有値の計算に失敗しました。")
	}
}

func GetStableEigenFunv(eig mat.Eigen) []float64 {
	var eigVec mat.CDense
	eig.VectorsTo(&eigVec)
	log.Print(eigVec.Dims())

	eigFloat64 := make([]float64, len(eig.Values(nil)))
	for i := 0; i < len(eigFloat64); i++ {
		eigFloat64[i] = real(eig.Values(nil)[i])
	}

	maxIdx := util.MaxIdx[float64](eigFloat64)

	eigf := make([]float64, len(eig.Values(nil)))
	for i := 0; i < len(eig.Values(nil)); i++ {
		eigf[i] = real(eigVec.At(i, maxIdx))
	}

	return eigf
}

func LDA(x, n []float64) float64 {
    engy := -3. / 4. * math.Pow(3. / math.Pi, 1./ 3.)
    engy *= util.Integral(x, n)
    return engy
}

func LDAPotential(n []float64) []float64 {
    c := -math.Pow(3 / math.Pi, 1. / 3.)
    p := make([]float64, len(n))
    for i, v := range n {
        p[i] = c * math.Pow(v, 1. / 3.)
    }
    return p
}

func main() {
	numGrid := 200
	x := util.Linspace(-5, 5, numGrid)
	p := potentials.ModifiedPoschlTeller{V0: 0}
	eig, err := SolveSchrodinger(x, p)
    if err != nil {
        return
    }
    phix := GetStableEigenFunv(eig)

    pltxy := make(plotter.XYs, len(phix))
    for i, v := range phix {
        pltxy[i].X = x[i]
		pltxy[i].Y = v
    }
	// インスタンスを生成
	plt := plot.New()

	// 表示項目の設定
	plt.Title.Text = "only english title"
	plt.X.Label.Text = "X axis"
	plt.Y.Label.Text = "Y axis"

	// グラフを描画
	err = plotutil.AddLinePoints(plt, pltxy)
	if err != nil {
		panic(err)
	}

	// 描画結果を保存
	// "5*vg.Inch" の数値を変更すれば，保存する画像のサイズを調整できます．
	if err := plt.Save(5*vg.Inch, 5*vg.Inch, "eigvec.png"); err != nil {
		panic(err)
	}
    log.Print("Success")
}
