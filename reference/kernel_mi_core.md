# カーネル法の相互情報量：本体

求める量は 2n x 2n 行列の logdet の差だが、ブロック構造と Schur
補行列により n x n の Cholesky 分解だけで等価に計算できる：
`MI = -1/2 * (logdet(tmp2^2 - K2 K1 tmp1^-2 K1 K2) - logdet(tmp2^2))`

## Usage

``` r
kernel_mi_core(E1, x2, kappa, sigma)
```

## Arguments

- E1:

  [`kernel_mi_prepare()`](https://morimotoosamu.github.io/lingamr/reference/kernel_mi_prepare.md)
  で前計算した変数1側の行列

- x2:

  変数2のベクトル

- kappa:

  正則化パラメータ

- sigma:

  ガウスカーネルの幅

## Value

相互情報量
