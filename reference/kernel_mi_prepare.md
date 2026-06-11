# カーネル法の相互情報量：変数1側の前計算

[`kernel_mi_core()`](https://morimotoosamu.github.io/lingamr/reference/kernel_mi_core.md)
で使う行列 `E1 = tmp1^-1 K1`（`tmp1 = K1 + n*kappa/2 * I`）
を計算する。候補変数ごとに1回だけ呼べばよく、ペアごとの再計算を避けられる。

## Usage

``` r
kernel_mi_prepare(x, kappa, sigma)
```

## Arguments

- x:

  変数1のベクトル

- kappa:

  正則化パラメータ

- sigma:

  ガウスカーネルの幅

## Value

行列 E1 (n x n)
