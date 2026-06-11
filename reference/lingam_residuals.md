# LiNGAM モデルの残差（誤差項）を計算する

入力の検証（`LingamResult` であること・X
が数値であること・次元の一致）を 行ったうえで `E = X - X B^T`
を返す。残差ベースの診断関数で共通利用する。

## Usage

``` r
lingam_residuals(X, lingam_result)
```

## Arguments

- X:

  元データ (matrix or data.frame)

- lingam_result:

  [`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md)
  の返り値

## Value

残差行列 (n_samples x n_features)。X の列名を保持する。
