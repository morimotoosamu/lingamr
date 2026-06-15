# Direct LiNGAM のブートストラップ

Direct LiNGAM のブートストラップ

## Usage

``` r
lingam_direct_bootstrap(
  X,
  n_sampling,
  prior_knowledge = NULL,
  apply_prior_knowledge_softly = FALSE,
  measure = "pwling",
  reg_method = "adaptive_lasso",
  lambda = "BIC",
  init_method = "ols",
  seed = NULL,
  verbose = TRUE,
  parallel = FALSE,
  n_cores = NULL
)
```

## Arguments

- X:

  数値行列 (n_samples x n_features)

- n_sampling:

  ブートストラップの反復回数

- prior_knowledge:

  事前知識行列 (NULL可)

- apply_prior_knowledge_softly:

  事前知識のソフト適用 (logical)

- measure:

  独立性の評価尺度 ("pwling" or "kernel")

- reg_method:

  回帰手法 ("ols", "lasso", "adaptive_lasso", "ridge")

- lambda:

  ラムダ選択 ("lambda.min", "lambda.1se", "AIC", "BIC","oracle")

- init_method:

  適応的LASSO回帰の初期重みの推定手法 ("ols" または "ridge")。
  [`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md)
  の同名引数と同じ。

- seed:

  乱数シード (NULL可)

- verbose:

  進捗を表示するか (logical)

- parallel:

  並列処理を行うか (logical)。`TRUE` の場合、各ブートストラップ
  反復を複数コアに分散して実行する。

- n_cores:

  使用するコア数 (整数, NULL可)。`NULL` の場合は安全のため最大 2
  コアに制限される。`parallel = FALSE` のときは無視される。

## Value

BootstrapResult (list)

## Details

`parallel = TRUE`
を指定すると、[`parallel::makePSOCKcluster()`](https://rdrr.io/r/parallel/makeCluster.html)
による
ソケットクラスターで反復を分散実行する。クラスターは処理終了時・エラー発生時
いずれの場合も [`on.exit()`](https://rdrr.io/r/base/on.exit.html)
により必ず解放される。

**再現性について:** 並列実行時は
[`parallel::clusterSetRNGStream()`](https://rdrr.io/r/parallel/RngStream.html)
による L'Ecuyer の並列乱数ストリームを用いる。同じ `seed`・同じ
`n_cores` であれば 結果は再現するが、逐次実行 (`parallel = FALSE`)
の結果とは数値的に一致しない。 厳密に逐次版と同じ結果が必要な場合は
`parallel = FALSE` を使用すること。

## Examples

``` r
LiNGAM_sample_1000 <- generate_lingam_sample_6()

# Fast example with OLS
bs <- lingam_direct_bootstrap(LiNGAM_sample_1000$data,
  n_sampling = 10L,
  reg_method = "ols",
  seed = 42
)
#> Bootstrap: 10 iterations, method=ols (sequential)
#>   iteration 1 / 10
#>   iteration 10 / 10
#> Completed in 0.1 seconds.
get_probabilities(bs)
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]  0.0  0.1  0.5  0.9  0.1  0.0
#> [2,]  0.9  0.0  0.9  0.9  0.4  0.5
#> [3,]  0.5  0.1  0.0  0.9  0.1  0.4
#> [4,]  0.1  0.1  0.1  0.0  0.1  0.1
#> [5,]  0.9  0.6  0.9  0.9  0.0  0.5
#> [6,]  1.0  0.5  0.6  0.9  0.5  0.0

# \donttest{
# With LASSO (requires glmnet)
bs_lasso <- lingam_direct_bootstrap(LiNGAM_sample_1000$data,
  n_sampling = 30L,
  seed = 42
)
#> Bootstrap: 30 iterations, method=adaptive_lasso (sequential)
#>   iteration 1 / 30
#>   iteration 10 / 30
#>   iteration 20 / 30
#>   iteration 30 / 30
#> Completed in 1.0 seconds.

# Parallel execution on 2 cores
bs_par <- lingam_direct_bootstrap(LiNGAM_sample_1000$data,
  n_sampling = 30L,
  seed = 42,
  parallel = TRUE,
  n_cores = 2L
)
#> Bootstrap: 30 iterations, method=adaptive_lasso (parallel, 2 cores)
#> Completed in 2.5 seconds.
# }
```
