# Estimate the total causal effects between all variables at once

Estimate the total causal effects between all variables at once

## Usage

``` r
estimate_all_total_effects(
  X,
  lingam_result,
  method = "adaptive_lasso",
  lambda = "BIC",
  init_method = "ols"
)
```

## Arguments

- X:

  Original data (n_samples x n_features)

- lingam_result:

  Return value of lingam_direct()

- method:

  Regression method ("ols", "lasso", "adaptive_lasso")

- lambda:

  Lambda selection ("lambda.min", "lambda.1se", "AIC", "BIC")

- init_method:

  Method for estimating the initial weights of adaptive LASSO regression
  ("ols" or "ridge")

## Value

Matrix of total causal effects (n_features x n_features). **Convention:
`TE[i, j]` is the total causal effect from variable j to variable i (j
-\> i).** Same index convention as the adjacency matrix
`adjacency_matrix`. The sum of direct and indirect effects.

## Examples

``` r
LiNGAM_sample_1000 <- generate_lingam_sample_6()

model <- LiNGAM_sample_1000$data |>
  lingam_direct()

LiNGAM_sample_1000$data |>
  estimate_all_total_effects(model)
#>          x0 x1        x2        x3 x4 x5
#> x0 0.000000  0  0.000000  3.033460  0  0
#> x1 2.896907  0  1.909712 21.058733  0  0
#> x2 0.000000  0  0.000000  5.992677  0  0
#> x3 0.000000  0  0.000000  0.000000  0  0
#> x4 8.001464  0 -1.308131 18.276121  0  0
#> x5 4.015103  0  0.000000 12.179395  0  0
```
