# Estimate the total causal effect between two specified variables

Estimate the total causal effect between two specified variables

## Usage

``` r
estimate_total_effect(
  X,
  lingam_result,
  from_index,
  to_index,
  method = "adaptive_lasso",
  lambda = "BIC",
  init_method = "ols"
)
```

## Arguments

- X:

  Original data (matrix or data.frame)

- lingam_result:

  Return value of lingam_direct()

- from_index:

  Cause variable (1-based index or variable name)

- to_index:

  Effect variable (1-based index or variable name)

- method:

  Regression method ("ols", "lasso", "adaptive_lasso"). Default is
  adaptive_lasso

- lambda:

  Lambda selection ("lambda.min", "lambda.1se", "AIC", "BIC", "oracle").
  Default is BIC

- init_method:

  Method for estimating the initial weights of adaptive LASSO regression
  ("ols" or "ridge")

## Value

Estimated total causal effect

## Examples

``` r
LiNGAM_sample_1000 <- generate_lingam_sample_6()

model <- LiNGAM_sample_1000$data |>
  lingam_direct()

LiNGAM_sample_1000$data |>
  estimate_total_effect(model, 4, 1)
#>      x3 
#> 3.03346 
```
