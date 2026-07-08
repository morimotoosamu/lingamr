# Estimate the total causal effect between two variables (ParceLiNGAM)

Analogous to
[`estimate_total_effect()`](https://morimotoosamu.github.io/lingamr/reference/estimate_total_effect.md),
but for
[`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md)
results, which may contain `NA` entries in the adjacency matrix.

## Usage

``` r
estimate_total_effect_parce(
  X,
  parce_result,
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

- parce_result:

  Return value of
  [`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md)

- from_index:

  Cause variable (1-based index or variable name)

- to_index:

  Effect variable (1-based index or variable name)

- method:

  Regression method ("ols", "lasso", "adaptive_lasso", "ridge"). Default
  is adaptive_lasso

- lambda:

  Lambda selection ("lambda.min", "lambda.1se", "AIC", "BIC", "oracle").
  Default is BIC

- init_method:

  Method for estimating the initial weights of adaptive LASSO regression
  ("ols" or "ridge")

## Value

Estimated total causal effect, or `NA` (with a warning) if `from_index`
is part of an unresolved block (its parents cannot be identified).

## Examples

``` r
confounded <- generate_parce_sample(n = 500, seed = 1)
result <- lingam_parce(confounded$data, reg_method = "ols")

# A well-identified pair returns a numeric estimate
estimate_total_effect_parce(confounded$data, result, from_index = 1, to_index = 5)
#> [1] 0.4965035
```
