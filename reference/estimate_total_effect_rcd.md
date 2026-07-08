# Estimate the total causal effect between two variables (RCD)

Analogous to
[`estimate_total_effect()`](https://morimotoosamu.github.io/lingamr/reference/estimate_total_effect.md),
but for
[`lingam_rcd()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd.md)
results, which may contain `NA` entries in the adjacency matrix.

## Usage

``` r
estimate_total_effect_rcd(
  X,
  rcd_result,
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

- rcd_result:

  Return value of
  [`lingam_rcd()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd.md)

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
is part of a suspected latent confounder pair (its parents cannot be
identified). Also warns (without altering the estimate) if `to_index` is
an ancestor of `from_index` according to `ancestors_list`, since that is
inconsistent with a `from -> to` effect.

## Examples

``` r
confounded <- generate_rcd_sample(n = 300, seed = 1)
result <- lingam_rcd(confounded$data)

# A well-identified pair returns a numeric estimate
estimate_total_effect_rcd(confounded$data, result, from_index = 6, to_index = 1)
#>      x5 
#> 1.05674 
```
