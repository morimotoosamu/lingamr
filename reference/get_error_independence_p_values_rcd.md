# Compute p-values for the independence of RCD residuals (HSIC-based)

Analogous to
[`get_error_independence_p_values_parce()`](https://morimotoosamu.github.io/lingamr/reference/get_error_independence_p_values_parce.md),
but for
[`lingam_rcd()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd.md)
results. Returns `NA` for any pair involving a variable whose row or
column in the adjacency matrix contains `NA` (residuals cannot be
computed for those variables).

## Usage

``` r
get_error_independence_p_values_rcd(X, rcd_result)
```

## Arguments

- X:

  Original data (matrix or data.frame)

- rcd_result:

  Return value of
  [`lingam_rcd()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd.md)

## Value

matrix of p-values (n_features x n_features)

## Examples

``` r
confounded <- generate_rcd_sample(n = 300, seed = 1)
result <- lingam_rcd(confounded$data)
round(get_error_independence_p_values_rcd(confounded$data, result), 3)
#>       x0    x1 x2    x3 x4    x5
#> x0    NA 0.000 NA 0.494 NA 0.088
#> x1 0.000    NA NA 0.808 NA 0.401
#> x2    NA    NA NA    NA NA    NA
#> x3 0.494 0.808 NA    NA NA 0.017
#> x4    NA    NA NA    NA NA    NA
#> x5 0.088 0.401 NA 0.017 NA    NA
```
