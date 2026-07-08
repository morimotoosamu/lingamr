# Compute p-values for the independence of ParceLiNGAM residuals (HSIC-based)

Analogous to
[`get_error_independence_p_values()`](https://morimotoosamu.github.io/lingamr/reference/get_error_independence_p_values.md),
but for
[`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md)
results. Uses the HSIC gamma-approximation test
([`hsic_test_gamma()`](https://morimotoosamu.github.io/lingamr/reference/hsic_test_gamma.md))
rather than a correlation test, and returns `NA` for any pair involving
a variable whose row or column in the adjacency matrix contains `NA`
(residuals cannot be computed for those variables).

## Usage

``` r
get_error_independence_p_values_parce(X, parce_result)
```

## Arguments

- X:

  Original data (matrix or data.frame)

- parce_result:

  Return value of
  [`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md)

## Value

matrix of p-values (n_features x n_features)

## Examples

``` r
confounded <- generate_parce_sample(n = 500, seed = 1)
result <- lingam_parce(confounded$data, reg_method = "ols")
round(get_error_independence_p_values_parce(confounded$data, result), 3)
#>       x0    x1 x2 x3    x4    x5
#> x0    NA 0.453 NA NA 0.795 0.367
#> x1 0.453    NA NA NA 0.924 0.471
#> x2    NA    NA NA NA    NA    NA
#> x3    NA    NA NA NA    NA    NA
#> x4 0.795 0.924 NA NA    NA 0.873
#> x5 0.367 0.471 NA NA 0.873    NA
```
