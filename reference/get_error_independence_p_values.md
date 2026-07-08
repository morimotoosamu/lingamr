# Compute p-values for the independence test of the errors

Compute p-values for the independence test of the errors

## Usage

``` r
get_error_independence_p_values(X, lingam_result, method = "spearman")
```

## Arguments

- X:

  original data (matrix or data.frame)

- lingam_result:

  return value of lingam_direct()

- method:

  type of correlation coefficient ("spearman", "pearson", "kendall").
  "kendall" uses the O(n^2)-per-pair algorithm in
  [`stats::cor.test()`](https://rdrr.io/r/stats/cor.test.html); for
  large `n` (beyond 5000) this warns, and "spearman" is a much faster
  alternative with similar rank-based semantics.

## Value

matrix of p-values (n_features x n_features)

## Examples

``` r
# Load the sample data
LiNGAM_sample_1000 <- generate_lingam_sample_6()

# Run Direct LiNGAM
result <- LiNGAM_sample_1000$data |>
  lingam_direct(reg_method = "ols")

# Compute p-values (default: Spearman)
p_vals <- get_error_independence_p_values(LiNGAM_sample_1000$data, result)
round(p_vals, 3)
#>       x0    x1    x2    x3    x4    x5
#> x0    NA 0.980 0.978 0.990 0.997 0.962
#> x1 0.980    NA 0.995 0.997 0.943 0.999
#> x2 0.978 0.995    NA 0.919 0.945 0.995
#> x3 0.990 0.997 0.919    NA 0.989 0.999
#> x4 0.997 0.943 0.945 0.989    NA 0.986
#> x5 0.962 0.999 0.995 0.999 0.986    NA

# Compute with Kendall
p_vals_k <- get_error_independence_p_values(LiNGAM_sample_1000$data, result, method = "kendall")
round(p_vals_k, 3)
#>       x0    x1    x2    x3    x4    x5
#> x0    NA 0.994 0.993 0.988 0.979 0.964
#> x1 0.994    NA 0.976 0.979 0.948 0.989
#> x2 0.993 0.976    NA 0.912 0.961 0.972
#> x3 0.988 0.979 0.912    NA 1.000 0.986
#> x4 0.979 0.948 0.961 1.000    NA 0.995
#> x5 0.964 0.989 0.972 0.986 0.995    NA
```
