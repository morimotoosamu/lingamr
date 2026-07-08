# Print method for var_stationarity

Print method for var_stationarity

## Usage

``` r
# S3 method for class 'var_stationarity'
print(x, ...)
```

## Arguments

- x:

  a `var_stationarity` object

- ...:

  additional arguments (unused)

## Value

The input object `x`, invisibly.

## Examples

``` r
s <- generate_varlingam_sample(n = 1000, seed = 42)
m <- lingam_var(s$data, lags = 1, reg_method = "ols", prune = FALSE)
print(check_var_stationarity(m))
#> === VAR Stationarity Check ===
#> Lag order:         1
#> Max |eigenvalue|:  0.5038  (threshold 1.00)
#> Stationary:        YES
```
