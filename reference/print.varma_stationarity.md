# Print method for varma_stationarity

Print method for varma_stationarity

## Usage

``` r
# S3 method for class 'varma_stationarity'
print(x, ...)
```

## Arguments

- x:

  a `varma_stationarity` object

- ...:

  additional arguments (unused)

## Value

The input object `x`, invisibly.

## Examples

``` r
s <- generate_varmalingam_sample(n = 1000, seed = 42)
m <- lingam_varma(s$data,
  order = c(1, 1), criterion = NULL,
  reg_method = "ols", prune = FALSE
)
print(check_varma_stationarity(m))
#> === VARMA Stationarity / Invertibility Check ===
#> Order (p, q):         (1, 1)
#> Max |AR eigenvalue|:  0.4763  (threshold 1.00)
#> Stationary:           YES
#> Max |MA eigenvalue|:  0.2974  (threshold 1.00)
#> Invertible:           YES
```
