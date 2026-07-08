# Print a VARBootstrapResult

Print a VARBootstrapResult

## Usage

``` r
# S3 method for class 'VARBootstrapResult'
print(x, ...)
```

## Arguments

- x:

  a VARBootstrapResult object

- ...:

  additional arguments (unused)

## Value

The input object `x`, invisibly.

## Examples

``` r
s <- generate_varlingam_sample(n = 500, seed = 42)
bs <- lingam_var_bootstrap(s$data,
  n_sampling = 10L, lags = 1, criterion = NULL,
  reg_method = "ols", prune = FALSE, seed = 1, verbose = FALSE
)
print(bs)
#> VARBootstrapResult: 10 samplings, 3 features, lag order 1
```
