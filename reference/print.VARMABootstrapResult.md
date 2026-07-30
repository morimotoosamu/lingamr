# Print a VARMABootstrapResult

Print a VARMABootstrapResult

## Usage

``` r
# S3 method for class 'VARMABootstrapResult'
print(x, ...)
```

## Arguments

- x:

  a VARMABootstrapResult object

- ...:

  additional arguments (unused)

## Value

The input object `x`, invisibly.

## Examples

``` r
s <- generate_varmalingam_sample(n = 300, seed = 42)
bs <- lingam_varma_bootstrap(s$data,
  n_sampling = 5L, order = c(1, 1), criterion = NULL,
  reg_method = "ols", prune = FALSE, seed = 1, verbose = FALSE
)
print(bs)
#> VARMABootstrapResult: 5 samplings, 3 features, order (1, 1)
```
