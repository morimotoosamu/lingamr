# Print method for ImputationBootstrapResult

Print method for ImputationBootstrapResult

## Usage

``` r
# S3 method for class 'ImputationBootstrapResult'
print(x, ...)
```

## Arguments

- x:

  ImputationBootstrapResult object

- ...:

  Additional arguments (unused)

## Value

The input object `x`, invisibly.

## Examples

``` r
set.seed(1)
sample6 <- generate_lingam_sample_6(n = 300, seed = 1)
X <- sample6$data
X$x5[sample.int(nrow(X), size = 30)] <- NA

if (requireNamespace("mice", quietly = TRUE)) {
  res <- bootstrap_with_imputation(X,
    n_sampling = 5L, n_repeats = 3L, seed = 42, verbose = FALSE
  )
  print(res)
}
#> ImputationBootstrapResult: 5 samplings x 3 repeats, 6 features, 30 missing cells (original data)
```
