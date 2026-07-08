# Print method for ParceLingamResult

Print method for ParceLingamResult

## Usage

``` r
# S3 method for class 'ParceLingamResult'
print(x, digits = 3, ...)
```

## Arguments

- x:

  ParceLingamResult object

- digits:

  Number of digits to display

- ...:

  Additional arguments (unused)

## Value

The input object `x`, invisibly.

## Examples

``` r
confounded <- generate_parce_sample(n = 300, seed = 42)
result <- lingam_parce(confounded$data, reg_method = "ols")
print(result)
#> Bottom-Up ParceLiNGAM Result
#>   Variables : 6
#>   Independence measure: hsic
#>   Causal order: x3 -> x2 -> x0 -> x5 -> x1 -> x4
#>   (NA entries in the adjacency matrix = unresolved order / suspected latent confounding)
#> 
#> Adjacency matrix (row = to, col = from):
#>       x0    x1     x2     x3 x4    x5
#> x0 0.000 0.000  0.027  0.463  0 0.000
#> x1 0.415 0.000  0.476  0.084  0 0.009
#> x2 0.000 0.000  0.000  0.856  0 0.000
#> x3 0.000 0.000  0.000  0.000  0 0.000
#> x4 0.405 0.023 -0.455 -0.035  0 0.086
#> x5 0.514 0.000  0.056 -0.024  0 0.000
```
