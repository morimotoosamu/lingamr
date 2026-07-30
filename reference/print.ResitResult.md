# Print method for ResitResult

Print method for ResitResult

## Usage

``` r
# S3 method for class 'ResitResult'
print(x, digits = 3, ...)
```

## Arguments

- x:

  ResitResult object

- digits:

  Number of digits to display

- ...:

  Additional arguments (unused)

## Value

The input object `x`, invisibly.

## Examples

``` r
# \donttest{
if (requireNamespace("mgcv", quietly = TRUE)) {
  nonlinear <- generate_resit_sample(n = 300, seed = 1)
  result <- lingam_resit(nonlinear$data)
  print(result)
}
#> RESIT Result
#>   Variables : 4
#>   Regressor : gam
#>   Causal order: x0 -> x1 -> x2 -> x3
#> 
#> Adjacency matrix (row = to, col = from):
#>   (entries are 0/1 edge indicators, not coefficients)
#>    x0 x1 x2 x3
#> x0  0  0  0  0
#> x1  1  0  0  0
#> x2  1  1  0  0
#> x3  0  0  1  0
# }
```
