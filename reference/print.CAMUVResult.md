# Print method for CAMUVResult

Print method for CAMUVResult

## Usage

``` r
# S3 method for class 'CAMUVResult'
print(x, digits = 3, ...)
```

## Arguments

- x:

  CAMUVResult object

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
  confounded <- generate_camuv_sample(n = 200, seed = 1)
  result <- lingam_camuv(confounded$data)
  print(result)
}
#> CAM-UV Result
#>   Variables : 6
#>   Regressor : gam
#> 
#> Parent sets:
#>   P(x0) = {}
#>   P(x1) = {x0}
#>   P(x2) = {}
#>   P(x3) = {x0}
#>   P(x4) = {x2}
#>   P(x5) = {}
#> 
#> Pairs with an unobserved causal/backdoor path (UCP/UBP):
#>   x2 -- x5
#>   x3 -- x4
#> 
#> Adjacency matrix (row = to, col = from):
#>   (entries are 0/1 edge indicators, not coefficients;
#>    NA = pair connected through an unobserved variable)
#>    x0 x1 x2 x3 x4 x5
#> x0  0  0  0  0  0  0
#> x1  1  0  0  0  0  0
#> x2  0  0  0  0  0 NA
#> x3  1  0  0  0 NA  0
#> x4  0  0  1 NA  0  0
#> x5  0  0 NA  0  0  0
# }
```
