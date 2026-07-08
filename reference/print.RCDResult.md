# Print method for RCDResult

Print method for RCDResult

## Usage

``` r
# S3 method for class 'RCDResult'
print(x, digits = 3, ...)
```

## Arguments

- x:

  RCDResult object

- digits:

  Number of digits to display

- ...:

  Additional arguments (unused)

## Value

The input object `x`, invisibly.

## Examples

``` r
confounded <- generate_rcd_sample(n = 300, seed = 1)
result <- lingam_rcd(confounded$data)
print(result)
#> RCD Result
#>   Variables : 6
#> 
#> Ancestor sets:
#>   M(x0) = {x1, x3, x5}
#>   M(x1) = {x5}
#>   M(x2) = {x0, x1, x3, x5}
#>   M(x3) = {x5}
#>   M(x4) = {x0, x1, x3, x5}
#>   M(x5) = {}
#> 
#>   (NA entries in the adjacency matrix = suspected shared latent confounder)
#> 
#> Adjacency matrix (row = to, col = from):
#>       x0    x1 x2    x3 x4    x5
#> x0 0.000 1.116  0 0.989  0 0.000
#> x1 0.000 0.000  0 0.000  0 0.588
#> x2 0.810 0.000  0 0.000 NA 0.000
#> x3 0.000 0.000  0 0.000  0 0.449
#> x4 1.015 0.000 NA 0.000  0 0.000
#> x5 0.000 0.000  0 0.000  0 0.000
```
