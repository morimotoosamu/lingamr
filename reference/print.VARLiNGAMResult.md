# Print method for VARLiNGAMResult

Print method for VARLiNGAMResult

## Usage

``` r
# S3 method for class 'VARLiNGAMResult'
print(x, digits = 3, ...)
```

## Arguments

- x:

  VARLiNGAMResult object

- digits:

  number of digits to display

- ...:

  additional arguments (unused)

## Value

The input object `x`, invisibly.

## Examples

``` r
sample <- generate_varlingam_sample(n = 500, seed = 42)
model <- lingam_var(sample$data, lags = 1, reg_method = "ols", prune = FALSE)
print(model)
#> VAR-LiNGAM Result
#>   Variables : 3
#>   Lag order : 1
#>   Causal order (instantaneous): x0 -> x1 -> x2
#> 
#> Instantaneous adjacency matrix B0 (row = to, col = from):
#>       x0     x1 x2
#> x0 0.000  0.000  0
#> x1 0.627  0.000  0
#> x2 0.020 -0.505  0
#> 
#> Lagged adjacency matrix B1 (row = to, col = from):
#>       x0     x1     x2
#> x0 0.341  0.051  0.314
#> x1 0.018  0.280 -0.026
#> x2 0.055 -0.075  0.423
```
