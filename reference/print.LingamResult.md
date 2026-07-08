# Print method for LingamResult

Print method for LingamResult

## Usage

``` r
# S3 method for class 'LingamResult'
print(x, digits = 3, ...)
```

## Arguments

- x:

  LingamResult object

- digits:

  Number of digits to display

- ...:

  Additional arguments (unused)

## Value

The input object `x`, invisibly.

## Examples

``` r
LiNGAM_sample_1000 <- generate_lingam_sample_6()
result <- lingam_direct(LiNGAM_sample_1000$data, reg_method = "ols")
print(result)
#> Direct LiNGAM Result
#>   Variables : 6
#>   Causal order: x3 -> x2 -> x0 -> x4 -> x5 -> x1
#> 
#> Adjacency matrix (row = to, col = from):
#>       x0 x1     x2     x3     x4    x5
#> x0 0.000  0 -0.040  3.274  0.000 0.000
#> x1 3.237  0  1.965  0.014 -0.034 0.006
#> x2 0.000  0  0.000  5.993  0.000 0.000
#> x3 0.000  0  0.000  0.000  0.000 0.000
#> x4 7.992  0 -1.062  0.394  0.000 0.000
#> x5 3.873  0  0.069 -0.315  0.018 0.000
```
