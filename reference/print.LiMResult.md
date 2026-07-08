# Print method for LiMResult

Print method for LiMResult

## Usage

``` r
# S3 method for class 'LiMResult'
print(x, digits = 3, ...)
```

## Arguments

- x:

  LiMResult object

- digits:

  Number of digits to display

- ...:

  Additional arguments (unused)

## Value

The input object `x`, invisibly.

## Examples

``` r
set.seed(1)
dat <- generate_lim_sample(n = 300)
result <- lingam_lim(dat$data, is_continuous = dat$is_continuous)
print(result)
#> LiM Result
#>   Variables : 3
#>   Variable types: continuous, discrete, continuous
#>   Causal order: x1 -> x2 -> x3
#> 
#> Adjacency matrix (row = to, col = from):
#>       x1 x2 x3
#> x1  0.00  0  0
#> x2 -0.22  0  0
#> x3  0.00  1  0
```
