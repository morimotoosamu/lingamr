# Print method for MultiGroupLingamResult

Print method for MultiGroupLingamResult

## Usage

``` r
# S3 method for class 'MultiGroupLingamResult'
print(x, digits = 3, ...)
```

## Arguments

- x:

  MultiGroupLingamResult object

- digits:

  Number of digits to display

- ...:

  Additional arguments (unused)

## Value

The input object `x`, invisibly.

## Examples

``` r
mg <- generate_multi_group_sample()
res <- lingam_multi_group(mg$data_list, reg_method = "ols")
print(res)
#> Multi-Group Direct LiNGAM Result
#>   Groups      : 2 (group1, group2)
#>   Variables   : 6
#>   Causal order (common): x3 -> x0 -> x5 -> x2 -> x4 -> x1
#> 
#> [group1] Adjacency matrix (row = to, col = from):
#>        x0 x1     x2     x3     x4    x5
#> x0  0.000  0  0.000  3.033  0.000 0.000
#> x1  3.237  0  1.965  0.014 -0.034 0.006
#> x2 -0.236  0  0.000  6.112  0.000 0.049
#> x3  0.000  0  0.000  0.000  0.000 0.000
#> x4  7.921  0 -1.063  0.399  0.000 0.018
#> x5  4.016  0  0.000 -0.003  0.000 0.000
#> 
#> [group2] Adjacency matrix (row = to, col = from):
#>       x0 x1     x2     x3    x4     x5
#> x0 0.000  0  0.000  3.504 0.000  0.000
#> x1 2.732  0  2.568  0.083 0.034  0.093
#> x2 0.154  0  0.000  6.322 0.000 -0.024
#> x3 0.000  0  0.000  0.000 0.000  0.000
#> x4 8.483  0 -1.487 -0.110 0.000  0.006
#> x5 4.515  0  0.000 -0.045 0.000  0.000
```
