# Convert a MultiGroupLingamResult to a tidy data.frame

Stacks the per-group edge lists into a single long-format data.frame
with a `group` column in front, one edge per row per group (the causal
order is shared across groups but the coefficients differ).

## Usage

``` r
# S3 method for class 'MultiGroupLingamResult'
tidy(x, threshold = 0, ...)
```

## Arguments

- x:

  The return value of
  [`lingam_multi_group()`](https://morimotoosamu.github.io/lingamr/reference/lingam_multi_group.md)
  (a `MultiGroupLingamResult` object)

- threshold:

  Coefficients with an absolute value at or below this are not treated
  as edges (default: 0)

- ...:

  Unused

## Value

data.frame(group, from, to, estimate)

## Examples

``` r
mg <- generate_multi_group_sample()
model <- lingam_multi_group(mg$data_list, reg_method = "ols")
tidy(model)
#>     group from to     estimate
#> 1  group1   x0 x1  3.236756411
#> 2  group1   x0 x2 -0.235903083
#> 3  group1   x0 x4  7.920618753
#> 4  group1   x0 x5  4.015868568
#> 5  group1   x2 x1  1.965485693
#> 6  group1   x2 x4 -1.062516157
#> 7  group1   x3 x0  3.033460095
#> 8  group1   x3 x1  0.013952441
#> 9  group1   x3 x2  6.112126791
#> 10 group1   x3 x4  0.399217297
#> 11 group1   x3 x5 -0.002581799
#> 12 group1   x4 x1 -0.033970558
#> 13 group1   x5 x1  0.005627409
#> 14 group1   x5 x2  0.048947659
#> 15 group1   x5 x4  0.017844824
#> 16 group2   x0 x1  2.732221630
#> 17 group2   x0 x2  0.154154505
#> 18 group2   x0 x4  8.482925725
#> 19 group2   x0 x5  4.514996602
#> 20 group2   x2 x1  2.568176885
#> 21 group2   x2 x4 -1.486976210
#> 22 group2   x3 x0  3.504076134
#> 23 group2   x3 x1  0.082655142
#> 24 group2   x3 x2  6.321827207
#> 25 group2   x3 x4 -0.110407739
#> 26 group2   x3 x5 -0.045156286
#> 27 group2   x4 x1  0.034024614
#> 28 group2   x5 x1  0.092726357
#> 29 group2   x5 x2 -0.023744698
#> 30 group2   x5 x4  0.005953567
```
