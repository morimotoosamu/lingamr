# Convert a LiMResult to a tidy data.frame

Converts the estimated adjacency matrix of a LiM model into a
long-format data.frame with one edge per row, exactly like
[`tidy.LingamResult()`](https://morimotoosamu.github.io/lingamr/reference/tidy.LingamResult.md).

## Usage

``` r
# S3 method for class 'LiMResult'
tidy(x, threshold = 0, ...)
```

## Arguments

- x:

  The return value of
  [`lingam_lim()`](https://morimotoosamu.github.io/lingamr/reference/lingam_lim.md)
  (a `LiMResult` object)

- threshold:

  Coefficients with an absolute value at or below this are not treated
  as edges (default: 0)

- ...:

  Unused

## Value

data.frame(from, to, estimate)

## Examples

``` r
set.seed(1)
dat <- generate_lim_sample(n = 300)
model <- lingam_lim(dat$data, is_continuous = dat$is_continuous)
tidy(model)
#>   from to   estimate
#> 1   x1 x2 -0.2197992
#> 2   x2 x3  1.0000000
```
