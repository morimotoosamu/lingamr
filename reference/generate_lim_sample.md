# Generate sample data for LiM (3 mixed variables)

Generates a small dataset with a known causal chain of continuous and
binary (0/1) discrete variables:
`x1 (continuous) -> x2 (discrete) -> x3 (continuous)`. Continuous
variables use Laplace-distributed noise (non-Gaussian, as required by
LiNGAM-family methods); the discrete variable is drawn from a Bernoulli
distribution whose logit is a linear function of its parent.

## Usage

``` r
generate_lim_sample(n = 1000L, seed = NULL)
```

## Arguments

- n:

  number of samples (default: 1000)

- seed:

  random seed. If `NULL` (default), no seed is set and results are not
  reproducible across calls.

## Value

A list with three elements:

- `data`: data frame with columns `x1`, `x2`, `x3`.

- `adjacency_matrix`: 3x3 true adjacency matrix, following the lingamr
  convention (`m[to, from]`, i.e. `adjacency_matrix["x2", "x1"]` is the
  coefficient of the x1 -\> x2 edge).

- `is_continuous`: logical vector `c(TRUE, FALSE, TRUE)` marking which
  columns of `data` are continuous.

## Examples

``` r
dat <- generate_lim_sample(n = 500, seed = 1)
result <- lingam_lim(dat$data, is_continuous = dat$is_continuous)
result$adjacency_matrix
#>    x1 x2        x3
#> x1  0  1 0.1358125
#> x2  0  0 0.0000000
#> x3  0  1 0.0000000
```
