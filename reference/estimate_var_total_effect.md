# Estimate a total causal effect in a VAR-LiNGAM model

Estimates the total causal effect from `from_index` (optionally at lag
`from_lag`) to `to_index` (at the current time) using the fitted
VAR-LiNGAM model. Port of the Python reference `estimate_total_effect`:
the destination variable is regressed on the source variable together
with the source's parents (a back-door adjustment), and the source's
coefficient is returned.

## Usage

``` r
estimate_var_total_effect(X, result, from_index, to_index, from_lag = 0)
```

## Arguments

- X:

  original data (matrix or data frame), rows ordered in time

- result:

  a `VARLiNGAMResult` from
  [`lingam_var()`](https://morimotoosamu.github.io/lingamr/reference/lingam_var.md)

- from_index:

  source variable (1-based index or variable name)

- to_index:

  destination variable (1-based index or variable name)

- from_lag:

  lag of the source variable (0 = current time, default)

## Value

the estimated total effect (scalar)

## Examples

``` r
sample <- generate_varlingam_sample(n = 1000, seed = 42)
model <- lingam_var(sample$data, lags = 1, reg_method = "ols", prune = FALSE)

# total effect of x0 (current) on x2 (current)
estimate_var_total_effect(sample$data, model, from_index = 1, to_index = 3)
#> [1] -0.2558022
```
