# Estimate a total causal effect in a VARMA-LiNGAM model

Estimates the total causal effect from `from_index` (optionally at lag
`from_lag`) to `to_index` (at the current time) using the fitted
VARMA-LiNGAM model. Port of the Python reference
`estimate_total_effect`: the destination variable is regressed on the
source variable together with the source's parents (a back-door
adjustment) over lagged X and lagged residual regressors, and the
source's coefficient is returned.

## Usage

``` r
estimate_varma_total_effect(X, result, from_index, to_index, from_lag = 0)
```

## Arguments

- X:

  original data (matrix or data frame), rows ordered in time; must be
  the data the model was fitted on (the residuals stored in `result` are
  aligned to it)

- result:

  a `VARMALiNGAMResult` from
  [`lingam_varma()`](https://morimotoosamu.github.io/lingamr/reference/lingam_varma.md)

- from_index:

  source variable (1-based index or variable name)

- to_index:

  destination variable (1-based index or variable name)

- from_lag:

  lag of the source variable (0 = current time, default)

## Value

the estimated total effect (scalar)

## Details

The Python reference requires the residual matrix `E` as an argument;
here it is reconstructed internally from `result$residuals` and the
fitted instantaneous matrix B0, so only the original data are needed.

## Examples

``` r
sample <- generate_varmalingam_sample(n = 1000, seed = 42)
model <- lingam_varma(sample$data,
  order = c(1, 1), criterion = NULL,
  reg_method = "ols", prune = FALSE
)

# total effect of x0 (current) on x2 (current)
estimate_varma_total_effect(sample$data, model, from_index = 1, to_index = 3)
#> [1] -0.2569613
```
