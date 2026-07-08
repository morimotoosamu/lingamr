# Default imputer: multiple imputation via mice::mice(method = "norm")

The closest standard R equivalent of the upstream Python default
(`sklearn.impute.IterativeImputer(sample_posterior = TRUE)`, Bayesian
linear regression chained equations with posterior sampling) is
`mice::mice(method = "norm")` (Bayesian linear regression). Numeric
values will not match the Python implementation; the multiple-imputation
design is equivalent.

## Usage

``` r
default_imputer(X_boot, n_repeats)
```

## Arguments

- X_boot:

  Bootstrap-resampled data (matrix, may contain NA)

- n_repeats:

  Number of imputed datasets to generate

## Value

A list of `n_repeats` complete numeric matrices
