# Evaluate the independence of a residual against a fixed set of predictors

Used for the `length(Uc) == 1` special case in
[`find_exo_vec()`](https://morimotoosamu.github.io/lingamr/reference/find_exo_vec.md),
where there is nothing left to compare against.

## Usage

``` r
parce_eval_independence(X, predictors, R, independence)
```

## Arguments

- X:

  data matrix

- predictors:

  predictor indices (may be empty)

- R:

  residual vector

- independence:

  "hsic" or "fcorr"

## Value

evaluation value (Fisher-combined p-value for hsic, max F-correlation
for fcorr)
