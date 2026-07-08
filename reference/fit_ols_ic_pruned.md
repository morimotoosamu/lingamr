# OLS fit for a single predictor, pruned by information criterion

glmnet requires at least two predictor columns, so the penalized methods
fall back to OLS when only one predictor remains. Plain OLS never yields
an exact zero, which would make single-predictor edges unprunable: the
second variable in the causal order always has exactly one predictor, so
a spurious edge would survive even for fully independent data. To
preserve the sparse behavior of the penalized methods, the OLS
coefficient is kept only when adding the predictor improves the
information criterion over the intercept-only model; otherwise it is set
to exactly zero.

## Usage

``` r
fit_ols_ic_pruned(y, Xp, lambda)
```

## Arguments

- y:

  response variable (numeric vector)

- Xp:

  single-column predictor matrix

- lambda:

  lambda selection method of the calling fit

## Value

length-1 coefficient vector (0 when the predictor is pruned)

## Details

The criterion is AIC for `lambda = "AIC"` and BIC otherwise (the CV /
oracle lambdas have no single-predictor counterpart, so the sparsest
criterion, BIC, is used for them as well).
