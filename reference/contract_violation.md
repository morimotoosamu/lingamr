# Signal a hook contract violation (imputer / cd_fit return value)

A plain error, but always called from outside the per-iteration
[`tryCatch()`](https://rdrr.io/r/base/conditions.html) in
[`bootstrap_with_imputation()`](https://morimotoosamu.github.io/lingamr/reference/bootstrap_with_imputation.md)
(see
[`validate_imputer_output()`](https://morimotoosamu.github.io/lingamr/reference/validate_imputer_output.md)
/
[`validate_cd_fit_output()`](https://morimotoosamu.github.io/lingamr/reference/validate_cd_fit_output.md)
call sites there), so that a programming error in a user-supplied hook
aborts the whole call immediately instead of being swallowed as a
per-iteration stochastic estimation failure.

## Usage

``` r
contract_violation(msg)
```

## Arguments

- msg:

  Error message
