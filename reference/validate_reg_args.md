# Validate the regression-method argument trio shared by bootstrap functions

Applies [`match.arg()`](https://rdrr.io/r/base/match.arg.html) to
`reg_method` / `lambda` / `init_method` and rejects the unsupported
`lambda = "oracle"` combination, with the same error text as the
algorithm entry points. Bootstrap functions call this before starting
the cluster so invalid arguments fail fast instead of erroring
confusingly inside workers.

## Usage

``` r
validate_reg_args(reg_method, lambda, init_method)
```

## Arguments

- reg_method:

  Regression method argument as received.

- lambda:

  Lambda selection argument as received.

- init_method:

  Adaptive-LASSO initial estimator argument as received.

## Value

`list(reg_method = , lambda = , init_method = )` with matched values.
