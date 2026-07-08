# Evaluate model fit of an estimated causal graph

Fits the causal graph implied by `adjacency_matrix` as a structural
equation model (SEM) via
[`lavaan::sem()`](https://rdrr.io/pkg/lavaan/man/sem.html) and returns
standard SEM fit measures (CFI, RMSEA, AIC/BIC, etc.). This is an R port
of the Python `lingam.utils.evaluate_model_fit()`, which delegates to
the Python package `semopy`; this R version delegates to `lavaan`
instead.

## Usage

``` r
evaluate_model_fit(adjacency_matrix, X, is_ordinal = NULL)
```

## Arguments

- adjacency_matrix:

  p x p numeric adjacency matrix (NA allowed for latent confounder
  pairs), or a lingamr result object (e.g. `LingamResult`,
  `ParceLingamResult`, `LiMResult`) with an `adjacency_matrix` element,
  from which the matrix is extracted automatically

- X:

  numeric matrix or data frame (n_samples x p) with no missing values

- is_ordinal:

  logical or 0/1 vector of length p. `TRUE` marks a variable as ordinal
  (categorical), fit with `lavaan`'s WLSMV-based estimator. `NULL`
  (default) treats all variables as continuous.

## Value

A one-row data.frame of fit measures: DoF, DoF Baseline, chi2, chi2
p-value, chi2 Baseline, CFI, GFI, AGFI, NFI, TLI, RMSEA, AIC, BIC,
LogLik. When `is_ordinal` is used, AIC/BIC/LogLik and some other
measures are not defined by the WLSMV estimator and are returned as
`NA`.

## Details

- **Optional dependency**: this function requires the lavaan package
  (listed in `Suggests`, not `Imports`). Install it with
  `install.packages("lavaan")`.

- **Latent confounders**: an `NA` element `B[i, j]` (as produced by e.g.
  [`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md)
  for a suspected latent confounder between variables i and j) is
  represented as a residual covariance `xi ~~ xj` in the lavaan model.
  This is algebraically equivalent to the two-indicator latent common
  cause (one loading fixed to 1) used by the Python `semopy`
  implementation, but is expressed with lavaan's standard
  residual-covariance idiom rather than an explicit latent variable.

- **Numerical values will not match `semopy` exactly**: `lavaan` and
  `semopy` use different default estimators/options. The fit measures
  returned are the same statistics, but exact numbers can differ.

- **Convention**: `adjacency_matrix` follows the lingamr convention
  `B[i, j]` = causal coefficient from variable j to variable i (j -\>
  i).

## References

Rosseel, Y. (2012). lavaan: An R Package for Structural Equation
Modeling. Journal of Statistical Software, 48(2), 1-36.
[doi:10.18637/jss.v048.i02](https://doi.org/10.18637/jss.v048.i02)

## Examples

``` r
if (requireNamespace("lavaan", quietly = TRUE)) {
  dat <- generate_lingam_sample_6()
  result <- lingam_direct(dat$data, reg_method = "ols")
  evaluate_model_fit(result, dat$data)
}
#>   DoF DoF Baseline chi2 chi2 p-value chi2 Baseline CFI GFI AGFI NFI TLI RMSEA
#> 1   0           15    0           NA       23023.7   1   1    1   1   1     0
#>        AIC      BIC    LogLik
#> 1 1860.598 1958.753 -910.2991
```
