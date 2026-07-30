# Filter VARMA residuals by deterministic recursion

Computes `n_t = x_t - c - sum Phi_tau x_{t-tau} - sum Theta_w n_{t-w}`
recursively. The first `max(p, q)` residuals are initialized with zeros
(the Python reference draws them from a standard normal, which makes the
fit non-deterministic; the transient effect dies out at the same rate
either way).

## Usage

``` r
filter_varma_residuals(X, phis, thetas, const = NULL)
```

## Arguments

- X:

  numeric matrix (n_samples x n_features), rows ordered in time

- phis:

  AR coefficient array (p, m, m)

- thetas:

  MA coefficient array (q, m, m)

- const:

  intercept vector (NULL = zero)

## Value

full-length residual matrix (n_samples x n_features); the first
`max(p, q)` rows are zero
