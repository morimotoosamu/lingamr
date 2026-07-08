# Local log-likelihood of continuous component i (super-Gaussian disturbance)

Direct port of `likelihood_i`/`log_p_super_gaussian`/`variance_i` from
`lingam/utils/__init__.py` (lines 866-939 of the Python source).

## Usage

``` r
lim_likelihood_i(X, i, b_full, b0)
```

## Details

A (near-)zero residual variance means component i is fit
deterministically by its parents; `-n * log(sqrt(var_i))` would then
diverge to `Inf` and make this candidate model look artificially
superior. Such candidates are instead scored with a large finite penalty
so they lose every BIC comparison in
[`lim_bic_loss()`](https://morimotoosamu.github.io/lingamr/reference/lim_bic_loss.md)/[`lim_local_search()`](https://morimotoosamu.github.io/lingamr/reference/lim_local_search.md)
without introducing `Inf`/`NaN` into the running score.
