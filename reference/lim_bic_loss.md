# BIC-penalized DAG score (sign-flipped, i.e. a "loss") used during the local search phase. `W` is interpreted as a 0/1 skeleton in the i -\> j orientation (parents of j are the nonzero rows of column j). With `is_poisson = TRUE`, discrete columns are scored as Poisson counts instead of Bernoulli variables (see the Details of [`lingam_lim()`](https://morimotoosamu.github.io/lingamr/reference/lingam_lim.md) for the documented deviations from the Python source).

BIC-penalized DAG score (sign-flipped, i.e. a "loss") used during the
local search phase. `W` is interpreted as a 0/1 skeleton in the i -\> j
orientation (parents of j are the nonzero rows of column j). With
`is_poisson = TRUE`, discrete columns are scored as Poisson counts
instead of Bernoulli variables (see the Details of
[`lingam_lim()`](https://morimotoosamu.github.io/lingamr/reference/lingam_lim.md)
for the documented deviations from the Python source).

## Usage

``` r
lim_bic_loss(W, X, is_continuous, is_poisson)
```
