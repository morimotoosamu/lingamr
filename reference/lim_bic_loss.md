# BIC-penalized DAG score (sign-flipped, i.e. a "loss") used during the local search phase. `W` is interpreted as a 0/1 skeleton in the i -\> j orientation (parents of j are the nonzero rows of column j).

BIC-penalized DAG score (sign-flipped, i.e. a "loss") used during the
local search phase. `W` is interpreted as a 0/1 skeleton in the i -\> j
orientation (parents of j are the nonzero rows of column j).

## Usage

``` r
lim_bic_loss(W, X, is_continuous)
```
