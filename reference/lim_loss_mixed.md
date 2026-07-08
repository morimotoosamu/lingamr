# Mixed loss (logistic for discrete columns, log-cosh/Laplace for continuous columns) and its gradient. `con` is a length-d vector, 1 = continuous, 0 = discrete. `W_dis_mask` / `W_con_mask` are precomputed d x d 0/1 masks.

Mixed loss (logistic for discrete columns, log-cosh/Laplace for
continuous columns) and its gradient. `con` is a length-d vector, 1 =
continuous, 0 = discrete. `W_dis_mask` / `W_con_mask` are precomputed d
x d 0/1 masks.

## Usage

``` r
lim_loss_mixed(W, X, con, W_dis_mask, W_con_mask)
```
