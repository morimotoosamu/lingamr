# Global (NOTEARS-style) optimization phase: augmented-Lagrangian outer loop with L-BFGS-B inner solves. Returns the thresholded W (i -\> j orientation).

Global (NOTEARS-style) optimization phase: augmented-Lagrangian outer
loop with L-BFGS-B inner solves. Returns the thresholded W (i -\> j
orientation).

## Usage

``` r
lim_global_optimize(X, con, lambda1, max_iter, h_tol, rho_max, w_threshold)
```
