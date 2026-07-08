# Kernel-based mutual information: low-rank core

Low-rank counterpart of
[`kernel_mi_core()`](https://morimotoosamu.github.io/lingamr/reference/kernel_mi_core.md).
With `K2 ~= G2 %*% t(G2)`, both `S = tmp2^2 - t(W) %*% W` and `tmp2^2`
reduce to `c0^2*I + G2 %*% C %*% t(G2)` for a d2 x d2 matrix `C`, via
the Woodbury identity on `W` and the matrix determinant lemma on the
resulting rank-d2 update. The `n*log(c0)` terms in `logdet(S)` and
`logdet(tmp2^2)` cancel in their difference, so only d2 x d2 matrices
remain.

## Usage

``` r
kernel_mi_core_lowrank(prep1, x2, kappa, sigma)
```

## Arguments

- prep1:

  Output of
  [`kernel_mi_prepare_lowrank()`](https://morimotoosamu.github.io/lingamr/reference/kernel_mi_prepare_lowrank.md)
  for variable 1

- x2:

  Vector of variable 2

- kappa:

  Regularization parameter

- sigma:

  Width of the Gaussian kernel

## Value

Mutual information
