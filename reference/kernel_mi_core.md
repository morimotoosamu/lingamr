# Kernel-based mutual information: core

The target quantity is the difference of logdets of 2n x 2n matrices,
but via the block structure and the Schur complement it can be computed
equivalently using only an n x n Cholesky decomposition:
`MI = -1/2 * (logdet(tmp2^2 - K2 K1 tmp1^-2 K1 K2) - logdet(tmp2^2))`

## Usage

``` r
kernel_mi_core(E1, x2, kappa, sigma)
```

## Arguments

- E1:

  Variable-1 matrix precomputed by
  [`kernel_mi_prepare()`](https://morimotoosamu.github.io/lingamr/reference/kernel_mi_prepare.md)

- x2:

  Vector of variable 2

- kappa:

  Regularization parameter

- sigma:

  Width of the Gaussian kernel

## Value

Mutual information
