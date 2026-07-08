# Kernel-based mutual information: low-rank precomputation for variable 1

Low-rank counterpart of
[`kernel_mi_prepare()`](https://morimotoosamu.github.io/lingamr/reference/kernel_mi_prepare.md).
Via the Woodbury identity, `E1 = tmp1^-1 K1` collapses to
`G1 %*% M1^-1 %*% t(G1)` where `M1 = c0*I + t(G1) %*% G1`, so the n x n
matrix `E1` never needs to be formed; only the d x d inverse `M1^-1`
does.

## Usage

``` r
kernel_mi_prepare_lowrank(x, kappa, sigma)
```

## Arguments

- x:

  Vector of variable 1

- kappa:

  Regularization parameter

- sigma:

  Width of the Gaussian kernel

## Value

list(G, Minv, A, c0) describing `E1` in factored form
