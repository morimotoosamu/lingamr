# Kernel-based mutual information: precomputation for variable 1

Computes the matrix `E1 = tmp1^-1 K1` (`tmp1 = K1 + n*kappa/2 * I`) used
in
[`kernel_mi_core()`](https://morimotoosamu.github.io/lingamr/reference/kernel_mi_core.md).
It only needs to be called once per candidate variable, avoiding
per-pair recomputation.

## Usage

``` r
kernel_mi_prepare(x, kappa, sigma)
```

## Arguments

- x:

  Vector of variable 1

- kappa:

  Regularization parameter

- sigma:

  Width of the Gaussian kernel

## Value

Matrix E1 (n x n)
