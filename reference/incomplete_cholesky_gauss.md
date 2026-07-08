# Gaussian kernel via pivoted incomplete Cholesky

Low-rank approximation `G` (n x d) with `tcrossprod(G) ~= K`, where
`K[i, j] = exp(-(x[i] - x[j])^2 / (2*sigma^2))`, via the greedy pivoted
incomplete Cholesky of Fine & Scheinberg (2001) / Bach & Jordan (2002,
kernel-ICA). Each step picks the index with the largest diagonal
residual and computes only the kernel column for that pivot, so the full
n x n Gram matrix is never formed. The Gaussian kernel's diagonal is 1,
so residuals start at 1 without evaluating `K`.

## Usage

``` r
incomplete_cholesky_gauss(
  x,
  sigma,
  tol = 1e-04,
  max_rank = min(length(x), 200)
)
```

## Arguments

- x:

  Input vector (length n)

- sigma:

  Width of the Gaussian kernel

- tol:

  Stop once the largest remaining diagonal residual is at or below this

- max_rank:

  Upper bound on the approximation rank

## Value

n x d matrix `G` with `tcrossprod(G) ~= K`
