# Gaussian Gram matrix and its double-centered version

Multivariate input is combined into a single Gaussian kernel over the
row-wise squared Euclidean distances (upstream behavior), not treated
column by column.

## Usage

``` r
hsic_gram_matrix(x, width)
```

## Arguments

- x:

  numeric vector or matrix (n x d)

- width:

  kernel width from
  [`hsic_kernel_width()`](https://morimotoosamu.github.io/lingamr/reference/hsic_kernel_width.md)

## Value

list(K = Gram matrix, Kc = centered Gram matrix)
