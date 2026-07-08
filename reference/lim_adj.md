# Convert the doubled parameter vector (length 2*d*d) into a d x d matrix W

Mirrors Python's `w[:d*d].reshape([d, d])` (row-major reshape), so `w`
and the reverse operation
[`lim_flatten_rowmajor()`](https://morimotoosamu.github.io/lingamr/reference/lim_flatten_rowmajor.md)
must stay consistent.

## Usage

``` r
lim_adj(w, d)
```
