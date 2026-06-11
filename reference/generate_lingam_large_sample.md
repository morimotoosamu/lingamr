# Generate large-scale sample data to benchmark Direct LiNGAM scalability

Generates a dataset with many variables to demonstrate the computational
scalability difference between Direct LiNGAM and ICA-LiNGAM.

## Usage

``` r
generate_lingam_large_sample(
  p = 20L,
  n = 1000L,
  max_parents = 3L,
  coef_min = 0.5,
  coef_max = 1.5,
  seed = 42L,
  noise_dist = "uniform"
)
```

## Arguments

- p:

  number of variables (default: 20)

- n:

  number of observations (default: 1000)

- max_parents:

  maximum number of parents per node (default: 3). Controls graph
  density. Each variable xi (i \>= 1) receives between 1 and
  `min(max_parents, i)` parents drawn from x0, ..., x(i-1).

- coef_min:

  minimum absolute value of edge coefficients (default: 0.5)

- coef_max:

  maximum absolute value of edge coefficients (default: 1.5)

- seed:

  random seed (default: 42)

- noise_dist:

  error term distribution. "uniform" : Uniform(0, 1) - default,
  non-Gaussian (LiNGAM works well) "gaussian" : Normal(0, 1) - LiNGAM
  may fail "lognormal" : Log-normal(0, 1) - skewed, non-Gaussian
  "exponential" : Exponential(1) - skewed, non-Gaussian "t3" :
  t-distribution (df=3) - heavy tails

## Value

A list with three elements:

- `data`: data.frame with `p` columns (x0, x1, ..., x(p-1)).

- `true_adjacency`: p x p matrix. `true_adjacency[i, j]` is the
  structural coefficient of the edge xj -\> xi (row = to, col = from).
  The matrix is strictly lower-triangular because variables are stored
  in causal order.

- `true_causal_order`: integer vector `0:(p-1)`. Variables are already
  in topological order, so the true causal order is always 0, 1, ...,
  p-1.

## Details

### Why Direct LiNGAM slows down with large p

At each of its `p` steps, Direct LiNGAM evaluates an independence
measure between every remaining candidate root and every other residual.
The total number of evaluations is:

\$\$\sum\_{k=1}^{p} k(k-1) \approx \frac{p^3}{3}\$\$

i.e., O(p^3). Each evaluation is itself O(n), giving O(p^3 n) overall.
For p = 10 this is about 330 evaluations; for p = 20 about 2,660; for p
= 40 about 21,320 — an 8x increase every time p doubles.

### Why ICA-LiNGAM scales better

ICA-LiNGAM applies FastICA once to the whole p x n data matrix. Each
FastICA iteration costs O(p^2 n), and the algorithm typically converges
in far fewer than p iterations. Additionally, these matrix operations
are fully vectorised (BLAS/LAPACK), whereas Direct LiNGAM iterates over
pairs in an R loop.

### Data-generating process

Variables are topologically ordered as x0, x1, ..., x(p-1). For each i
\>= 1, the number of parents is sampled uniformly from 1 to
`min(max_parents, i)`, and the parents are drawn without replacement
from x0, ..., x(i-1). Edge coefficients are drawn uniformly from
\[-coef_max, -coef_min\] U \[coef_min, coef_max\]. The resulting
adjacency matrix is strictly lower-triangular.

## Examples

``` r
# 20変数のデータを生成してスパース性を確認
dataset <- generate_lingam_large_sample(p = 20, n = 500)
dim(dataset$data)                    # 500 x 20
#> [1] 500  20
sum(dataset$true_adjacency != 0)     # 辺の本数
#> [1] 32
dataset$true_causal_order            # 0, 1, ..., 19
#>  [1]  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19

# \donttest{
# 変数数が増えると Direct LiNGAM の実行時間が急増する
t10 <- system.time(lingam_direct(generate_lingam_large_sample(p = 10)$data))
t20 <- system.time(lingam_direct(generate_lingam_large_sample(p = 20)$data))
cat(sprintf("p=10: %.1f sec,  p=20: %.1f sec\n", t10["elapsed"], t20["elapsed"]))
#> p=10: 0.0 sec,  p=20: 0.2 sec
# }
```
