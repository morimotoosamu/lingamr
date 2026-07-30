# LiM: LiNGAM for Mixed Data

Estimates a causal structure from data containing a mixture of
continuous and discrete variables, following Zeng et al. (2022). The
method combines a NOTEARS-style continuous optimization (global phase)
with a combinatorial local search over edge directions, pruning, and
edge addition. Discrete variables are binary (0/1) by default; set
`is_poisson = TRUE` to treat them as Poisson-distributed counts instead.

## Usage

``` r
lingam_lim(
  X,
  is_continuous,
  lambda1 = 0.1,
  max_iter = 150L,
  h_tol = 1e-08,
  rho_max = 1e+16,
  w_threshold = 0.1,
  only_global = FALSE,
  is_poisson = FALSE
)
```

## Arguments

- X:

  Numeric matrix (n_samples x n_features) or data frame

- is_continuous:

  Logical vector of length `ncol(X)`. `TRUE` marks a continuous
  variable, `FALSE` marks a discrete variable (binary 0/1 by default;
  non-negative integer counts when `is_poisson = TRUE`).

- lambda1:

  L1 penalty parameter (default: 0.1)

- max_iter:

  Maximum number of dual ascent (outer loop) steps (default: 150)

- h_tol:

  Tolerance for the acyclicity constraint h(W) (default: 1e-8)

- rho_max:

  Maximum value of the augmented-Lagrangian penalty rho (default: 1e16)

- w_threshold:

  Edges with \|weight\| below this value are dropped after the global
  optimization phase (default: 0.1)

- only_global:

  If `TRUE`, skip the combinatorial local search phase and return the
  thresholded global-optimization result directly (default: FALSE)

- is_poisson:

  If `TRUE`, all discrete variables (`is_continuous = FALSE`) are
  treated as Poisson-distributed counts (non-negative integers) rather
  than binary variables, and the local search phase scores them with
  Poisson regression log-likelihoods (default: FALSE). See Details.

## Value

A `LiMResult` object (list) containing the following elements:

- `adjacency_matrix`: adjacency matrix B (n_features x n_features).
  **Convention: `B[i, j]` is the causal coefficient from variable j to
  variable i (j -\> i)**, i.e. the same convention as
  [`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md).
  Zero elements indicate no causal relationship.

- `causal_order`: estimated causal order (integer vector of 1-based
  indices), derived from `adjacency_matrix` via a topological sort. `NA`
  (with a warning) if the estimated matrix is not acyclic.

- `is_continuous`: the input `is_continuous` vector, stored for
  reference.

- `is_poisson`: the input `is_poisson` flag, stored for reference.

## Details

By default, discrete variables must be binary (0/1). With
`is_poisson = TRUE`, discrete variables are instead treated as
Poisson-distributed counts: the local search phase scores each discrete
variable with a Poisson regression log-likelihood (an unregularized
`glm(family = poisson())` fit on all parents, with an intercept), or
with the closed-form intercept-only Poisson maximum likelihood when the
variable has no parents. The global optimization phase is unchanged in
both modes and keeps the logistic surrogate loss for discrete columns;
this matches the behavior of the Python implementation's
`fit(is_poisson=True)`, whose `loss_type = "poisson"` option (a
whole-matrix Poisson loss that cannot be combined with mixed data) is
intentionally not ported. Consequently, with `only_global = TRUE` the
local search is skipped and `is_poisson` only affects input validation.

The Poisson scoring deliberately deviates from the Python
implementation, which fits a separate univariate regression per parent
(accumulating only the last parent's likelihood due to a loop bug),
ignores the fitted intercept, uses scikit-learn's `PoissonRegressor`
with its default L2 regularization, and scores parentless count
variables with the Bernoulli frequency-table code (a multinomial
likelihood on a different scale). This implementation uses the full
multivariate maximum-likelihood fit instead; numeric results will
therefore not match the Python implementation.

The Python implementation's `adjacency_matrix_` uses the opposite
convention (`W[i, j]` = i -\> j). This R implementation transposes the
internal result so that `adjacency_matrix` follows the lingamr
convention (`B[i, j]` = j -\> i), consistent with
[`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md).

In the local search phase, edges that are reversed or newly added are
assigned a weight of exactly 1 rather than a re-estimated coefficient
(this matches the original Python implementation). As a result, non-zero
entries of `adjacency_matrix` are a mix of global-phase estimated
coefficients and local-phase placeholder weights of 1.

The local search's BIC score for discrete variables with parents uses
R's `glm(family = binomial())`, an unregularized maximum-likelihood fit.
The Python implementation uses scikit-learn's `LogisticRegression`,
which applies L2 regularization by default; consequently, numeric
results will not exactly match the Python implementation.

## References

Zeng Y, Shimizu S, Matsui H, Sun F. Causal discovery for linear mixed
data. In: Proceedings of the First Conference on Causal Learning and
Reasoning (CLeaR 2022). PMLR 177, pp. 994-1009, 2022.

## Examples

``` r
# Reproducibility requires set.seed(), since the optimization starts from
# a random initial point.
set.seed(1)
dat <- generate_lim_sample(n = 300)
result <- lingam_lim(dat$data, is_continuous = dat$is_continuous)
print(result)
#> LiM Result
#>   Variables : 3
#>   Variable types: continuous, discrete, continuous
#>   Causal order: x1 -> x2 -> x3
#> 
#> Adjacency matrix (row = to, col = from):
#>       x1 x2 x3
#> x1  0.00  0  0
#> x2 -0.22  0  0
#> x3  0.00  1  0
```
