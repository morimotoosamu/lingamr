# Plot the causal graph of a ParceLingamResult with ggplot2

Draws the estimated causal structure as a ggplot2-based directed graph,
like
[`autoplot.LingamResult()`](https://morimotoosamu.github.io/lingamr/reference/autoplot.LingamResult.md).
Variable pairs whose adjacency-matrix entries are `NA` (unresolved order
/ suspected latent confounding) are drawn as dashed, arrowless segments.

## Usage

``` r
# S3 method for class 'ParceLingamResult'
autoplot(
  object,
  threshold = 0,
  node_size = 16,
  node_color = "lightblue",
  label_edges = TRUE,
  label_pos = 0.35,
  ...
)
```

## Arguments

- object:

  Return value of
  [`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md)
  (a `ParceLingamResult` object)

- threshold:

  Coefficients with an absolute value at or below this are not treated
  as edges (default: 0)

- node_size:

  Node size (default: 16)

- node_color:

  Node fill color (default: "lightblue")

- label_edges:

  Whether to display coefficient labels on edges (default: TRUE)

- label_pos:

  Position of each coefficient label along its edge, as a fraction from
  the source (0) to the target (1). The default 0.35 places labels
  off-center (toward the source) so labels on crossing edges do not
  overlap near the midpoint.

- ...:

  Unused

## Value

A ggplot object

## Examples

``` r
# \donttest{
if (requireNamespace("ggplot2", quietly = TRUE) &&
    requireNamespace("igraph", quietly = TRUE)) {
  library(ggplot2)
  dat <- generate_parce_sample(n = 500, seed = 42)
  model <- lingam_parce(dat$data)
  autoplot(model)
}

# }
```
