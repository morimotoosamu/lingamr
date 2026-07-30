# Plot the causal graph of a ResitResult with ggplot2

Draws the estimated causal structure as a ggplot2-based directed graph,
like
[`autoplot.LingamResult()`](https://morimotoosamu.github.io/lingamr/reference/autoplot.LingamResult.md).
Because RESIT returns 0/1 edge indicators rather than coefficients, edge
labels are hidden by default (`label_edges = FALSE`); a constant "1" on
every edge carries no information.

## Usage

``` r
# S3 method for class 'ResitResult'
autoplot(
  object,
  threshold = 0,
  node_size = 16,
  node_color = "lightblue",
  label_edges = FALSE,
  label_pos = 0.35,
  ...
)
```

## Arguments

- object:

  Return value of
  [`lingam_resit()`](https://morimotoosamu.github.io/lingamr/reference/lingam_resit.md)
  (a `ResitResult` object)

- threshold:

  Coefficients with an absolute value at or below this are not treated
  as edges (default: 0)

- node_size:

  Node size (default: 16)

- node_color:

  Node fill color (default: "lightblue")

- label_edges:

  Whether to display edge labels (default: FALSE, since every edge would
  be labeled "1")

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
    requireNamespace("igraph", quietly = TRUE) &&
    requireNamespace("mgcv", quietly = TRUE)) {
  library(ggplot2)
  nonlinear <- generate_resit_sample(n = 300, seed = 1)
  model <- lingam_resit(nonlinear$data)
  autoplot(model)
}

# }
```
