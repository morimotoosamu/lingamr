# Plot the causal graph of a LingamResult with ggplot2

Draws the estimated causal structure as a ggplot2-based directed graph.
Node positions are computed with igraph's hierarchical layout
(sugiyama), so the causal flow is generally arranged from top to bottom.
Because the output is a static image, it is stable in RMarkdown /
Quarto. If you need an interactive HTML figure, use
[`plot_adjacency()`](https://morimotoosamu.github.io/lingamr/reference/plot_adjacency.md)
(DiagrammeR-based).

## Usage

``` r
# S3 method for class 'LingamResult'
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
  [`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md)
  (a `LingamResult` object)

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

## Details

[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html) is
a ggplot2 generic, so you must load ggplot2 with
[`library(ggplot2)`](https://ggplot2.tidyverse.org) before using it.
Plotting requires ggplot2 and igraph.

## Examples

``` r
# \donttest{
if (requireNamespace("ggplot2", quietly = TRUE) &&
    requireNamespace("igraph", quietly = TRUE)) {
  library(ggplot2)
  dat <- generate_lingam_sample_6()
  model <- lingam_direct(dat$data, reg_method = "ols")
  autoplot(model)
}

# }
```
