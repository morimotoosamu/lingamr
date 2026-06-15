# Plot a causal graph from an adjacency matrix with DiagrammeR

Plot a causal graph from an adjacency matrix with DiagrammeR

## Usage

``` r
plot_adjacency(
  B,
  labels = NULL,
  threshold = 0,
  rankdir = "TB",
  title = "Estimated Causal Structure",
  shape = "circle",
  fillcolor = "lightyellow",
  bordercolor = "black",
  fontsize_node = 14,
  fontsize_edge = 10,
  edge_color = "gray40",
  edge_label_color = "red",
  true_B = NULL,
  color_tp = "forestgreen",
  color_fp = "firebrick",
  color_fn = "darkorange",
  debug = FALSE
)
```

## Arguments

- B:

  Adjacency matrix (n_features x n_features). **Convention: `B[i, j]` is
  the causal coefficient from variable j to variable i (j -\> i).** The
  `adjacency_matrix` from
  [`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md)
  can be passed directly.

- labels:

  Vector of variable names (if NULL, x0, x1, ... are generated
  automatically)

- threshold:

  Minimum absolute coefficient value to display (default: 0)

- rankdir:

  Layout direction (default: "LR") "LR" = left -\> right, "RL" = right
  -\> left, "TB" = top -\> bottom, "BT" = bottom -\> top

- title:

  Graph title (default: "Estimated Causal Structure")

- shape:

  Node shape (default: "circle") "circle", "box", "ellipse", "diamond",
  "plaintext", "square", "triangle", "hexagon", "octagon", etc.

- fillcolor:

  Node fill color (default: "lightyellow")

- bordercolor:

  Border color

- fontsize_node:

  Node font size (default: 14)

- fontsize_edge:

  Edge label font size (default: 10)

- edge_color:

  Edge color (default: "gray40"). Unused when `true_B` is specified.

- edge_label_color:

  Edge label color (default: "red"). Unused when `true_B` is specified.

- true_B:

  True adjacency matrix (may be NULL). When specified, edges are
  classified into three colors:

  - Correct edges (estimated and true): solid line in `color_tp`

  - False positives (estimated but not true): solid line in `color_fp`

  - Missed edges (not estimated but true): dashed line in `color_fn`
    (showing the true coefficient)

- color_tp:

  Color for correct edges (default: "forestgreen")

- color_fp:

  Color for false-positive edges (default: "firebrick")

- color_fn:

  Color for missed edges (default: "darkorange")

- debug:

  Enable debug mode (logical)

## Value

A grViz object (when DiagrammeR is available)

## Examples

``` r
LiNGAM_sample_1000 <- generate_lingam_sample_6()

LiNGAM_sample_1000$true_adjacency |>
  plot_adjacency(title = "True Causal Structure")

{"x":{"diagram":"digraph estimated_structure {\n\n  graph [rankdir = \"TB\",\n         label = \"True Causal Structure\",\n         labelloc = \"t\",\n         fontname = \"Helvetica-Bold\",\n         fontsize = 14]\n\n  node [shape = \"circle\",\n        style = \"solid,filled\",\n        fillcolor = \"#FFFFE0\",\n        color = \"#000000\",\n        fontname = \"Helvetica\",\n        fontsize = 14,\n        width = 0.6]\n\n  edge [fontname = \"Helvetica\",\n        fontsize = 10,\n        fontcolor = \"#FF0000\",\n        color = \"#666666\"]\n\n  x3 -> x0 [label = \" 3.00\"]\n  x0 -> x1 [label = \" 3.00\"]\n  x2 -> x1 [label = \" 2.00\"]\n  x3 -> x2 [label = \" 6.00\"]\n  x0 -> x4 [label = \" 8.00\"]\n  x2 -> x4 [label = \" -1.00\"]\n  x0 -> x5 [label = \" 4.00\"]\n}\n","config":{"engine":"dot","options":null}},"evals":[],"jsHooks":[]}
model <- LiNGAM_sample_1000$data |>
  lingam_direct()

model$adjacency_matrix |>
  plot_adjacency()

{"x":{"diagram":"digraph estimated_structure {\n\n  graph [rankdir = \"TB\",\n         label = \"Estimated Causal Structure\",\n         labelloc = \"t\",\n         fontname = \"Helvetica-Bold\",\n         fontsize = 14]\n\n  node [shape = \"circle\",\n        style = \"solid,filled\",\n        fillcolor = \"#FFFFE0\",\n        color = \"#000000\",\n        fontname = \"Helvetica\",\n        fontsize = 14,\n        width = 0.6]\n\n  edge [fontname = \"Helvetica\",\n        fontsize = 10,\n        fontcolor = \"#FF0000\",\n        color = \"#666666\"]\n\n  x3 -> x0 [label = \" 3.03\"]\n  x0 -> x1 [label = \" 2.99\"]\n  x2 -> x1 [label = \" 2.00\"]\n  x3 -> x2 [label = \" 5.99\"]\n  x0 -> x4 [label = \" 8.02\"]\n  x2 -> x4 [label = \" -1.01\"]\n  x0 -> x5 [label = \" 4.02\"]\n}\n","config":{"engine":"dot","options":null}},"evals":[],"jsHooks":[]}
# \donttest{
# Compare with the true structure
# (correct = green, false positive = red, missed = orange dashed)
model$adjacency_matrix |>
  plot_adjacency(true_B = LiNGAM_sample_1000$true_adjacency)

{"x":{"diagram":"digraph estimated_structure {\n\n  graph [rankdir = \"TB\",\n         label = \"Estimated Causal Structure\",\n         labelloc = \"t\",\n         fontname = \"Helvetica-Bold\",\n         fontsize = 14]\n\n  node [shape = \"circle\",\n        style = \"solid,filled\",\n        fillcolor = \"#FFFFE0\",\n        color = \"#000000\",\n        fontname = \"Helvetica\",\n        fontsize = 14,\n        width = 0.6]\n\n  edge [fontname = \"Helvetica\",\n        fontsize = 10,\n        fontcolor = \"#888888\",\n        color = \"#888888\"]\n\n  x3 -> x0 [label = \" 3.03\", color = \"#228B22\", fontcolor = \"#228B22\"]\n  x0 -> x1 [label = \" 2.99\", color = \"#228B22\", fontcolor = \"#228B22\"]\n  x2 -> x1 [label = \" 2.00\", color = \"#228B22\", fontcolor = \"#228B22\"]\n  x3 -> x2 [label = \" 5.99\", color = \"#228B22\", fontcolor = \"#228B22\"]\n  x0 -> x4 [label = \" 8.02\", color = \"#228B22\", fontcolor = \"#228B22\"]\n  x2 -> x4 [label = \" -1.01\", color = \"#228B22\", fontcolor = \"#228B22\"]\n  x0 -> x5 [label = \" 4.02\", color = \"#228B22\", fontcolor = \"#228B22\"]\n}\n","config":{"engine":"dot","options":null}},"evals":[],"jsHooks":[]}# }
```
