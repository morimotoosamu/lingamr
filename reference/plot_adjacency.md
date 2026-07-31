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

  Layout direction (default: "TB") "LR" = left -\> right, "RL" = right
  -\> left, "TB" = top -\> bottom, "BT" = bottom -\> top

- title:

  Graph title (default: "Estimated Causal Structure")

- shape:

  Node shape (default: "circle") "circle", "box", "ellipse", "diamond",
  "plaintext", "square", "triangle", "hexagon", "octagon", etc.

- fillcolor:

  Node fill color (default: "lightyellow")

- bordercolor:

  Border color (default: "black")

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

  Enable debug mode (default: FALSE)

## Value

A grViz object (when DiagrammeR is available)

## Examples

``` r
if (requireNamespace("DiagrammeR", quietly = TRUE)) {
  LiNGAM_sample_1000 <- generate_lingam_sample_6()

  LiNGAM_sample_1000$true_adjacency |>
    plot_adjacency(title = "True Causal Structure")

  model <- LiNGAM_sample_1000$data |>
    lingam_direct(reg_method = "ols")

  model$adjacency_matrix |>
    plot_adjacency()

  # \donttest{
  # Compare with the true structure
  # (correct = green, false positive = red, missed = orange dashed)
  model$adjacency_matrix |>
    plot_adjacency(true_B = LiNGAM_sample_1000$true_adjacency)
  # }
}

{"x":{"diagram":"digraph estimated_structure {\n\n  graph [rankdir = \"TB\",\n         label = \"Estimated Causal Structure\",\n         labelloc = \"t\",\n         fontname = \"Helvetica-Bold\",\n         fontsize = 14]\n\n  node [shape = \"circle\",\n        style = \"solid,filled\",\n        fillcolor = \"#FFFFE0\",\n        color = \"#000000\",\n        fontname = \"Helvetica\",\n        fontsize = 14,\n        width = 0.6]\n\n  edge [fontname = \"Helvetica\",\n        fontsize = 10,\n        fontcolor = \"#888888\",\n        color = \"#888888\"]\n\n  x0 -> x1 [label = \" 3.24\", color = \"#228B22\", fontcolor = \"#228B22\"]\n  x0 -> x4 [label = \" 7.99\", color = \"#228B22\", fontcolor = \"#228B22\"]\n  x0 -> x5 [label = \" 3.87\", color = \"#228B22\", fontcolor = \"#228B22\"]\n  x2 -> x1 [label = \" 1.97\", color = \"#228B22\", fontcolor = \"#228B22\"]\n  x2 -> x4 [label = \" -1.06\", color = \"#228B22\", fontcolor = \"#228B22\"]\n  x3 -> x0 [label = \" 3.27\", color = \"#228B22\", fontcolor = \"#228B22\"]\n  x3 -> x2 [label = \" 5.99\", color = \"#228B22\", fontcolor = \"#228B22\"]\n  x2 -> x0 [label = \" -0.04\", color = \"#B22222\", fontcolor = \"#B22222\"]\n  x2 -> x5 [label = \" 0.07\", color = \"#B22222\", fontcolor = \"#B22222\"]\n  x3 -> x1 [label = \" 0.01\", color = \"#B22222\", fontcolor = \"#B22222\"]\n  x3 -> x4 [label = \" 0.39\", color = \"#B22222\", fontcolor = \"#B22222\"]\n  x3 -> x5 [label = \" -0.31\", color = \"#B22222\", fontcolor = \"#B22222\"]\n  x4 -> x1 [label = \" -0.03\", color = \"#B22222\", fontcolor = \"#B22222\"]\n  x4 -> x5 [label = \" 0.02\", color = \"#B22222\", fontcolor = \"#B22222\"]\n  x5 -> x1 [label = \" 0.01\", color = \"#B22222\", fontcolor = \"#B22222\"]\n}\n","config":{"engine":"dot","options":null}},"evals":[],"jsHooks":[]}
```
