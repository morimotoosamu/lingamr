# Draw bootstrap probabilities with DiagrammeR

Draw bootstrap probabilities with DiagrammeR

## Usage

``` r
plot_bootstrap_probabilities(
  result,
  labels = NULL,
  min_causal_effect = NULL,
  min_probability = 0.5,
  rankdir = "TB",
  shape = "circle"
)
```

## Arguments

- result:

  BootstrapResult object

- labels:

  Vector of variable names (NULL allowed)

- min_causal_effect:

  Minimum causal effect to display

- min_probability:

  Minimum probability to display

- rankdir:

  Layout direction

- shape:

  Node shape

## Value

grViz object

## Examples

``` r
LiNGAM_sample_1000 <- generate_lingam_sample_6()

bs_model <- lingam_direct_bootstrap(LiNGAM_sample_1000$data, n_sampling = 30L, seed = 42)
#> Bootstrap: 30 iterations, method=adaptive_lasso (sequential)
#>   iteration 1 / 30
#>   iteration 10 / 30
#>   iteration 20 / 30
#>   iteration 30 / 30
#> Completed in 0.7 seconds.
plot_bootstrap_probabilities(bs_model)

{"x":{"diagram":"digraph bootstrap_result {\n  graph [rankdir = TB, fontsize = 14,\n         label = \"Bootstrap Probabilities\",\n         labelloc = t, fontname = \"Helvetica-Bold\"]\n  node [shape = circle, style = filled, fillcolor = lightyellow,\n        fontname = Helvetica, fontsize = 14, width = 0.6]\n  edge [fontname = Helvetica, fontsize = 10, fontcolor = blue, color = gray40]\n\n  x3 -> x0 [label = \" 0.97\", penwidth = 3.4]\n  x0 -> x1 [label = \" 0.97\", penwidth = 3.4]\n  x2 -> x1 [label = \" 0.97\", penwidth = 3.4]\n  x3 -> x2 [label = \" 0.97\", penwidth = 3.4]\n  x0 -> x4 [label = \" 0.97\", penwidth = 3.4]\n  x2 -> x4 [label = \" 0.97\", penwidth = 3.4]\n  x0 -> x5 [label = \" 1.00\", penwidth = 3.5]\n}\n","config":{"engine":"dot","options":null}},"evals":[],"jsHooks":[]}
```
