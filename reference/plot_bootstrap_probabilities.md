# ブートストラップ確率を DiagrammeR で描画

ブートストラップ確率を DiagrammeR で描画

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

  BootstrapResult オブジェクト

- labels:

  変数名ベクトル (NULL可)

- min_causal_effect:

  表示する最小因果効果

- min_probability:

  表示する最小確率

- rankdir:

  レイアウト方向

- shape:

  ノード形状

## Value

grViz オブジェクト

## Examples

``` r
LiNGAM_sample_1000 <- generate_lingam_sample_6()

bs_model <- lingam_direct_bootstrap(LiNGAM_sample_1000$data, n_sampling = 30L, seed = 42)
#> Bootstrap: 30 iterations, method=adaptive_lasso (sequential)
#>   iteration 1 / 30
#>   iteration 10 / 30
#>   iteration 20 / 30
#>   iteration 30 / 30
#> Completed in 2.0 seconds.
plot_bootstrap_probabilities(bs_model)

{"x":{"diagram":"digraph bootstrap_result {\n  graph [rankdir = TB, fontsize = 14,\n         label = \"Bootstrap Probabilities\",\n         labelloc = t, fontname = \"Helvetica-Bold\"]\n  node [shape = circle, style = filled, fillcolor = lightyellow,\n        fontname = Helvetica, fontsize = 14, width = 0.6]\n  edge [fontname = Helvetica, fontsize = 10, fontcolor = blue, color = gray40]\n\n  x3 -> x0 [label = \" 0.97\", penwidth = 3.4]\n  x0 -> x1 [label = \" 0.97\", penwidth = 3.4]\n  x2 -> x1 [label = \" 0.97\", penwidth = 3.4]\n  x3 -> x2 [label = \" 0.97\", penwidth = 3.4]\n  x0 -> x4 [label = \" 0.97\", penwidth = 3.4]\n  x2 -> x4 [label = \" 0.97\", penwidth = 3.4]\n  x0 -> x5 [label = \" 1.00\", penwidth = 3.5]\n}\n","config":{"engine":"dot","options":null}},"evals":[],"jsHooks":[]}
```
