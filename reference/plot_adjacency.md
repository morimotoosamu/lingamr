# 隣接行列から DiagrammeR で因果グラフを描画

隣接行列から DiagrammeR で因果グラフを描画

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

  隣接行列 (n_features x n_features)。 **規則: `B[i, j]` は変数 j
  から変数 i への因果係数（j → i）。**
  [`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md)
  の `adjacency_matrix` をそのまま渡せる。

- labels:

  変数名ベクトル (NULL の場合は x0, x1, ... を自動生成)

- threshold:

  表示する最小係数の絶対値 (default: 0)

- rankdir:

  レイアウト方向 (default: "LR") "LR" = 左→右, "RL" = 右→左, "TB" =
  上→下, "BT" = 下→上

- title:

  グラフのタイトル (default: "Estimated Causal Structure")

- shape:

  ノードの形状 (default: "circle") "circle", "box", "ellipse",
  "diamond", "plaintext", "square", "triangle", "hexagon", "octagon"
  など

- fillcolor:

  ノードの塗りつぶし色 (default: "lightyellow")

- bordercolor:

  枠線の色

- fontsize_node:

  ノードのフォントサイズ (default: 14)

- fontsize_edge:

  エッジラベルのフォントサイズ (default: 10)

- edge_color:

  エッジの色 (default: "gray40")。`true_B` 指定時は未使用。

- edge_label_color:

  エッジラベルの色 (default: "red")。`true_B` 指定時は未使用。

- true_B:

  真の隣接行列 (NULL 可)。指定するとエッジを3色で分類する：

  - 正解エッジ（推定あり・真あり）: `color_tp` の実線

  - 過検出（推定あり・真なし）: `color_fp` の実線

  - 見逃し（推定なし・真あり）: `color_fn` の破線（真の係数を表示）

- color_tp:

  正解エッジの色 (default: "forestgreen")

- color_fp:

  過検出エッジの色 (default: "firebrick")

- color_fn:

  見逃しエッジの色 (default: "darkorange")

- debug:

  デバッグモードの有効化 (logical)

## Value

grViz オブジェクト（DiagrammeR が利用可能な場合）

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
# 真の構造と比較（正解=緑, 過検出=赤, 見逃し=オレンジ破線）
model$adjacency_matrix |>
  plot_adjacency(true_B = LiNGAM_sample_1000$true_adjacency)

{"x":{"diagram":"digraph estimated_structure {\n\n  graph [rankdir = \"TB\",\n         label = \"Estimated Causal Structure\",\n         labelloc = \"t\",\n         fontname = \"Helvetica-Bold\",\n         fontsize = 14]\n\n  node [shape = \"circle\",\n        style = \"solid,filled\",\n        fillcolor = \"#FFFFE0\",\n        color = \"#000000\",\n        fontname = \"Helvetica\",\n        fontsize = 14,\n        width = 0.6]\n\n  edge [fontname = \"Helvetica\",\n        fontsize = 10,\n        fontcolor = \"#888888\",\n        color = \"#888888\"]\n\n  x3 -> x0 [label = \" 3.03\", color = \"#228B22\", fontcolor = \"#228B22\"]\n  x0 -> x1 [label = \" 2.99\", color = \"#228B22\", fontcolor = \"#228B22\"]\n  x2 -> x1 [label = \" 2.00\", color = \"#228B22\", fontcolor = \"#228B22\"]\n  x3 -> x2 [label = \" 5.99\", color = \"#228B22\", fontcolor = \"#228B22\"]\n  x0 -> x4 [label = \" 8.02\", color = \"#228B22\", fontcolor = \"#228B22\"]\n  x2 -> x4 [label = \" -1.01\", color = \"#228B22\", fontcolor = \"#228B22\"]\n  x0 -> x5 [label = \" 4.02\", color = \"#228B22\", fontcolor = \"#228B22\"]\n}\n","config":{"engine":"dot","options":null}},"evals":[],"jsHooks":[]}# }
```
