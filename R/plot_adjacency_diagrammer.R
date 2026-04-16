# =============================================================================
# Direct LiNGAM - R Implementation
# Based on the Python implementation from the LiNGAM Project
# https://sites.google.com/view/sshimizu06/lingam
# https://github.com/cdt15/lingam
# =============================================================================


#' 隣接行列から DiagrammeR で因果グラフを描画
#'
#' @param B 隣接行列
#' @param labels 変数名ベクトル (NULL の場合は x0, x1, ... を自動生成)
#' @param threshold 表示する最小係数の絶対値 (default: 0.01)
#' @param rankdir レイアウト方向 (default: "LR")
#'   "LR" = 左→右, "RL" = 右→左, "TB" = 上→下, "BT" = 下→上
#' @param graph_label グラフのタイトル (default: "Estimated Causal Structure")
#' @param shape ノードの形状 (default: "circle")
#'   "circle", "box", "ellipse", "diamond", "plaintext",
#'   "square", "triangle", "hexagon", "octagon" など
#' @param fillcolor ノードの塗りつぶし色 (default: "lightyellow")
#' @param fontsize_node ノードのフォントサイズ (default: 14)
#' @param fontsize_edge エッジラベルのフォントサイズ (default: 10)
#' @param edge_color エッジの色 (default: "gray40")
#' @param edge_label_color エッジラベルの色 (default: "red")
#' @return grViz オブジェクト（DiagrammeR が利用可能な場合）
plot_adjacency_diagrammer <- function(B,
                                      labels = NULL,
                                      threshold = 0.01,
                                      rankdir = "LR",
                                      graph_label = "Estimated Causal Structure",
                                      shape = "circle",
                                      fillcolor = "lightyellow",
                                      fontsize_node = 14,
                                      fontsize_edge = 10,
                                      edge_color = "gray40",
                                      edge_label_color = "red") {
  # --- DiagrammeR パッケージの確認 ---
  if (!requireNamespace("DiagrammeR", quietly = TRUE)) {
    message("============================================================")
    message("DiagrammeR パッケージがインストールされていません。")
    message("以下のコマンドでインストールしてください：")
    message("")
    message('  install.packages("DiagrammeR")')
    message("")
    message("インストール後、再度この関数を実行してください。")
    message("============================================================")
    return(invisible(NULL))
  }

  # --- 引数のバリデーション ---
  valid_rankdir <- c("LR", "RL", "TB", "BT")
  if (!(rankdir %in% valid_rankdir)) {
    stop(sprintf(
      "rankdir は %s のいずれかを指定してください。",
      paste(valid_rankdir, collapse = ", ")
    ))
  }

  valid_shapes <- c(
    "circle", "box", "ellipse", "diamond", "plaintext",
    "square", "triangle", "hexagon", "octagon",
    "doublecircle", "rect", "oval", "egg",
    "pentagon", "star", "none"
  )
  if (!(shape %in% valid_shapes)) {
    stop(sprintf(
      "shape は %s のいずれかを指定してください。",
      paste(valid_shapes, collapse = ", ")
    ))
  }

  # --- 引数の処理 ---
  p <- ncol(B)
  if (is.null(labels)) labels <- paste0("x", seq_len(p) - 1)

  # --- エッジ記述の生成 ---
  edge_lines <- c()
  for (i in 1:p) {
    for (j in 1:p) {
      if (abs(B[i, j]) > threshold) {
        coef_label <- sprintf("%.2f", B[i, j])
        edge_lines <- c(
          edge_lines,
          sprintf("  %s -> %s [label = ' %s']", labels[j], labels[i], coef_label)
        )
      }
    }
  }

  if (length(edge_lines) == 0) {
    message("閾値を超えるエッジが見つかりませんでした。threshold を小さくしてください。")
    return(invisible(NULL))
  }

  # --- DOT 記述の生成 ---
  dot <- paste0(
    "digraph estimated_structure {\n",
    sprintf("  graph [rankdir = %s, fontsize = 14,\n", rankdir),
    sprintf("         label = '%s',\n", graph_label),
    "         labelloc = t, fontname = 'Helvetica-Bold']\n",
    sprintf(
      "  node [shape = %s, style = filled, fillcolor = %s,\n",
      shape, fillcolor
    ),
    sprintf(
      "        fontname = Helvetica, fontsize = %d, width = 0.6]\n",
      fontsize_node
    ),
    sprintf(
      "  edge [fontname = Helvetica, fontsize = %d, fontcolor = %s, color = %s]\n\n",
      fontsize_edge, edge_label_color, edge_color
    ),
    paste(edge_lines, collapse = "\n"), "\n",
    "}\n"
  )

  # --- 描画 ---
  DiagrammeR::grViz(dot)
}
