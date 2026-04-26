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
#' @param title グラフのタイトル (default: "Estimated Causal Structure")
#' @param shape ノードの形状 (default: "circle")
#'   "circle", "box", "ellipse", "diamond", "plaintext",
#'   "square", "triangle", "hexagon", "octagon" など
#' @param fillcolor ノードの塗りつぶし色 (default: "lightyellow")
#' @param fontsize_node ノードのフォントサイズ (default: 14)
#' @param fontsize_edge エッジラベルのフォントサイズ (default: 10)
#' @param edge_color エッジの色 (default: "gray40")
#' @param edge_label_color エッジラベルの色 (default: "red")
#' @param bordercolor 枠線の色
#' @param debug デバッグモードの有効化 (logical)
#' @return grViz オブジェクト（DiagrammeR が利用可能な場合）
#' @importFrom grDevices col2rgb
#' @export
#' @examples
#' data(LiNGAM_sample_1000)
#'
#' model <- LiNGAM_sample_1000 |>
#'   direct_lingam()
#'
#' model$adjacency_matrix |>
#'   plot_adjacency_diagrammer()
plot_adjacency_diagrammer <- function(B,
                                      labels = NULL,
                                      threshold = 0.01,
                                      rankdir = "TB",
                                      title = "Estimated Causal Structure",
                                      shape = "circle",
                                      fillcolor = "lightyellow",
                                      bordercolor = "black",
                                      fontsize_node = 14,
                                      fontsize_edge = 10,
                                      edge_color = "gray40",
                                      edge_label_color = "red",
                                      debug = FALSE) {
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

  # --- 色名をカラーコードに変換するヘルパー ---
  to_hex <- function(color_str) {
    if (grepl("^#", color_str)) {
      return(color_str)
    } # 既にカラーコードならそのまま
    tryCatch(
      {
        rgb_val <- grDevices::col2rgb(color_str)
        sprintf("#%02X%02X%02X", rgb_val[1], rgb_val[2], rgb_val[3])
      },
      error = function(e) {
        warning(sprintf("Unknown color '%s', using black.", color_str))
        "#000000"
      }
    )
  }

  # --- 全ての色をカラーコードに変換 ---
  fillcolor <- to_hex(fillcolor)
  edge_color <- to_hex(edge_color)
  edge_label_color <- to_hex(edge_label_color)

  if (is.null(bordercolor)) {
    bordercolor <- fillcolor
  } else {
    bordercolor <- to_hex(bordercolor)
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

  # --- title のエスケープ ---
  safe_label <- gsub("'", "\\\\'", title)

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

  # --- DOT 記述の構築 ---
  dot <- paste0(
    "digraph estimated_structure {\n",
    "\n",
    "  graph [rankdir = '", rankdir, "',\n",
    "         label = '", safe_label, "',\n",
    "         labelloc = 't',\n",
    "         fontname = 'Helvetica-Bold',\n",
    "         fontsize = 14]\n",
    "\n",
    "  node [shape = '", shape, "',\n",
    "        style = 'solid,filled',\n",
    "        fillcolor = '", fillcolor, "',\n",
    "        color = '", bordercolor, "',\n",
    "        fontname = 'Helvetica',\n",
    "        fontsize = ", fontsize_node, ",\n",
    "        width = 0.6]\n",
    "\n",
    "  edge [fontname = 'Helvetica',\n",
    "        fontsize = ", fontsize_edge, ",\n",
    "        fontcolor = '", edge_label_color, "',\n",
    "        color = '", edge_color, "']\n",
    "\n",
    paste(edge_lines, collapse = "\n"), "\n",
    "}\n"
  )

  if (debug) {
    cat("=== Generated DOT ===\n")
    cat(dot)
    cat("=====================\n")
  }

  DiagrammeR::grViz(dot)
}
