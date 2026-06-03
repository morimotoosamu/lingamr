# =============================================================================
# Direct LiNGAM - R Implementation
# Based on the Python implementation from the LiNGAM Project
# https://sites.google.com/view/sshimizu06/lingam
# https://github.com/cdt15/lingam
# =============================================================================


#' 隣接行列から DiagrammeR で因果グラフを描画
#'
#' @param B 隣接行列 (n_features x n_features)。
#'   **規則: `B[i, j]` は変数 j から変数 i への因果係数（j → i）。**
#'   `lingam_direct()` の `adjacency_matrix` をそのまま渡せる。
#' @param labels 変数名ベクトル (NULL の場合は x0, x1, ... を自動生成)
#' @param threshold 表示する最小係数の絶対値 (default: 0)
#' @param rankdir レイアウト方向 (default: "LR")
#'   "LR" = 左→右, "RL" = 右→左, "TB" = 上→下, "BT" = 下→上
#' @param title グラフのタイトル (default: "Estimated Causal Structure")
#' @param shape ノードの形状 (default: "circle")
#'   "circle", "box", "ellipse", "diamond", "plaintext",
#'   "square", "triangle", "hexagon", "octagon" など
#' @param fillcolor ノードの塗りつぶし色 (default: "lightyellow")
#' @param bordercolor 枠線の色
#' @param fontsize_node ノードのフォントサイズ (default: 14)
#' @param fontsize_edge エッジラベルのフォントサイズ (default: 10)
#' @param edge_color エッジの色 (default: "gray40")。`true_B` 指定時は未使用。
#' @param edge_label_color エッジラベルの色 (default: "red")。`true_B` 指定時は未使用。
#' @param true_B 真の隣接行列 (NULL 可)。指定するとエッジを3色で分類する：
#'   * 正解エッジ（推定あり・真あり）: `color_tp` の実線
#'   * 過検出（推定あり・真なし）: `color_fp` の実線
#'   * 見逃し（推定なし・真あり）: `color_fn` の破線（真の係数を表示）
#' @param color_tp 正解エッジの色 (default: "forestgreen")
#' @param color_fp 過検出エッジの色 (default: "firebrick")
#' @param color_fn 見逃しエッジの色 (default: "darkorange")
#' @param debug デバッグモードの有効化 (logical)
#' @return grViz オブジェクト（DiagrammeR が利用可能な場合）
#' @importFrom grDevices col2rgb
#' @export
#' @examples
#' LiNGAM_sample_1000 <- generate_lingam_sample_6()
#'
#' LiNGAM_sample_1000$true_adjacency |>
#'   plot_adjacency(title = "True Causal Structure")
#'
#' model <- LiNGAM_sample_1000$data |>
#'   lingam_direct()
#'
#' model$adjacency_matrix |>
#'   plot_adjacency()
#'
#' \donttest{
#' # 真の構造と比較（正解=緑, 過検出=赤, 見逃し=オレンジ破線）
#' model$adjacency_matrix |>
#'   plot_adjacency(true_B = LiNGAM_sample_1000$true_adjacency)
#' }
plot_adjacency <- function(B,
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
                           debug = FALSE) {
  # --- DiagrammeR パッケージの確認 ---
  if (!requireNamespace("DiagrammeR", quietly = TRUE)) {
    stop("Package 'DiagrammeR' is required. Please install it.", call. = FALSE)
  }

  # --- 色名をカラーコードに変換するヘルパー ---
  to_hex <- function(color_str) {
    if (grepl("^#", color_str)) return(color_str)
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
  fillcolor       <- to_hex(fillcolor)
  edge_color      <- to_hex(edge_color)
  edge_label_color <- to_hex(edge_label_color)
  color_tp        <- to_hex(color_tp)
  color_fp        <- to_hex(color_fp)
  color_fn        <- to_hex(color_fn)

  if (is.null(bordercolor)) {
    bordercolor <- fillcolor
  } else {
    bordercolor <- to_hex(bordercolor)
  }

  # --- 引数のバリデーション ---
  rankdir <- match.arg(rankdir, c("LR", "RL", "TB", "BT"))

  valid_shapes <- c(
    "circle", "box", "ellipse", "diamond", "plaintext",
    "square", "triangle", "hexagon", "octagon",
    "doublecircle", "rect", "oval", "egg",
    "pentagon", "star", "none"
  )
  if (!(shape %in% valid_shapes)) {
    stop(sprintf("'shape' must be one of: %s.", paste(valid_shapes, collapse = ", ")))
  }

  p <- ncol(B)

  if (!is.null(true_B)) {
    if (!is.matrix(true_B) || !all(dim(true_B) == c(p, p))) {
      stop("'true_B' must be a matrix with the same dimensions as 'B'.", call. = FALSE)
    }
  }

  if (is.null(labels)) labels <- paste0("x", seq_len(p) - 1)

  # --- title のエスケープ ---
  safe_label <- gsub("'", "\\\\'", title)

  # --- エッジ記述の生成 ---
  edge_lines <- c()

  if (is.null(true_B)) {
    # 通常モード: 全エッジを同色で描画
    for (i in seq_len(p)) {
      for (j in seq_len(p)) {
        if (abs(B[i, j]) > threshold) {
          coef_label <- sprintf("%.2f", B[i, j])
          edge_lines <- c(edge_lines,
            sprintf("  %s -> %s [label = ' %s']", labels[j], labels[i], coef_label)
          )
        }
      }
    }
  } else {
    # 比較モード: TP / FP / FN を色分け
    estimated <- abs(B) > threshold
    true_exist <- abs(true_B) > 0

    for (i in seq_len(p)) {
      for (j in seq_len(p)) {
        if (i == j) next

        is_est  <- estimated[i, j]
        is_true <- true_exist[i, j]

        if (!is_est && !is_true) next

        if (is_est && is_true) {
          # TP: 正解エッジ（緑・実線）
          coef_label <- sprintf("%.2f", B[i, j])
          edge_lines <- c(edge_lines, sprintf(
            "  %s -> %s [label = ' %s', color = '%s', fontcolor = '%s']",
            labels[j], labels[i], coef_label, color_tp, color_tp
          ))
        } else if (is_est && !is_true) {
          # FP: 過検出（赤・実線）
          coef_label <- sprintf("%.2f", B[i, j])
          edge_lines <- c(edge_lines, sprintf(
            "  %s -> %s [label = ' %s', color = '%s', fontcolor = '%s']",
            labels[j], labels[i], coef_label, color_fp, color_fp
          ))
        } else {
          # FN: 見逃し（オレンジ・破線、真の係数を表示）
          coef_label <- sprintf("%.2f", true_B[i, j])
          edge_lines <- c(edge_lines, sprintf(
            "  %s -> %s [label = ' %s', color = '%s', fontcolor = '%s', style = 'dashed']",
            labels[j], labels[i], coef_label, color_fn, color_fn
          ))
        }
      }
    }
  }

  if (length(edge_lines) == 0) {
    message("No edges found above the threshold. Please use a smaller 'threshold' value.")
    return(invisible(NULL))
  }

  # --- DOT 記述の構築 ---
  # true_B 指定時はエッジごとに色を上書きするので、グローバル edge 色は中立色にする
  global_edge_color       <- if (is.null(true_B)) edge_color       else "#888888"
  global_edge_label_color <- if (is.null(true_B)) edge_label_color else "#888888"

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
    "        fontcolor = '", global_edge_label_color, "',\n",
    "        color = '", global_edge_color, "']\n",
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
