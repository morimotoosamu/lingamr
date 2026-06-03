# =============================================================================
# ggplot2 autoplot method for LingamResult
# =============================================================================

#' LingamResult の因果グラフを ggplot2 で描画
#'
#' 推定された因果構造を ggplot2 ベースの有向グラフとして描画する。ノード配置は
#' igraph の階層レイアウト (sugiyama) で計算し、因果の流れが概ね上から下へ並ぶ。
#' 静的画像として出力されるため RMarkdown / Quarto で安定する。対話的な HTML
#' 図が必要な場合は [plot_adjacency()]（DiagrammeR ベース）を使う。
#'
#' `autoplot()` は ggplot2 のジェネリックなので、利用前に `library(ggplot2)` で
#' ggplot2 を読み込む必要がある。描画には ggplot2 と igraph が必要。
#'
#' @param object [lingam_direct()] の返り値 (`LingamResult` オブジェクト)
#' @param threshold この絶対値以下の係数はエッジとみなさない (default: 0)
#' @param node_size ノードの大きさ (default: 16)
#' @param node_color ノードの塗り色 (default: "lightblue")
#' @param label_edges エッジに係数ラベルを表示するか (default: TRUE)
#' @param ... 未使用
#' @return ggplot オブジェクト
#' @exportS3Method ggplot2::autoplot
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("igraph", quietly = TRUE)) {
#'   library(ggplot2)
#'   dat <- generate_lingam_sample_6()
#'   model <- lingam_direct(dat$data, reg_method = "ols")
#'   autoplot(model)
#' }
#' }
autoplot.LingamResult <- function(object, threshold = 0,
                                  node_size = 16, node_color = "lightblue",
                                  label_edges = TRUE, ...) {
  for (pkg in c("ggplot2", "igraph")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required for autoplot(). Please install it.",
           call. = FALSE)
    }
  }

  B <- object$adjacency_matrix
  p <- ncol(B)
  var_names <- colnames(B)
  if (is.null(var_names)) var_names <- paste0("x", seq_len(p) - 1)

  edges <- tidy(object, threshold = threshold)

  # --- igraph グラフ作成（孤立ノードも含める）---
  g <- igraph::graph_from_data_frame(
    d = edges[, c("from", "to"), drop = FALSE],
    directed = TRUE,
    vertices = data.frame(name = var_names, stringsAsFactors = FALSE)
  )

  # --- 階層レイアウト座標（y を反転して上流を上に）---
  lay <- igraph::layout_with_sugiyama(g)$layout
  node_df <- data.frame(
    name = var_names,
    x = lay[, 1],
    y = -lay[, 2],
    stringsAsFactors = FALSE
  )

  # --- エッジ座標をノード座標から構成 ---
  if (nrow(edges) > 0) {
    fi <- match(edges$from, node_df$name)
    ti <- match(edges$to, node_df$name)
    edge_df <- data.frame(
      x     = node_df$x[fi],
      y     = node_df$y[fi],
      xend  = node_df$x[ti],
      yend  = node_df$y[ti],
      label = sprintf("%.2f", edges$estimate),
      stringsAsFactors = FALSE
    )
  } else {
    edge_df <- data.frame(x = numeric(0), y = numeric(0),
                          xend = numeric(0), yend = numeric(0),
                          label = character(0), stringsAsFactors = FALSE)
  }

  # --- 描画 ---
  pl <- ggplot2::ggplot()

  if (nrow(edge_df) > 0) {
    pl <- pl + ggplot2::geom_segment(
      data = edge_df,
      mapping = ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      arrow = ggplot2::arrow(length = ggplot2::unit(3, "mm"), type = "closed"),
      color = "gray40"
    )
    if (label_edges) {
      pl <- pl + ggplot2::geom_text(
        data = edge_df,
        mapping = ggplot2::aes(x = (x + xend) / 2, y = (y + yend) / 2, label = label),
        color = "firebrick", size = 3
      )
    }
  }

  pl +
    ggplot2::geom_point(
      data = node_df,
      mapping = ggplot2::aes(x = x, y = y),
      size = node_size, color = node_color
    ) +
    ggplot2::geom_text(
      data = node_df,
      mapping = ggplot2::aes(x = x, y = y, label = name)
    ) +
    ggplot2::theme_void() +
    ggplot2::labs(title = "Estimated Causal Structure")
}
