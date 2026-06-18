# =============================================================================
# ggplot2 autoplot method for LingamResult
# =============================================================================

#' Plot the causal graph of a LingamResult with ggplot2
#'
#' Draws the estimated causal structure as a ggplot2-based directed graph. Node
#' positions are computed with igraph's hierarchical layout (sugiyama), so the
#' causal flow is generally arranged from top to bottom. Because the output is a
#' static image, it is stable in RMarkdown / Quarto. If you need an interactive
#' HTML figure, use [plot_adjacency()] (DiagrammeR-based).
#'
#' `autoplot()` is a ggplot2 generic, so you must load ggplot2 with
#' `library(ggplot2)` before using it. Plotting requires ggplot2 and igraph.
#'
#' @param object Return value of [lingam_direct()] (a `LingamResult` object)
#' @param threshold Coefficients with an absolute value at or below this are not treated as edges (default: 0)
#' @param node_size Node size (default: 16)
#' @param node_color Node fill color (default: "lightblue")
#' @param label_edges Whether to display coefficient labels on edges (default: TRUE)
#' @param ... Unused
#' @return A ggplot object
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
  var_names <- get_var_names(B)

  edges <- tidy(object, threshold = threshold)

  # --- Build the igraph graph (including isolated nodes) ---
  g <- igraph::graph_from_data_frame(
    d = edges[, c("from", "to"), drop = FALSE],
    directed = TRUE,
    vertices = data.frame(name = var_names)
  )

  # --- Hierarchical layout coordinates ---
  # sugiyama places upstream (source) nodes at the largest y, so the raw y is
  # used directly to draw the causal flow from top (upstream) to bottom
  # (downstream). (Negating it would put the flow upside down.)
  lay <- igraph::layout_with_sugiyama(g)$layout
  node_df <- data.frame(
    name = var_names,
    x = lay[, 1],
    y = lay[, 2]
  )

  # --- Build edge coordinates from the node coordinates ---
  if (nrow(edges) > 0) {
    fi <- match(edges$from, node_df$name)
    ti <- match(edges$to, node_df$name)
    edge_df <- data.frame(
      x     = node_df$x[fi],
      y     = node_df$y[fi],
      xend  = node_df$x[ti],
      yend  = node_df$y[ti],
      label = sprintf("%.2f", edges$estimate)
    )
  } else {
    edge_df <- data.frame(x = numeric(0), y = numeric(0),
                          xend = numeric(0), yend = numeric(0),
                          label = character(0))
  }

  # Pull each edge's endpoints inward so the closed arrowhead sits in a gap
  # before the target node instead of being hidden under the node circle.
  if (nrow(edge_df) > 0) {
    dx <- edge_df$xend - edge_df$x
    dy <- edge_df$yend - edge_df$y
    edge_df$x    <- edge_df$x    + 0.10 * dx
    edge_df$y    <- edge_df$y    + 0.10 * dy
    edge_df$xend <- edge_df$xend - 0.18 * dx
    edge_df$yend <- edge_df$yend - 0.18 * dy
  }

  # --- Plotting ---
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
