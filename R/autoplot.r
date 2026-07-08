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
#' @param label_pos Position of each coefficient label along its edge, as a
#'   fraction from the source (0) to the target (1). The default 0.35 places
#'   labels off-center (toward the source) so labels on crossing edges do not
#'   overlap near the midpoint.
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
                                  label_edges = TRUE, label_pos = 0.35, ...) {
  autoplot_causal_graph(object$adjacency_matrix,
    threshold = threshold, node_size = node_size, node_color = node_color,
    label_edges = label_edges, label_pos = label_pos
  )
}


#' Shared drawing backend for the autoplot methods
#'
#' Draws an adjacency matrix as a directed graph. `NA` entries (unresolved
#' order / suspected latent confounding in ParceLiNGAM and RCD results) are
#' drawn as dashed, arrowless segments when `dashed_na = TRUE`; they never
#' contribute to the hierarchical layout.
#'
#' @param B adjacency matrix (`B[i, j]` is j -> i)
#' @param dashed_na whether to draw `NA` entries as dashed segments
#' @inheritParams autoplot.LingamResult
#' @return A ggplot object
#' @keywords internal
#' @noRd
autoplot_causal_graph <- function(B, threshold = 0,
                                  node_size = 16, node_color = "lightblue",
                                  label_edges = TRUE, label_pos = 0.35,
                                  dashed_na = FALSE) {
  for (pkg in c("ggplot2", "igraph")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required for autoplot(). Please install it.",
           call. = FALSE)
    }
  }

  var_names <- get_var_names(B)

  all_edges <- adjacency_edges(B, threshold = threshold, include_na = dashed_na)
  edges <- all_edges[!is.na(all_edges$estimate), , drop = FALSE]
  na_edges <- all_edges[is.na(all_edges$estimate), , drop = FALSE]
  # NA entries are usually marked in both directions; draw one dashed segment
  # per unordered pair.
  if (nrow(na_edges) > 0) {
    key <- paste(
      pmin(na_edges$from, na_edges$to),
      pmax(na_edges$from, na_edges$to)
    )
    na_edges <- na_edges[!duplicated(key), , drop = FALSE]
  }

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
    # place labels off-center (toward the source) so labels on crossing edges
    # do not collide near the midpoint
    edge_df$lx <- edge_df$x + label_pos * (edge_df$xend - edge_df$x)
    edge_df$ly <- edge_df$y + label_pos * (edge_df$yend - edge_df$y)
  }

  # --- Dashed segments for NA (confounded / unresolved) pairs ---
  if (nrow(na_edges) > 0) {
    fi <- match(na_edges$from, node_df$name)
    ti <- match(na_edges$to, node_df$name)
    na_df <- data.frame(
      x    = node_df$x[fi],
      y    = node_df$y[fi],
      xend = node_df$x[ti],
      yend = node_df$y[ti]
    )
    # Pull the endpoints inward so the segments stop at the node circles.
    dx <- na_df$xend - na_df$x
    dy <- na_df$yend - na_df$y
    na_df$x    <- na_df$x    + 0.10 * dx
    na_df$y    <- na_df$y    + 0.10 * dy
    na_df$xend <- na_df$xend - 0.10 * dx
    na_df$yend <- na_df$yend - 0.10 * dy
  } else {
    na_df <- data.frame(x = numeric(0), y = numeric(0),
                        xend = numeric(0), yend = numeric(0))
  }

  # --- Plotting ---
  pl <- ggplot2::ggplot()

  if (nrow(na_df) > 0) {
    pl <- pl + ggplot2::geom_segment(
      data = na_df,
      mapping = ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      linetype = "dashed", color = "gray60"
    )
  }

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
        mapping = ggplot2::aes(x = lx, y = ly, label = label),
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


#' Plot the causal graph of a LiMResult with ggplot2
#'
#' Draws the estimated causal structure of a LiM model as a ggplot2-based
#' directed graph, exactly like [autoplot.LingamResult()].
#'
#' @param object Return value of [lingam_lim()] (a `LiMResult` object)
#' @inheritParams autoplot.LingamResult
#' @return A ggplot object
#' @exportS3Method ggplot2::autoplot
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("igraph", quietly = TRUE)) {
#'   library(ggplot2)
#'   set.seed(1)
#'   dat <- generate_lim_sample(n = 300)
#'   model <- lingam_lim(dat$data, is_continuous = dat$is_continuous)
#'   autoplot(model)
#' }
#' }
autoplot.LiMResult <- function(object, threshold = 0,
                               node_size = 16, node_color = "lightblue",
                               label_edges = TRUE, label_pos = 0.35, ...) {
  autoplot_causal_graph(object$adjacency_matrix,
    threshold = threshold, node_size = node_size, node_color = node_color,
    label_edges = label_edges, label_pos = label_pos
  )
}


#' Plot the causal graph of a ParceLingamResult with ggplot2
#'
#' Draws the estimated causal structure as a ggplot2-based directed graph,
#' like [autoplot.LingamResult()]. Variable pairs whose adjacency-matrix
#' entries are `NA` (unresolved order / suspected latent confounding) are
#' drawn as dashed, arrowless segments.
#'
#' @param object Return value of [lingam_parce()] (a `ParceLingamResult` object)
#' @inheritParams autoplot.LingamResult
#' @return A ggplot object
#' @exportS3Method ggplot2::autoplot
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("igraph", quietly = TRUE)) {
#'   library(ggplot2)
#'   dat <- generate_parce_sample(n = 500, seed = 42)
#'   model <- lingam_parce(dat$data)
#'   autoplot(model)
#' }
#' }
autoplot.ParceLingamResult <- function(object, threshold = 0,
                                       node_size = 16, node_color = "lightblue",
                                       label_edges = TRUE, label_pos = 0.35, ...) {
  autoplot_causal_graph(object$adjacency_matrix,
    threshold = threshold, node_size = node_size, node_color = node_color,
    label_edges = label_edges, label_pos = label_pos, dashed_na = TRUE
  )
}


#' Plot the causal graph of an RCDResult with ggplot2
#'
#' Draws the estimated causal structure as a ggplot2-based directed graph,
#' like [autoplot.LingamResult()]. Variable pairs suspected to share a latent
#' confounder (`NA` entries in the adjacency matrix) are drawn as dashed,
#' arrowless segments.
#'
#' @param object Return value of [lingam_rcd()] (an `RCDResult` object)
#' @inheritParams autoplot.LingamResult
#' @return A ggplot object
#' @exportS3Method ggplot2::autoplot
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("igraph", quietly = TRUE)) {
#'   library(ggplot2)
#'   confounded <- generate_rcd_sample(n = 300, seed = 1)
#'   model <- lingam_rcd(confounded$data)
#'   autoplot(model)
#' }
#' }
autoplot.RCDResult <- function(object, threshold = 0,
                               node_size = 16, node_color = "lightblue",
                               label_edges = TRUE, label_pos = 0.35, ...) {
  autoplot_causal_graph(object$adjacency_matrix,
    threshold = threshold, node_size = node_size, node_color = node_color,
    label_edges = label_edges, label_pos = label_pos, dashed_na = TRUE
  )
}


#' Plot one group of a MultiGroupLingamResult with ggplot2
#'
#' Extracts a single group's result with [get_group_result()] and draws it
#' like [autoplot.LingamResult()], with the group shown in the subtitle. The
#' causal order is shared across groups, but the coefficients (and therefore
#' the plotted edges) are group specific.
#'
#' @param object Return value of [lingam_multi_group()]
#'   (a `MultiGroupLingamResult` object)
#' @param group Which group to plot: a group name (character) or a 1-based
#'   index (default: 1)
#' @inheritParams autoplot.LingamResult
#' @return A ggplot object
#' @exportS3Method ggplot2::autoplot
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("igraph", quietly = TRUE)) {
#'   library(ggplot2)
#'   mg <- generate_multi_group_sample()
#'   model <- lingam_multi_group(mg$data_list, reg_method = "ols")
#'   autoplot(model, group = 2)
#' }
#' }
autoplot.MultiGroupLingamResult <- function(object, group = 1,
                                            threshold = 0,
                                            node_size = 16, node_color = "lightblue",
                                            label_edges = TRUE, label_pos = 0.35, ...) {
  gr <- get_group_result(object, group)
  group_label <- if (is.character(group)) {
    group
  } else {
    nm <- names(object$adjacency_matrices)
    if (!is.null(nm)) nm[group] else paste("group", group)
  }
  autoplot_causal_graph(gr$adjacency_matrix,
    threshold = threshold, node_size = node_size, node_color = node_color,
    label_edges = label_edges, label_pos = label_pos
  ) +
    ggplot2::labs(subtitle = group_label)
}
