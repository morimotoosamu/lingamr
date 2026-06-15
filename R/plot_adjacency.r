# =============================================================================
# Direct LiNGAM - R Implementation
# Based on the Python implementation from the LiNGAM Project
# https://sites.google.com/view/sshimizu06/lingam
# https://github.com/cdt15/lingam
# =============================================================================


#' Plot a causal graph from an adjacency matrix with DiagrammeR
#'
#' @param B Adjacency matrix (n_features x n_features).
#'   **Convention: `B[i, j]` is the causal coefficient from variable j to
#'   variable i (j -> i).** The `adjacency_matrix` from `lingam_direct()` can
#'   be passed directly.
#' @param labels Vector of variable names (if NULL, x0, x1, ... are generated
#'   automatically)
#' @param threshold Minimum absolute coefficient value to display (default: 0)
#' @param rankdir Layout direction (default: "LR")
#'   "LR" = left -> right, "RL" = right -> left, "TB" = top -> bottom,
#'   "BT" = bottom -> top
#' @param title Graph title (default: "Estimated Causal Structure")
#' @param shape Node shape (default: "circle")
#'   "circle", "box", "ellipse", "diamond", "plaintext",
#'   "square", "triangle", "hexagon", "octagon", etc.
#' @param fillcolor Node fill color (default: "lightyellow")
#' @param bordercolor Border color
#' @param fontsize_node Node font size (default: 14)
#' @param fontsize_edge Edge label font size (default: 10)
#' @param edge_color Edge color (default: "gray40"). Unused when `true_B` is
#'   specified.
#' @param edge_label_color Edge label color (default: "red"). Unused when
#'   `true_B` is specified.
#' @param true_B True adjacency matrix (may be NULL). When specified, edges are
#'   classified into three colors:
#'   * Correct edges (estimated and true): solid line in `color_tp`
#'   * False positives (estimated but not true): solid line in `color_fp`
#'   * Missed edges (not estimated but true): dashed line in `color_fn`
#'     (showing the true coefficient)
#' @param color_tp Color for correct edges (default: "forestgreen")
#' @param color_fp Color for false-positive edges (default: "firebrick")
#' @param color_fn Color for missed edges (default: "darkorange")
#' @param debug Enable debug mode (logical)
#' @return A grViz object (when DiagrammeR is available)
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
#' # Compare with the true structure
#' # (correct = green, false positive = red, missed = orange dashed)
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
  # --- Check for the DiagrammeR package ---
  if (!requireNamespace("DiagrammeR", quietly = TRUE)) {
    stop("Package 'DiagrammeR' is required. Please install it.", call. = FALSE)
  }

  # --- Helper to convert a color name to a color code ---
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

  # --- Convert all colors to color codes ---
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

  # --- Validate arguments ---
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

  # --- Escape the title ---
  safe_label <- gsub("'", "\\\\'", title)

  # --- Generate edge descriptions ---
  edge_lines <- c()

  if (is.null(true_B)) {
    # Normal mode: draw all edges in the same color
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
    # Comparison mode: color-code TP / FP / FN
    estimated <- abs(B) > threshold
    true_exist <- abs(true_B) > 0

    for (i in seq_len(p)) {
      for (j in seq_len(p)) {
        if (i == j) next

        is_est  <- estimated[i, j]
        is_true <- true_exist[i, j]

        if (!is_est && !is_true) next

        if (is_est && is_true) {
          # TP: correct edge (green, solid line)
          coef_label <- sprintf("%.2f", B[i, j])
          edge_lines <- c(edge_lines, sprintf(
            "  %s -> %s [label = ' %s', color = '%s', fontcolor = '%s']",
            labels[j], labels[i], coef_label, color_tp, color_tp
          ))
        } else if (is_est && !is_true) {
          # FP: false positive (red, solid line)
          coef_label <- sprintf("%.2f", B[i, j])
          edge_lines <- c(edge_lines, sprintf(
            "  %s -> %s [label = ' %s', color = '%s', fontcolor = '%s']",
            labels[j], labels[i], coef_label, color_fp, color_fp
          ))
        } else {
          # FN: missed edge (orange, dashed line, showing the true coefficient)
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

  # --- Build the DOT description ---
  # When true_B is specified, per-edge colors are overridden, so set the
  # global edge color to a neutral color
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
