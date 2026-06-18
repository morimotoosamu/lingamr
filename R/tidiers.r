# =============================================================================
# broom (generics) tidiers for lingamr
# =============================================================================

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance


#' Convert a LingamResult to a tidy data.frame
#'
#' Converts the estimated adjacency matrix into a long-format data.frame with
#' one edge per row. Following the `B[i, j]` convention (the coefficient for
#' j -> i), the `from` column is the cause and the `to` column is the effect.
#' Convenient for visualization with ggplot2 or ggraph and for filtering with dplyr.
#'
#' @param x The return value of [lingam_direct()] (a `LingamResult` object)
#' @param threshold Coefficients with an absolute value at or below this are not
#'   treated as edges (default: 0)
#' @param ... Unused
#' @return data.frame(from, to, estimate). `from`/`to` are variable names
#'   (strings) and `estimate` is the causal coefficient. Returns a 0-row
#'   data.frame if there are no edges.
#' @export
#' @examples
#' dat <- generate_lingam_sample_6()
#' model <- lingam_direct(dat$data, reg_method = "ols")
#' tidy(model)
tidy.LingamResult <- function(x, threshold = 0, ...) {
  B <- x$adjacency_matrix
  var_names <- get_var_names(B)

  idx <- which(abs(B) > threshold, arr.ind = TRUE)
  if (nrow(idx) == 0) {
    return(data.frame(
      from = character(0), to = character(0),
      estimate = numeric(0)
    ))
  }

  # B[i, j] is j -> i. Row i is the "to", column j is the "from".
  ord <- order(idx[, 2], idx[, 1])
  idx <- idx[ord, , drop = FALSE]
  data.frame(
    from     = var_names[idx[, 2]],
    to       = var_names[idx[, 1]],
    estimate = B[idx]
  )
}


#' Get a one-row summary of a LingamResult
#'
#' Summarizes the entire model in a single row. The data `X` is not required
#' because no residuals are computed. If residual-based diagnostics are needed,
#' use [summary_lingam()] instead.
#'
#' @param x The return value of [lingam_direct()] (a `LingamResult` object)
#' @param ... Unused
#' @return A one-row data.frame(n_variables, n_edges, causal_order)
#' @export
#' @examples
#' dat <- generate_lingam_sample_6()
#' model <- lingam_direct(dat$data, reg_method = "ols")
#' glance(model)
glance.LingamResult <- function(x, ...) {
  B <- x$adjacency_matrix
  p <- ncol(B)
  var_names <- get_var_names(B)
  data.frame(
    n_variables  = p,
    n_edges      = sum(abs(B) > 0),
    causal_order = paste(var_names[x$causal_order], collapse = " -> ")
  )
}


#' Convert a BootstrapResult to a tidy data.frame
#'
#' Returns a summary of the occurrence count, proportion, and effect size for
#' each causal direction. Internally it calls [get_causal_direction_counts()],
#' so that function's arguments can be passed through `...`.
#'
#' @param x The return value of [lingam_direct_bootstrap()] (a `BootstrapResult` object)
#' @param ... Arguments passed to [get_causal_direction_counts()]
#'   (such as `n_directions`, `min_causal_effect`, `split_by_causal_effect_sign`, `labels`)
#' @return data.frame (from, to, count, proportion, ...)
#' @export
#' @examples
#' dat <- generate_lingam_sample_6()
#' bs <- lingam_direct_bootstrap(dat$data, n_sampling = 30L, seed = 42)
#' tidy(bs)
tidy.BootstrapResult <- function(x, ...) {
  get_causal_direction_counts(x, ...)
}
