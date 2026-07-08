# =============================================================================
# broom (generics) tidiers for lingamr
# =============================================================================

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance


#' Long-format edges from an adjacency matrix
#'
#' Shared backend for the `tidy()` methods. Follows the `B[i, j]` = j -> i
#' convention: column j is the cause (`from`), row i is the effect (`to`).
#'
#' @param B adjacency matrix
#' @param threshold coefficients with `abs(B) <= threshold` are dropped
#' @param include_na if `TRUE`, `NA` entries (unresolved order / suspected
#'   latent confounding in ParceLiNGAM and RCD results) are kept as rows with
#'   `estimate = NA` instead of being dropped
#' @return data.frame(from, to, estimate)
#' @keywords internal
#' @noRd
adjacency_edges <- function(B, threshold = 0, include_na = FALSE) {
  var_names <- get_var_names(B)

  keep <- !is.na(B) & abs(B) > threshold
  if (include_na) keep <- keep | is.na(B)
  idx <- which(keep, arr.ind = TRUE)
  if (nrow(idx) == 0) {
    return(data.frame(
      from = character(0), to = character(0),
      estimate = numeric(0)
    ))
  }

  ord <- order(idx[, 2], idx[, 1])
  idx <- idx[ord, , drop = FALSE]
  data.frame(
    from     = var_names[idx[, 2]],
    to       = var_names[idx[, 1]],
    estimate = B[idx]
  )
}


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
  adjacency_edges(x$adjacency_matrix, threshold)
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
#' bs <- lingam_direct_bootstrap(dat$data, n_sampling = 30L, reg_method = "ols", seed = 42)
#' tidy(bs)
tidy.BootstrapResult <- function(x, ...) {
  get_causal_direction_counts(x, ...)
}


#' Convert a LiMResult to a tidy data.frame
#'
#' Converts the estimated adjacency matrix of a LiM model into a long-format
#' data.frame with one edge per row, exactly like [tidy.LingamResult()].
#'
#' @param x The return value of [lingam_lim()] (a `LiMResult` object)
#' @param threshold Coefficients with an absolute value at or below this are not
#'   treated as edges (default: 0)
#' @param ... Unused
#' @return data.frame(from, to, estimate)
#' @export
#' @examples
#' set.seed(1)
#' dat <- generate_lim_sample(n = 300)
#' model <- lingam_lim(dat$data, is_continuous = dat$is_continuous)
#' tidy(model)
tidy.LiMResult <- function(x, threshold = 0, ...) {
  adjacency_edges(x$adjacency_matrix, threshold)
}


#' Convert a ParceLingamResult to a tidy data.frame
#'
#' Converts the estimated adjacency matrix into a long-format data.frame with
#' one edge per row, like [tidy.LingamResult()]. `NA` entries of the adjacency
#' matrix (variable pairs whose order could not be resolved / suspected latent
#' confounding) are kept as rows with `estimate = NA` so they remain visible;
#' drop them with e.g. `subset(tidy(x), !is.na(estimate))` if not needed.
#'
#' @param x The return value of [lingam_parce()] (a `ParceLingamResult` object)
#' @param threshold Coefficients with an absolute value at or below this are not
#'   treated as edges (default: 0). `NA` entries are always kept.
#' @param ... Unused
#' @return data.frame(from, to, estimate)
#' @export
#' @examples
#' dat <- generate_parce_sample(n = 500, seed = 42)
#' model <- lingam_parce(dat$data)
#' tidy(model)
tidy.ParceLingamResult <- function(x, threshold = 0, ...) {
  adjacency_edges(x$adjacency_matrix, threshold, include_na = TRUE)
}


#' Convert an RCDResult to a tidy data.frame
#'
#' Converts the estimated adjacency matrix into a long-format data.frame with
#' one edge per row, like [tidy.LingamResult()]. `NA` entries of the adjacency
#' matrix (variable pairs suspected to share a latent confounder; marked in
#' both directions) are kept as rows with `estimate = NA` so they remain
#' visible; drop them with e.g. `subset(tidy(x), !is.na(estimate))` if not
#' needed.
#'
#' @param x The return value of [lingam_rcd()] (an `RCDResult` object)
#' @param threshold Coefficients with an absolute value at or below this are not
#'   treated as edges (default: 0). `NA` entries are always kept.
#' @param ... Unused
#' @return data.frame(from, to, estimate)
#' @export
#' @examples
#' confounded <- generate_rcd_sample(n = 300, seed = 1)
#' model <- lingam_rcd(confounded$data)
#' tidy(model)
tidy.RCDResult <- function(x, threshold = 0, ...) {
  adjacency_edges(x$adjacency_matrix, threshold, include_na = TRUE)
}


#' Convert a MultiGroupLingamResult to a tidy data.frame
#'
#' Stacks the per-group edge lists into a single long-format data.frame with a
#' `group` column in front, one edge per row per group (the causal order is
#' shared across groups but the coefficients differ).
#'
#' @param x The return value of [lingam_multi_group()]
#'   (a `MultiGroupLingamResult` object)
#' @param threshold Coefficients with an absolute value at or below this are not
#'   treated as edges (default: 0)
#' @param ... Unused
#' @return data.frame(group, from, to, estimate)
#' @export
#' @examples
#' mg <- generate_multi_group_sample()
#' model <- lingam_multi_group(mg$data_list, reg_method = "ols")
#' tidy(model)
tidy.MultiGroupLingamResult <- function(x, threshold = 0, ...) {
  groups <- names(x$adjacency_matrices)
  if (is.null(groups)) groups <- paste0("group", seq_along(x$adjacency_matrices))

  out <- lapply(seq_along(x$adjacency_matrices), function(i) {
    edges <- adjacency_edges(x$adjacency_matrices[[i]], threshold)
    data.frame(group = rep(groups[i], nrow(edges)), edges)
  })
  do.call(rbind, c(out, list(make.row.names = FALSE)))
}


#' Convert a MultiGroupBootstrapResult to a tidy data.frame
#'
#' Stacks each group's causal direction counts (via
#' [get_causal_direction_counts()]) into a single data.frame with a `group`
#' column in front. Arguments for [get_causal_direction_counts()] can be
#' passed through `...`.
#'
#' @param x The return value of [lingam_multi_group_bootstrap()]
#'   (a `MultiGroupBootstrapResult` object)
#' @param ... Arguments passed to [get_causal_direction_counts()]
#'   (such as `n_directions`, `min_causal_effect`, `split_by_causal_effect_sign`)
#' @return data.frame (group, from, to, count, proportion, ...)
#' @export
#' @examples
#' mg <- generate_multi_group_sample()
#' bs <- lingam_multi_group_bootstrap(mg$data_list,
#'   n_sampling = 10L, reg_method = "ols", seed = 42
#' )
#' tidy(bs)
tidy.MultiGroupBootstrapResult <- function(x, ...) {
  groups <- names(x)
  if (is.null(groups)) groups <- paste0("group", seq_along(x))

  out <- lapply(seq_along(x), function(i) {
    counts <- get_causal_direction_counts(x[[i]], ...)
    data.frame(group = rep(groups[i], nrow(counts)), counts)
  })
  do.call(rbind, c(out, list(make.row.names = FALSE)))
}


#' Convert an ImputationBootstrapResult to a tidy data.frame
#'
#' Collapses the imputation dimension with [as_bootstrap_result()] and then
#' summarizes the causal direction counts like [tidy.BootstrapResult()].
#'
#' @param x The return value of [bootstrap_with_imputation()]
#'   (an `ImputationBootstrapResult` object)
#' @param aggregate How to collapse the `n_repeats` imputation dimension,
#'   passed to [as_bootstrap_result()] ("median" or "mean")
#' @param ... Arguments passed to [get_causal_direction_counts()]
#' @return data.frame (from, to, count, proportion, ...)
#' @export
#' @examplesIf requireNamespace("mice", quietly = TRUE)
#' dat <- generate_lingam_sample_6(n = 200, seed = 1)$data
#' dat[sample(nrow(dat), 20), 1] <- NA
#' bs <- bootstrap_with_imputation(dat, n_sampling = 5L, n_repeats = 2L, seed = 42)
#' tidy(bs)
tidy.ImputationBootstrapResult <- function(x, aggregate = c("median", "mean"), ...) {
  tidy(as_bootstrap_result(x, aggregate = match.arg(aggregate)), ...)
}


#' Get a one-row summary of a LiMResult
#'
#' Like [glance.LingamResult()], with an additional `n_discrete` column giving
#' the number of discrete variables in the model.
#'
#' @param x The return value of [lingam_lim()] (a `LiMResult` object)
#' @param ... Unused
#' @return A one-row data.frame(n_variables, n_edges, n_discrete, causal_order)
#' @export
#' @examples
#' set.seed(1)
#' dat <- generate_lim_sample(n = 300)
#' model <- lingam_lim(dat$data, is_continuous = dat$is_continuous)
#' glance(model)
glance.LiMResult <- function(x, ...) {
  B <- x$adjacency_matrix
  var_names <- get_var_names(B)
  data.frame(
    n_variables  = ncol(B),
    n_edges      = sum(abs(B) > 0),
    n_discrete   = sum(!x$is_continuous),
    causal_order = paste(var_names[x$causal_order], collapse = " -> ")
  )
}


#' Get a one-row summary of a ParceLingamResult
#'
#' Like [glance.LingamResult()]. `n_edges` counts non-`NA` edges only, and
#' `n_na_entries` counts the adjacency-matrix entries left `NA` (unresolved
#' order / suspected latent confounding). Unresolved blocks in the causal
#' order are shown in parentheses, as in the print method.
#'
#' @param x The return value of [lingam_parce()] (a `ParceLingamResult` object)
#' @param ... Unused
#' @return A one-row data.frame(n_variables, n_edges, n_na_entries, causal_order)
#' @export
#' @examples
#' dat <- generate_parce_sample(n = 500, seed = 42)
#' model <- lingam_parce(dat$data)
#' glance(model)
glance.ParceLingamResult <- function(x, ...) {
  B <- x$adjacency_matrix
  var_names <- get_var_names(B)
  order_labels <- vapply(x$causal_order, function(blk) {
    if (length(blk) > 1) {
      paste0("(", paste(var_names[blk], collapse = ", "), ")")
    } else {
      var_names[blk]
    }
  }, character(1))
  data.frame(
    n_variables  = ncol(B),
    n_edges      = sum(abs(B) > 0, na.rm = TRUE),
    n_na_entries = sum(is.na(B)),
    causal_order = paste(order_labels, collapse = " -> ")
  )
}


#' Get a one-row summary of an RCDResult
#'
#' Like [glance.LingamResult()], but without a causal order (RCD does not
#' estimate one). `n_edges` counts non-`NA` edges only, and
#' `n_confounded_pairs` counts the variable pairs whose adjacency-matrix
#' entries are `NA` (suspected shared latent confounder).
#'
#' @param x The return value of [lingam_rcd()] (an `RCDResult` object)
#' @param ... Unused
#' @return A one-row data.frame(n_variables, n_edges, n_confounded_pairs)
#' @export
#' @examples
#' confounded <- generate_rcd_sample(n = 300, seed = 1)
#' model <- lingam_rcd(confounded$data)
#' glance(model)
glance.RCDResult <- function(x, ...) {
  B <- x$adjacency_matrix
  idx <- which(is.na(B), arr.ind = TRUE)
  # NA entries are marked in both directions; count unordered pairs.
  n_pairs <- if (nrow(idx) == 0) {
    0L
  } else {
    nrow(unique(t(apply(idx, 1, sort))))
  }
  data.frame(
    n_variables        = ncol(B),
    n_edges            = sum(abs(B) > 0, na.rm = TRUE),
    n_confounded_pairs = n_pairs
  )
}


#' Get a one-row summary of a MultiGroupLingamResult
#'
#' Summarizes the joint model in a single row. The causal order is shared
#' across groups; per-group edge counts are available via
#' `glance(get_group_result(x, i))`.
#'
#' @param x The return value of [lingam_multi_group()]
#'   (a `MultiGroupLingamResult` object)
#' @param ... Unused
#' @return A one-row data.frame(n_groups, n_variables, causal_order)
#' @export
#' @examples
#' mg <- generate_multi_group_sample()
#' model <- lingam_multi_group(mg$data_list, reg_method = "ols")
#' glance(model)
glance.MultiGroupLingamResult <- function(x, ...) {
  B1 <- x$adjacency_matrices[[1]]
  var_names <- get_var_names(B1)
  data.frame(
    n_groups     = length(x$adjacency_matrices),
    n_variables  = ncol(B1),
    causal_order = paste(var_names[x$causal_order], collapse = " -> ")
  )
}
