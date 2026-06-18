# =============================================================================
# VAR-LiNGAM - Total causal effects
# Based on the Python implementation from the LiNGAM Project
# https://github.com/cdt15/lingam  (lingam/var_lingam.py)
#
# License: MIT + file LICENSE
#
# Original work:
#   Copyright (c) 2019 T.Ikeuchi, G.Haraoka, M.Ide, W.Kurebayashi, S.Shimizu
#
# Portions of this work:
#   Copyright (c) 2026 O.Morimoto
# =============================================================================


#' Roll matrix rows (numpy np.roll equivalent, axis = 0)
#'
#' Shifts the rows of `M` downward by `shift`, wrapping the last `shift` rows
#' around to the top. Used to build the lagged design for total-effect
#' regression. The wrap-around contaminates the first `shift` rows, matching
#' the Python reference (the effect is negligible for long series).
#'
#' @param M numeric matrix
#' @param shift non-negative integer number of rows to shift down
#' @return matrix with rolled rows
#' @keywords internal
roll_rows <- function(M, shift) {
  n <- nrow(M)
  if (shift == 0L) return(M)
  idx <- ((seq_len(n) - 1L - shift) %% n) + 1L
  M[idx, , drop = FALSE]
}


#' Estimate a total causal effect in a VAR-LiNGAM model
#'
#' Estimates the total causal effect from `from_index` (optionally at lag
#' `from_lag`) to `to_index` (at the current time) using the fitted VAR-LiNGAM
#' model. Port of the Python reference `estimate_total_effect`: the destination
#' variable is regressed on the source variable together with the source's
#' parents (a back-door adjustment), and the source's coefficient is returned.
#'
#' @param X original data (matrix or data frame), rows ordered in time
#' @param result a `VARLiNGAMResult` from [lingam_var()]
#' @param from_index source variable (1-based index or variable name)
#' @param to_index destination variable (1-based index or variable name)
#' @param from_lag lag of the source variable (0 = current time, default)
#' @return the estimated total effect (scalar)
#' @export
#' @examples
#' sample <- generate_varlingam_sample(n = 1000, seed = 42)
#' model <- lingam_var(sample$data, lags = 1, reg_method = "ols", prune = FALSE)
#'
#' # total effect of x0 (current) on x2 (current)
#' estimate_var_total_effect(sample$data, model, from_index = 1, to_index = 3)
estimate_var_total_effect <- function(X, result, from_index, to_index, from_lag = 0) {
  if (!inherits(result, "VARLiNGAMResult")) {
    stop("result must be a VARLiNGAMResult (output of lingam_var()).", call. = FALSE)
  }

  if (is.data.frame(X)) {
    col_names <- colnames(X)
    X <- as.matrix(X)
  } else {
    X <- as.matrix(X)
    col_names <- colnames(X)
  }
  if (!is.numeric(X)) stop("X must be a numeric matrix or data frame.", call. = FALSE)

  lags <- result$lags
  am <- result$adjacency_matrices       # (lags + 1, p, p)
  n_features <- dim(am)[2]
  if (ncol(X) != n_features) {
    stop(sprintf(
      "X has %d variables but result was estimated from %d.",
      ncol(X), n_features
    ), call. = FALSE)
  }

  from_lag <- suppressWarnings(as.integer(from_lag))
  if (length(from_lag) != 1 || is.na(from_lag) || from_lag < 0) {
    stop("from_lag must be a non-negative integer.", call. = FALSE)
  }

  # --- resolve variable name -> 1-based index ---
  resolve_index <- function(idx, arg_name) {
    if (is.character(idx)) {
      if (is.null(col_names)) {
        stop(sprintf("'%s' was given as a name, but X has no column names.", arg_name),
             call. = FALSE)
      }
      pos <- match(idx, col_names)
      if (is.na(pos)) {
        stop(sprintf("Variable '%s' not found in X.", idx), call. = FALSE)
      }
      return(pos)
    }
    idx <- as.integer(idx)
    if (length(idx) != 1 || is.na(idx) || idx < 1 || idx > n_features) {
      stop(sprintf("'%s' must be between 1 and %d.", arg_name, n_features), call. = FALSE)
    }
    idx
  }
  from_index <- resolve_index(from_index, "from_index")
  to_index <- resolve_index(to_index, "to_index")

  # --- warn if the instantaneous causal order is reversed (from after to) ---
  if (from_lag == 0L) {
    causal_order <- result$causal_order
    from_order <- which(causal_order == from_index)
    to_order <- which(causal_order == to_index)
    if (length(from_order) && length(to_order) && from_order > to_order) {
      from_label <- if (!is.null(col_names)) col_names[from_index] else paste0("x", from_index - 1L)
      to_label <- if (!is.null(col_names)) col_names[to_index] else paste0("x", to_index - 1L)
      warning(sprintf(
        "Causal order of %s (to) is earlier than %s (from). Result may be incorrect.",
        to_label, from_label
      ), call. = FALSE)
    }
  }

  # --- joined adjacency matrix: cbind(B0, B1, ..., Bp), p x p*(1 + lags) ---
  # am_joined[i, j] keeps the j -> i convention across all lag blocks.
  am_joined <- do.call(cbind, lapply(seq_len(lags + 1L), function(k) am[k, , ]))

  # --- joined data: blocks [X_t, X_{t-1}, ..., X_{t-(lags+from_lag)}] ---
  n_blocks <- 1L + lags + from_lag
  X_joined <- matrix(0, nrow = nrow(X), ncol = n_features * n_blocks)
  for (b in seq_len(n_blocks)) {
    pos <- (b - 1L) * n_features
    # block b holds X shifted down by (b - 1) rows (b = 1 is contemporaneous).
    X_joined[, (pos + 1L):(pos + n_features)] <- roll_rows(X, b - 1L)
  }

  # --- predictors: the source plus its parents, shifted into the from_lag block ---
  parents <- which(abs(am_joined[from_index, ]) > 0)
  from_col <- from_index + n_features * from_lag
  parents_col <- if (from_lag == 0L) parents else parents + n_features * from_lag
  predictors <- unique(c(from_col, parents_col))
  from_pos <- which(predictors == from_col)

  # destination is always at the current time (to_lag = 0 -> block 1).
  y <- X_joined[, to_index]
  Xp <- X_joined[, predictors, drop = FALSE]
  coefs <- fit_ols(y, Xp)
  unname(coefs[from_pos])
}


#' Total causal effect from a joined adjacency matrix (graph-based)
#'
#' Computes the total effect by summing path products over the time-expanded
#' graph, reusing [calculate_total_effect()]. Port of the Python reference
#' `estimate_total_effect2`; used internally by the VAR-LiNGAM bootstrap.
#'
#' @param am_joined joined adjacency matrix (n_features x n_features*(1 + lags)),
#'   `B[i, j]` is the coefficient of j -> i
#' @param from_index source column in the joined index space (1-based)
#' @param to_index destination column in the joined index space (1-based)
#' @return the total effect (scalar)
#' @keywords internal
var_total_effect_graph <- function(am_joined, from_index, to_index) {
  # Pad with zero rows so the joined matrix becomes square (lagged blocks have
  # no outgoing-from-themselves rows of their own), then enumerate paths.
  ncol_j <- ncol(am_joined)
  nrow_j <- nrow(am_joined)
  am_sq <- matrix(0, nrow = ncol_j, ncol = ncol_j)
  am_sq[seq_len(nrow_j), ] <- am_joined
  calculate_total_effect(am_sq, from_index, to_index)
}
