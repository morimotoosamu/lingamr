# =============================================================================
# VARMA-LiNGAM - Total causal effects
# Based on the Python implementation from the LiNGAM Project
# https://github.com/cdt15/lingam  (lingam/varma_lingam.py)
#
# License: MIT + file LICENSE
#
# Original work:
#   Copyright (c) 2019 T.Ikeuchi, G.Haraoka, M.Ide, W.Kurebayashi, S.Shimizu
#
# Portions of this work:
#   Copyright (c) 2026 O.Morimoto
# =============================================================================


#' Total-effect regression core for VARMA-LiNGAM
#'
#' Regresses the destination variable on the source variable plus the source's
#' parents (a back-door adjustment) over a joined design of lagged X and
#' lagged LiNGAM residuals, and returns the source's coefficient. Shared by
#' [estimate_varma_total_effect()] and the VARMA-LiNGAM bootstrap.
#'
#' The X region of the design holds lags `0..(p + from_lag)` and the residual
#' region holds lags `1..(q + from_lag)`, so parents shifted into the
#' `from_lag` block always reference filled columns. (The Python reference
#' allocates the extra `from_lag` blocks but never fills them, so its shifted
#' predictors can point at zero or misaligned columns when `from_lag > 0`.)
#'
#' @param X numeric matrix (n x m), rows ordered in time
#' @param ee_full LiNGAM residuals `e_t = (I - B0) n_t`, full length (n x m)
#' @param am_joined joined causal matrix `cbind(psi_0..psi_p, omega_1..omega_q)`
#'   (m x m(1 + p + q))
#' @param order VARMA order c(p, q)
#' @param from_index source variable (1-based integer)
#' @param to_index destination variable (1-based integer)
#' @param from_lag lag of the source variable (non-negative integer)
#' @return the estimated total effect (scalar)
#' @keywords internal
varma_total_effect_core <- function(X, ee_full, am_joined, order,
                                    from_index, to_index, from_lag) {
  p_order <- order[1]
  q_order <- order[2]
  m <- ncol(X)

  # --- joined design: [X_t, ..., X_{t-(p+from_lag)}, e_{t-1}, ..., e_{t-(q+from_lag)}] ---
  n_x_blocks <- 1L + p_order + from_lag
  n_e_blocks <- q_order + from_lag
  X_joined <- matrix(0, nrow = nrow(X), ncol = m * (n_x_blocks + n_e_blocks))
  for (b in seq_len(n_x_blocks)) {
    pos <- (b - 1L) * m
    # block b holds X shifted down by (b - 1) rows (b = 1 is contemporaneous).
    X_joined[, pos + seq_len(m)] <- roll_rows(X, b - 1L)
  }
  for (b in seq_len(n_e_blocks)) {
    pos <- (n_x_blocks + b - 1L) * m
    X_joined[, pos + seq_len(m)] <- roll_rows(ee_full, b)
  }

  # --- predictors: the source plus its parents, shifted into the from_lag block ---
  parents <- which(abs(am_joined[from_index, ]) > 0)
  x_region_width <- m * (1L + p_order)
  shift_parent <- function(c) {
    if (c <= x_region_width) {
      # X-region parent: same lag offset within the (extended) X region.
      c + m * from_lag
    } else {
      # Omega-region parent e_{t-b}: shift to e_{t-b-from_lag} in the design.
      m * n_x_blocks + (c - x_region_width) + m * from_lag
    }
  }
  parents_col <- vapply(parents, shift_parent, integer(1))
  from_col <- from_index + m * from_lag
  predictors <- unique(c(from_col, parents_col))
  from_pos <- which(predictors == from_col)

  # destination is always at the current time (block 1).
  y <- X_joined[, to_index]
  Xp <- X_joined[, predictors, drop = FALSE]
  coefs <- fit_ols(y, Xp)
  unname(coefs[from_pos])
}


#' Estimate a total causal effect in a VARMA-LiNGAM model
#'
#' Estimates the total causal effect from `from_index` (optionally at lag
#' `from_lag`) to `to_index` (at the current time) using the fitted
#' VARMA-LiNGAM model. Port of the Python reference `estimate_total_effect`:
#' the destination variable is regressed on the source variable together with
#' the source's parents (a back-door adjustment) over lagged X and lagged
#' residual regressors, and the source's coefficient is returned.
#'
#' @param X original data (matrix or data frame), rows ordered in time; must
#'   be the data the model was fitted on (the residuals stored in `result`
#'   are aligned to it)
#' @param result a `VARMALiNGAMResult` from [lingam_varma()]
#' @param from_index source variable (1-based index or variable name)
#' @param to_index destination variable (1-based index or variable name)
#' @param from_lag lag of the source variable (0 = current time, default)
#' @return the estimated total effect (scalar)
#' @details
#' The Python reference requires the residual matrix `E` as an argument; here
#' it is reconstructed internally from `result$residuals` and the fitted
#' instantaneous matrix B0, so only the original data are needed.
#' @export
#' @examples
#' sample <- generate_varmalingam_sample(n = 1000, seed = 42)
#' model <- lingam_varma(sample$data,
#'   order = c(1, 1), criterion = NULL,
#'   reg_method = "ols", prune = FALSE
#' )
#'
#' # total effect of x0 (current) on x2 (current)
#' estimate_varma_total_effect(sample$data, model, from_index = 1, to_index = 3)
estimate_varma_total_effect <- function(X, result, from_index, to_index, from_lag = 0) {
  if (!inherits(result, "VARMALiNGAMResult")) {
    stop("result must be a VARMALiNGAMResult (output of lingam_varma()).", call. = FALSE)
  }

  X <- as.matrix(X)
  col_names <- colnames(X)
  if (!is.numeric(X)) stop("X must be a numeric matrix or data frame.", call. = FALSE)

  order <- result$order
  psis <- result$adjacency_matrices$psis
  omegas <- result$adjacency_matrices$omegas
  n_features <- dim(psis)[2]
  if (ncol(X) != n_features) {
    stop(sprintf(
      "X has %d variables but result was estimated from %d.",
      ncol(X), n_features
    ), call. = FALSE)
  }
  k0 <- max(order[1], order[2])
  if (nrow(X) != nrow(result$residuals) + k0) {
    stop("X must be the data the model was fitted on (row count does not ",
      "match the stored residuals).",
      call. = FALSE
    )
  }

  from_lag <- suppressWarnings(as.integer(from_lag))
  if (length(from_lag) != 1 || is.na(from_lag) || from_lag < 0) {
    stop("from_lag must be a non-negative integer.", call. = FALSE)
  }

  # --- resolve variable name -> 1-based index ---
  from_index <- resolve_var_index(from_index, "from_index", col_names, n_features)
  to_index <- resolve_var_index(to_index, "to_index", col_names, n_features)

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

  # --- reconstruct the LiNGAM residuals e_t = (I - B0) n_t, full length ---
  E_full <- matrix(0, nrow(X), n_features)
  E_full[(k0 + 1L):nrow(X), ] <- result$residuals
  B0 <- psis[1, , ]
  ee_full <- E_full %*% t(diag(n_features) - B0)

  # --- joined causal matrix: cbind(psi_0..psi_p, omega_1..omega_q) ---
  am_joined <- joined_varma_matrix(psis, omegas)

  varma_total_effect_core(X, ee_full, am_joined, order,
    from_index, to_index, from_lag
  )
}


#' Join psi and omega arrays into a single wide matrix
#'
#' Returns `cbind(psi_0, ..., psi_p, omega_1, ..., omega_q)` with shape
#' (m x m(1 + p + q)), keeping the `[i, j]` = j -> i convention in each block.
#'
#' @param psis array (1 + p, m, m)
#' @param omegas array (q, m, m)
#' @return matrix (m x m(1 + p + q))
#' @keywords internal
joined_varma_matrix <- function(psis, omegas) {
  blocks <- lapply(seq_len(dim(psis)[1]), function(k) psis[k, , ])
  if (dim(omegas)[1] > 0L) {
    blocks <- c(blocks, lapply(seq_len(dim(omegas)[1]), function(w) omegas[w, , ]))
  }
  do.call(cbind, blocks)
}
