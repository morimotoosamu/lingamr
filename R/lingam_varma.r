# =============================================================================
# VARMA-LiNGAM - R Implementation
# Based on the Python implementation from the LiNGAM Project
# https://sites.google.com/view/sshimizu06/lingam
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


#' VARMA-LiNGAM for time series causal discovery
#'
#' Fits a vector autoregressive moving-average (VARMA) model to time series
#' data and applies Direct LiNGAM to the residuals to recover the
#' instantaneous (lag-0) causal structure. The lagged causal matrices (psi)
#' and the moving-average causal matrices (omega) are then derived from the
#' VARMA coefficients and the instantaneous structure.
#'
#' @param X numeric matrix or data frame (n_samples x n_features). Rows are
#'   ordered in time (earliest first).
#' @param order VARMA order `c(p, q)` (AR and MA lags, non-negative integers,
#'   not both zero). When `criterion` is not NULL, the best order in
#'   `0:p x 0:q` (excluding `c(0, 0)`) is selected by the information
#'   criterion; otherwise `order` is used directly.
#' @param criterion order-selection criterion ("bic", "aic", or "hqic"), or
#'   NULL to use `order` directly without selection.
#' @param measure independence measure passed to [lingam_direct()]
#'   ("pwling" or "kernel").
#' @param reg_method regression method for the instantaneous adjacency matrix:
#'   "adaptive_lasso" (default), "lasso", "ols", or "ridge" (see [lingam_direct()]).
#' @param lambda penalty (lambda) selection for the instantaneous matrix:
#'   "BIC" (default), "AIC", "lambda.min", "lambda.1se", or "oracle"
#'   (see [lingam_direct()]).
#' @param init_method initial-weight method for adaptive LASSO
#'   (see [lingam_direct()]).
#' @param prune logical; if `TRUE` (default, matching the Python reference),
#'   all adjacency matrices (psi and omega) are refined together by adaptive
#'   LASSO so weak edges are shrunk toward zero. Requires the glmnet package.
#'   Set `FALSE` to keep the raw transformed matrices (no glmnet needed when
#'   `reg_method = "ols"`).
#' @param ar_coefs optional known AR coefficient array (p, n_features,
#'   n_features) in reduced form. Must be supplied together with `ma_coefs`;
#'   when both are given the VARMA estimation step is skipped, `order` is
#'   derived from the arrays, and `criterion` is ignored.
#' @param ma_coefs optional known MA coefficient array (q, n_features,
#'   n_features) in reduced form. See `ar_coefs`.
#' @return A `VARMALiNGAMResult` object (list) containing:
#' * `adjacency_matrices`: list with two arrays.
#'   `psis` (1 + p, n_features, n_features): slice `[1, , ]` is the
#'   instantaneous matrix B0 (= psi_0) and slice `[tau + 1, , ]` is the lagged
#'   matrix psi_tau. `omegas` (q, n_features, n_features): slice `[w, , ]` is
#'   the MA causal matrix omega_w. Convention: entry `[i, j]` is the effect
#'   from variable j to variable i.
#' * `causal_order`: estimated causal order of the instantaneous structure
#'   (1-based indices).
#' * `residuals`: filtered VARMA residuals n_t
#'   (n_samples - max(p, q), n_features).
#' * `order`: the VARMA order `c(p, q)` actually used.
#' * `ar_coefs`, `ma_coefs`: reduced-form coefficient arrays Phi (p, m, m) and
#'   Theta (q, m, m) used to derive the causal matrices.
#' * `const`: estimated intercept vector of the reduced-form VARMA (zero when
#'   `ar_coefs`/`ma_coefs` are supplied).
#' @details
#' The structural model is
#' `x_t = B0 x_t + sum_{tau=1}^{p} psi_tau x_{t-tau} + e_t +
#' sum_{w=1}^{q} omega_w e_{t-w}`, where B0 is the instantaneous effect matrix
#' (strictly acyclic) and e_t are mutually independent non-Gaussian errors.
#' The reduced-form VARMA
#' `x_t = c + sum Phi_tau x_{t-tau} + n_t + sum Theta_w n_{t-w}` is estimated
#' first; Direct LiNGAM applied to the residuals n_t yields B0, and the causal
#' matrices follow `psi_tau = (I - B0) Phi_tau` and
#' `omega_w = (I - B0) Theta_w (I - B0)^{-1}`.
#'
#' Differences from the Python reference implementation:
#' * VARMA coefficients are estimated by the two-stage Hannan-Rissanen
#'   procedure (a long VAR yields residual estimates, then X is regressed on
#'   its own lags and the lagged residuals, both stages with an intercept)
#'   instead of state-space maximum likelihood (statsmodels `VARMAX`).
#'   Estimates are consistent but not identical to the Python output; the
#'   `max_iter` parameter of the Python class has no counterpart here because
#'   the procedure is non-iterative.
#' * The MA causal matrices use the full similarity transform
#'   `omega_w = (I - B0) Theta_w (I - B0)^{-1}`. The Python implementation
#'   drops the trailing inverse factor due to a three-argument `np.dot` call
#'   (the third argument is treated as an output buffer), so its unpruned
#'   omegas equal `(I - B0) Theta_w`.
#' * Residual filtering initializes the first `max(p, q)` residuals with
#'   zeros instead of random draws, making the fit deterministic.
#' @references
#' Kawahara, Y., Shimizu, S., & Washio, T. (2011). Analyzing relationships
#' among ARMA processes based on non-Gaussianity of external influences.
#' *Neurocomputing*, 74(12-13), 2212-2221. Ported from the Python
#' implementation cdt15/lingam (<https://github.com/cdt15/lingam>).
#' @export
#' @examples
#' sample <- generate_varmalingam_sample(n = 300, seed = 1)
#'
#' # OLS instantaneous structure without pruning (no extra packages required)
#' model <- lingam_varma(sample$data,
#'   order = c(1, 1), criterion = NULL,
#'   reg_method = "ols", prune = FALSE
#' )
#' model$causal_order
#' round(model$adjacency_matrices$psis[1, , ], 2) # instantaneous B0
lingam_varma <- function(X,
                         order = c(1L, 1L),
                         criterion = "bic",
                         measure = "pwling",
                         reg_method = "adaptive_lasso",
                         lambda = "BIC",
                         init_method = "ols",
                         prune = TRUE,
                         ar_coefs = NULL,
                         ma_coefs = NULL) {
  col_names <- if (is.data.frame(X)) names(X) else colnames(X)
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be a numeric matrix or data frame.", call. = FALSE)
  # Missing values would otherwise propagate silently through crossprod/.lm.fit.
  if (anyNA(X)) stop("X must not contain missing values (NA).", call. = FALSE)
  if (ncol(X) < 2) stop("X must have at least 2 variables (columns).", call. = FALSE)
  if (nrow(X) < 3) stop("X must have at least 3 observations (rows).", call. = FALSE)
  if (!is.null(col_names)) colnames(X) <- col_names
  n_features <- ncol(X)

  order <- suppressWarnings(as.integer(order))
  if (length(order) != 2 || anyNA(order) || any(order < 0) || sum(order) < 1) {
    stop("order must be two non-negative integers c(p, q), not both zero.", call. = FALSE)
  }
  if (!is.null(criterion)) {
    criterion <- match.arg(criterion, c("bic", "aic", "hqic"))
  }
  if (!is.logical(prune) || length(prune) != 1 || is.na(prune)) {
    stop("prune must be a single logical (TRUE or FALSE).", call. = FALSE)
  }
  if (is.null(ar_coefs) != is.null(ma_coefs)) {
    stop("ar_coefs and ma_coefs must be supplied together (or both NULL).", call. = FALSE)
  }

  if (!is.null(ar_coefs)) {
    ar_coefs <- validate_varma_coef_array(ar_coefs, n_features, "ar_coefs")
    ma_coefs <- validate_varma_coef_array(ma_coefs, n_features, "ma_coefs")
    order <- c(dim(ar_coefs)[1], dim(ma_coefs)[1])
    if (sum(order) < 1) {
      stop("ar_coefs and ma_coefs must contain at least one lag in total.", call. = FALSE)
    }
    phis <- ar_coefs
    thetas <- ma_coefs
    const <- numeric(n_features)
  } else {
    # --- VARMA order selection ---
    if (!is.null(criterion)) {
      order <- select_varma_order(X,
        max_p = order[1], max_q = order[2],
        criterion = criterion
      )
    }
    # --- VARMA(p, q) fit by two-stage Hannan-Rissanen ---
    hr <- fit_varma_hr(X, order[1], order[2])
    phis <- hr$phis
    thetas <- hr$thetas
    const <- hr$const
  }
  p_order <- order[1]
  q_order <- order[2]

  # A non-invertible MA polynomial makes the residual recursion in
  # filter_varma_residuals explosive; warn but keep going (the Python
  # reference does not enforce invertibility either).
  if (q_order >= 1L && max(companion_moduli(thetas)) >= 1) {
    warning("Estimated MA coefficients are not invertible (largest eigenvalue ",
      "modulus >= 1); residual filtering may be unstable. Consider a smaller ",
      "MA order.",
      call. = FALSE
    )
  }

  # --- residuals n_t by deterministic recursive filtering ---
  E_full <- filter_varma_residuals(X, phis, thetas, const)
  k0 <- max(p_order, q_order)
  resid <- E_full[(k0 + 1L):nrow(X), , drop = FALSE]
  colnames(resid) <- colnames(X)

  # --- instantaneous structure via Direct LiNGAM on residuals ---
  dl <- lingam_direct(resid,
    measure = measure, reg_method = reg_method,
    lambda = lambda, init_method = init_method
  )
  B0 <- dl$adjacency_matrix
  causal_order <- dl$causal_order

  # --- causal matrices: psi_tau = (I - B0) Phi_tau,
  #     omega_w = (I - B0) Theta_w (I - B0)^{-1} ---
  # B0 is a permuted strictly lower-triangular matrix, so det(I - B0) = 1 and
  # the inverse always exists.
  I_m <- diag(n_features)
  ib0 <- I_m - B0
  ib0_inv <- solve(ib0)
  psis <- array(0, dim = c(1L + p_order, n_features, n_features))
  psis[1, , ] <- B0
  for (tau in seq_len(p_order)) {
    psis[tau + 1L, , ] <- ib0 %*% phis[tau, , ]
  }
  omegas <- array(0, dim = c(q_order, n_features, n_features))
  for (w in seq_len(q_order)) {
    omegas[w, , ] <- ib0 %*% thetas[w, , ] %*% ib0_inv
  }

  # --- optional pruning: refit psi and omega jointly by adaptive LASSO ---
  # Matches the Python reference (prune = TRUE by default). Uses the causal
  # order to pick contemporaneous ancestors and refits each target on its
  # ancestors plus all lagged variables and lagged LiNGAM residuals.
  if (prune) {
    ee_full <- E_full %*% t(ib0) # LiNGAM residuals e_t = (I - B0) n_t
    pruned <- prune_varma_lingam(X, ee_full, causal_order, order,
      lambda = lambda, init_method = init_method
    )
    psis <- pruned$psis
    omegas <- pruned$omegas
  }
  # Name the slices ("lag0" = instantaneous B0, "lag1".. = psi_tau, "ma1".. =
  # omega_w) so they can be referenced by name.
  dimnames(psis) <- list(
    c("lag0", if (p_order > 0) paste0("lag", seq_len(p_order))),
    colnames(X), colnames(X)
  )
  dimnames(omegas) <- list(
    if (q_order > 0) paste0("ma", seq_len(q_order)) else character(0),
    colnames(X), colnames(X)
  )

  result <- list(
    adjacency_matrices = list(psis = psis, omegas = omegas),
    causal_order = causal_order,
    residuals = resid,
    order = order,
    ar_coefs = phis,
    ma_coefs = thetas,
    const = const
  )
  class(result) <- "VARMALiNGAMResult"
  result
}


#' Validate a user-supplied VARMA coefficient array
#'
#' @param coefs candidate array
#' @param n_features number of variables in X
#' @param arg_name argument name for error messages
#' @return the array, with dimensions checked
#' @keywords internal
validate_varma_coef_array <- function(coefs, n_features, arg_name) {
  if (!is.array(coefs) || length(dim(coefs)) != 3 || !is.numeric(coefs)) {
    stop(sprintf(
      "%s must be a numeric array of shape (lags, n_features, n_features).",
      arg_name
    ), call. = FALSE)
  }
  d <- dim(coefs)
  if (d[2] != n_features || d[3] != n_features) {
    stop(sprintf(
      "%s must have shape (lags, %d, %d) to match X.",
      arg_name, n_features, n_features
    ), call. = FALSE)
  }
  if (anyNA(coefs)) {
    stop(sprintf("%s must not contain missing values (NA).", arg_name), call. = FALSE)
  }
  coefs
}


#' Choose the long autoregression order for Hannan-Rissanen
#'
#' The first Hannan-Rissanen stage fits a long VAR(h) whose residuals proxy
#' the unobserved innovations. `h` follows the Box-Jenkins-style growth rule
#' `max(p + q, ceiling(1.5 (p + q)), floor(10 log10(n)))`, which satisfies the
#' consistency requirement that h grows with n while h^3/n -> 0.
#'
#' @param n number of observations
#' @param p_order AR order
#' @param q_order MA order
#' @param n_features number of variables
#' @return the long AR order h (integer)
#' @keywords internal
hr_long_ar_order <- function(n, p_order, q_order, n_features) {
  pq <- p_order + q_order
  min_h <- max(1L, pq)
  h <- max(min_h, as.integer(ceiling(1.5 * pq)), as.integer(floor(10 * log10(n))))
  # Saturation guard (same rationale as select_var_lag): the stage-1 VAR(h)
  # with intercept must leave enough residual degrees of freedom for the
  # residuals to be usable innovation proxies.
  df_ok <- function(h) (n - h) - (h * n_features + 1L) > 2L * n_features
  h_initial <- h
  while (h > min_h && !df_ok(h)) h <- h - 1L
  if (!df_ok(h)) {
    stop("Not enough observations to estimate the long autoregression for ",
      "Hannan-Rissanen.",
      call. = FALSE
    )
  }
  if (h < h_initial) {
    warning("Long autoregression order reduced from ", h_initial, " to ", h,
      " because of the small sample size.",
      call. = FALSE
    )
  }
  h
}


#' Fit a VARMA(p, q) model by two-stage Hannan-Rissanen
#'
#' Stage 1 fits a long VAR(h) with intercept whose residuals estimate the
#' innovations; stage 2 regresses X on its own lags 1..p and the stage-1
#' residual lags 1..q (with intercept) by OLS.
#'
#' @param X numeric matrix (n_samples x n_features), rows ordered in time
#' @param p_order AR order (non-negative integer)
#' @param q_order MA order (non-negative integer)
#' @param h long AR order for stage 1 (NULL = use [hr_long_ar_order()];
#'   0 when q = 0, where stage 1 is skipped)
#' @return list with `phis` (array (p, m, m)), `thetas` (array (q, m, m)),
#'   `const` (numeric m), `residuals` (stage-2 OLS residuals, one row per
#'   observation in the estimation window t = h + q + 1 .. n), and `h`
#' @keywords internal
fit_varma_hr <- function(X, p_order, q_order, h = NULL) {
  n <- nrow(X)
  m <- ncol(X)

  if (q_order == 0L) {
    # Pure AR: stage 2 needs no innovation proxies, so skip the long
    # autoregression entirely and use the full VAR(p) window.
    h <- 0L
    ehat_full <- matrix(0, n, m) # never read (the q loop below is empty)
  } else {
    if (is.null(h)) h <- hr_long_ar_order(n, p_order, q_order, m)

    # --- stage 1: long VAR(h) with intercept -> innovation proxies ---
    Y1 <- X[(h + 1L):n, , drop = FALSE]
    Z1 <- cbind(1, build_lag_matrix(X, h))
    ehat_full <- matrix(0, n, m) # aligned by absolute time t
    ehat_full[(h + 1L):n, ] <- stats::.lm.fit(Z1, Y1)$residuals
  }

  # --- stage 2: OLS of X_t on X lags and residual lags over the window where
  #     all regressors exist ---
  t0 <- max(h + q_order, p_order) + 1L
  if (t0 > n) {
    stop("Not enough observations for the requested VARMA order.", call. = FALSE)
  }
  Y2 <- X[t0:n, , drop = FALSE]
  # Column 1 is all ones and is never overwritten below: the intercept.
  Z2 <- matrix(1, nrow = n - t0 + 1L, ncol = 1L + (p_order + q_order) * m)
  for (tau in seq_len(p_order)) {
    Z2[, 1L + (tau - 1L) * m + seq_len(m)] <- X[(t0 - tau):(n - tau), , drop = FALSE]
  }
  for (w in seq_len(q_order)) {
    Z2[, 1L + p_order * m + (w - 1L) * m + seq_len(m)] <-
      ehat_full[(t0 - w):(n - w), , drop = FALSE]
  }
  fit <- stats::.lm.fit(Z2, Y2)
  coef <- fit$coefficients # (1 + (p+q)*m) x m
  if (anyNA(coef)) {
    stop("VARMA design matrix is rank deficient; reduce the order or provide ",
      "more data.",
      call. = FALSE
    )
  }

  const <- as.numeric(coef[1L, ])
  # split into Phi_tau / Theta_w: coef rows follow the Z2 block layout, and
  # coef[block_row j, i] is the effect of regressor j on variable i, so each
  # coefficient block transposes into the (m, m) matrix (same as fit_var_ols).
  phis <- array(0, dim = c(p_order, m, m))
  for (tau in seq_len(p_order)) {
    phis[tau, , ] <- t(coef[1L + (tau - 1L) * m + seq_len(m), , drop = FALSE])
  }
  thetas <- array(0, dim = c(q_order, m, m))
  for (w in seq_len(q_order)) {
    thetas[w, , ] <- t(coef[1L + p_order * m + (w - 1L) * m + seq_len(m), , drop = FALSE])
  }
  list(phis = phis, thetas = thetas, const = const, residuals = fit$residuals, h = h)
}


#' Filter VARMA residuals by deterministic recursion
#'
#' Computes `n_t = x_t - c - sum Phi_tau x_{t-tau} - sum Theta_w n_{t-w}`
#' recursively. The first `max(p, q)` residuals are initialized with zeros
#' (the Python reference draws them from a standard normal, which makes the
#' fit non-deterministic; the transient effect dies out at the same rate
#' either way).
#'
#' @param X numeric matrix (n_samples x n_features), rows ordered in time
#' @param phis AR coefficient array (p, m, m)
#' @param thetas MA coefficient array (q, m, m)
#' @param const intercept vector (NULL = zero)
#' @return full-length residual matrix (n_samples x n_features); the first
#'   `max(p, q)` rows are zero
#' @keywords internal
filter_varma_residuals <- function(X, phis, thetas, const = NULL) {
  n <- nrow(X)
  m <- ncol(X)
  p_order <- dim(phis)[1]
  q_order <- dim(thetas)[1]
  k0 <- max(p_order, q_order)
  if (is.null(const)) const <- numeric(m)

  E <- matrix(0, n, m)
  for (t in (k0 + 1L):n) {
    pred <- const
    for (tau in seq_len(p_order)) {
      pred <- pred + phis[tau, , ] %*% X[t - tau, ]
    }
    for (w in seq_len(q_order)) {
      pred <- pred + thetas[w, , ] %*% E[t - w, ]
    }
    E[t, ] <- X[t, ] - pred
  }
  E
}


#' Select the VARMA order by information criterion
#'
#' All candidate orders `(p, q)` in `0:max_p x 0:max_q` (excluding `(0, 0)`)
#' are compared on a **common sample**: the stage-1 long VAR is estimated once
#' with `h` derived from `(max_p, max_q)`, and every stage-2 candidate
#' regression uses the same window t = h + max_q + 1 .. n, so the criteria are
#' comparable across orders (same rationale as [select_var_lag()]).
#'
#' @param X numeric matrix (n_samples x n_features)
#' @param max_p maximum AR order to consider
#' @param max_q maximum MA order to consider
#' @param criterion "bic", "aic", or "hqic"
#' @return the selected order `c(p, q)` (integer vector)
#' @keywords internal
select_varma_order <- function(X, max_p, max_q, criterion = "bic") {
  n <- nrow(X)
  m <- ncol(X)
  h <- hr_long_ar_order(n, max_p, max_q, m)

  # Stage 1 once, shared by every candidate.
  Y1 <- X[(h + 1L):n, , drop = FALSE]
  Z1 <- cbind(1, build_lag_matrix(X, h))
  ehat_full <- matrix(0, n, m)
  ehat_full[(h + 1L):n, ] <- stats::.lm.fit(Z1, Y1)$residuals

  # Common estimation window for all candidates.
  t0 <- h + max_q + 1L
  n_obs <- n - t0 + 1L
  if (n_obs < 1) {
    stop("Not enough observations for the requested maximum VARMA order.", call. = FALSE)
  }
  Y <- X[t0:n, , drop = FALSE]
  # Full design blocks at the maximum orders; candidates use column subsets.
  Zx <- matrix(0, n_obs, max_p * m)
  for (tau in seq_len(max_p)) {
    Zx[, (tau - 1L) * m + seq_len(m)] <- X[(t0 - tau):(n - tau), , drop = FALSE]
  }
  Ze <- matrix(0, n_obs, max_q * m)
  for (w in seq_len(max_q)) {
    Ze[, (w - 1L) * m + seq_len(m)] <- ehat_full[(t0 - w):(n - w), , drop = FALSE]
  }

  # Fallback stays inside the requested search space (e.g. a pure-MA search
  # with max_p = 0 must not fall back to an AR term).
  best_order <- c(min(1L, max_p), min(1L, max_q))
  best_ic <- Inf
  any_evaluated <- FALSE
  for (p_cand in 0:max_p) {
    for (q_cand in 0:max_q) {
      if (p_cand == 0L && q_cand == 0L) next
      # Per-equation parameter count including the intercept.
      n_params_eq <- (p_cand + q_cand) * m + 1L
      # Saturation guard: see select_var_lag for the rationale.
      if (n_obs - n_params_eq <= 2L * m) next

      Z <- cbind(
        1,
        Zx[, seq_len(p_cand * m), drop = FALSE],
        Ze[, seq_len(q_cand * m), drop = FALSE]
      )
      fit <- stats::.lm.fit(Z, Y)
      if (anyNA(fit$coefficients)) next
      any_evaluated <- TRUE

      sigma <- crossprod(fit$residuals) / n_obs
      log_det <- as.numeric(determinant(sigma, logarithm = TRUE)$modulus)
      # Total free parameters; the intercept is shared by all candidates and
      # therefore dropped from the count.
      n_params <- (p_cand + q_cand) * m * m
      ic <- switch(criterion,
        bic  = log_det + (log(n_obs) / n_obs) * n_params,
        aic  = log_det + (2 / n_obs) * n_params,
        hqic = log_det + (2 * log(log(n_obs)) / n_obs) * n_params
      )
      if (is.finite(ic) && ic < best_ic) {
        best_ic <- ic
        best_order <- c(p_cand, q_cand)
      }
    }
  }
  if (!any_evaluated) {
    warning("Too few observations relative to the number of variables to ",
      "reliably compare VARMA orders up to (", max_p, ", ", max_q, "); ",
      "falling back to order = c(", best_order[1], ", ", best_order[2], ").",
      call. = FALSE
    )
  }
  as.integer(best_order)
}


#' Eigenvalue moduli of the companion matrix
#'
#' Builds the companion matrix of a lag-coefficient array (k, m, m) and
#' returns all eigenvalue moduli. Used for AR stationarity and MA
#' invertibility checks (both require every modulus to be below 1).
#'
#' @param coefs coefficient array (k, m, m)
#' @return numeric vector of eigenvalue moduli (empty when k = 0)
#' @keywords internal
companion_moduli <- function(coefs) {
  k <- dim(coefs)[1]
  if (k == 0L) {
    return(numeric(0))
  }
  m <- dim(coefs)[2]
  companion <- matrix(0, k * m, k * m)
  for (j in seq_len(k)) {
    companion[seq_len(m), (j - 1L) * m + seq_len(m)] <- coefs[j, , ]
  }
  if (k > 1L) {
    companion[(m + 1L):(k * m), seq_len((k - 1L) * m)] <- diag((k - 1L) * m)
  }
  Mod(eigen(companion, only.values = TRUE)$values)
}


#' Prune VARMA-LiNGAM adjacency matrices by adaptive LASSO
#'
#' Re-estimates the instantaneous matrix B0, every lagged matrix psi_tau, and
#' every MA matrix omega_w jointly, shrinking weak edges to zero. Port of the
#' Python reference `_pruning`. For each target variable, the predictors are
#' its contemporaneous ancestors (those preceding it in `causal_order`) plus
#' all variables at lags 1..p and the LiNGAM residuals e_t at lags 1..q.
#' Unlike the Python reference (which wraps rows around with `np.roll`), the
#' regression window is restricted to rows where every regressor is observed.
#'
#' @param X numeric matrix (n_samples x n_features), rows ordered in time
#' @param ee_full LiNGAM residuals `e_t = (I - B0) n_t`, full length
#'   (n_samples rows; the first max(p, q) rows are zero)
#' @param causal_order instantaneous causal order (1-based indices)
#' @param order VARMA order c(p, q)
#' @param lambda lambda selection passed to [fit_adaptive_lasso()]
#' @param init_method initial-weight method for adaptive LASSO
#' @return list with `psis` (array (1 + p, m, m); slice 1 is B0) and `omegas`
#'   (array (q, m, m))
#' @keywords internal
prune_varma_lingam <- function(X, ee_full, causal_order, order,
                               lambda = "BIC", init_method = "ols") {
  p_order <- order[1]
  q_order <- order[2]
  n <- nrow(X)
  m <- ncol(X)
  k0 <- max(p_order, q_order)
  # Window where X lags are observed and the ee lags come from filtered rows
  # (rows 1..k0 of ee_full are the zero-initialized transient).
  t0 <- 2L * k0 + 1L
  n_obs <- n - t0 + 1L
  if (n_obs < 1) {
    stop("Not enough observations to prune the requested VARMA order.", call. = FALSE)
  }

  # Contemporaneous block X_t (the regression targets).
  Y_full <- X[t0:n, , drop = FALSE]
  # Lagged design [X_{t-1}, ..., X_{t-p}, e_{t-1}, ..., e_{t-q}].
  Z <- matrix(0, n_obs, (p_order + q_order) * m)
  for (tau in seq_len(p_order)) {
    Z[, (tau - 1L) * m + seq_len(m)] <- X[(t0 - tau):(n - tau), , drop = FALSE]
  }
  for (w in seq_len(q_order)) {
    Z[, p_order * m + (w - 1L) * m + seq_len(m)] <-
      ee_full[(t0 - w):(n - w), , drop = FALSE]
  }

  psis <- array(0, dim = c(1L + p_order, m, m))
  omegas <- array(0, dim = c(q_order, m, m))
  for (i in seq_len(m)) {
    # Number of variables ahead of i in the causal order = its ancestors.
    co_no <- which(causal_order == i) - 1L
    ancestors <- if (co_no >= 1L) causal_order[seq_len(co_no)] else integer(0)

    obj <- Y_full[, i]
    # Predictors: contemporaneous ancestors first, then all lagged blocks.
    exp_mat <- cbind(Y_full[, ancestors, drop = FALSE], Z)
    coefs <- fit_adaptive_lasso(obj, exp_mat, lambda = lambda, init_method = init_method)

    # First co_no coefficients map back to B0[i, ancestors].
    if (co_no >= 1L) psis[1, i, ancestors] <- coefs[seq_len(co_no)]
    # The remaining coefficients map back to psi_tau[i, ] and omega_w[i, ].
    for (tau in seq_len(p_order)) {
      psis[tau + 1L, i, ] <- coefs[co_no + (tau - 1L) * m + seq_len(m)]
    }
    for (w in seq_len(q_order)) {
      omegas[w, i, ] <- coefs[co_no + p_order * m + (w - 1L) * m + seq_len(m)]
    }
  }
  list(psis = psis, omegas = omegas)
}


#' Print method for VARMALiNGAMResult
#'
#' @param x VARMALiNGAMResult object
#' @param digits number of digits to display
#' @param ... additional arguments (unused)
#' @return The input object `x`, invisibly.
#' @export
#' @examples
#' sample <- generate_varmalingam_sample(n = 300, seed = 42)
#' model <- lingam_varma(sample$data,
#'   order = c(1, 1), criterion = NULL,
#'   reg_method = "ols", prune = FALSE
#' )
#' print(model)
print.VARMALiNGAMResult <- function(x, digits = 3, ...) {
  n <- length(x$causal_order)
  psis <- x$adjacency_matrices$psis
  omegas <- x$adjacency_matrices$omegas
  var_names <- dimnames(psis)[[2]]
  order_labels <- if (!is.null(var_names)) {
    var_names[x$causal_order]
  } else {
    paste0("x", x$causal_order - 1L)
  }
  cat("VARMA-LiNGAM Result\n")
  cat(sprintf("  Variables : %d\n", n))
  cat(sprintf("  Order (p, q) : (%d, %d)\n", x$order[1], x$order[2]))
  cat(sprintf("  Causal order (instantaneous): %s\n", paste(order_labels, collapse = " -> ")))
  cat("\nInstantaneous adjacency matrix B0 (row = to, col = from):\n")
  print(round(psis[1, , ], digits = digits))
  # Also show each lagged matrix so the full model is visible at a glance.
  for (tau in seq_len(x$order[1])) {
    cat(sprintf("\nLagged adjacency matrix psi%d (row = to, col = from):\n", tau))
    print(round(psis[tau + 1L, , ], digits = digits))
  }
  for (w in seq_len(x$order[2])) {
    cat(sprintf("\nMA adjacency matrix omega%d (row = to, col = from):\n", w))
    print(round(omegas[w, , ], digits = digits))
  }
  invisible(x)
}
