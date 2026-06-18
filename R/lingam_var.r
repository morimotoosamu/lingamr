# =============================================================================
# VAR-LiNGAM - R Implementation
# Based on the Python implementation from the LiNGAM Project
# https://sites.google.com/view/sshimizu06/lingam
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


#' VAR-LiNGAM for time series causal discovery
#'
#' Fits a vector autoregressive (VAR) model to time series data and applies
#' Direct LiNGAM to the residuals to recover the instantaneous (lag-0) causal
#' structure. The lagged causal matrices are then derived from the VAR
#' coefficients and the instantaneous structure.
#'
#' @param X numeric matrix or data frame (n_samples x n_features). Rows are
#'   ordered in time (earliest first).
#' @param lags maximum lag order. When `criterion` is not NULL, the best lag
#'   in `1:lags` is selected by the information criterion; otherwise `lags` is
#'   used directly.
#' @param criterion lag-selection criterion ("bic", "aic", "hqic", or "fpe"),
#'   or NULL to use `lags` directly without selection.
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
#'   all adjacency matrices (instantaneous B0 and the lagged B_k) are refined
#'   together by adaptive LASSO so weak edges are shrunk toward zero. Requires
#'   the glmnet package. Set `FALSE` to keep the raw `B_k = (I - B0) M_k`
#'   matrices (no glmnet needed when `reg_method = "ols"`).
#' @return A `VARLiNGAMResult` object (list) containing:
#' * `adjacency_matrices`: array (1 + lags, n_features, n_features). The first
#'   slice `[1, , ]` is the instantaneous matrix B0; slice `[k + 1, , ]` is the
#'   lagged matrix B_k for lag k (k = 1..lags). Convention: `B[i, j]` is the
#'   effect from variable j to variable i.
#' * `causal_order`: estimated causal order of the instantaneous structure
#'   (1-based indices).
#' * `residuals`: VAR residuals (n_samples - lags, n_features).
#' * `lags`: the lag order actually used.
#' @details
#' The model is `X_t = B0 X_t + sum_{k=1}^{p} B_k X_{t-k} + e_t`, where B0 is
#' the instantaneous effect matrix (strictly acyclic) and e_t are mutually
#' independent non-Gaussian errors. VAR coefficients `M_k` are estimated by
#' ordinary least squares (no intercept); residuals `e_t = X_t - sum M_k X_{t-k}`
#' are passed to [lingam_direct()] to obtain B0, and the lagged matrices follow
#' `B_k = (I - B0) M_k`.
#' @references
#' Hyvärinen, A., Zhang, K., Shimizu, S., & Hoyer, P. O. (2010). Estimation of a
#' structural vector autoregression model using non-Gaussianity. *Journal of
#' Machine Learning Research*, 11, 1709-1731. Ported from the Python
#' implementation cdt15/lingam (<https://github.com/cdt15/lingam>). See also the
#' VARLiNGAM R code of Moneta et al.
#' (<https://sites.google.com/site/dorisentner/publications/VARLiNGAM>).
#' @importFrom stats sd
#' @export
#' @examples
#' sample <- generate_varlingam_sample(n = 500, seed = 42)
#'
#' # OLS instantaneous structure without pruning (no extra packages required)
#' model <- lingam_var(sample$data, lags = 1, reg_method = "ols", prune = FALSE)
#' round(model$adjacency_matrices[1, , ], 2)  # instantaneous B0
lingam_var <- function(X,
                       lags = 1L,
                       criterion = "bic",
                       measure = "pwling",
                       reg_method = "adaptive_lasso",
                       lambda = "BIC",
                       init_method = "ols",
                       prune = TRUE) {
  col_names <- if (is.data.frame(X)) names(X) else colnames(X)
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be a numeric matrix or data frame.", call. = FALSE)
  # Missing values would otherwise propagate silently through crossprod/.lm.fit.
  if (anyNA(X)) stop("X must not contain missing values (NA).", call. = FALSE)
  if (ncol(X) < 2) stop("X must have at least 2 variables (columns).", call. = FALSE)
  if (nrow(X) < 3) stop("X must have at least 3 observations (rows).", call. = FALSE)
  if (!is.null(col_names)) colnames(X) <- col_names

  lags <- suppressWarnings(as.integer(lags))
  if (length(lags) != 1 || is.na(lags) || lags < 1) {
    stop("lags must be a positive integer.", call. = FALSE)
  }
  if (!is.null(criterion)) {
    criterion <- match.arg(criterion, c("bic", "aic", "hqic", "fpe"))
  }
  if (!is.logical(prune) || length(prune) != 1 || is.na(prune)) {
    stop("prune must be a single logical (TRUE or FALSE).", call. = FALSE)
  }
  n_features <- ncol(X)

  # --- lag order selection ---
  if (!is.null(criterion)) {
    lags <- select_var_lag(X, max_lag = lags, criterion = criterion)
  }

  # --- VAR(p) fit by OLS (no intercept) ---
  vf <- fit_var_ols(X, lags)
  M <- vf$coefs          # array (lags, n_features, n_features)
  resid <- vf$residuals  # (n_samples - lags, n_features)
  colnames(resid) <- colnames(X)

  # --- instantaneous structure via Direct LiNGAM on residuals ---
  dl <- lingam_direct(resid,
    measure = measure, reg_method = reg_method,
    lambda = lambda, init_method = init_method
  )
  B0 <- dl$adjacency_matrix
  causal_order <- dl$causal_order

  # --- lagged matrices: B_k = (I - B0) M_k ---
  I_p <- diag(n_features)
  B <- array(0, dim = c(lags + 1L, n_features, n_features))
  B[1, , ] <- B0
  for (k in seq_len(lags)) {
    B[k + 1L, , ] <- (I_p - B0) %*% M[k, , ]
  }

  # --- optional pruning: refit B0 and all B_k jointly by adaptive LASSO ---
  # Matches the Python reference (prune = TRUE by default). Uses the causal
  # order to pick contemporaneous ancestors and refits each target on its
  # ancestors plus all lagged variables, shrinking weak edges to zero.
  if (prune) {
    B <- prune_var_lingam(X, causal_order, lags,
      lambda = lambda, init_method = init_method
    )
  }
  # Name the lag dimension ("lag0" = instantaneous B0, "lag1".. = B_k) so slices
  # can be referenced by name, e.g. adjacency_matrices["lag1", , ].
  lag_labels <- c("lag0", paste0("lag", seq_len(lags)))
  dimnames(B) <- list(lag_labels, colnames(X), colnames(X))

  result <- list(
    adjacency_matrices = B,
    causal_order = causal_order,
    residuals = resid,
    lags = lags
  )
  class(result) <- "VARLiNGAMResult"
  result
}


#' Build the lagged design matrix for VAR models
#'
#' Constructs Z such that column block k (columns `(k-1)*p+1` to `k*p`)
#' contains `X_{t-k}` for `t = lags+1, ..., n`.
#'
#' @param X numeric matrix (n_samples x n_features)
#' @param lags lag order
#' @return matrix of shape `(n - lags, lags * p)`
#' @keywords internal
build_lag_matrix <- function(X, lags) {
  n <- nrow(X)
  p <- ncol(X)
  n_obs <- n - lags
  Z <- matrix(0, nrow = n_obs, ncol = lags * p)
  for (k in seq_len(lags)) {
    Z[, ((k - 1L) * p + 1L):(k * p)] <- X[(lags + 1L - k):(n - k), , drop = FALSE]
  }
  Z
}


#' Fit a VAR(p) model by OLS (no intercept)
#'
#' @param X numeric matrix (n_samples x n_features), rows ordered in time
#' @param lags lag order (positive integer)
#' @return list with `coefs` (array (lags, n_features, n_features); `coefs[k, , ]`
#'   is M_k such that `X_t = sum_k M_k X_{t-k} + e_t`) and `residuals`
#'   (n_samples - lags, n_features)
#' @keywords internal
fit_var_ols <- function(X, lags) {
  n <- nrow(X)
  p <- ncol(X)
  n_obs <- n - lags
  if (n_obs < 1) {
    stop("Not enough observations for the requested lag order.", call. = FALSE)
  }
  # response: X_t for t = (lags + 1) .. n
  Y <- X[(lags + 1L):n, , drop = FALSE]
  # design: [X_{t-1}, ..., X_{t-lags}] stacked column-wise
  Z <- build_lag_matrix(X, lags)
  # OLS via QR decomposition (numerically more stable than solving the
  # normal equations crossprod(Z), which squares the condition number).
  # No intercept column is added, matching the Python reference (trend = "n").
  fit <- stats::.lm.fit(Z, Y)
  coef <- fit$coefficients   # (lags*p) x p
  resid <- fit$residuals     # n_obs x p
  # split into M_k: coef[(k-1)*p + j, i] = M_k[i, j], so M_k = t(block_k)
  M <- array(0, dim = c(lags, p, p))
  for (k in seq_len(lags)) {
    M[k, , ] <- t(coef[((k - 1L) * p + 1L):(k * p), , drop = FALSE])
  }
  list(coefs = M, residuals = resid)
}


#' Prune VAR-LiNGAM adjacency matrices by adaptive LASSO
#'
#' Re-estimates the instantaneous matrix B0 and every lagged matrix B_k jointly,
#' shrinking weak edges to zero. Port of the Python reference `_pruning`.
#' For each target variable, the predictors are its contemporaneous ancestors
#' (those preceding it in `causal_order`) plus all variables at lags 1..lags;
#' the coefficients are fitted by adaptive LASSO and written back into B.
#'
#' @param X numeric matrix (n_samples x n_features), rows ordered in time
#' @param causal_order instantaneous causal order (1-based indices)
#' @param lags lag order
#' @param lambda lambda selection passed to [fit_adaptive_lasso()]
#' @param init_method initial-weight method for adaptive LASSO
#' @return array (lags + 1, n_features, n_features); slice 1 is B0, slice k+1 is B_k
#' @keywords internal
prune_var_lingam <- function(X, causal_order, lags, lambda = "BIC", init_method = "ols") {
  n <- nrow(X)
  p <- ncol(X)
  n_obs <- n - lags

  # Contemporaneous block X_t for t = (lags + 1) .. n (the regression targets).
  Y_full <- X[(lags + 1L):n, , drop = FALSE]
  # Lagged design [X_{t-1}, ..., X_{t-lags}] stacked column-wise (same layout
  # as fit_var_ols), so lag-k variables occupy columns ((k-1)*p+1):(k*p).
  Z <- build_lag_matrix(X, lags)

  B <- array(0, dim = c(lags + 1L, p, p))
  for (i in seq_len(p)) {
    # Number of variables ahead of i in the causal order = its ancestors.
    co_no <- which(causal_order == i) - 1L
    ancestors <- if (co_no >= 1L) causal_order[seq_len(co_no)] else integer(0)

    obj <- Y_full[, i]
    # Predictors: contemporaneous ancestors first, then all lagged variables.
    exp_mat <- cbind(Y_full[, ancestors, drop = FALSE], Z)
    coefs <- fit_adaptive_lasso(obj, exp_mat, lambda = lambda, init_method = init_method)

    # First co_no coefficients map back to B0[i, ancestors].
    if (co_no >= 1L) B[1, i, ancestors] <- coefs[seq_len(co_no)]
    # The remaining coefficients map back to each lagged matrix B_k[i, ].
    for (k in seq_len(lags)) {
      idx <- co_no + (k - 1L) * p + seq_len(p)
      B[k + 1L, i, ] <- coefs[idx]
    }
  }
  B
}


#' Select the VAR lag order by information criterion
#'
#' All candidate lag orders are compared on a **common sample**: the first
#' `max_lag` observations are dropped for every candidate so that each VAR(lag)
#' is estimated over the same time window (t = max_lag + 1 .. n). This mirrors
#' statsmodels' `VAR.select_order` and makes the criteria comparable across
#' lags (otherwise a longer lag would be scored on fewer observations).
#'
#' @param X numeric matrix (n_samples x n_features)
#' @param max_lag maximum lag order to consider
#' @param criterion "bic", "aic", "hqic", or "fpe"
#' @return the selected lag order (integer)
#' @keywords internal
select_var_lag <- function(X, max_lag, criterion = "bic") {
  n <- nrow(X)
  p <- ncol(X)
  # Common sample size shared by all candidate lags.
  n_obs <- n - max_lag
  if (n_obs < 1) {
    stop("Not enough observations for the requested maximum lag order.", call. = FALSE)
  }
  # Response is fixed to the common window for every candidate.
  Y <- X[(max_lag + 1L):n, , drop = FALSE]

  best_lag <- 1L
  best_ic <- Inf
  for (lag in seq_len(max_lag)) {
    # Build the lagged design over the same window: [X_{t-1}, ..., X_{t-lag}].
    Z <- matrix(0, nrow = n_obs, ncol = lag * p)
    for (k in seq_len(lag)) {
      Z[, ((k - 1L) * p + 1L):(k * p)] <- X[(max_lag + 1L - k):(n - k), , drop = FALSE]
    }
    resid <- stats::.lm.fit(Z, Y)$residuals
    # Residual covariance (MLE scaling, divide by n_obs) as used by the criteria.
    sigma <- crossprod(resid) / n_obs
    log_det <- as.numeric(determinant(sigma, logarithm = TRUE)$modulus)
    # Total free parameters; per-equation count is lag * p (no intercept).
    n_params <- lag * p * p
    n_params_eq <- lag * p
    ic <- switch(criterion,
      bic  = log_det + (log(n_obs) / n_obs) * n_params,
      aic  = log_det + (2 / n_obs) * n_params,
      hqic = log_det + (2 * log(log(n_obs)) / n_obs) * n_params,
      # Final Prediction Error: det(sigma) * ((n + k_eq) / (n - k_eq))^p
      fpe  = exp(log_det) * ((n_obs + n_params_eq) / (n_obs - n_params_eq))^p
    )
    if (is.finite(ic) && ic < best_ic) {
      best_ic <- ic
      best_lag <- lag
    }
  }
  best_lag
}


#' Print method for VARLiNGAMResult
#'
#' @param x VARLiNGAMResult object
#' @param digits number of digits to display
#' @param ... additional arguments (unused)
#' @export
print.VARLiNGAMResult <- function(x, digits = 3, ...) {
  n <- length(x$causal_order)
  var_names <- dimnames(x$adjacency_matrices)[[2]]
  order_labels <- if (!is.null(var_names)) {
    var_names[x$causal_order]
  } else {
    paste0("x", x$causal_order - 1L)
  }
  cat("VAR-LiNGAM Result\n")
  cat(sprintf("  Variables : %d\n", n))
  cat(sprintf("  Lag order : %d\n", x$lags))
  cat(sprintf("  Causal order (instantaneous): %s\n", paste(order_labels, collapse = " -> ")))
  cat("\nInstantaneous adjacency matrix B0 (row = to, col = from):\n")
  print(round(x$adjacency_matrices[1, , ], digits = digits))
  # Also show each lagged matrix B_k so the full model is visible at a glance.
  for (k in seq_len(x$lags)) {
    cat(sprintf("\nLagged adjacency matrix B%d (row = to, col = from):\n", k))
    print(round(x$adjacency_matrices[k + 1L, , ], digits = digits))
  }
  invisible(x)
}


#' Generate sample data from a VAR-LiNGAM model
#'
#' Generates a 3-variable time series following a VAR-LiNGAM model with a
#' strictly acyclic instantaneous structure B0, a lag-1 coefficient matrix M1,
#' and non-Gaussian (uniform) errors.
#'
#' @param n number of time points to return (after burn-in)
#' @param seed random seed (NULL allowed)
#' @return list with `data` (data frame, n x 3), `true_B0` (instantaneous
#'   matrix), and `true_M1` (lag-1 coefficient matrix)
#' @importFrom stats runif
#' @export
#' @examples
#' sample <- generate_varlingam_sample(n = 500, seed = 1)
#' head(sample$data)
generate_varlingam_sample <- function(n = 1000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  p <- 3L

  # instantaneous structure: x0 -> x1 -> x2 (strictly acyclic, lower-triangular)
  B0 <- matrix(0, p, p)
  B0[2, 1] <- 0.6
  B0[3, 2] <- -0.5

  # lag-1 coefficients
  M1 <- matrix(0, p, p)
  diag(M1) <- c(0.4, 0.3, 0.5)
  M1[1, 3] <- 0.3

  burn_in <- 100L
  total <- n + burn_in
  # non-Gaussian errors (uniform)
  e <- matrix(stats::runif(total * p, min = -1, max = 1), nrow = total, ncol = p)

  X <- matrix(0, nrow = total, ncol = p)
  ib0_inv <- solve(diag(p) - B0)
  for (t in 2:total) {
    X[t, ] <- ib0_inv %*% (M1 %*% X[t - 1L, ] + e[t, ])
  }
  X <- X[(burn_in + 1L):total, , drop = FALSE]
  colnames(X) <- paste0("x", seq_len(p) - 1L)

  list(
    data = as.data.frame(X),
    true_B0 = B0,
    true_M1 = M1
  )
}
