# =============================================================================
# Direct LiNGAM - Adjacency matrix estimation and regression backends (OLS / LASSO / Adaptive LASSO)
# Based on the Python implementation from the LiNGAM Project
# https://sites.google.com/view/sshimizu06/lingam
# https://github.com/cdt15/lingam
#
# License: MIT + file LICENSE
#
# Original work:
#   Copyright (c) 2019 T.Ikeuchi, G.Haraoka, M.Ide, W.Kurebayashi, S.Shimizu
#
# Portions of this work:
#   Copyright (c) 2026 O.Morimoto
# =============================================================================


# Search range for the lambda values passed to glmnet.
# Specified explicitly because the default automatic generation can stop early.
# Shared by fit_lasso() and fit_adaptive_lasso().
lasso_lambda_seq <- exp(seq(2, -7, length.out = 80))

# Lambda grid for Ridge. Set wider because larger values than LASSO can be optimal.
ridge_lambda_seq <- exp(seq(6, -7, length.out = 100))


#' Check whether glmnet is available
#'
#' If it is not available, raise an error indicating which regression method
#' required it.
#'
#' @param method name of the regression method that requires glmnet (for the
#'   error message)
#' @keywords internal
check_glmnet_available <- function(method) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop(sprintf(
      "Package 'glmnet' is required for reg_method = \"%s\". Please install it.",
      method
    ), call. = FALSE)
  }
}


#' Estimate the adjacency matrix from a causal order
#'
#' @param X original data
#' @param causal_order causal order (vector of 1-based indices)
#' @param prior_knowledge prior-knowledge matrix (NULL allowed)
#' @param method regression method
#'   "ols"           : ordinary least squares (default)
#'   "lasso"         : LASSO regression (glmnet)
#'   "adaptive_lasso": Adaptive LASSO (two-stage)
#'   "ridge"         : Ridge regression (glmnet)
#' @param init_method estimation method for the initial weights of Adaptive LASSO
#'   "ols"   : ordinary least squares (default)
#'   "ridge" : Ridge regression
#' @param lambda LASSO penalty (NULL = automatic selection by cross-validation)
#'   "lambda.min" : minimum prediction error
#'   "lambda.1se" : 1SE rule (sparser)
#'   "AIC"       : minimum AIC (no CV, fast)
#'   "BIC"        : minimum BIC (no CV, fast, sparsest), default
#'   "oracle"     : Adaptive LASSO only. Not usable with Ridge.
#' @return adjacency matrix B (n_features x n_features)
#' @keywords internal
estimate_adjacency_matrix <- function(X,
                                      causal_order,
                                      prior_knowledge = NULL,
                                      method = "adaptive_lasso",
                                      lambda = "BIC",
                                      init_method = "ols") {
  valid_methods <- c("ols", "lasso", "adaptive_lasso", "ridge")
  if (!(method %in% valid_methods)) {
    stop(sprintf(
      "'method' must be one of: %s.",
      paste(valid_methods, collapse = ", ")
    ))
  }

  # Check glmnet (not needed for OLS, which runs with base R only)
  if (method != "ols") {
    check_glmnet_available(method)
  }

  n_features <- ncol(X)
  B <- matrix(0, nrow = n_features, ncol = n_features)

  for (idx in seq_along(causal_order)) {
    target <- causal_order[idx]
    if (idx == 1) next

    # variables that precede this variable in the causal order
    predictors <- causal_order[1:(idx - 1)]

    # constrain by prior knowledge
    if (!is.null(prior_knowledge)) {
      keep <- sapply(predictors, function(p) {
        val <- prior_knowledge[target, p]
        is.na(val) || val != 0
      })
      predictors <- predictors[keep]
    }

    if (length(predictors) == 0) next

    y <- X[, target]
    Xp <- X[, predictors, drop = FALSE]

    B[target, predictors] <- fit_coef_by_method(y, Xp, method, lambda, init_method)
  }

  return(B)
}


#' Dispatch a single regression to the backend selected by `method`
#'
#' Central dispatcher shared by every place that fits "y on Xp with the
#' user-chosen regression method" (adjacency estimation, total effects,
#' Parce/RCD variants). Callers are expected to have validated `method`,
#' `lambda`, and `init_method` already; no validation happens here so that
#' error behavior stays with the caller.
#'
#' @param y response variable (numeric vector)
#' @param Xp predictor matrix
#' @param method one of "ols", "lasso", "adaptive_lasso", "ridge"
#' @param lambda lambda selection rule (ignored for OLS)
#' @param init_method initial estimator for adaptive LASSO (ignored otherwise)
#' @return coefficient vector (excluding intercept)
#' @keywords internal
fit_coef_by_method <- function(y, Xp, method, lambda, init_method) {
  switch(method,
    "ols"            = fit_ols(y, Xp),
    "lasso"          = fit_lasso(y, Xp, lambda = lambda),
    "adaptive_lasso" = fit_adaptive_lasso(y, Xp,
                                          lambda = lambda,
                                          init_method = init_method),
    "ridge"          = fit_ridge_reg(y, Xp, lambda = lambda)
  )
}


#' OLS regression
#'
#' @param y response variable (numeric vector)
#' @param Xp predictor matrix
#' @return coefficient vector (excluding intercept)
#' @keywords internal
fit_ols <- function(y, Xp) {
  fit <- stats::lm.fit(x = cbind(1, as.matrix(Xp)), y = y)
  fit$coefficients[-1]
}


#' OLS fit for a single predictor, pruned by information criterion
#'
#' glmnet requires at least two predictor columns, so the penalized methods
#' fall back to OLS when only one predictor remains. Plain OLS never yields an
#' exact zero, which would make single-predictor edges unprunable: the second
#' variable in the causal order always has exactly one predictor, so a
#' spurious edge would survive even for fully independent data. To preserve
#' the sparse behavior of the penalized methods, the OLS coefficient is kept
#' only when adding the predictor improves the information criterion over the
#' intercept-only model; otherwise it is set to exactly zero.
#'
#' The criterion is AIC for `lambda = "AIC"` and BIC otherwise (the CV /
#' oracle lambdas have no single-predictor counterpart, so the sparsest
#' criterion, BIC, is used for them as well).
#'
#' @param y response variable (numeric vector)
#' @param Xp single-column predictor matrix
#' @param lambda lambda selection method of the calling fit
#' @return length-1 coefficient vector (0 when the predictor is pruned)
#' @keywords internal
fit_ols_ic_pruned <- function(y, Xp, lambda) {
  fit <- stats::lm.fit(x = cbind(1, as.matrix(Xp)), y = y)
  coefs <- fit$coefficients[-1]
  n <- length(y)
  rss_full <- sum(fit$residuals^2)
  rss_null <- sum((y - mean(y))^2)
  # Degenerate fits (zero residual variance on either side) carry no usable
  # IC information; keep the OLS estimate as-is.
  if (rss_full <= 0 || rss_null <= 0) return(coefs)
  penalty <- if (identical(lambda, "AIC")) 2 else log(n)
  ic_full <- n * log(rss_full / n) + penalty
  ic_null <- n * log(rss_null / n)
  if (ic_full < ic_null) coefs else 0
}


#' Select lambda by information criterion
#'
#' @param glmnet_model a glmnet model object
#' @return list with lambda_AIC_best, lambda_BIC_best, idx_AIC_best,
#'   idx_BIC_best, ic_table
#' @keywords internal
ic_glmnet <- function(glmnet_model) {
  tLL <- glmnet_model$nulldev - deviance(glmnet_model)
  k <- glmnet_model$df
  n <- glmnet_model$nobs
  AIC <- -tLL + 2 * k + 2 * k * (k + 1) / pmax(n - k - 1, 1)
  BIC <- log(n) * k - tLL
  ic_table <- data.frame(
    lambda = glmnet_model$lambda,
    df     = k,
    AIC   = AIC,
    BIC    = BIC
  )
  idx_AIC_best <- which.min(ic_table$AIC)
  idx_BIC_best <- which.min(ic_table$BIC)
  list(
    lambda_AIC_best = ic_table$lambda[idx_AIC_best],
    lambda_BIC_best  = ic_table$lambda[idx_BIC_best],
    idx_AIC_best     = idx_AIC_best,
    idx_BIC_best     = idx_BIC_best,
    ic_table         = ic_table
  )
}


#' Scale factor used to make the (otherwise fixed, absolute) lambda search
#' grids adapt to the response's natural scale.
#'
#' `lasso_lambda_seq` / `ridge_lambda_seq` are fixed absolute grids. Because
#' glmnet's `standardize = TRUE` only standardizes the predictors (not the
#' response `y`), the penalty strength needed for meaningful shrinkage scales
#' with the magnitude of `y`. Without this scaling, multiplying the whole
#' input data by a constant changes which edges the default
#' `adaptive_lasso` + `lambda = "BIC"`/`"AIC"` pipeline selects, even though
#' the underlying relationships are identical up to that constant.
#' @param y response variable (numeric vector)
#' @return a positive scale factor
#' @keywords internal
lambda_scale_factor <- function(y) {
  s <- sd_pop(y)
  if (s < 1e-10) 1e-10 else s
}


#' Penalized regression via glmnet (IC or CV lambda selection)
#'
#' Internal helper shared by [fit_lasso()] and [fit_ridge_reg()]. Both
#' functions differ only in `alpha` and `lambda_seq`; this function
#' encapsulates the duplicated IC / CV branches.
#'
#' @param y response variable (numeric vector)
#' @param Xp_mat predictor matrix (already coerced to matrix)
#' @param alpha glmnet mixing parameter: 1 = LASSO, 0 = Ridge
#' @param lambda lambda selection method ("AIC", "BIC", "lambda.min", "lambda.1se")
#' @param lambda_seq numeric vector of (relative) lambda values, scaled internally
#'   by [lambda_scale_factor()] to the response's natural scale before use
#' @return coefficient vector (excluding intercept)
#' @keywords internal
fit_penalized_regression <- function(y, Xp_mat, alpha, lambda, lambda_seq) {
  lambda_seq <- lambda_seq * lambda_scale_factor(y)
  if (lambda %in% c("AIC", "BIC")) {
    fit <- glmnet::glmnet(
      x = Xp_mat, y = y,
      alpha = alpha, intercept = TRUE, standardize = TRUE,
      lambda = lambda_seq
    )
    ic <- ic_glmnet(fit)
    # The selected lambda is a single point on fit$lambda, so extract it
    # directly by column index instead of going through the interpolation in
    # coef(fit, s = ...) (the result is identical and faster).
    k_best <- if (lambda == "AIC") ic$idx_AIC_best else ic$idx_BIC_best
    return(as.numeric(fit$beta[, k_best]))
  }

  cv_fit <- glmnet::cv.glmnet(
    x = Xp_mat, y = y,
    alpha = alpha, intercept = TRUE, standardize = TRUE,
    lambda = lambda_seq
  )
  lambda_val <- cv_fit[[lambda]]
  return(as.numeric(stats::coef(cv_fit, s = lambda_val))[-1])
}


#' LASSO regression (lambda selection by information criterion or CV)
#'
#' @param y response variable
#' @param Xp predictor matrix
#' @param lambda lambda selection method
#'   "lambda.min" : minimum CV prediction error
#'   "lambda.1se" : CV 1SE rule
#'   "AIC"       : minimum AIC
#'   "BIC"        : minimum BIC, default
#' @return coefficient vector
#' @keywords internal
fit_lasso <- function(y, Xp, lambda = "BIC") {
  # glmnet needs >= 2 columns; use the IC-pruned OLS fallback so that a lone
  # predictor can still be shrunk to exactly zero (see fit_ols_ic_pruned).
  if (ncol(Xp) == 1) return(fit_ols_ic_pruned(y, Xp, lambda))
  check_glmnet_available("lasso")
  Xp_mat <- as.matrix(Xp)
  fit_penalized_regression(y, Xp_mat, alpha = 1, lambda = lambda, lambda_seq = lasso_lambda_seq)
}


#' Ridge regression (lambda selection by information criterion or CV)
#'
#' @param y response variable
#' @param Xp predictor matrix
#' @param lambda lambda selection method
#'   "lambda.min" : minimum CV prediction error
#'   "lambda.1se" : CV 1SE rule
#'   "AIC"       : minimum AIC
#'   "BIC"        : minimum BIC, default
#'   "oracle" is not usable (Adaptive LASSO only).
#' @return coefficient vector
#' @keywords internal
fit_ridge_reg <- function(y, Xp, lambda = "BIC") {
  # Ridge performs no sparse estimation, so the plain OLS fallback (without
  # the IC pruning used by the lasso-family fits) is the consistent choice.
  if (ncol(Xp) == 1) return(fit_ols(y, Xp))
  if (lambda == "oracle") {
    stop("lambda = \"oracle\" is only supported for reg_method = \"adaptive_lasso\".",
         call. = FALSE)
  }
  check_glmnet_available("ridge")
  Xp_mat <- as.matrix(Xp)
  fit_penalized_regression(y, Xp_mat, alpha = 0, lambda = lambda, lambda_seq = ridge_lambda_seq)
}


#' Adaptive LASSO
#' @param y response variable
#' @param Xp predictor matrix
#' @param lambda lambda selection method ("lambda.min", "lambda.1se", "AIC", "BIC", "oracle")
#' @param gamma_weight exponent of the weights
#' @param init_method estimation method for the initial weights ("ols" or "ridge")
#' @return coefficient vector
#' @keywords internal
fit_adaptive_lasso <- function(y, Xp,
                               lambda = "BIC",
                               gamma_weight = 1.0,
                               init_method = "ols") {
  # glmnet needs >= 2 columns; use the IC-pruned OLS fallback so that a lone
  # predictor can still be shrunk to exactly zero (see fit_ols_ic_pruned).
  if (ncol(Xp) == 1) return(fit_ols_ic_pruned(y, Xp, lambda))
  check_glmnet_available("adaptive_lasso")

  Xp_mat <- as.matrix(Xp)
  n <- nrow(Xp_mat) # obtain the sample size n

  # --- Step 1: compute the initial estimator (in the original scale) ---
  if (init_method == "ols") {
    init_fit <- stats::lm.fit(x = cbind(1, Xp_mat), y = y)
    init_coefs <- as.numeric(init_fit$coefficients[-1])
  } else {
    ridge_cv <- glmnet::cv.glmnet(
      x = Xp_mat, y = y, alpha = 0,
      intercept = TRUE, standardize = TRUE
    )
    init_coefs <- as.numeric(stats::coef(ridge_cv, s = "lambda.min"))[-1]
  }
  init_coefs[is.na(init_coefs)] <- 0

  # --- Step 2: compute penalty.factor ---
  x_sds <- apply(Xp_mat, 2, sd_pop)
  x_sds[x_sds < 1e-10] <- 1e-10

  init_coefs_std <- init_coefs * x_sds

  pf <- 1 / (abs(init_coefs_std)^gamma_weight)
  pf[is.infinite(pf) | is.na(pf)] <- 1e10

  # Scale the (otherwise fixed, absolute) lambda grid to the response's
  # natural scale; see lambda_scale_factor() for why this is needed.
  y_scale <- lambda_scale_factor(y)
  lambda_seq_scaled <- lasso_lambda_seq * y_scale

  # --- Step 3: run Adaptive LASSO ---
  fit <- glmnet::glmnet(
    x = Xp_mat, y = y, alpha = 1,
    intercept = TRUE, standardize = TRUE,
    penalty.factor = pf,
    lambda = lambda_seq_scaled
  )

  if (lambda %in% c("AIC", "BIC")) {
    ic <- ic_glmnet(fit)
    # A single point on fit$lambda, so extract it directly by column index
    # without interpolation.
    k_best <- if (lambda == "AIC") ic$idx_AIC_best else ic$idx_BIC_best
    coef_vec <- as.numeric(fit$beta[, k_best])
    return(coef_vec)
  }

  if (lambda == "oracle") {
    # The oracle lambda is not on the search grid, so interpolate with coef().
    # Scaled by y_scale for the same reason as the grid itself, so it stays
    # in the same units as `fit$lambda` for coef()'s interpolation.
    lambda_val <- (5 / (n^(1.75))) * y_scale
  } else {
    cv_fit <- glmnet::cv.glmnet(
      x = Xp_mat, y = y, alpha = 1,
      intercept = TRUE, standardize = TRUE,
      penalty.factor = pf,
      lambda = lambda_seq_scaled
    )

    lambda_val <- cv_fit[[lambda]]
  }

  coef_vec <- as.numeric(stats::coef(fit, s = lambda_val))[-1]

  return(coef_vec)
}
