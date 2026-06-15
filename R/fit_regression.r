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

    # --- branch by regression method ---
    coefs <- switch(method,
      "ols"            = fit_ols(y, Xp),
      "lasso"          = fit_lasso(y, Xp, lambda = lambda),
      "adaptive_lasso" = fit_adaptive_lasso(y, Xp,
                                            lambda = lambda,
                                            init_method = init_method),
      "ridge"          = fit_ridge_reg(y, Xp, lambda = lambda)
    )

    B[target, predictors] <- coefs
  }

  return(B)
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
  if (ncol(Xp) == 1) {
    return(fit_ols(y, Xp))
  }

  check_glmnet_available("lasso")

  Xp_mat <- as.matrix(Xp)

  if (lambda %in% c("AIC", "BIC")) {
    fit <- glmnet::glmnet(
      x = Xp_mat, y = y,
      alpha = 1, intercept = TRUE, standardize = TRUE,
      lambda = lasso_lambda_seq
    )

    ic <- ic_glmnet(fit)
    # The selected lambda is a single point on fit$lambda, so extract it
    # directly by column index instead of going through the interpolation in
    # coef(fit, s = ...) (the result is identical and faster).
    k_best <- if (lambda == "AIC") ic$idx_AIC_best else ic$idx_BIC_best
    coef_vec <- as.numeric(fit$beta[, k_best])

  } else {
    cv_fit <- glmnet::cv.glmnet(
      x = Xp_mat, y = y,
      alpha = 1, intercept = TRUE, standardize = TRUE,
      lambda = lasso_lambda_seq
    )

    lambda_val <- cv_fit[[lambda]]
    coef_vec <- as.numeric(stats::coef(cv_fit, s = lambda_val))[-1]
  }

  return(coef_vec)
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
  if (ncol(Xp) == 1) {
    return(fit_ols(y, Xp))
  }

  if (lambda == "oracle") {
    stop("lambda = \"oracle\" is only supported for reg_method = \"adaptive_lasso\".",
         call. = FALSE)
  }

  check_glmnet_available("ridge")

  Xp_mat <- as.matrix(Xp)

  if (lambda %in% c("AIC", "BIC")) {
    fit <- glmnet::glmnet(
      x = Xp_mat, y = y,
      alpha = 0, intercept = TRUE, standardize = TRUE,
      lambda = ridge_lambda_seq
    )

    ic <- ic_glmnet(fit)
    k_best <- if (lambda == "AIC") ic$idx_AIC_best else ic$idx_BIC_best
    coef_vec <- as.numeric(fit$beta[, k_best])

  } else {
    cv_fit <- glmnet::cv.glmnet(
      x = Xp_mat, y = y,
      alpha = 0, intercept = TRUE, standardize = TRUE,
      lambda = ridge_lambda_seq
    )

    lambda_val <- cv_fit[[lambda]]
    coef_vec <- as.numeric(stats::coef(cv_fit, s = lambda_val))[-1]
  }

  return(coef_vec)
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
  if (ncol(Xp) == 1) return(fit_ols(y, Xp))
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

  # --- Step 3: run Adaptive LASSO ---
  fit <- glmnet::glmnet(
    x = Xp_mat, y = y, alpha = 1,
    intercept = TRUE, standardize = TRUE,
    penalty.factor = pf,
    lambda = lasso_lambda_seq
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
    lambda_val <- 5 / (n^(1.75))
  } else {
    cv_fit <- glmnet::cv.glmnet(
      x = Xp_mat, y = y, alpha = 1,
      intercept = TRUE, standardize = TRUE,
      penalty.factor = pf,
      lambda = lasso_lambda_seq
    )

    lambda_val <- cv_fit[[lambda]]
  }

  coef_vec <- as.numeric(stats::coef(fit, s = lambda_val))[-1]

  return(coef_vec)
}
