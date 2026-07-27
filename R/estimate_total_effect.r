#' Estimate the total causal effect between two specified variables
#'
#' @param X Original data (matrix or data.frame)
#' @param lingam_result Return value of lingam_direct()
#' @param from_index Cause variable (1-based index or variable name)
#' @param to_index Effect variable (1-based index or variable name)
#' @param method Regression method ("ols", "lasso", "adaptive_lasso", "ridge"). Default is adaptive_lasso
#' @param lambda Lambda selection ("lambda.min", "lambda.1se", "AIC", "BIC", "oracle"). Default is BIC
#' @param init_method Method for estimating the initial weights of adaptive LASSO regression ("ols" or "ridge")
#' @return Estimated total causal effect
#' @export
#' @examples
#' LiNGAM_sample_1000 <- generate_lingam_sample_6()
#'
#' model <- LiNGAM_sample_1000$data |>
#'   lingam_direct(reg_method = "ols")
#'
#' LiNGAM_sample_1000$data |>
#'   estimate_total_effect(model, 4, 1)
estimate_total_effect <- function(X, lingam_result, from_index, to_index,
                                  method = "adaptive_lasso", lambda = "BIC",
                                  init_method = "ols") {
  validate_lingam_result(lingam_result)
  method <- match.arg(method, c("adaptive_lasso", "lasso", "ols", "ridge"))
  lambda <- match.arg(lambda, c("BIC", "AIC", "lambda.min", "lambda.1se", "oracle"))
  init_method <- match.arg(init_method, c("ols", "ridge"))

  X <- as.matrix(X)
  col_names <- colnames(X)
  if (!is.numeric(X)) stop("X must be a numeric matrix or data frame.", call. = FALSE)
  if (anyNA(X)) stop("X must not contain missing values (NA).", call. = FALSE)

  n_features <- ncol(X)
  if (ncol(lingam_result$adjacency_matrix) != n_features) {
    stop(
      "X has ", n_features, " variables but lingam_result was estimated from ",
      ncol(lingam_result$adjacency_matrix), " variables.",
      call. = FALSE
    )
  }

  # --- Convert variable name -> index (shared helper) ---
  from_index <- resolve_var_index(from_index, "from_index", col_names, n_features)
  to_index <- resolve_var_index(to_index, "to_index", col_names, n_features)
  if (from_index == to_index) stop("from_index and to_index must differ.")

  adjacency_matrix <- lingam_result$adjacency_matrix
  causal_order <- lingam_result$causal_order

  # --- Check causal order ---
  from_order <- which(causal_order == from_index)
  to_order <- which(causal_order == to_index)

  from_label <- if (!is.null(col_names)) col_names[from_index] else paste0("x", from_index)
  to_label <- if (!is.null(col_names)) col_names[to_index] else paste0("x", to_index)

  if (from_order > to_order) {
    warning(sprintf(
      "Causal order of %s (to) is earlier than %s (from). Result may be incorrect.",
      to_label, from_label
    ))
  }

  # --- Identify parent variables and run regression ---
  parents <- which(abs(adjacency_matrix[from_index, ]) > 0)
  predictors <- unique(c(from_index, parents))
  from_pos <- which(predictors == from_index)

  y <- X[, to_index]
  Xp <- X[, predictors, drop = FALSE]

  coefs <- fit_coef_by_method(y, Xp, method, lambda, init_method)

  return(coefs[from_pos])
}


#' Estimate the total causal effects between all variables at once
#'
#' @param X Original data (n_samples x n_features)
#' @param lingam_result Return value of lingam_direct()
#' @param method Regression method ("ols", "lasso", "adaptive_lasso", "ridge")
#' @param lambda Lambda selection ("lambda.min", "lambda.1se", "AIC", "BIC")
#' @param init_method Method for estimating the initial weights of adaptive LASSO regression ("ols" or "ridge")
#' @return Matrix of total causal effects (n_features x n_features).
#'   **Convention: `TE[i, j]` is the total causal effect from variable j to variable i (j -> i).**
#'   Same index convention as the adjacency matrix `adjacency_matrix`. The sum of direct and indirect effects.
#' @export
#' @examples
#' LiNGAM_sample_1000 <- generate_lingam_sample_6()
#'
#' model <- LiNGAM_sample_1000$data |>
#'   lingam_direct(reg_method = "ols")
#'
#' LiNGAM_sample_1000$data |>
#'   estimate_all_total_effects(model)
estimate_all_total_effects <- function(X,
                                       lingam_result,
                                       method = "adaptive_lasso",
                                       lambda = "BIC",
                                       init_method = "ols") {
  validate_lingam_result(lingam_result)
  method <- match.arg(method, c("adaptive_lasso", "lasso", "ols", "ridge"))
  lambda <- match.arg(lambda, c("BIC", "AIC", "lambda.min", "lambda.1se", "oracle"))
  init_method <- match.arg(init_method, c("ols", "ridge"))

  X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be a numeric matrix or data frame.", call. = FALSE)
  if (anyNA(X)) stop("X must not contain missing values (NA).", call. = FALSE)
  n_features <- ncol(X)
  if (ncol(lingam_result$adjacency_matrix) != n_features) {
    stop(
      "X has ", n_features, " variables but lingam_result was estimated from ",
      ncol(lingam_result$adjacency_matrix), " variables.",
      call. = FALSE
    )
  }

  causal_order <- lingam_result$causal_order
  adj_matrix <- lingam_result$adjacency_matrix

  TE <- matrix(0, nrow = n_features, ncol = n_features)
  if (!is.null(colnames(X))) {
    rownames(TE) <- colnames(X)
    colnames(TE) <- colnames(X)
  }

  # The covariance matrix is loop-invariant, so compute it only once (used only in the OLS path)
  cov_mat <- if (method == "ols") cov(X) else NULL

  for (i in 1:(n_features - 1)) {
    from_idx <- causal_order[i]

    parents <- which(abs(adj_matrix[from_idx, ]) > 0)
    predictors <- unique(c(from_idx, parents))
    from_pos <- which(predictors == from_idx)

    downstream <- causal_order[(i + 1):n_features]

    if (method == "ols") {
      # --- OLS: batch computation based on the covariance matrix (fastest) ---
      cov_xx <- cov_mat[predictors, predictors, drop = FALSE]
      cov_xy <- cov_mat[predictors, downstream, drop = FALSE]
      beta_mat <- solve(cov_xx, cov_xy)
      TE[downstream, from_idx] <- beta_mat[from_pos, ]
    } else {
      # --- LASSO / Adaptive LASSO / Ridge ---
      Xp <- X[, predictors, drop = FALSE]
      for (to_idx in downstream) {
        y <- X[, to_idx]
        # method is guaranteed non-OLS here (OLS takes the cov-batch fast
        # path above), so the dispatcher's "ols" branch is never reached.
        coefs <- fit_coef_by_method(y, Xp, method, lambda, init_method)
        TE[to_idx, from_idx] <- coefs[from_pos]
      }
    }
  }

  return(TE)
}
