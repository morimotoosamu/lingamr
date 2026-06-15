#' Summarize the goodness-of-fit of a Direct LiNGAM model at once
#'
#' For a fitted Direct LiNGAM model, this verifies how well the two main
#' assumptions on which LiNGAM relies (mutual independence of residuals and
#' non-Gaussianity of residuals) hold, all at once, and displays the results
#' together. Internally it calls [get_error_independence_p_values()] and
#' [test_residual_normality()].
#'
#' Gaussian-likelihood-based criteria such as BIC/AIC are not included because
#' they are theoretically inconsistent with LiNGAM's assumption that "the errors
#' are non-Gaussian". Instead, the verification results of the assumptions
#' themselves are summarized.
#'
#' @param X The original data (matrix or data.frame), the one used to estimate
#'   `lingam_result`.
#' @param lingam_result The return value of [lingam_direct()] (a `LingamResult`
#'   object)
#' @param independence_method The type of correlation coefficient used in the
#'   residual independence test
#'   ("spearman", "pearson", "kendall"). Passed to
#'   [get_error_independence_p_values()].
#' @param normality_method The method for the residual normality test
#'   ("shapiro", "ks", "ad", "lillie", "jb"). Passed to
#'   [test_residual_normality()].
#' @param alpha Significance level (default: 0.05)
#' @return A list of class `lingam_summary`, containing the following elements:
#' * `n_variables`, `n_samples`: Number of variables / number of observations
#' * `causal_order`: Causal order (variable-name labels)
#' * `n_edges`: Number of nonzero elements in the adjacency matrix (number of
#'   estimated edges)
#' * `independence_p_values`: Matrix of p-values from the independence test
#'   between residuals
#' * `n_dependent_pairs`, `n_pairs`: Number of pairs with p < alpha / total
#'   number of pairs
#' * `min_independence_p`: Minimum p-value of the independence test
#' * `normality`: Result of the normality test (a `lingam_normality_test`
#'   object)
#' * `n_non_gaussian`: Number of variables judged to be non-Gaussian
#' * `alpha`, `independence_method`, `normality_method`: The settings used
#' @export
#' @examples
#' LiNGAM_sample_1000 <- generate_lingam_sample_6()
#'
#' model <- lingam_direct(LiNGAM_sample_1000$data, reg_method = "ols")
#'
#' summary_lingam(LiNGAM_sample_1000$data, model)
summary_lingam <- function(X, lingam_result,
                           independence_method = "spearman",
                           normality_method = "shapiro",
                           alpha = 0.05) {
  validate_lingam_result(lingam_result)
  independence_method <- match.arg(independence_method,
                                   c("spearman", "pearson", "kendall"))
  normality_method <- match.arg(normality_method,
                                c("shapiro", "ks", "ad", "lillie", "jb"))
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single number in (0, 1).", call. = FALSE)
  }

  X <- as.matrix(X)
  n_features <- ncol(X)
  if (ncol(lingam_result$adjacency_matrix) != n_features) {
    stop(
      "X has ", n_features, " variables but lingam_result was estimated from ",
      ncol(lingam_result$adjacency_matrix), " variables.",
      call. = FALSE
    )
  }

  B <- lingam_result$adjacency_matrix
  var_names <- colnames(B)
  if (is.null(var_names)) var_names <- paste0("x", seq_len(n_features) - 1)

  # --- Assumption 1: Independence of residuals ---
  p_mat <- get_error_independence_p_values(X, lingam_result,
                                           method = independence_method)
  upper_p <- p_mat[upper.tri(p_mat)]
  n_pairs <- length(upper_p)
  n_dependent_pairs <- sum(upper_p < alpha, na.rm = TRUE)
  min_independence_p <- if (all(is.na(upper_p))) NA_real_ else min(upper_p, na.rm = TRUE)

  # --- Assumption 2: Non-Gaussianity of residuals ---
  norm_df <- test_residual_normality(X, lingam_result,
                                     method = normality_method, alpha = alpha)
  n_non_gaussian <- attr(norm_df, "n_non_gaussian")

  out <- list(
    n_variables           = n_features,
    n_samples             = nrow(X),
    causal_order          = var_names[lingam_result$causal_order],
    n_edges               = sum(abs(B) > 0),
    independence_p_values = p_mat,
    n_pairs               = n_pairs,
    n_dependent_pairs     = n_dependent_pairs,
    min_independence_p    = min_independence_p,
    normality             = norm_df,
    n_non_gaussian        = n_non_gaussian,
    alpha                 = alpha,
    independence_method   = independence_method,
    normality_method      = normality_method
  )
  class(out) <- "lingam_summary"
  out
}


#' print method for lingam_summary
#'
#' @param x A `lingam_summary` object
#' @param ... Additional arguments (unused)
#' @export
print.lingam_summary <- function(x, ...) {
  cat("=== Direct LiNGAM Model Summary ===\n")
  cat(sprintf("Variables:    %d\n", x$n_variables))
  cat(sprintf("Observations: %d\n", x$n_samples))
  cat(sprintf("Edges:        %d\n", x$n_edges))
  cat(sprintf("Causal order: %s\n", paste(x$causal_order, collapse = " -> ")))

  # --- Assumption 1: Independence ---
  cat("\n--- Assumption 1: Independence of residuals ---\n")
  cat(sprintf("Method:           %s\n", x$independence_method))
  cat(sprintf("Dependent pairs:  %d / %d  (p < %.3f)\n",
              x$n_dependent_pairs, x$n_pairs, x$alpha))
  cat(sprintf("Min p-value:      %s\n",
              if (is.na(x$min_independence_p)) "NA" else sprintf("%.4f", x$min_independence_p)))
  if (x$n_dependent_pairs == 0) {
    cat("=> Residuals appear mutually independent (assumption supported).\n")
  } else {
    cat(sprintf("=> WARNING: %d residual pair(s) appear dependent. Model may be misspecified.\n",
                x$n_dependent_pairs))
  }

  # --- Assumption 2: Non-Gaussianity ---
  cat("\n--- Assumption 2: Non-Gaussianity of residuals ---\n")
  cat(sprintf("Method:           %s\n", x$normality_method))
  cat(sprintf("Non-Gaussian:     %d / %d  (p <= %.3f)\n",
              x$n_non_gaussian, x$n_variables, x$alpha))
  if (x$n_non_gaussian == x$n_variables) {
    cat("=> All residuals are non-Gaussian (assumption supported).\n")
  } else {
    cat(sprintf("=> WARNING: %d variable(s) appear Gaussian. LiNGAM may be unreliable.\n",
                x$n_variables - x$n_non_gaussian))
  }

  invisible(x)
}
