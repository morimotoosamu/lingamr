#' Compute residuals (error terms) of a LiNGAM model
#'
#' After validating the inputs (that `lingam_result` is a `LingamResult`,
#' that X is numeric, and that the dimensions match), returns `E = X - X B^T`.
#' Shared by the residual-based diagnostic functions.
#'
#' @param X original data (matrix or data.frame)
#' @param lingam_result return value of [lingam_direct()]
#' @return residual matrix (n_samples x n_features). Preserves the column names of X.
#' @keywords internal
lingam_residuals <- function(X, lingam_result) {
  validate_lingam_result(lingam_result)
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be a numeric matrix or data.frame.", call. = FALSE)
  B <- lingam_result$adjacency_matrix
  if (ncol(B) != ncol(X)) {
    stop(
      "X has ", ncol(X), " variables but lingam_result was estimated from ",
      ncol(B), " variables.",
      call. = FALSE
    )
  }
  X - X %*% t(B)
}


#' Skewness (divided by n)
#' @param x numeric vector
#' @keywords internal
skewness_pop <- function(x) {
  n <- length(x)
  x_c <- x - mean(x)
  (sum(x_c^3) / n) / (sum(x_c^2) / n)^(3 / 2)
}


#' Kurtosis (divided by n; excess kurtosis, which is 0 for a normal distribution)
#' @param x numeric vector
#' @keywords internal
kurtosis_pop <- function(x) {
  n <- length(x)
  x_c <- x - mean(x)
  (sum(x_c^4) / n) / (sum(x_c^2) / n)^2 - 3
}


#' Compute p-values for the independence test of the errors
#'
#' @param X original data (matrix or data.frame)
#' @param lingam_result return value of lingam_direct()
#' @param method type of correlation coefficient ("spearman", "pearson", "kendall")
#' @return matrix of p-values (n_features x n_features)
#' @importFrom stats cor.test
#' @export
#' @examples
#' # Load the sample data
#' LiNGAM_sample_1000 <- generate_lingam_sample_6()
#'
#' # Run Direct LiNGAM
#' result <- LiNGAM_sample_1000$data |>
#'   lingam_direct()
#'
#' # Compute p-values (default: Spearman)
#' p_vals <- get_error_independence_p_values(LiNGAM_sample_1000$data, result)
#' round(p_vals, 3)
#'
#' # Compute with Kendall
#' p_vals_k <- get_error_independence_p_values(LiNGAM_sample_1000$data, result, method = "kendall")
#' round(p_vals_k, 3)
get_error_independence_p_values <- function(X, lingam_result, method = "spearman") {
  method <- match.arg(method, c("spearman", "pearson", "kendall"))

  # Compute residuals (error terms) (with input validation)
  E <- lingam_residuals(X, lingam_result)
  n_features <- ncol(E)

  # Generate indices for all pairs (excluding the diagonal)
  pairs <- which(upper.tri(matrix(TRUE, n_features, n_features)), arr.ind = TRUE)

  # Compute all upper-triangular pairs at once
  p_upper <- apply(pairs, 1, function(idx) {
    stats::cor.test(E[, idx[1]], E[, idx[2]], method = method)$p.value
  })

  # Store in a symmetric matrix
  p_values <- matrix(NA, nrow = n_features, ncol = n_features)
  p_values[pairs] <- p_upper
  p_values[pairs[, 2:1]] <- p_upper # copy to the lower triangle

  colnames(p_values) <- rownames(p_values) <- colnames(E)

  return(p_values)
}


#' Test normality of residuals from Direct LiNGAM
#'
#' Calculate residuals (error terms) from the estimated adjacency matrix
#' and test their normality. Since LiNGAM assumes non-Gaussian errors,
#' rejecting normality (small p-value) supports the LiNGAM model assumption.
#'
#' @param X original data matrix or data.frame
#' @param lingam_result result from lingam_direct()
#' @param method normality test method
#'   "shapiro"      : Shapiro-Wilk test (default, n <= 5000)
#'   "ks"           : Kolmogorov-Smirnov test (n > 5000)
#'   "ad"           : Anderson-Darling test (requires nortest package)
#'   "lillie"       : Lilliefors test (requires nortest package)
#'   "jb"           : Jarque-Bera test (requires tseries package)
#' @param alpha significance level (default: 0.05)
#' @return data.frame with test results for each variable
#' @export
#' @examples
#' # Load the sample data
#' LiNGAM_sample_1000 <- generate_lingam_sample_6()
#'
#' # Run Direct LiNGAM
#' result <- lingam_direct(LiNGAM_sample_1000$data)
#'
#' # Shapiro-Wilk (default)
#' test_residual_normality(LiNGAM_sample_1000$data, result)
test_residual_normality <- function(X, lingam_result,
                                    method = "shapiro",
                                    alpha = 0.05) {
  # --- Calculate residuals: e = X - X %*% t(B) (with input validation) ---
  E <- lingam_residuals(X, lingam_result)
  n_features <- ncol(E)
  n_samples  <- nrow(E)

  # --- Variable names ---
  var_names <- colnames(E)
  if (is.null(var_names)) var_names <- paste0("x", seq_len(n_features) - 1)

  # --- Method validation ---
  valid_methods <- c("shapiro", "ks", "ad", "lillie", "jb")
  if (!(method %in% valid_methods)) {
    stop(sprintf("method must be one of: %s",
                 paste(valid_methods, collapse = ", ")))
  }

  # --- Sample size warnings ---
  if (method == "shapiro" && n_samples > 5000) {
    warning("Shapiro-Wilk test may not work well with n > 5000. ",
            "Consider using 'ad' or 'lillie' instead.")
  }

  # --- Package check ---
  if (method %in% c("ad", "lillie")) {
    if (!requireNamespace("nortest", quietly = TRUE)) {
      stop("Package 'nortest' is required for method '", method, "'.\n",
           "Install it with: install.packages('nortest')")
    }
  }
  if (method == "jb") {
    if (!requireNamespace("tseries", quietly = TRUE)) {
      stop("Package 'tseries' is required for method 'jb'.\n",
           "Install it with: install.packages('tseries')")
    }
  }

  # --- Test function selector ---
  test_fn <- switch(method,
                    "shapiro" = function(x) {
                      if (length(x) > 5000) {
                        # Subsample for shapiro
                        x <- sample(x, 5000)
                      }
                      stats::shapiro.test(x)
                    },
                    "ks"      = function(x) stats::ks.test(x, "pnorm", mean(x), sd(x)),
                    "ad"      = function(x) nortest::ad.test(x),
                    "lillie"  = function(x) nortest::lillie.test(x),
                    "jb"      = function(x) tseries::jarque.bera.test(x)
  )

  # --- Run tests for each variable ---
  results <- data.frame(
    variable     = var_names,
    method       = method,
    statistic    = NA_real_,
    p_value      = NA_real_,
    is_normal    = NA,
    is_non_gauss = NA,
    mean         = NA_real_,
    sd           = NA_real_,
    skewness     = NA_real_,
    kurtosis     = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(n_features)) {
    e_i <- E[, i]

    # Basic stats
    results$mean[i]     <- mean(e_i)
    results$sd[i]       <- sd(e_i)
    results$skewness[i] <- skewness_pop(e_i)
    results$kurtosis[i] <- kurtosis_pop(e_i)

    # Normality test
    test_result <- tryCatch(
      test_fn(e_i),
      error = function(err) {
        warning(sprintf("Test failed for %s: %s", var_names[i], err$message))
        list(statistic = NA, p.value = NA)
      }
    )

    results$statistic[i] <- as.numeric(test_result$statistic)
    results$p_value[i]   <- test_result$p.value
    results$is_normal[i] <- test_result$p.value > alpha
    results$is_non_gauss[i] <- test_result$p.value <= alpha
  }

  # --- Add summary ---
  attr(results, "alpha") <- alpha
  attr(results, "n_samples") <- n_samples
  attr(results, "n_features") <- n_features
  attr(results, "n_non_gaussian") <- sum(results$is_non_gauss, na.rm = TRUE)
  class(results) <- c("lingam_normality_test", "data.frame")

  return(results)
}


#' Print method for lingam_normality_test
#' @param x lingam_normality_test object
#' @param ... additional arguments
#' @export
print.lingam_normality_test <- function(x, ...) {
  alpha <- attr(x, "alpha")
  n_features <- attr(x, "n_features")
  n_non_gauss <- attr(x, "n_non_gaussian")

  cat("=== Residual Normality Test ===\n")
  cat(sprintf("Method:         %s\n", x$method[1]))
  cat(sprintf("Sample size:    %d\n", attr(x, "n_samples")))
  cat(sprintf("Significance:   %.3f\n", alpha))
  cat(sprintf("Non-Gaussian:   %d / %d variables\n\n", n_non_gauss, n_features))

  # --- Build a data.frame for display (with types explicitly set to character) ---
  display_df <- data.frame(
    variable     = x$variable,
    statistic    = sprintf("%.4f", x$statistic),
    p_value      = sapply(x$p_value, function(pv) {
      if (is.na(pv)) return("NA")
      if (pv < 1e-16) return("< 2.2e-16")
      if (pv < 0.0001) return(sprintf("%.2e", pv))
      sprintf("%.4f", pv)
    }),
    is_non_gauss = x$is_non_gauss,
    skewness     = sprintf("%.3f", x$skewness),
    kurtosis     = sprintf("%.3f", x$kurtosis),
    stringsAsFactors = FALSE
  )

  print(display_df, row.names = FALSE)

  cat("\nInterpretation:\n")
  cat("  is_non_gauss = TRUE  -> rejects normality (supports LiNGAM assumption)\n")
  cat("  is_non_gauss = FALSE -> cannot reject normality (LiNGAM may not fit)\n")

  if (n_non_gauss < n_features) {
    cat(sprintf("\nWARNING: %d variable(s) appear Gaussian.\n",
                n_features - n_non_gauss))
    cat("  LiNGAM assumes non-Gaussian errors. Results may be unreliable.\n")
  } else {
    cat("\nAll residuals are non-Gaussian. LiNGAM assumption is supported.\n")
  }

  invisible(x)
}

#' plot QQ
#' @param X original data (matrix or data.frame)
#' @param lingam_result return value of lingam_direct()
#' @param ncol Number of columns.
#' @param nrow Number of rows.
#' @export
#' @examples
#' # Load the sample data
#' LiNGAM_sample_1000 <- generate_lingam_sample_6()
#'
#' # Run Direct LiNGAM
#' result <- lingam_direct(LiNGAM_sample_1000$data)
#'
#' plot_residual_qq(LiNGAM_sample_1000$data, result)
plot_residual_qq <- function(X, lingam_result, ncol = 3, nrow = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please install it.", call. = FALSE)
  }

  E <- lingam_residuals(X, lingam_result)

  var_names <- colnames(E)
  if (is.null(var_names)) var_names <- paste0("x", seq_len(ncol(E)) - 1)

  df <- data.frame(
    variable = rep(var_names, each = nrow(E)),
    residual = as.vector(E)
  )

  ggplot2::ggplot(df, ggplot2::aes(sample = residual)) +
    ggplot2::stat_qq() +
    ggplot2::stat_qq_line(color = "red") +
    ggplot2::facet_wrap(~ variable, ncol = ncol, nrow = nrow, scales = "free") +
    ggplot2::labs(
      title = "Q-Q Plot of LiNGAM Residuals",
      subtitle = "Deviation from line indicates non-Gaussianity (good for LiNGAM)",
      x = "Theoretical Quantiles",
      y = "Sample Quantiles"
    ) +
    ggplot2::theme_minimal()
}


