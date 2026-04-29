#' 誤差の独立性検定の p 値を計算
#'
#' @param X 元データ (matrix or data.frame)
#' @param lingam_result direct_lingam() の返り値
#' @param method 相関係数の種類 ("spearman", "pearson", "kendall")
#' @return p 値の行列 (n_features x n_features)
#' @importFrom stats cor.test
#' @export
#' @examples
#' # サンプルデータの呼び出し
#' LiNGAM_sample_1000 <- generate_lingam_sample_6()
#'
#' # Direct LiNGAM の実行
#' result <- LiNGAM_sample_1000$data |>
#'   direct_lingam()
#'
#' # p 値の計算（デフォルト: Spearman）
#' p_vals <- get_error_independence_p_values(LiNGAM_sample_1000$data, result)
#' round(p_vals, 3)
#'
#' # Kendall で計算
#' p_vals_k <- get_error_independence_p_values(LiNGAM_sample_1000$data, result, method = "kendall")
#' round(p_vals_k, 3)
get_error_independence_p_values <- function(X, lingam_result, method = "spearman") {
  X <- as.matrix(X)
  B <- lingam_result$adjacency_matrix
  n_features <- ncol(X)

  # 残差（誤差項）の計算
  E <- X - X %*% t(B)

  # 全ペアのインデックスを生成（対角を除く）
  pairs <- which(upper.tri(matrix(TRUE, n_features, n_features)), arr.ind = TRUE)

  # 上三角の全ペアを一括計算
  p_upper <- apply(pairs, 1, function(idx) {
    stats::cor.test(E[, idx[1]], E[, idx[2]], method = method)$p.value
  })

  # 対称行列に格納
  p_values <- matrix(NA, nrow = n_features, ncol = n_features)
  p_values[pairs] <- p_upper
  p_values[pairs[, 2:1]] <- p_upper # 下三角にコピー

  colnames(p_values) <- rownames(p_values) <- colnames(X)

  return(p_values)
}


#' Test normality of residuals from Direct LiNGAM
#'
#' Calculate residuals (error terms) from the estimated adjacency matrix
#' and test their normality. Since LiNGAM assumes non-Gaussian errors,
#' rejecting normality (small p-value) supports the LiNGAM model assumption.
#'
#' @param X original data matrix or data.frame
#' @param lingam_result result from direct_lingam()
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
#' # サンプルデータの呼び出し
#' LiNGAM_sample_1000 <- generate_lingam_sample_6()
#'
#' # Direct LiNGAM の実行
#' result <- direct_lingam(LiNGAM_sample_1000$data)
#'
#' # Shapiro-Wilk (default)
#' test_residual_normality(LiNGAM_sample_1000$data, result)
test_residual_normality <- function(X, lingam_result,
                                    method = "shapiro",
                                    alpha = 0.05) {
  # --- Input validation ---
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be a numeric matrix or data.frame.")

  B <- lingam_result$adjacency_matrix
  n_features <- ncol(X)
  n_samples  <- nrow(X)

  # --- Variable names ---
  var_names <- colnames(X)
  if (is.null(var_names)) var_names <- paste0("x", seq_len(n_features) - 1)

  # --- Calculate residuals: e = X - X %*% t(B) ---
  E <- X - X %*% t(B)

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

  # --- Calculate skewness and kurtosis ---
  skewness <- function(x) {
    n <- length(x)
    x_c <- x - mean(x)
    (sum(x_c^3) / n) / (sum(x_c^2) / n)^(3/2)
  }

  kurtosis <- function(x) {
    n <- length(x)
    x_c <- x - mean(x)
    (sum(x_c^4) / n) / (sum(x_c^2) / n)^2 - 3
  }

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
    results$skewness[i] <- skewness(e_i)
    results$kurtosis[i] <- kurtosis(e_i)

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

  # --- Display用のdata.frameを作成（型を明示的にcharacterに） ---
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
#' @param X 元データ (matrix or data.frame)
#' @param lingam_result direct_lingam() の返り値
#' @param ncol Number of columns.
#' @param nrow Number of rows.
#' @export
#' @examples
#' # サンプルデータの呼び出し
#' LiNGAM_sample_1000 <- generate_lingam_sample_6()
#'
#' # Direct LiNGAM の実行
#' result <- direct_lingam(LiNGAM_sample_1000$data)
#'
#' plot_residual_qq(LiNGAM_sample_1000$data, result)
plot_residual_qq <- function(X, lingam_result, ncol = 3, nrow = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please install it.", call. = FALSE)
  }

  X <- as.matrix(X)
  B <- lingam_result$adjacency_matrix
  E <- X - X %*% t(B)

  var_names <- colnames(X)
  if (is.null(var_names)) var_names <- paste0("x", seq_len(ncol(X)) - 1)

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


