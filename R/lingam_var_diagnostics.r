# =============================================================================
# VAR-LiNGAM - Diagnostics (stationarity & residual non-Gaussianity)
#
# License: MIT + file LICENSE
#
# The residual normality tests and Q-Q plots are inspired by the diagnostics
# (Gauss_Tests / Gauss_Stats) in the VARLiNGAM R code of Moneta, Entner, Hoyer
# & Coad: https://sites.google.com/site/dorisentner/publications/VARLiNGAM
#
# Portions of this work:
#   Copyright (c) 2026 O.Morimoto
# =============================================================================


#' Check the stationarity of a fitted VAR-LiNGAM model
#'
#' Recovers the reduced-form VAR coefficients `M_k = (I - B0)^{-1} B_k` from the
#' structural matrices and inspects the eigenvalues of the VAR companion matrix.
#' The process is stationary when every eigenvalue lies strictly inside the unit
#' circle (all moduli < 1); a modulus on or outside it signals a (near-)unit root,
#' under which the VAR-LiNGAM estimates are unreliable.
#'
#' @param result a `VARLiNGAMResult` from [lingam_var()]
#' @param tol stationarity threshold for the eigenvalue moduli (default 1)
#' @return a `var_stationarity` object (list) with `moduli` (sorted descending),
#'   `max_modulus`, `is_stationary` (logical), `lags`, and `tol`.
#' @references
#' Stationarity diagnostics in the spirit of the VARLiNGAM R code of Moneta, A.,
#' Entner, D., Hoyer, P. O., & Coad, A. (2013), *Oxford Bulletin of Economics
#' and Statistics*, 75(5), 705-730.
#' <https://sites.google.com/site/dorisentner/publications/VARLiNGAM>
#' @export
#' @examples
#' s <- generate_varlingam_sample(n = 1000, seed = 42)
#' m <- lingam_var(s$data, lags = 1, reg_method = "ols", prune = FALSE)
#' check_var_stationarity(m)
check_var_stationarity <- function(result, tol = 1) {
  if (!inherits(result, "VARLiNGAMResult")) {
    stop("result must be a VARLiNGAMResult (output of lingam_var()).", call. = FALSE)
  }
  am <- result$adjacency_matrices
  p <- dim(am)[2]
  lags <- result$lags

  B0 <- am[1, , ]
  IB0_inv <- solve(diag(p) - B0)
  # reduced-form coefficients M_k = (I - B0)^{-1} B_k for k = 1..lags
  M <- lapply(seq_len(lags), function(k) IB0_inv %*% am[k + 1L, , ])

  # VAR(p) companion matrix (p*lags x p*lags):
  #   [ M_1 M_2 ... M_p ]
  #   [  I   0  ...  0  ]
  #   [  0   I  ...  0  ]   (sub-diagonal identity blocks)
  d <- p * lags
  companion <- matrix(0, nrow = d, ncol = d)
  companion[seq_len(p), ] <- do.call(cbind, M)
  if (lags > 1L) {
    companion[(p + 1L):d, seq_len(p * (lags - 1L))] <- diag(p * (lags - 1L))
  }

  moduli <- Mod(eigen(companion, only.values = TRUE)$values)
  max_mod <- max(moduli)

  obj <- list(
    moduli        = sort(moduli, decreasing = TRUE),
    max_modulus   = max_mod,
    is_stationary = max_mod < tol,
    lags          = lags,
    tol           = tol
  )
  class(obj) <- "var_stationarity"
  obj
}


#' Print method for var_stationarity
#'
#' @param x a `var_stationarity` object
#' @param ... additional arguments (unused)
#' @method print var_stationarity
#' @export
print.var_stationarity <- function(x, ...) {
  cat("=== VAR Stationarity Check ===\n")
  cat(sprintf("Lag order:         %d\n", x$lags))
  cat(sprintf("Max |eigenvalue|:  %.4f  (threshold %.2f)\n", x$max_modulus, x$tol))
  cat(sprintf("Stationary:        %s\n", if (x$is_stationary) "YES" else "NO"))
  if (!x$is_stationary) {
    cat("\nWARNING: the estimated VAR is non-stationary (a root lies on or\n")
    cat("  outside the unit circle). VAR-LiNGAM estimates may be unreliable.\n")
  }
  invisible(x)
}


#' Residual matrix to diagnose for a VAR-LiNGAM model
#'
#' Returns the series targeted by the residual diagnostics: either the LiNGAM
#' innovations `e_t = (I - B0) n_t` (the independent errors) or the reduced-form
#' VAR residuals `n_t`. Shared by the normality tests and the QQ plot.
#'
#' @param result a `VARLiNGAMResult`
#' @param on "innovations" or "var"
#' @return residual matrix (n_obs x n_features), column names preserved
#' @keywords internal
compute_varlingam_residuals <- function(result, on = c("innovations", "var")) {
  on <- match.arg(on)
  N <- result$residuals   # VAR residuals n_t (rows = time, cols = variables)
  if (on == "var") {
    return(N)
  }
  # e_t = (I - B0) n_t; in row-per-observation form: E = N %*% t(I - B0)
  p <- dim(result$adjacency_matrices)[2]
  B0 <- result$adjacency_matrices[1, , ]
  E <- N %*% t(diag(p) - B0)
  colnames(E) <- colnames(N)
  E
}


#' Build a zero-adjacency LingamResult
#'
#' A stand-in `LingamResult` whose adjacency matrix is all zeros, so that
#' [lingam_residuals()] returns its input unchanged. This lets the VAR
#' diagnostics reuse the Direct LiNGAM residual routines on an already-computed
#' residual matrix.
#'
#' @param p number of features
#' @return a `LingamResult` with a p x p zero adjacency matrix
#' @keywords internal
zero_lingam_result <- function(p) {
  structure(
    list(adjacency_matrix = matrix(0, nrow = p, ncol = p), causal_order = seq_len(p)),
    class = "LingamResult"
  )
}


#' Test the non-Gaussianity of VAR-LiNGAM residuals
#'
#' LiNGAM assumes the error terms are non-Gaussian, so rejecting normality
#' (small p-value) supports the model assumption. By default the test is run on
#' the LiNGAM innovations `e_t = (I - B0) n_t` (the independent errors the model
#' assumes), where `n_t` are the stored VAR residuals; set `on = "var"` to test
#' the reduced-form VAR residuals `n_t` directly instead.
#'
#' @param result a `VARLiNGAMResult` from [lingam_var()]
#' @param method normality test ("shapiro", "ks", "ad", "lillie", "jb");
#'   see [test_residual_normality()] for package requirements
#' @param alpha significance level (default 0.05)
#' @param on which series to test: "innovations" (default, `e_t = (I - B0) n_t`)
#'   or "var" (the reduced-form VAR residuals `n_t`)
#' @return a `lingam_normality_test` data frame (one row per variable), printed
#'   via [print.lingam_normality_test()].
#' @references
#' Residual non-Gaussianity diagnostics inspired by the VARLiNGAM R code
#' (Gauss_Tests) of Moneta, A., Entner, D., Hoyer, P. O., & Coad, A. (2013),
#' *Oxford Bulletin of Economics and Statistics*, 75(5), 705-730.
#' <https://sites.google.com/site/dorisentner/publications/VARLiNGAM>
#' @export
#' @examples
#' s <- generate_varlingam_sample(n = 1000, seed = 42)
#' m <- lingam_var(s$data, lags = 1, reg_method = "ols", prune = FALSE)
#' test_varlingam_residual_normality(m)
test_varlingam_residual_normality <- function(result,
                                              method = "shapiro",
                                              alpha = 0.05,
                                              on = c("innovations", "var")) {
  if (!inherits(result, "VARLiNGAMResult")) {
    stop("result must be a VARLiNGAMResult (output of lingam_var()).", call. = FALSE)
  }
  on <- match.arg(on)
  E <- compute_varlingam_residuals(result, on)

  # Reuse the Direct LiNGAM normality routine via a zero adjacency matrix (see
  # zero_lingam_result): lingam_residuals() returns E unchanged, so the same
  # tests, summary statistics, result class, and print method all apply.
  test_residual_normality(E, zero_lingam_result(ncol(E)), method = method, alpha = alpha)
}


#' Run several normality tests on VAR-LiNGAM residuals at once
#'
#' Convenience wrapper (analogous to the Moneta `Gauss_Tests`) that applies
#' multiple normality tests to the residuals and returns a single table with one
#' p-value column per method plus per-variable skewness and excess kurtosis.
#' Methods whose optional package is unavailable are skipped with a warning.
#'
#' @param result a `VARLiNGAMResult` from [lingam_var()]
#' @param methods character vector of tests to run; any of "shapiro", "ks",
#'   "ad", "lillie", "jb" (default runs shapiro/ad/lillie/jb)
#' @param alpha significance level (default 0.05)
#' @param on which series to test: "innovations" (default) or "var"
#' @return a data frame with columns `variable`, `skewness`, `kurtosis`, one
#'   `p_<method>` column per method, and `all_non_gauss` (TRUE when every run
#'   test rejects normality for that variable).
#' @references
#' Analogous to the multi-test residual check (Gauss_Tests) in the VARLiNGAM R
#' code of Moneta, A., Entner, D., Hoyer, P. O., & Coad, A. (2013), *Oxford
#' Bulletin of Economics and Statistics*, 75(5), 705-730.
#' <https://sites.google.com/site/dorisentner/publications/VARLiNGAM>
#' @export
#' @examples
#' s <- generate_varlingam_sample(n = 1000, seed = 42)
#' m <- lingam_var(s$data, lags = 1, reg_method = "ols", prune = FALSE)
#' test_varlingam_residual_normality_all(m, methods = c("shapiro", "jb"))
test_varlingam_residual_normality_all <- function(result,
                                                  methods = c("shapiro", "ad", "lillie", "jb"),
                                                  alpha = 0.05,
                                                  on = c("innovations", "var")) {
  if (!inherits(result, "VARLiNGAMResult")) {
    stop("result must be a VARLiNGAMResult (output of lingam_var()).", call. = FALSE)
  }
  on <- match.arg(on)
  valid <- c("shapiro", "ks", "ad", "lillie", "jb")
  methods <- unique(match.arg(methods, valid, several.ok = TRUE))

  # Skip methods whose optional package is missing, but run the remaining ones.
  # list (not a named vector) so that pkg_for[[m]] returns NULL for methods
  # that need no extra package (shapiro/ks) instead of erroring.
  pkg_for <- list(ad = "nortest", lillie = "nortest", jb = "tseries")
  available <- vapply(methods, function(m) {
    pkg <- pkg_for[[m]]
    is.null(pkg) || requireNamespace(pkg, quietly = TRUE)
  }, logical(1))
  if (any(!available)) {
    warning("Skipping methods with missing packages: ",
            paste(methods[!available], collapse = ", "), call. = FALSE)
    methods <- methods[available]
  }
  if (length(methods) == 0L) {
    stop("No usable normality-test method (required packages not installed).", call. = FALSE)
  }

  # Skewness/kurtosis are method-independent, so take them from the first run.
  first <- test_varlingam_residual_normality(result, method = methods[1], alpha = alpha, on = on)
  out <- data.frame(
    variable = first$variable,
    skewness = first$skewness,
    kurtosis = first$kurtosis,
    stringsAsFactors = FALSE
  )
  for (m in methods) {
    r <- test_varlingam_residual_normality(result, method = m, alpha = alpha, on = on)
    out[[paste0("p_", m)]] <- r$p_value
  }
  # A variable is flagged non-Gaussian only when *every* run test rejects.
  p_cols <- grep("^p_", names(out), value = TRUE)
  out$all_non_gauss <- apply(out[, p_cols, drop = FALSE], 1L, function(p) all(p <= alpha, na.rm = TRUE))

  attr(out, "alpha") <- alpha
  attr(out, "on") <- on
  out
}


#' Q-Q plots of VAR-LiNGAM residuals
#'
#' Draws per-variable normal Q-Q plots of the residuals (analogous to the Moneta
#' `Gauss_Stats` visual check). Deviations from the reference line indicate
#' non-Gaussianity, which supports the LiNGAM assumption. Requires ggplot2.
#'
#' @param result a `VARLiNGAMResult` from [lingam_var()]
#' @param on which series to plot: "innovations" (default) or "var"
#' @param ncol number of facet columns
#' @param nrow number of facet rows (NULL = automatic)
#' @return a ggplot object
#' @references
#' Analogous to the residual visual check (Gauss_Stats) in the VARLiNGAM R code
#' of Moneta, A., Entner, D., Hoyer, P. O., & Coad, A. (2013), *Oxford Bulletin
#' of Economics and Statistics*, 75(5), 705-730.
#' <https://sites.google.com/site/dorisentner/publications/VARLiNGAM>
#' @export
#' @examples
#' s <- generate_varlingam_sample(n = 1000, seed = 42)
#' m <- lingam_var(s$data, lags = 1, reg_method = "ols", prune = FALSE)
#' \donttest{
#' plot_varlingam_residual_qq(m)
#' }
plot_varlingam_residual_qq <- function(result, on = c("innovations", "var"),
                                       ncol = 3, nrow = NULL) {
  if (!inherits(result, "VARLiNGAMResult")) {
    stop("result must be a VARLiNGAMResult (output of lingam_var()).", call. = FALSE)
  }
  on <- match.arg(on)
  E <- compute_varlingam_residuals(result, on)
  # Reuse the Direct LiNGAM QQ plot via a zero adjacency matrix (see
  # zero_lingam_result), which makes it plot E directly.
  plot_residual_qq(E, zero_lingam_result(ncol(E)), ncol = ncol, nrow = nrow)
}
