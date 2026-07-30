# =============================================================================
# VARMA-LiNGAM - Diagnostics (stationarity/invertibility & residual
# non-Gaussianity)
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


#' @keywords internal
validate_varmalingam_result <- function(result) {
  if (!inherits(result, "VARMALiNGAMResult")) {
    stop("result must be a VARMALiNGAMResult (output of lingam_varma()).", call. = FALSE)
  }
}


#' Check the stationarity and invertibility of a fitted VARMA-LiNGAM model
#'
#' Inspects the eigenvalues of the companion matrices of the reduced-form AR
#' coefficients (Phi, stationarity) and MA coefficients (Theta, invertibility)
#' stored in the result. The process is stationary when every AR eigenvalue
#' lies strictly inside the unit circle, and invertible when every MA
#' eigenvalue does; a modulus on or outside the circle signals a (near-)unit
#' root or a non-invertible MA polynomial, under which the VARMA-LiNGAM
#' estimates (and the residual filtering) are unreliable. Invertibility is
#' worth checking here because the Hannan-Rissanen estimator does not enforce
#' it.
#'
#' @param result a `VARMALiNGAMResult` from [lingam_varma()]
#' @param tol threshold for the eigenvalue moduli (default 1)
#' @return a `varma_stationarity` object (list) with `ar_moduli` /
#'   `ma_moduli` (sorted descending; empty when p = 0 / q = 0),
#'   `max_ar_modulus`, `max_ma_modulus`, `is_stationary`, `is_invertible`,
#'   `order`, and `tol`.
#' @references
#' Stationarity diagnostics in the spirit of the VARLiNGAM R code of Moneta, A.,
#' Entner, D., Hoyer, P. O., & Coad, A. (2013), *Oxford Bulletin of Economics
#' and Statistics*, 75(5), 705-730.
#' <https://sites.google.com/site/dorisentner/publications/VARLiNGAM>
#' @export
#' @examples
#' s <- generate_varmalingam_sample(n = 1000, seed = 42)
#' m <- lingam_varma(s$data,
#'   order = c(1, 1), criterion = NULL,
#'   reg_method = "ols", prune = FALSE
#' )
#' check_varma_stationarity(m)
check_varma_stationarity <- function(result, tol = 1) {
  validate_varmalingam_result(result)

  ar_moduli <- companion_moduli(result$ar_coefs)
  ma_moduli <- companion_moduli(result$ma_coefs)
  max_ar <- if (length(ar_moduli)) max(ar_moduli) else 0
  max_ma <- if (length(ma_moduli)) max(ma_moduli) else 0

  obj <- list(
    ar_moduli      = sort(ar_moduli, decreasing = TRUE),
    max_ar_modulus = max_ar,
    is_stationary  = max_ar < tol,
    ma_moduli      = sort(ma_moduli, decreasing = TRUE),
    max_ma_modulus = max_ma,
    is_invertible  = max_ma < tol,
    order          = result$order,
    tol            = tol
  )
  class(obj) <- "varma_stationarity"
  obj
}


#' Print method for varma_stationarity
#'
#' @param x a `varma_stationarity` object
#' @param ... additional arguments (unused)
#' @return The input object `x`, invisibly.
#' @method print varma_stationarity
#' @export
#' @examples
#' s <- generate_varmalingam_sample(n = 1000, seed = 42)
#' m <- lingam_varma(s$data,
#'   order = c(1, 1), criterion = NULL,
#'   reg_method = "ols", prune = FALSE
#' )
#' print(check_varma_stationarity(m))
print.varma_stationarity <- function(x, ...) {
  cat("=== VARMA Stationarity / Invertibility Check ===\n")
  cat(sprintf("Order (p, q):         (%d, %d)\n", x$order[1], x$order[2]))
  cat(sprintf(
    "Max |AR eigenvalue|:  %.4f  (threshold %.2f)\n",
    x$max_ar_modulus, x$tol
  ))
  cat(sprintf("Stationary:           %s\n", if (x$is_stationary) "YES" else "NO"))
  cat(sprintf(
    "Max |MA eigenvalue|:  %.4f  (threshold %.2f)\n",
    x$max_ma_modulus, x$tol
  ))
  cat(sprintf("Invertible:           %s\n", if (x$is_invertible) "YES" else "NO"))
  if (!x$is_stationary) {
    cat("\nWARNING: the estimated VARMA is non-stationary (an AR root lies on\n")
    cat("  or outside the unit circle). VARMA-LiNGAM estimates may be unreliable.\n")
  }
  if (!x$is_invertible) {
    cat("\nWARNING: the estimated MA polynomial is not invertible; the residual\n")
    cat("  filtering underlying the fit may be unstable.\n")
  }
  invisible(x)
}


#' Residual matrix to diagnose for a VARMA-LiNGAM model
#'
#' Returns the series targeted by the residual diagnostics: either the LiNGAM
#' innovations `e_t = (I - B0) n_t` (the independent errors) or the
#' reduced-form VARMA residuals `n_t`. Shared by the normality tests and the
#' QQ plot.
#'
#' @param result a `VARMALiNGAMResult`
#' @param on "innovations" or "varma"
#' @return residual matrix (n_obs x n_features), column names preserved
#' @keywords internal
compute_varmalingam_residuals <- function(result, on = c("innovations", "varma")) {
  on <- match.arg(on)
  N <- result$residuals # VARMA residuals n_t (rows = time, cols = variables)
  if (on == "varma") {
    return(N)
  }
  # e_t = (I - B0) n_t; in row-per-observation form: E = N %*% t(I - B0)
  p <- dim(result$adjacency_matrices$psis)[2]
  B0 <- result$adjacency_matrices$psis[1, , ]
  E <- N %*% t(diag(p) - B0)
  colnames(E) <- colnames(N)
  E
}


#' Test the non-Gaussianity of VARMA-LiNGAM residuals
#'
#' LiNGAM assumes the error terms are non-Gaussian, so rejecting normality
#' (small p-value) supports the model assumption. By default the test is run on
#' the LiNGAM innovations `e_t = (I - B0) n_t` (the independent errors the model
#' assumes), where `n_t` are the stored VARMA residuals; set `on = "varma"` to
#' test the reduced-form VARMA residuals `n_t` directly instead.
#'
#' @param result a `VARMALiNGAMResult` from [lingam_varma()]
#' @param method normality test ("shapiro", "ks", "ad", "lillie", "jb");
#'   see [test_residual_normality()] for package requirements
#' @param alpha significance level (default 0.05)
#' @param on which series to test: "innovations" (default, `e_t = (I - B0) n_t`)
#'   or "varma" (the reduced-form VARMA residuals `n_t`)
#' @return a `lingam_normality_test` data frame (one row per variable), printed
#'   via [print.lingam_normality_test()].
#' @references
#' Residual non-Gaussianity diagnostics inspired by the VARLiNGAM R code
#' (Gauss_Tests) of Moneta, A., Entner, D., Hoyer, P. O., & Coad, A. (2013),
#' *Oxford Bulletin of Economics and Statistics*, 75(5), 705-730.
#' <https://sites.google.com/site/dorisentner/publications/VARLiNGAM>
#' @export
#' @examples
#' s <- generate_varmalingam_sample(n = 1000, seed = 42)
#' m <- lingam_varma(s$data,
#'   order = c(1, 1), criterion = NULL,
#'   reg_method = "ols", prune = FALSE
#' )
#' test_varmalingam_residual_normality(m)
test_varmalingam_residual_normality <- function(result,
                                                method = "shapiro",
                                                alpha = 0.05,
                                                on = c("innovations", "varma")) {
  validate_varmalingam_result(result)
  on <- match.arg(on)
  E <- compute_varmalingam_residuals(result, on)

  # Reuse the Direct LiNGAM normality routine via a zero adjacency matrix (see
  # zero_lingam_result): lingam_residuals() returns E unchanged, so the same
  # tests, summary statistics, result class, and print method all apply.
  test_residual_normality(E, zero_lingam_result(ncol(E)), method = method, alpha = alpha)
}


#' Run several normality tests on VARMA-LiNGAM residuals at once
#'
#' Convenience wrapper (analogous to the Moneta `Gauss_Tests`) that applies
#' multiple normality tests to the residuals and returns a single table with one
#' p-value column per method plus per-variable skewness and excess kurtosis.
#' Methods whose optional package is unavailable are skipped with a warning.
#'
#' @param result a `VARMALiNGAMResult` from [lingam_varma()]
#' @param methods character vector of tests to run; any of "shapiro", "ks",
#'   "ad", "lillie", "jb" (default runs shapiro/ad/lillie/jb)
#' @param alpha significance level (default 0.05)
#' @param on which series to test: "innovations" (default) or "varma"
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
#' s <- generate_varmalingam_sample(n = 1000, seed = 42)
#' m <- lingam_varma(s$data,
#'   order = c(1, 1), criterion = NULL,
#'   reg_method = "ols", prune = FALSE
#' )
#' test_varmalingam_residual_normality_all(m, methods = c("shapiro", "jb"))
test_varmalingam_residual_normality_all <- function(result,
                                                    methods = c("shapiro", "ad", "lillie", "jb"),
                                                    alpha = 0.05,
                                                    on = c("innovations", "varma")) {
  validate_varmalingam_result(result)
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
      paste(methods[!available], collapse = ", "),
      call. = FALSE
    )
    methods <- methods[available]
  }
  if (length(methods) == 0L) {
    stop("No usable normality-test method (required packages not installed).", call. = FALSE)
  }

  # Skewness/kurtosis are method-independent, so take them from the first run.
  first <- test_varmalingam_residual_normality(result, method = methods[1], alpha = alpha, on = on)
  out <- data.frame(
    variable = first$variable,
    skewness = first$skewness,
    kurtosis = first$kurtosis
  )
  for (m in methods) {
    r <- test_varmalingam_residual_normality(result, method = m, alpha = alpha, on = on)
    out[[paste0("p_", m)]] <- r$p_value
  }
  # A variable is flagged non-Gaussian only when *every* run test rejects.
  p_cols <- grep("^p_", names(out), value = TRUE)
  out$all_non_gauss <- apply(out[, p_cols, drop = FALSE], 1L, function(p) all(p <= alpha, na.rm = TRUE))

  attr(out, "alpha") <- alpha
  attr(out, "on") <- on
  out
}


#' Q-Q plots of VARMA-LiNGAM residuals
#'
#' Draws per-variable normal Q-Q plots of the residuals (analogous to the Moneta
#' `Gauss_Stats` visual check). Deviations from the reference line indicate
#' non-Gaussianity, which supports the LiNGAM assumption. Requires ggplot2.
#'
#' @param result a `VARMALiNGAMResult` from [lingam_varma()]
#' @param on which series to plot: "innovations" (default) or "varma"
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
#' s <- generate_varmalingam_sample(n = 1000, seed = 42)
#' m <- lingam_varma(s$data,
#'   order = c(1, 1), criterion = NULL,
#'   reg_method = "ols", prune = FALSE
#' )
#' \donttest{
#' plot_varmalingam_residual_qq(m)
#' }
plot_varmalingam_residual_qq <- function(result, on = c("innovations", "varma"),
                                         ncol = 3, nrow = NULL) {
  validate_varmalingam_result(result)
  on <- match.arg(on)
  E <- compute_varmalingam_residuals(result, on)
  # Reuse the Direct LiNGAM QQ plot via a zero adjacency matrix (see
  # zero_lingam_result), which makes it plot E directly.
  plot_residual_qq(E, zero_lingam_result(ncol(E)), ncol = ncol, nrow = nrow)
}
