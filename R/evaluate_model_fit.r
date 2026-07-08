# =============================================================================
# Evaluate model fit - SEM-based fit measures for an estimated adjacency matrix
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


#' Check whether lavaan is available
#' @keywords internal
check_lavaan_available <- function() {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop(
      "Package 'lavaan' is required for evaluate_model_fit(). ",
      "Install it with install.packages(\"lavaan\").",
      call. = FALSE
    )
  }
}


#' Extract an adjacency matrix from a result object or matrix
#' @keywords internal
extract_adjacency_matrix <- function(adjacency_matrix) {
  if (is.list(adjacency_matrix) && !is.data.frame(adjacency_matrix) &&
      !is.null(adjacency_matrix$adjacency_matrix)) {
    return(adjacency_matrix$adjacency_matrix)
  }
  adjacency_matrix
}


#' Build a lavaan model string from an adjacency matrix
#'
#' Converts a lingamr-convention adjacency matrix (`B[i, j]` = causal
#' coefficient from j to i) into a lavaan model syntax string. Non-zero
#' elements become regression paths (`xi ~ xj`); `NA` elements (used by
#' e.g. [lingam_parce()] to mark a latent confounder between two variables)
#' become a residual covariance (`xi ~~ xj`) between the two machine-named
#' variables, which is the standard lavaan idiom equivalent to a two-indicator
#' latent common cause with one loading fixed (as used by the Python
#' `semopy`-based original).
#'
#' @param B adjacency matrix (machine names `x0, x1, ...` expected as
#'   row/column indices; the caller supplies names via `var_names`)
#' @param var_names machine variable names, length `ncol(B)`
#' @return character scalar, lavaan model syntax (possibly empty string)
#' @keywords internal
build_lavaan_model <- function(B, var_names) {
  p <- ncol(B)
  reg_lines <- character(0)
  cov_pairs <- character(0)
  seen_pairs <- character(0)

  for (i in seq_len(p)) {
    row <- B[i, ]
    if (!anyNA(row) && all(row == 0)) next # exogenous variable: no equation

    rhs <- character(0)
    for (j in seq_len(p)) {
      if (i == j) next
      elem <- unname(row[j])
      if (is.na(elem)) {
        pair_idx <- sort(c(i, j))
        pair_key <- paste(pair_idx, collapse = "_")
        if (!(pair_key %in% seen_pairs)) {
          seen_pairs <- c(seen_pairs, pair_key)
          cov_pairs <- c(cov_pairs, sprintf(
            "%s ~~ %s", var_names[pair_idx[1]], var_names[pair_idx[2]]
          ))
        }
      } else if (!isTRUE(all.equal(elem, 0))) {
        rhs <- c(rhs, var_names[j])
      }
    }
    if (length(rhs) > 0) {
      reg_lines <- c(reg_lines, sprintf("%s ~ %s", var_names[i], paste(rhs, collapse = " + ")))
    }
  }

  paste(c(reg_lines, cov_pairs), collapse = "\n")
}


#' Map lavaan fitMeasures() to the semopy-style column set
#' @keywords internal
lavaan_fit_measures_to_df <- function(fit) {
  fm <- lavaan::fitMeasures(fit)
  get_fm <- function(name) {
    if (name %in% names(fm)) unname(fm[[name]]) else NA_real_
  }
  data.frame(
    DoF             = get_fm("df"),
    `DoF Baseline`  = get_fm("baseline.df"),
    chi2            = get_fm("chisq"),
    `chi2 p-value`  = get_fm("pvalue"),
    `chi2 Baseline` = get_fm("baseline.chisq"),
    CFI             = get_fm("cfi"),
    GFI             = get_fm("gfi"),
    AGFI            = get_fm("agfi"),
    NFI             = get_fm("nfi"),
    TLI             = get_fm("tli"),
    RMSEA           = get_fm("rmsea"),
    AIC             = get_fm("aic"),
    BIC             = get_fm("bic"),
    LogLik          = get_fm("logl"),
    check.names = FALSE
  )
}


#' Evaluate model fit of an estimated causal graph
#'
#' Fits the causal graph implied by `adjacency_matrix` as a structural
#' equation model (SEM) via `lavaan::sem()` and returns standard SEM fit
#' measures (CFI, RMSEA, AIC/BIC, etc.). This is an R port of the Python
#' `lingam.utils.evaluate_model_fit()`, which delegates to the Python
#' package `semopy`; this R version delegates to `lavaan` instead.
#'
#' @details
#' * **Optional dependency**: this function requires the \pkg{lavaan}
#'   package (listed in `Suggests`, not `Imports`). Install it with
#'   `install.packages("lavaan")`.
#' * **Latent confounders**: an `NA` element `B[i, j]` (as produced by
#'   e.g. [lingam_parce()] for a suspected latent confounder between
#'   variables i and j) is represented as a residual covariance
#'   `xi ~~ xj` in the lavaan model. This is algebraically equivalent to
#'   the two-indicator latent common cause (one loading fixed to 1) used
#'   by the Python `semopy` implementation, but is expressed with lavaan's
#'   standard residual-covariance idiom rather than an explicit latent
#'   variable.
#' * **Numerical values will not match `semopy` exactly**: `lavaan` and
#'   `semopy` use different default estimators/options. The fit measures
#'   returned are the same statistics, but exact numbers can differ.
#' * **Convention**: `adjacency_matrix` follows the lingamr convention
#'   `B[i, j]` = causal coefficient from variable j to variable i (j -> i).
#'
#' @param adjacency_matrix p x p numeric adjacency matrix (NA allowed for
#'   latent confounder pairs), or a lingamr result object (e.g.
#'   `LingamResult`, `ParceLingamResult`, `LiMResult`) with an
#'   `adjacency_matrix` element, from which the matrix is extracted
#'   automatically
#' @param X numeric matrix or data frame (n_samples x p) with no missing
#'   values
#' @param is_ordinal logical or 0/1 vector of length p. `TRUE` marks a
#'   variable as ordinal (categorical), fit with `lavaan`'s WLSMV-based
#'   estimator. `NULL` (default) treats all variables as continuous.
#' @return A one-row data.frame of fit measures: DoF, DoF Baseline, chi2,
#'   chi2 p-value, chi2 Baseline, CFI, GFI, AGFI, NFI, TLI, RMSEA, AIC,
#'   BIC, LogLik. When `is_ordinal` is used, AIC/BIC/LogLik and some other
#'   measures are not defined by the WLSMV estimator and are returned as
#'   `NA`.
#' @references
#' Rosseel, Y. (2012). lavaan: An R Package for Structural Equation
#' Modeling. Journal of Statistical Software, 48(2), 1-36.
#' \doi{10.18637/jss.v048.i02}
#' @export
#' @examples
#' if (requireNamespace("lavaan", quietly = TRUE)) {
#'   dat <- generate_lingam_sample_6()
#'   result <- lingam_direct(dat$data, reg_method = "ols")
#'   evaluate_model_fit(result, dat$data)
#' }
evaluate_model_fit <- function(adjacency_matrix, X, is_ordinal = NULL) {
  check_lavaan_available()

  B <- extract_adjacency_matrix(adjacency_matrix)
  B <- as.matrix(B)
  if (nrow(B) != ncol(B)) stop("adjacency_matrix must be a square matrix.", call. = FALSE)

  p <- ncol(B)
  X <- as.matrix(as.data.frame(X))
  if (!is.numeric(X)) stop("X must be a numeric matrix or data frame.", call. = FALSE)
  if (ncol(X) != p) stop("ncol(X) must equal the number of columns in adjacency_matrix.", call. = FALSE)
  if (anyNA(X)) stop("X must not contain missing values (NA).", call. = FALSE)

  if (!is.null(is_ordinal)) {
    if (length(is_ordinal) != p) {
      stop("is_ordinal must have length equal to the number of columns in adjacency_matrix.", call. = FALSE)
    }
    is_ordinal <- as.logical(is_ordinal)
  } else {
    is_ordinal <- rep(FALSE, p)
  }

  var_names <- paste0("x", seq_len(p) - 1L)
  colnames(X) <- var_names
  ordered_vars <- var_names[is_ordinal]

  model_str <- build_lavaan_model(B, var_names)
  if (nchar(trimws(model_str)) == 0) {
    stop("evaluate_model_fit() requires at least one edge: the adjacency matrix has no edges.", call. = FALSE)
  }

  fit <- tryCatch(
    lavaan::sem(
      model_str,
      data = as.data.frame(X),
      ordered = if (length(ordered_vars) > 0) ordered_vars else NULL
    ),
    error = function(e) {
      stop(
        "evaluate_model_fit(): lavaan::sem() failed to fit the model: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  if (!isTRUE(lavaan::lavInspect(fit, "converged"))) {
    stop(
      "evaluate_model_fit(): the SEM did not converge; fit measures would ",
      "not be meaningful. This can happen with near-collinear predictors, ",
      "too few observations relative to the number of edges, or an ",
      "over-identified model.",
      call. = FALSE
    )
  }

  lavaan_fit_measures_to_df(fit)
}
