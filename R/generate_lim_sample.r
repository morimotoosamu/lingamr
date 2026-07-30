# =============================================================================
# Sample data generation for LiM (LiNGAM for Mixed data)
# =============================================================================


#' Generate sample data for LiM (3 mixed variables)
#'
#' Generates a small dataset with a known causal chain of continuous and
#' discrete variables:
#' \code{x1 (continuous) -> x2 (discrete) -> x3 (continuous)}.
#' Continuous variables use Laplace-distributed noise (non-Gaussian, as
#' required by LiNGAM-family methods). By default the discrete variable is
#' drawn from a Bernoulli distribution whose logit is a linear function of
#' its parent; with `is_poisson = TRUE` it is instead drawn from a Poisson
#' distribution whose log-mean is a linear function of its parent.
#'
#' @param n number of samples (default: 1000)
#' @param seed random seed. If `NULL` (default), no seed is set and results
#'   are not reproducible across calls.
#' @param is_poisson if `TRUE`, generate the discrete variable `x2` as
#'   Poisson-distributed counts instead of binary 0/1 values (default: FALSE)
#' @return A list with four elements:
#' * `data`: data frame with columns `x1`, `x2`, `x3`.
#' * `adjacency_matrix`: 3x3 true adjacency matrix, following the lingamr
#'   convention (`m[to, from]`, i.e. `adjacency_matrix["x2", "x1"]` is the
#'   coefficient of the x1 -> x2 edge). The x1 -> x2 coefficient is on the
#'   linear-predictor scale (logit for binary, log for Poisson).
#' * `is_continuous`: logical vector `c(TRUE, FALSE, TRUE)` marking which
#'   columns of `data` are continuous.
#' * `is_poisson`: the input `is_poisson` flag, stored for reference.
#'
#' @examples
#' dat <- generate_lim_sample(n = 500, seed = 1)
#' result <- lingam_lim(dat$data, is_continuous = dat$is_continuous)
#' result$adjacency_matrix
#'
#' @export
generate_lim_sample <- function(n = 1000L, seed = NULL, is_poisson = FALSE) {
  if (!is.numeric(n) || length(n) != 1 || n < 2) {
    stop("n must be an integer >= 2.", call. = FALSE)
  }
  if (!is.null(seed) && !is.numeric(seed)) {
    stop("seed must be numeric or NULL.", call. = FALSE)
  }
  if (!is.logical(is_poisson) || length(is_poisson) != 1 || is.na(is_poisson)) {
    stop("is_poisson must be a single logical value (TRUE or FALSE).", call. = FALSE)
  }
  n <- as.integer(n)
  if (!is.null(seed)) set.seed(as.integer(seed))

  # smaller coefficients in the Poisson case: x2 ranges over counts rather
  # than 0/1, so a logit-scale 1.3 would let x2 (and hence x3) blow up
  if (is_poisson) {
    coef_12 <- 0.7
    coef_23 <- 0.7
  } else {
    coef_12 <- 1.3
    coef_23 <- 1.3
  }

  # Laplace(0, 1) noise: difference of two independent Exponential(1) variables
  laplace_noise <- function(n) stats::rexp(n) - stats::rexp(n)

  x1 <- laplace_noise(n)
  if (is_poisson) {
    eta2 <- 0.5 + coef_12 * x1
    # cap the Poisson mean at 50: the Laplace tail of x1 would otherwise
    # occasionally produce huge counts that dominate the scale of x3
    x2 <- stats::rpois(n, pmin(exp(eta2), 50))
  } else {
    eta2 <- coef_12 * x1
    x2 <- stats::rbinom(n, 1, stats::plogis(eta2))
  }
  x3 <- coef_23 * x2 + laplace_noise(n)

  X <- data.frame(x1 = x1, x2 = x2, x3 = x3)

  m_true <- matrix(0, 3, 3, dimnames = list(names(X), names(X)))
  m_true["x2", "x1"] <- coef_12
  m_true["x3", "x2"] <- coef_23

  list(
    data = X,
    adjacency_matrix = m_true,
    is_continuous = c(TRUE, FALSE, TRUE),
    is_poisson = is_poisson
  )
}
