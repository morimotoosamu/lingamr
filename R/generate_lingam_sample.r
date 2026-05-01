#' Create noise generation function
#'
#' Internal helper to create a noise function for the specified distribution.
#'
#' @param noise_dist distribution name
#'   "uniform"     : Uniform(0, 1) - non-Gaussian (LiNGAM works well)
#'   "gaussian"    : Normal(0, 1)  - LiNGAM may fail
#'   "lognormal"   : Log-normal(0, 1) - skewed, non-Gaussian
#'   "exponential" : Exponential(1) - skewed, non-Gaussian
#'   "t3"          : t-distribution (df=3) - heavy tails
#' @importFrom stats rexp rlnorm rt runif rnorm
#' @return function(n) that generates n random numbers
#' @keywords internal
make_noise_fn <- function(noise_dist) {
  switch(noise_dist,
         "uniform" = function(n) {
           runif(n)
         },
         "gaussian" = function(n) {
           rnorm(n, mean = 0, sd = 1)
         },
         "lognormal" = function(n) {
           x <- rlnorm(n, meanlog = 0, sdlog = 1)
           # 平均0になるよう中心化
           x - mean(x)
         },
         "exponential" = function(n) {
           rexp(n, rate = 1) - 1  # 平均0になるよう中心化
         },
         "t3" = function(n) {
           rt(n, df = 3)
         },
         stop(sprintf(
           "noise_dist must be one of: uniform, gaussian, lognormal, exponential, t3"
         ))
  )
}


#' Generate sample data for Direct LiNGAM (6 variables)
#'
#' Generates a sample dataset with a known causal structure.
#' The true causal structure is:
#'   x3 -> x0 (coef = 3.0)
#'   x3 -> x2 (coef = 6.0)
#'   x0 -> x1 (coef = 3.0)
#'   x2 -> x1 (coef = 2.0)
#'   x0 -> x5 (coef = 4.0)
#'   x0 -> x4 (coef = 8.0)
#'   x2 -> x4 (coef = -1.0)
#'
#' @param n number of samples (default: 1000)
#' @param seed random seed (default: 42)
#' @param noise_dist error term distribution
#'   "uniform"     : Uniform(0, 1) - default, non-Gaussian (LiNGAM works well)
#'   "gaussian"    : Normal(0, 1)  - LiNGAM may fail
#'   "lognormal"   : Log-normal(0, 1) - skewed, non-Gaussian
#'   "exponential" : Exponential(1) - skewed, non-Gaussian
#'   "t3"          : t-distribution (df=3) - heavy tails
#' @return list(data, true_adjacency)
#'
#' @examples
#' # Non-Gaussian (LiNGAM works well)
#' X_nongauss <- generate_lingam_sample_6(noise_dist = "uniform")
#' result <- lingam_direct(X_nongauss$data)
#' result$causal_order
#'
#' # Gaussian (LiNGAM may fail)
#' X_gauss <- generate_lingam_sample_6(noise_dist = "gaussian")
#' result <- lingam_direct(X_gauss$data)
#' result$causal_order
#'
#' @export
generate_lingam_sample_6 <- function(n = 1000L,
                                     seed = 42L,
                                     noise_dist = "uniform") {

  # --- Input validation ---
  if (!is.numeric(n) || n < 2) stop("n must be an integer >= 2.")
  if (!is.numeric(seed))       stop("seed must be numeric.")

  valid_dists <- c("uniform", "gaussian", "lognormal", "exponential", "t3")
  if (!(noise_dist %in% valid_dists)) {
    stop(sprintf("noise_dist must be one of: %s",
                 paste(valid_dists, collapse = ", ")))
  }

  n    <- as.integer(n)
  seed <- as.integer(seed)

  # --- ノイズ生成関数 ---
  noise_fn <- make_noise_fn(noise_dist)

  # --- Error terms（各変数に独立したシード）---
  set.seed(seed)      ; e0 <- noise_fn(n)
  set.seed(seed + 1L) ; e1 <- noise_fn(n)
  set.seed(seed + 2L) ; e2 <- noise_fn(n)
  set.seed(seed + 3L) ; e3 <- noise_fn(n)
  set.seed(seed + 4L) ; e4 <- noise_fn(n)
  set.seed(seed + 5L) ; e5 <- noise_fn(n)

  # --- Data generation（因果順序に従う）---
  x3 <- e3
  x0 <- 3.0 * x3 + e0
  x2 <- 6.0 * x3 + e2
  x1 <- 3.0 * x0 + 2.0 * x2 + e1
  x5 <- 4.0 * x0 + e5
  x4 <- 8.0 * x0 - 1.0 * x2 + e4

  X <- data.frame(x0 = x0, x1 = x1, x2 = x2,
             x3 = x3, x4 = x4, x5 = x5)

  # --- 真の隣接行列 ---
  var_names <- names(X)
  p <- length(var_names)
  m_true <- matrix(0, p, p, dimnames = list(var_names, var_names))

  # [row = to, col = from] = coefficient
  m_true["x0", "x3"] <-  3.0
  m_true["x2", "x3"] <-  6.0
  m_true["x1", "x0"] <-  3.0
  m_true["x1", "x2"] <-  2.0
  m_true["x5", "x0"] <-  4.0
  m_true["x4", "x0"] <-  8.0
  m_true["x4", "x2"] <- -1.0

  return(list(
    data           = X,
    true_adjacency = m_true
  ))
}


#' Generate 10-variable sample data for Direct LiNGAM
#'
#' Generates a sample dataset with a known causal structure.
#' The true causal structure is:
#'   x3 -> x0 (coef =  3.0)
#'   x3 -> x2 (coef =  6.0)
#'   x3 -> x9 (coef =  7.0)
#'   x0 -> x1 (coef =  3.0)
#'   x0 -> x5 (coef =  4.0)
#'   x0 -> x4 (coef =  8.0)
#'   x0 -> x7 (coef =  3.0)
#'   x2 -> x1 (coef =  2.0)
#'   x2 -> x4 (coef = -1.0)
#'   x2 -> x8 (coef =  0.5)
#'   x1 -> x6 (coef =  2.0)
#'   x5 -> x8 (coef =  2.0)
#'   x4 -> x7 (coef =  1.5)
#'   x6 -> x9 (coef =  1.0)
#'
#' @param n number of samples (default: 1000)
#' @param seed random seed (default: 42)
#' @param noise_dist error term distribution
#'   "uniform"     : Uniform(0, 1) - default, non-Gaussian (LiNGAM works well)
#'   "gaussian"    : Normal(0, 1)  - LiNGAM may fail
#'   "lognormal"   : Log-normal(0, 1) - skewed, non-Gaussian
#'   "exponential" : Exponential(1) - skewed, non-Gaussian
#'   "t3"          : t-distribution (df=3) - heavy tails
#' @return list(data, true_adjacency)
#'
#' @examples
#' # Non-Gaussian (LiNGAM works well)
#' X_nongauss <- generate_lingam_sample_10(noise_dist = "uniform")
#' result <- lingam_direct(X_nongauss$data)
#' result$causal_order
#'
#' # Gaussian (LiNGAM may fail)
#' X_gauss <- generate_lingam_sample_10(noise_dist = "gaussian")
#' result <- lingam_direct(X_gauss$data)
#' result$causal_order
#'
#' @export
generate_lingam_sample_10 <- function(n = 1000L,
                                      seed = 42L,
                                      noise_dist = "uniform") {

  # --- Input validation ---
  if (!is.numeric(n) || n < 2) stop("n must be an integer >= 2.")
  if (!is.numeric(seed))       stop("seed must be numeric.")

  valid_dists <- c("uniform", "gaussian", "lognormal", "exponential", "t3")
  if (!(noise_dist %in% valid_dists)) {
    stop(sprintf("noise_dist must be one of: %s",
                 paste(valid_dists, collapse = ", ")))
  }

  n    <- as.integer(n)
  seed <- as.integer(seed)

  # --- ノイズ生成関数 ---
  noise_fn <- make_noise_fn(noise_dist)

  # --- Error terms（各変数に独立したシード）---
  set.seed(seed)      ; e0 <- noise_fn(n)
  set.seed(seed + 1L) ; e1 <- noise_fn(n)
  set.seed(seed + 2L) ; e2 <- noise_fn(n)
  set.seed(seed + 3L) ; e3 <- noise_fn(n)
  set.seed(seed + 4L) ; e4 <- noise_fn(n)
  set.seed(seed + 5L) ; e5 <- noise_fn(n)
  set.seed(seed + 6L) ; e6 <- noise_fn(n)
  set.seed(seed + 7L) ; e7 <- noise_fn(n)
  set.seed(seed + 8L) ; e8 <- noise_fn(n)
  set.seed(seed + 9L) ; e9 <- noise_fn(n)

  # --- Data generation（因果順序に従う）---
  x3 <- e3
  x0 <- 3.0 * x3 + e0
  x2 <- 6.0 * x3 + e2
  x1 <- 3.0 * x0 + 2.0 * x2 + e1
  x5 <- 4.0 * x0 + e5
  x4 <- 8.0 * x0 - 1.0 * x2 + e4
  x6 <- 2.0 * x1 + e6
  x7 <- 1.5 * x4 + 3.0 * x0 + e7
  x8 <- 0.5 * x2 + 2.0 * x5 + e8
  x9 <- 7.0 * x3 + 1.0 * x6 + e9

  X <- data.frame(x0 = x0, x1 = x1, x2 = x2, x3 = x3, x4 = x4,
             x5 = x5, x6 = x6, x7 = x7, x8 = x8, x9 = x9)

  # --- 真の隣接行列 ---
  var_names <- names(X)
  p <- length(var_names)
  m_true <- matrix(0, p, p, dimnames = list(var_names, var_names))

  # [row=to, col=from] = coefficient
  m_true["x0", "x3"] <-  3.0
  m_true["x2", "x3"] <-  6.0
  m_true["x9", "x3"] <-  7.0
  m_true["x1", "x0"] <-  3.0
  m_true["x5", "x0"] <-  4.0
  m_true["x4", "x0"] <-  8.0
  m_true["x7", "x0"] <-  3.0
  m_true["x1", "x2"] <-  2.0
  m_true["x4", "x2"] <- -1.0
  m_true["x8", "x2"] <-  0.5
  m_true["x6", "x1"] <-  2.0
  m_true["x8", "x5"] <-  2.0
  m_true["x7", "x4"] <-  1.5
  m_true["x9", "x6"] <-  1.0

  return(list(
    data           = X,
    true_adjacency = m_true
  ))
}


#' Generate a challenging sample data for Direct LiNGAM
#'
#' Generates a dataset with conditions that make causal estimation difficult:
#'   1. High multicollinearity among predictors
#'   2. Moderate sample size relative to variables
#'   3. True coefficients of similar magnitude
#'
#' These conditions destabilize OLS initial estimates in Adaptive LASSO,
#' making Ridge-initialized Adaptive LASSO preferable.
#'
#' @param n number of samples (default: 200)
#' @param seed random seed (default: 42)
#' @param collinearity strength of multicollinearity (0 to 1, default: 0.95)
#' @return list(data, true_adjacency)
#'
#' @examples
#' result <- generate_lingam_hard_sample()
#' result$true_adjacency
#' head(result$data)
#'
#' @export
generate_lingam_hard_sample <- function(n = 200L,
                                        seed = 42L,
                                        collinearity = 0.95) {

  # --- Input validation ---
  if (!is.numeric(n) || n < 2) stop("n must be an integer >= 2.")
  if (!is.numeric(seed))       stop("seed must be numeric.")
  if (collinearity < 0 || collinearity >= 1) {
    stop("collinearity must be in [0, 1).")
  }

  n    <- as.integer(n)
  seed <- as.integer(seed)

  # --- 多重共線性の強さパラメータ ---
  noise_scale <- sqrt(1 - collinearity^2)

  # --- Error terms ---
  set.seed(seed)      ; e0  <- runif(n)
  set.seed(seed + 1L) ; e1  <- runif(n)
  set.seed(seed + 2L) ; e2  <- runif(n)
  set.seed(seed + 3L) ; e3  <- runif(n)
  set.seed(seed + 4L) ; e4  <- runif(n)
  set.seed(seed + 5L) ; e5  <- runif(n)
  set.seed(seed + 6L) ; e6  <- runif(n)
  set.seed(seed + 7L) ; e7  <- runif(n)
  set.seed(seed + 8L) ; e8  <- runif(n)
  set.seed(seed + 9L) ; e_c <- runif(n)

  # --- 強い共通因子（多重共線性の源泉）---
  common <- e_c

  x0 <- collinearity * common + noise_scale * e0
  x1 <- collinearity * common + noise_scale * e1
  x2 <- collinearity * common + noise_scale * e2
  x3 <- collinearity * common + noise_scale * e3
  x4 <- collinearity * common + noise_scale * e4

  x5 <- 1.5 * x0 + 1.5 * x1 + 1.5 * x2 + e5
  x6 <- 1.0 * x1 + 1.0 * x2 + 1.0 * x3 + 1.0 * x4 + e6
  x7 <- 2.0 * x0 + 2.0 * x3 + e7
  x8 <- 1.0 * x5 + 1.0 * x6 + e8

  X <- data.frame(
    x0 = x0, x1 = x1, x2 = x2, x3 = x3, x4 = x4,
    x5 = x5, x6 = x6, x7 = x7, x8 = x8
  )

  # --- 真の隣接行列 ---
  var_names <- names(X)
  p <- length(var_names)
  m_true <- matrix(0, p, p, dimnames = list(var_names, var_names))

  m_true["x5", "x0"] <- 1.5
  m_true["x5", "x1"] <- 1.5
  m_true["x5", "x2"] <- 1.5
  m_true["x6", "x1"] <- 1.0
  m_true["x6", "x2"] <- 1.0
  m_true["x6", "x3"] <- 1.0
  m_true["x6", "x4"] <- 1.0
  m_true["x7", "x0"] <- 2.0
  m_true["x7", "x3"] <- 2.0
  m_true["x8", "x5"] <- 1.0
  m_true["x8", "x6"] <- 1.0

  return(list(
    data           = X,
    true_adjacency = m_true
  ))
}
