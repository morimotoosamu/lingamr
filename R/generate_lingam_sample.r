# =============================================================================
# サンプルデータ生成関数
#
# 共通ヘルパー:
#   make_noise_fn()         : ノイズ分布名 → 乱数生成関数（分布名の検証も担う）
#   validate_sample_args()  : n / seed の共通バリデーション
#   generate_noise_matrix() : 変数ごとに独立シードでノイズ列を生成
#   build_true_adjacency()  : エッジ指定から真の隣接行列を構築
# =============================================================================


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


#' サンプル生成関数の n / seed を検証して整数化する
#'
#' @param n サンプルサイズ
#' @param seed 乱数シード
#' @return list(n, seed)（いずれも integer）
#' @keywords internal
validate_sample_args <- function(n, seed) {
  if (!is.numeric(n) || n < 2) stop("n must be an integer >= 2.")
  if (!is.numeric(seed))       stop("seed must be numeric.")
  list(n = as.integer(n), seed = as.integer(seed))
}


#' 変数ごとに独立シードでノイズ行列を生成する
#'
#' 列 k は `set.seed(seed + k - 1)` の直後に `noise_fn(n)` で生成される。
#' 変数ごとにシードを固定することで、同じ seed なら常に同じノイズ列が得られる。
#'
#' @param n サンプルサイズ
#' @param n_vars 変数（列）の数
#' @param seed 基準シード。列 k には seed + k - 1 を用いる
#' @param noise_fn `function(n)` 形式のノイズ生成関数
#' @return ノイズ行列 (n x n_vars)
#' @keywords internal
generate_noise_matrix <- function(n, n_vars, seed, noise_fn) {
  E <- matrix(0, nrow = n, ncol = n_vars)
  for (k in seq_len(n_vars)) {
    set.seed(seed + k - 1L)
    E[, k] <- noise_fn(n)
  }
  E
}


#' エッジ指定から真の隣接行列を構築する
#'
#' @param var_names 変数名ベクトル
#' @param from エッジの原因変数名ベクトル
#' @param to エッジの結果変数名ベクトル（from と同じ長さ）
#' @param coef エッジ係数ベクトル（from と同じ長さ）
#' @return 隣接行列 (p x p)。`m[to, from] = coef`（行 = to, 列 = from）
#' @keywords internal
build_true_adjacency <- function(var_names, from, to, coef) {
  p <- length(var_names)
  m <- matrix(0, p, p, dimnames = list(var_names, var_names))
  m[cbind(to, from)] <- coef
  m
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
  args <- validate_sample_args(n, seed)
  n <- args$n
  seed <- args$seed

  # ノイズ生成関数（不正な noise_dist はここでエラーになる）
  noise_fn <- make_noise_fn(noise_dist)

  # --- Error terms（各変数に独立したシード）---
  E <- generate_noise_matrix(n, 6L, seed, noise_fn)

  # --- Data generation（因果順序に従う）---
  x3 <- E[, 4]
  x0 <- 3.0 * x3 + E[, 1]
  x2 <- 6.0 * x3 + E[, 3]
  x1 <- 3.0 * x0 + 2.0 * x2 + E[, 2]
  x5 <- 4.0 * x0 + E[, 6]
  x4 <- 8.0 * x0 - 1.0 * x2 + E[, 5]

  X <- data.frame(x0 = x0, x1 = x1, x2 = x2,
             x3 = x3, x4 = x4, x5 = x5)

  m_true <- build_true_adjacency(
    names(X),
    from = c("x3", "x3", "x0", "x2", "x0", "x0", "x2"),
    to   = c("x0", "x2", "x1", "x1", "x5", "x4", "x4"),
    coef = c(3.0, 6.0, 3.0, 2.0, 4.0, 8.0, -1.0)
  )

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
  args <- validate_sample_args(n, seed)
  n <- args$n
  seed <- args$seed

  # ノイズ生成関数（不正な noise_dist はここでエラーになる）
  noise_fn <- make_noise_fn(noise_dist)

  # --- Error terms（各変数に独立したシード）---
  E <- generate_noise_matrix(n, 10L, seed, noise_fn)

  # --- Data generation（因果順序に従う）---
  x3 <- E[, 4]
  x0 <- 3.0 * x3 + E[, 1]
  x2 <- 6.0 * x3 + E[, 3]
  x1 <- 3.0 * x0 + 2.0 * x2 + E[, 2]
  x5 <- 4.0 * x0 + E[, 6]
  x4 <- 8.0 * x0 - 1.0 * x2 + E[, 5]
  x6 <- 2.0 * x1 + E[, 7]
  x7 <- 1.5 * x4 + 3.0 * x0 + E[, 8]
  x8 <- 0.5 * x2 + 2.0 * x5 + E[, 9]
  x9 <- 7.0 * x3 + 1.0 * x6 + E[, 10]

  X <- data.frame(x0 = x0, x1 = x1, x2 = x2, x3 = x3, x4 = x4,
             x5 = x5, x6 = x6, x7 = x7, x8 = x8, x9 = x9)

  m_true <- build_true_adjacency(
    names(X),
    from = c("x3", "x3", "x3", "x0", "x0", "x0", "x0",
             "x2", "x2", "x2", "x1", "x5", "x4", "x6"),
    to   = c("x0", "x2", "x9", "x1", "x5", "x4", "x7",
             "x1", "x4", "x8", "x6", "x8", "x7", "x9"),
    coef = c(3.0, 6.0, 7.0, 3.0, 4.0, 8.0, 3.0,
             2.0, -1.0, 0.5, 2.0, 2.0, 1.5, 1.0)
  )

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
  args <- validate_sample_args(n, seed)
  n <- args$n
  seed <- args$seed
  if (collinearity < 0 || collinearity >= 1) {
    stop("collinearity must be in [0, 1).")
  }

  # --- 多重共線性の強さパラメータ ---
  noise_scale <- sqrt(1 - collinearity^2)

  # --- Error terms（各変数に独立したシード。10列目は共通因子）---
  E <- generate_noise_matrix(n, 10L, seed, runif)

  # --- 強い共通因子（多重共線性の源泉）---
  common <- E[, 10]

  x0 <- collinearity * common + noise_scale * E[, 1]
  x1 <- collinearity * common + noise_scale * E[, 2]
  x2 <- collinearity * common + noise_scale * E[, 3]
  x3 <- collinearity * common + noise_scale * E[, 4]
  x4 <- collinearity * common + noise_scale * E[, 5]

  x5 <- 1.5 * x0 + 1.5 * x1 + 1.5 * x2 + E[, 6]
  x6 <- 1.0 * x1 + 1.0 * x2 + 1.0 * x3 + 1.0 * x4 + E[, 7]
  x7 <- 2.0 * x0 + 2.0 * x3 + E[, 8]
  x8 <- 1.0 * x5 + 1.0 * x6 + E[, 9]

  X <- data.frame(
    x0 = x0, x1 = x1, x2 = x2, x3 = x3, x4 = x4,
    x5 = x5, x6 = x6, x7 = x7, x8 = x8
  )

  m_true <- build_true_adjacency(
    names(X),
    from = c("x0", "x1", "x2", "x1", "x2", "x3", "x4", "x0", "x3", "x5", "x6"),
    to   = c("x5", "x5", "x5", "x6", "x6", "x6", "x6", "x7", "x7", "x8", "x8"),
    coef = c(1.5, 1.5, 1.5, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 1.0, 1.0)
  )

  return(list(
    data           = X,
    true_adjacency = m_true
  ))
}

#' Generate Paradoxical Data Where DirectLiNGAM Struggles
#'
#' Generates a synthetic dataset designed to favor ICA-LiNGAM (due to standardized scales)
#' while challenging DirectLiNGAM (due to heavy measurement noise on the root variable,
#' which triggers error propagation). The true causal structure is a serial chain:
#' \code{x0 -> x1 -> x2 -> x3} (each coefficient 0.8).
#'
#' @param n number of samples (default: 2000)
#' @param seed random seed (default: 42)
#'
#' @return list(data, true_adjacency)
#' * \code{data}: a data frame with 4 standardized variables (\code{x0}, \code{x1},
#'   \code{x2}, \code{x3}); each column has a mean of 0 and a standard deviation of 1.
#' * \code{true_adjacency}: the 4x4 true adjacency matrix of the data-generating
#'   chain, following the \code{m[row = to, col = from]} convention and holding the
#'   structural coefficients (0.8) on the latent, pre-standardization scale.
#'
#' @details
#' This function intentionally injects strong measurement error into the root (causal upstream)
#' variable \code{x0}. This noise corrupts the independence tests performed at the initial step
#' of DirectLiNGAM, frequently causing it to misidentify the root variable and leading to
#' a cascading failure (error propagation) throughout the causal ordering.
#'
#' On the other hand, the output data is completely standardized using the \code{scale()} function.
#' This eliminates any differences in scale among the variables, thereby neutralizing the major
#' weakness of ICA-LiNGAM (scale-dependence) and allowing it to perform relatively better.
#'
#' Because the data are standardized and the root carries measurement error, the coefficients
#' estimated by \code{lingam_direct()} will not exactly match the 0.8 values stored in
#' \code{true_adjacency}.
#'
#' @examples
#' # Generate the dataset
#' paradox <- generate_lingam_paradox_data(n = 1000, seed = 123)
#'
#' # Verify the dataset
#' head(paradox$data)
#' sapply(paradox$data, sd)
#'
#' # True data-generating structure
#' paradox$true_adjacency
#'
#' @export
generate_lingam_paradox_data <- function(n = 2000L, seed = 42L) {
  args <- validate_sample_args(n, seed)
  n <- args$n
  seed <- args$seed

  # 1. Set seed for reproducibility
  set.seed(seed)

  # 2. Generate intrinsic noise for each variable (Non-Gaussian uniform distribution)
  e0 <- runif(n, min = -1.0, max = 1.0)
  e1 <- runif(n, min = -0.5, max = 0.5)
  e2 <- runif(n, min = -0.5, max = 0.5)
  e3 <- runif(n, min = -0.5, max = 0.5)

  # 3. Build the true causal chain (x0 -> x1 -> x2 -> x3)
  x0_true <- e0
  x1      <- 0.8 * x0_true + e1
  x2      <- 0.8 * x1 + e2
  x3      <- 0.8 * x2 + e3

  # 4. Inject heavy measurement noise into the root variable x0
  measurement_error_x0 <- runif(n, min = -1.2, max = 1.2)
  x0_observed          <- x0_true + measurement_error_x0

  # 5. Create the data frame
  df <- data.frame(x0 = x0_observed, x1 = x1, x2 = x2, x3 = x3)

  # 6. Standardize all variables to unify their scales
  df_scaled <- as.data.frame(scale(df))

  m_true <- build_true_adjacency(
    names(df_scaled),
    from = c("x0", "x1", "x2"),
    to   = c("x1", "x2", "x3"),
    coef = c(0.8, 0.8, 0.8)
  )

  return(list(
    data           = df_scaled,
    true_adjacency = m_true
  ))
}


#' Generate large-scale sample data to benchmark Direct LiNGAM scalability
#'
#' Generates a dataset with many variables to demonstrate the computational
#' scalability difference between Direct LiNGAM and ICA-LiNGAM.
#'
#' @param p number of variables (default: 20)
#' @param n number of observations (default: 1000)
#' @param max_parents maximum number of parents per node (default: 3).
#'   Controls graph density. Each variable xi (i >= 1) receives between 1 and
#'   `min(max_parents, i)` parents drawn from x0, ..., x(i-1).
#' @param coef_min minimum absolute value of edge coefficients (default: 0.5)
#' @param coef_max maximum absolute value of edge coefficients (default: 1.5)
#' @param seed random seed (default: 42)
#' @param noise_dist error term distribution.
#'   "uniform"     : Uniform(0, 1) - default, non-Gaussian (LiNGAM works well)
#'   "gaussian"    : Normal(0, 1)  - LiNGAM may fail
#'   "lognormal"   : Log-normal(0, 1) - skewed, non-Gaussian
#'   "exponential" : Exponential(1) - skewed, non-Gaussian
#'   "t3"          : t-distribution (df=3) - heavy tails
#'
#' @return A list with three elements:
#'   * `data`: data.frame with `p` columns (x0, x1, ..., x(p-1)).
#'   * `true_adjacency`: p x p matrix. `true_adjacency[i, j]` is the structural
#'     coefficient of the edge xj -> xi (row = to, col = from). The matrix is
#'     strictly lower-triangular because variables are stored in causal order.
#'   * `true_causal_order`: integer vector `0:(p-1)`. Variables are already in
#'     topological order, so the true causal order is always 0, 1, ..., p-1.
#'
#' @details
#' ## Why Direct LiNGAM slows down with large p
#'
#' At each of its `p` steps, Direct LiNGAM evaluates an independence measure
#' between every remaining candidate root and every other residual.
#' The total number of evaluations is:
#'
#' \deqn{\sum_{k=1}^{p} k(k-1) \approx \frac{p^3}{3}}
#'
#' i.e., O(p^3). Each evaluation is itself O(n), giving O(p^3 n) overall.
#' For p = 10 this is about 330 evaluations; for p = 20 about 2,660;
#' for p = 40 about 21,320 --- an 8x increase every time p doubles.
#'
#' ## Why ICA-LiNGAM scales better
#'
#' ICA-LiNGAM applies FastICA once to the whole p x n data matrix.
#' Each FastICA iteration costs O(p^2 n), and the algorithm typically
#' converges in far fewer than p iterations.  Additionally, these matrix
#' operations are fully vectorised (BLAS/LAPACK), whereas Direct LiNGAM
#' iterates over pairs in an R loop.
#'
#' ## Data-generating process
#'
#' Variables are topologically ordered as x0, x1, ..., x(p-1).
#' For each i >= 1, the number of parents is sampled uniformly from
#' 1 to `min(max_parents, i)`, and the parents are drawn without replacement
#' from x0, ..., x(i-1).  Edge coefficients are drawn uniformly from
#' \[-coef_max, -coef_min\] U \[coef_min, coef_max\].
#' The resulting adjacency matrix is strictly lower-triangular.
#'
#' @examples
#' # 20変数のデータを生成してスパース性を確認
#' dataset <- generate_lingam_large_sample(p = 20, n = 500)
#' dim(dataset$data)                    # 500 x 20
#' sum(dataset$true_adjacency != 0)     # 辺の本数
#' dataset$true_causal_order            # 0, 1, ..., 19
#'
#' \donttest{
#' # 変数数が増えると Direct LiNGAM の実行時間が急増する
#' t10 <- system.time(lingam_direct(generate_lingam_large_sample(p = 10)$data))
#' t20 <- system.time(lingam_direct(generate_lingam_large_sample(p = 20)$data))
#' cat(sprintf("p=10: %.1f sec,  p=20: %.1f sec\n", t10["elapsed"], t20["elapsed"]))
#' }
#'
#' @export
generate_lingam_large_sample <- function(
    p           = 20L,
    n           = 1000L,
    max_parents = 3L,
    coef_min    = 0.5,
    coef_max    = 1.5,
    seed        = 42L,
    noise_dist  = "uniform") {

  # --- Input validation ---
  if (!is.numeric(p) || length(p) != 1L || p < 2)
    stop("p must be an integer >= 2.", call. = FALSE)
  if (!is.numeric(n) || length(n) != 1L || n < 2)
    stop("n must be an integer >= 2.", call. = FALSE)
  if (!is.numeric(max_parents) || length(max_parents) != 1L || max_parents < 1)
    stop("max_parents must be an integer >= 1.", call. = FALSE)
  if (!is.numeric(coef_min) || !is.numeric(coef_max) ||
      coef_min <= 0 || coef_max <= coef_min)
    stop("coef_min and coef_max must satisfy 0 < coef_min < coef_max.", call. = FALSE)
  if (!is.numeric(seed))
    stop("seed must be numeric.", call. = FALSE)

  valid_dists <- c("uniform", "gaussian", "lognormal", "exponential", "t3")
  if (!(noise_dist %in% valid_dists)) {
    stop(sprintf("noise_dist must be one of: %s",
                 paste(valid_dists, collapse = ", ")), call. = FALSE)
  }

  p           <- as.integer(p)
  n           <- as.integer(n)
  max_parents <- as.integer(max_parents)
  seed        <- as.integer(seed)

  var_names <- paste0("x", seq_len(p) - 1L)
  noise_fn  <- make_noise_fn(noise_dist)

  # --- Step 1: Build random DAG adjacency matrix ---
  # Single set.seed covers both structure and data generation deterministically.
  set.seed(seed)

  B <- matrix(0.0, p, p, dimnames = list(var_names, var_names))

  for (i in seq(2L, p)) {
    # Candidate parents: R indices 1..(i-1) => variables x0..x(i-2)
    candidates  <- seq_len(i - 1L)
    n_pa_max    <- min(max_parents, length(candidates))
    n_pa_chosen <- sample.int(n_pa_max, 1L)               # how many parents
    pa_indices  <- sort(sample(candidates, n_pa_chosen))   # which parents
    for (pa in pa_indices) {
      sign_coef <- sample(c(-1.0, 1.0), 1L)
      B[i, pa]  <- sign_coef * runif(1L, coef_min, coef_max)
    }
  }

  # --- Step 2: Generate observations in causal order ---
  X_mat <- matrix(0.0, n, p, dimnames = list(NULL, var_names))
  for (i in seq_len(p)) {
    xi_val <- noise_fn(n)                       # intrinsic noise
    if (i > 1L) {
      pa_idx <- which(B[i, ] != 0.0)
      for (pa in pa_idx) {
        xi_val <- xi_val + B[i, pa] * X_mat[, pa]
      }
    }
    X_mat[, i] <- xi_val
  }

  return(list(
    data              = as.data.frame(X_mat),
    true_adjacency    = B,
    true_causal_order = seq_len(p) - 1L
  ))
}
