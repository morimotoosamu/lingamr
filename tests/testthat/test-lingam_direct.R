test_that("lingam_direct returns LingamResult with correct structure", {
  dat <- sample6_500_s42()
  res <- fit_direct_500()

  expect_s3_class(res, "LingamResult")
  expect_named(res, c("adjacency_matrix", "causal_order"))
  expect_true(is.matrix(res$adjacency_matrix))
  expect_equal(dim(res$adjacency_matrix), c(6L, 6L))
  expect_equal(length(res$causal_order), 6L)
  # column names are preserved
  expect_equal(colnames(res$adjacency_matrix), names(dat$data))
})

test_that("lingam_direct identifies x3 as root (first in causal order)", {
  # x3 is an exogenous variable with no parents, so it should be first in the causal order
  dat <- sample6_2000_s42()
  res <- fit_direct_2000()

  expect_equal(res$causal_order[1], which(names(dat$data) == "x3"))
})

test_that("lingam_direct accepts data.frame and preserves colnames", {
  dat <- sample6_300()
  res <- fit_direct_300()

  expect_equal(colnames(res$adjacency_matrix), names(dat$data))
  expect_equal(rownames(res$adjacency_matrix), names(dat$data))
})

test_that("lingam_direct accepts matrix input", {
  dat <- sample6_300()
  mat <- as.matrix(dat$data)
  res <- lingam_direct(mat, reg_method = "ols")

  expect_s3_class(res, "LingamResult")
})

test_that("lingam_direct errors on invalid inputs", {
  dat <- sample6_100()

  expect_error(lingam_direct(dat$data, measure = "bad_measure"), "should be one of")
  expect_error(lingam_direct(dat$data, reg_method = "bad_method"), "should be one of")
  expect_error(lingam_direct(dat$data, lambda = "bad_lambda"), "should be one of")
  expect_error(lingam_direct(matrix(letters[1:4], nrow = 2)), "numeric matrix")  # non-numeric
  expect_error(
    lingam_direct(as.data.frame(matrix(1, nrow = 5, ncol = 1))),
    "at least 2 variables"
  )  # 1 variable
  expect_error(
    lingam_direct(matrix(numeric(0), nrow = 0, ncol = 3)),
    "at least 2 observations"
  )  # 0 rows
})

test_that("lingam_direct rejects constant and perfectly collinear columns", {
  set.seed(1)
  X_const <- cbind(x0 = runif(100), x1 = runif(100), x2 = rep(1, 100))
  # the offending column is named in the message
  expect_error(lingam_direct(X_const, reg_method = "ols"), "constant columns: x2")

  x <- runif(200)
  X_dup <- cbind(a = x, b = x, c = runif(200))
  expect_error(lingam_direct(X_dup, reg_method = "ols"), "linearly dependent")
  # a column equal to another plus a constant offset is also caught
  X_offset <- cbind(a = x, b = x + 5, c = runif(200))
  expect_error(lingam_direct(X_offset, reg_method = "ols"), "linearly dependent")
})

test_that("lingam_direct rejects prior knowledge values outside {-1, 0, 1}", {
  dat <- sample6_100()
  pk <- matrix(-1, 6, 6)
  pk[1, 2] <- 0.5
  expect_error(
    lingam_direct(dat$data, prior_knowledge = pk, reg_method = "ols"),
    "must contain only -1"
  )
  pk[1, 2] <- 2
  expect_error(
    lingam_direct(dat$data, prior_knowledge = pk, reg_method = "ols"),
    "must contain only -1"
  )
})

test_that("single-predictor edges are pruned for independent variables", {
  skip_if_not_installed("glmnet")
  # The second variable in the causal order has exactly one predictor, which
  # is fitted by the IC-pruned OLS fallback (glmnet needs >= 2 columns). For
  # fully independent data every edge, including that one, must be zero.
  set.seed(7)
  X <- cbind(a = runif(2000), b = runif(2000), c = runif(2000))
  res <- lingam_direct(X)  # default adaptive_lasso + BIC
  expect_equal(sum(abs(res$adjacency_matrix) > 0), 0L)
})

test_that("print.LingamResult runs without error", {
  res <- fit_direct_200()

  expect_output(print(res), "Direct LiNGAM Result")
  expect_output(print(res), "Causal order")
  expect_output(print(res), "Adjacency matrix")
})

test_that("lingam_direct with prior_knowledge runs without error", {
  dat <- sample6_300()
  pk <- make_prior_knowledge(6,
    exogenous_variables = 4,
    labels = names(dat$data)
  )
  res <- lingam_direct(dat$data, prior_knowledge = pk, reg_method = "ols")

  expect_s3_class(res, "LingamResult")
})

test_that("soft prior knowledge with -1 (unknown) entries runs without error", {
  # regression test: -1 becomes NA during preprocessing, and
  # if (sum(...) == 0) inside search_candidate() returned NA and crashed
  dat <- sample6_200()
  pk <- matrix(-1L, 6, 6)
  pk[2, 1] <- 1L

  res <- lingam_direct(dat$data,
    prior_knowledge = pk,
    apply_prior_knowledge_softly = TRUE,
    reg_method = "ols"
  )

  expect_s3_class(res, "LingamResult")
  expect_setequal(res$causal_order, 1:6)
})

test_that("search_candidate handles NA entries like Python's NaN semantics", {
  # In the Python version, for a row containing NaN, sum() == 0 is False (not an error)
  U <- 1:3

  # row 1: all 0 (known exogenous variable) -> enters Uc
  # row 2: contains NA -> not an exogenous candidate
  # row 3: contains 1 -> not an exogenous candidate
  Aknw <- matrix(NA_real_, 3, 3)
  Aknw[1, 2:3] <- 0
  Aknw[2, 3] <- 0
  Aknw[3, 2] <- 1

  cand <- search_candidate(U, Aknw,
    apply_prior_knowledge_softly = TRUE,
    partial_orders = NULL
  )
  expect_equal(cand$Uc, 1L)

  # even when all are NA (entirely unknown), it does not crash, and Uc becomes all of U
  Aknw_all_na <- matrix(NA_real_, 3, 3)
  cand_na <- search_candidate(U, Aknw_all_na,
    apply_prior_knowledge_softly = TRUE,
    partial_orders = NULL
  )
  expect_equal(cand_na$Uc, U)
  expect_equal(cand_na$Vj, integer(0))
})

test_that("measure = 'kernel' returns a valid LingamResult", {
  dat <- generate_lingam_sample_6(n = 100, seed = 42)
  res <- lingam_direct(dat$data, measure = "kernel", reg_method = "ols")

  expect_s3_class(res, "LingamResult")
  expect_equal(sort(res$causal_order), 1:6)
  expect_equal(dim(res$adjacency_matrix), c(6L, 6L))
})

test_that("kernel mutual information (prepare + core) returns finite values", {
  set.seed(1)
  x <- rnorm(100)
  y <- 0.8 * x + runif(100)
  z <- rnorm(100)
  kappa <- 2e-2
  sigma <- 1.0

  E1 <- kernel_mi_prepare(x, kappa, sigma)
  mi_dep <- kernel_mi_core(E1, y, kappa, sigma)
  mi_ind <- kernel_mi_core(E1, z, kappa, sigma)

  expect_true(is.finite(mi_dep) && is.finite(mi_ind))
  # the MI of the dependent pair is greater than that of the independent pair
  expect_gt(mi_dep, mi_ind)
})

test_that("incomplete_cholesky_gauss reconstructs the Gaussian Gram matrix", {
  set.seed(3)
  n <- 300
  sigma <- 1.0
  x <- rnorm(n)

  G <- incomplete_cholesky_gauss(x, sigma)
  K_true <- exp(-1 / (2 * sigma^2) * outer(x, x, "-")^2)
  K_approx <- tcrossprod(G)

  expect_lt(max(abs(K_true - K_approx)), 1e-3)
})

test_that("low-rank kernel MI agrees with the exact computation", {
  set.seed(4)
  n <- 300
  kappa <- 2e-2
  sigma <- 1.0
  x <- rnorm(n)
  y_dep <- 0.7 * x + rnorm(n, sd = 0.5)
  y_ind <- rnorm(n)

  E1 <- kernel_mi_prepare(x, kappa, sigma)
  mi_dep_exact <- kernel_mi_core(E1, y_dep, kappa, sigma)
  mi_ind_exact <- kernel_mi_core(E1, y_ind, kappa, sigma)

  prep1 <- kernel_mi_prepare_lowrank(x, kappa, sigma)
  mi_dep_lr <- kernel_mi_core_lowrank(prep1, y_dep, kappa, sigma)
  mi_ind_lr <- kernel_mi_core_lowrank(prep1, y_ind, kappa, sigma)

  expect_lt(abs(mi_dep_exact - mi_dep_lr), 1e-3)
  expect_lt(abs(mi_ind_exact - mi_ind_lr), 1e-3)
  # the low-rank path must preserve which pair has the larger MI
  expect_equal(mi_dep_exact > mi_ind_exact, mi_dep_lr > mi_ind_lr)
})

test_that("search_causal_order_kernel uses the low-rank path above n = 1000 and still finds a valid order", {
  set.seed(5)
  n <- 1500
  e1 <- rnorm(n)
  e2 <- rnorm(n)
  e3 <- rnorm(n)
  x1 <- e1
  x2 <- 0.8 * x1 + e2
  x3 <- 0.5 * x1 + 0.6 * x2 + e3
  X <- cbind(x1, x2, x3)

  result <- search_causal_order_kernel(X, U = 1:3, Uc = 1:3, Vj = integer(0))
  expect_true(result %in% 1:3)
})

test_that("reg_method = 'ols' does not require glmnet", {
  dat <- sample6_200()

  # simulate an environment where glmnet is not available
  local_mocked_bindings(
    check_glmnet_available = function(method) {
      stop(sprintf(
        "Package 'glmnet' is required for reg_method = \"%s\". Please install it.",
        method
      ), call. = FALSE)
    }
  )

  expect_no_error(lingam_direct(dat$data, reg_method = "ols"))
  expect_error(lingam_direct(dat$data, reg_method = "lasso"), "glmnet")
  expect_error(lingam_direct(dat$data, reg_method = "adaptive_lasso"), "glmnet")
  expect_error(lingam_direct(dat$data, reg_method = "ridge"), "glmnet")
})

test_that("reg_method = 'ridge' returns a valid LingamResult", {
  skip_if_not_installed("glmnet")
  dat <- generate_lingam_sample_6(n = 300, seed = 42)
  res <- lingam_direct(dat$data, reg_method = "ridge")

  expect_s3_class(res, "LingamResult")
  expect_equal(dim(res$adjacency_matrix), c(6L, 6L))
  expect_equal(length(res$causal_order), 6L)
  # Ridge is not sparse, so there are many non-zero coefficients
  expect_gt(sum(res$adjacency_matrix != 0), 0L)
})

test_that("reg_method = 'ridge' with lambda = 'oracle' errors", {
  dat <- sample6_200()

  expect_error(
    lingam_direct(dat$data, reg_method = "ridge", lambda = "oracle"),
    "oracle"
  )
})
