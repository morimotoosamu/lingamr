test_that("evaluate_model_fit() errors clearly when lavaan is not installed", {
  local_mocked_bindings(
    check_lavaan_available = function() {
      stop(
        "Package 'lavaan' is required for evaluate_model_fit(). ",
        "Install it with install.packages(\"lavaan\").",
        call. = FALSE
      )
    }
  )
  B <- matrix(c(0, 0, 1, 0), nrow = 2)
  X <- matrix(rnorm(20), ncol = 2)
  expect_error(evaluate_model_fit(B, X), "lavaan")
})

skip_if_not_installed("lavaan")
skip_if_not_installed("glmnet") # default reg_method = "adaptive_lasso" needs a sparse matrix

test_that("evaluate_model_fit() returns a good fit for the correct model", {
  dat <- generate_lingam_sample_6(n = 1000, seed = 1)
  res <- lingam_direct(dat$data)
  fit_df <- evaluate_model_fit(res, dat$data)

  expect_s3_class(fit_df, "data.frame")
  expect_equal(nrow(fit_df), 1)
  expect_true(fit_df$CFI > 0.95)
  expect_true(fit_df$RMSEA < 0.1)
})

test_that("evaluate_model_fit() fit worsens when edges are reversed", {
  dat <- generate_lingam_sample_6(n = 1000, seed = 1)
  res <- lingam_direct(dat$data)

  fit_correct <- evaluate_model_fit(res, dat$data)

  B_reversed <- t(res$adjacency_matrix)
  fit_reversed <- evaluate_model_fit(B_reversed, dat$data)

  expect_true(fit_reversed$CFI < fit_correct$CFI)
})

test_that("evaluate_model_fit() accepts a result object directly", {
  dat <- generate_lingam_sample_6(n = 300, seed = 2)
  res <- lingam_direct(dat$data)

  fit_from_object <- evaluate_model_fit(res, dat$data)
  fit_from_matrix <- evaluate_model_fit(res$adjacency_matrix, dat$data)

  expect_equal(fit_from_object, fit_from_matrix)
})

test_that("evaluate_model_fit() handles NA (latent confounder) entries", {
  dat <- generate_lingam_sample_6(n = 500, seed = 3)
  res <- lingam_direct(dat$data)

  B <- res$adjacency_matrix
  B[1, 2] <- NA
  B[2, 1] <- NA

  fit_df <- expect_no_error(evaluate_model_fit(B, dat$data))
  expect_s3_class(fit_df, "data.frame")
  expect_equal(nrow(fit_df), 1)
})

test_that("evaluate_model_fit() validates its arguments", {
  X <- matrix(rnorm(30), ncol = 3)
  B_ok <- matrix(0, nrow = 3, ncol = 3)
  B_ok[1, 2] <- 1

  B_not_square <- matrix(0, nrow = 2, ncol = 3)
  expect_error(evaluate_model_fit(B_not_square, X), "square")

  B_wrong_dim <- matrix(0, nrow = 2, ncol = 2)
  expect_error(evaluate_model_fit(B_wrong_dim, X), "ncol")

  X_na <- X
  X_na[1, 1] <- NA
  expect_error(evaluate_model_fit(B_ok, X_na), "missing")

  expect_error(evaluate_model_fit(B_ok, X, is_ordinal = c(TRUE, FALSE)), "is_ordinal")
})

test_that("evaluate_model_fit() supports is_ordinal", {
  # moderate (not near-deterministic) signal so dichotomizing x0 keeps the
  # model estimable; B[i, j] = j -> i, so x1 -> x0
  set.seed(42)
  n <- 500
  x1 <- rnorm(n)
  x0 <- 0.5 * x1 + rnorm(n)
  X <- cbind(x0, x1)
  B <- matrix(c(0, 0.5, 0, 0), nrow = 2, byrow = TRUE)

  X_ord <- X
  X_ord[, 1] <- as.integer(X[, 1] > median(X[, 1]))
  is_ordinal <- c(TRUE, FALSE)

  expect_no_error(evaluate_model_fit(B, X_ord, is_ordinal = is_ordinal))
})

test_that("evaluate_model_fit() errors clearly when the SEM does not converge", {
  local_mocked_bindings(lavInspect = function(...) FALSE, .package = "lavaan")

  dat <- generate_lingam_sample_6(n = 300, seed = 4)
  res <- lingam_direct(dat$data)

  expect_error(evaluate_model_fit(res, dat$data), "did not converge")
})

test_that("evaluate_model_fit() errors clearly for an edgeless adjacency matrix", {
  X <- matrix(rnorm(30), ncol = 3)
  B_empty <- matrix(0, nrow = 3, ncol = 3)
  expect_error(evaluate_model_fit(B_empty, X), "no edges")
})
