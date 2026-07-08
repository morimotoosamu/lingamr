test_that("lingam_rcd returns RCDResult with correct structure", {
  dat <- generate_rcd_sample(n = 300, seed = 42)
  res <- lingam_rcd(dat$data)

  expect_s3_class(res, "RCDResult")
  expect_named(res, c("adjacency_matrix", "ancestors_list"))
  expect_true(is.matrix(res$adjacency_matrix))
  expect_equal(dim(res$adjacency_matrix), c(6L, 6L))
  expect_equal(colnames(res$adjacency_matrix), names(dat$data))
  expect_true(is.list(res$ancestors_list))
  expect_length(res$ancestors_list, 6L)
  expect_equal(names(res$ancestors_list), names(dat$data))
})

test_that("lingam_rcd validates its inputs", {
  dat <- generate_rcd_sample(n = 300, seed = 42)

  expect_error(lingam_rcd(matrix(letters[1:4], nrow = 2)), "numeric matrix")
  expect_error(lingam_rcd(dat$data[1:2, ]), "at least 3 observations")
  expect_error(lingam_rcd(dat$data, max_explanatory_num = 0L), "max_explanatory_num")
  expect_error(lingam_rcd(dat$data, cor_alpha = -1), "cor_alpha")
  expect_error(lingam_rcd(dat$data, ind_alpha = -1), "ind_alpha")
  expect_error(lingam_rcd(dat$data, shapiro_alpha = -1), "shapiro_alpha")
  expect_error(lingam_rcd(dat$data, independence = "foo"), "should be one of")
  expect_error(lingam_rcd(dat$data, ind_corr = -1), "ind_corr")
  expect_error(lingam_rcd(dat$data, MLHSICR = "yes"), "MLHSICR")
})

test_that("lingam_rcd recovers the true ancestor sets", {
  dat <- generate_rcd_sample(n = 300, seed = 42)
  res <- lingam_rcd(dat$data)

  for (i in seq_along(dat$ancestors_list)) {
    expect_true(all(dat$ancestors_list[[i]] %in% res$ancestors_list[[i]]))
  }

  exact_match <- vapply(seq_along(dat$ancestors_list), function(i) {
    setequal(dat$ancestors_list[[i]], res$ancestors_list[[i]])
  }, logical(1))
  expect_true(any(exact_match))
})

test_that("lingam_rcd detects the confounded pair as NA and recovers true edges", {
  dat <- generate_rcd_sample(n = 300, seed = 42)
  res <- lingam_rcd(dat$data)

  x2 <- dat$confounded_pair[1]
  x4 <- dat$confounded_pair[2]
  expect_true(is.na(res$adjacency_matrix[x2, x4]))
  expect_true(is.na(res$adjacency_matrix[x4, x2]))

  x0 <- 1L
  x1 <- 2L
  x3 <- 4L
  x5 <- 6L
  # true edges (all positive coefficients): x1 -> x0, x3 -> x0, x5 -> x1,
  # x0 -> x2, x0 -> x4
  expect_gt(res$adjacency_matrix[x0, x1], 0)
  expect_gt(res$adjacency_matrix[x0, x3], 0)
  expect_gt(res$adjacency_matrix[x1, x5], 0)
  expect_gt(res$adjacency_matrix[x2, x0], 0)
  expect_gt(res$adjacency_matrix[x4, x0], 0)
})

test_that("MLHSICR = TRUE runs without error and returns valid structure", {
  skip_on_cran()
  dat <- generate_rcd_sample(n = 300, seed = 42)
  res <- lingam_rcd(dat$data, MLHSICR = TRUE, max_explanatory_num = 2L)

  expect_s3_class(res, "RCDResult")
  expect_equal(dim(res$adjacency_matrix), c(6L, 6L))
  expect_length(res$ancestors_list, 6L)
})

test_that("independence = 'fcorr' returns a valid RCDResult", {
  dat <- generate_rcd_sample(n = 300, seed = 42)
  res <- lingam_rcd(dat$data, independence = "fcorr")

  expect_s3_class(res, "RCDResult")
  expect_equal(dim(res$adjacency_matrix), c(6L, 6L))
})

test_that("mlhsicr_regression returns finite coef/resid on well-behaved input", {
  set.seed(50)
  n <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 0.5 * x1 - 0.3 * x2 + rnorm(n)
  Y <- cbind(y = y, x1 = x1, x2 = x2)

  out <- mlhsicr_regression(Y, xi = 1L, xj_list = c(2L, 3L))
  expect_true(all(is.finite(out$coef)))
  expect_true(all(is.finite(out$resid)))
})

test_that("mlhsicr_regression falls back to the OLS solution when the HSIC objective cannot be minimized", {
  # a near-constant target collapses width_xi to ~0, which combined with any
  # nonzero coefficient immediately pushes width <= 0 (the objective's
  # invalid-region sentinel) everywhere optim() explores; the guard should
  # fall back to the OLS coefficients/residual instead of returning a
  # coefficient vector that never actually minimized HSIC.
  set.seed(51)
  n <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- rep(1e-12, n)
  Y <- cbind(y = y, x1 = x1, x2 = x2)

  ols <- rcd_ols_resid_coef(y, cbind(x1, x2))
  out <- mlhsicr_regression(Y, xi = 1L, xj_list = c(2L, 3L))

  expect_true(all(is.finite(out$coef)))
  expect_true(all(is.finite(out$resid)))
  expect_equal(out$coef, as.vector(ols$coef))
})

test_that("estimate_total_effect_rcd warns and returns NA for a confounded 'from', numeric otherwise", {
  dat <- generate_rcd_sample(n = 300, seed = 42)
  res <- lingam_rcd(dat$data)

  x2 <- dat$confounded_pair[1]
  x4 <- dat$confounded_pair[2]
  expect_warning(
    eff <- estimate_total_effect_rcd(dat$data, res, from_index = x2, to_index = x4),
    "latent confounder"
  )
  expect_true(is.na(eff))

  # a well-identified pair (x5 -> x0, via x1/x3) returns a numeric estimate
  x0 <- 1L
  x5 <- 6L
  eff2 <- estimate_total_effect_rcd(dat$data, res, from_index = x5, to_index = x0)
  expect_true(is.numeric(eff2) && !is.na(eff2))
})

test_that("get_error_independence_p_values_rcd gives NA for NA-linked pairs and [0,1] otherwise", {
  dat <- generate_rcd_sample(n = 300, seed = 42)
  res <- lingam_rcd(dat$data)

  pvals <- get_error_independence_p_values_rcd(dat$data, res)
  expect_equal(dim(pvals), c(6L, 6L))

  x2 <- dat$confounded_pair[1]
  x4 <- dat$confounded_pair[2]
  expect_true(is.na(pvals[x2, x4]))

  finite_vals <- pvals[!is.na(pvals)]
  expect_true(all(finite_vals >= 0 & finite_vals <= 1))
})

test_that("print.RCDResult runs without error and shows ancestor sets", {
  dat <- generate_rcd_sample(n = 300, seed = 42)
  res <- lingam_rcd(dat$data)

  expect_output(print(res), "RCD")
  expect_output(print(res), "M\\(")
})

test_that("lingam_rcd_bootstrap returns a usable BootstrapResult", {
  dat <- generate_rcd_sample(n = 300, seed = 42)

  bs <- lingam_rcd_bootstrap(dat$data, n_sampling = 4L, seed = 42, verbose = FALSE)

  expect_s3_class(bs, "BootstrapResult")
  expect_null(bs$causal_orders)
  expect_false(anyNA(bs$adjacency_matrices))
  expect_false(anyNA(bs$total_effects))

  p <- get_probabilities(bs)
  expect_equal(dim(p), c(6L, 6L))
  expect_false(anyNA(p))
})

test_that("lingam_rcd_bootstrap is reproducible with the same seed", {
  dat <- generate_rcd_sample(n = 300, seed = 42)

  bs1 <- lingam_rcd_bootstrap(dat$data, n_sampling = 3L, seed = 99L, verbose = FALSE)
  bs2 <- lingam_rcd_bootstrap(dat$data, n_sampling = 3L, seed = 99L, verbose = FALSE)

  expect_equal(bs1$adjacency_matrices, bs2$adjacency_matrices)
})

test_that("parallel lingam_rcd_bootstrap is reproducible and NA-free", {
  skip_on_cran()
  skip_if_not(parallel::detectCores() >= 2L, "requires >= 2 cores")

  dat <- generate_rcd_sample(n = 300, seed = 42)

  bs1 <- lingam_rcd_bootstrap(dat$data, n_sampling = 4L, seed = 7L, verbose = FALSE,
                              parallel = TRUE, n_cores = 2L)
  bs2 <- lingam_rcd_bootstrap(dat$data, n_sampling = 4L, seed = 7L, verbose = FALSE,
                              parallel = TRUE, n_cores = 2L)

  expect_equal(bs1$adjacency_matrices, bs2$adjacency_matrices)
  expect_false(anyNA(bs1$adjacency_matrices))
})

test_that("lingam_rcd does not error with n > SHAPIRO_MAX_N", {
  skip_on_cran()
  set.seed(1)
  n <- 6000L
  x1 <- rnorm(n, 0, 0.5)^3
  x0 <- 1.0 * x1 + rnorm(n, 0, 0.5)^3
  dat <- data.frame(x0 = x0, x1 = x1)

  expect_no_error(res <- lingam_rcd(dat, max_explanatory_num = 1L))
  expect_s3_class(res, "RCDResult")
})
