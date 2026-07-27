test_that("lingam_parce returns ParceLingamResult with correct structure", {
  dat <- generate_parce_sample(n = 500, seed = 42)
  res <- lingam_parce(dat$data, reg_method = "ols")

  expect_s3_class(res, "ParceLingamResult")
  expect_named(res, c("adjacency_matrix", "causal_order", "p_values", "independence"))
  expect_true(is.matrix(res$adjacency_matrix))
  expect_equal(dim(res$adjacency_matrix), c(6L, 6L))
  expect_equal(colnames(res$adjacency_matrix), names(dat$data))
  expect_true(is.list(res$causal_order))
  expect_setequal(unlist(res$causal_order), 1:6)
  # every variable appears exactly once across the causal order
  expect_equal(length(unlist(res$causal_order)), 6L)
})

test_that("lingam_parce validates its inputs", {
  dat <- generate_parce_sample(n = 200, seed = 1)

  expect_error(lingam_parce(matrix(letters[1:4], nrow = 2)), "numeric matrix")
  expect_error(
    lingam_parce(as.data.frame(matrix(1, nrow = 5, ncol = 1))),
    "at least 2 variables"
  )
  expect_error(lingam_parce(dat$data, alpha = -1), "alpha")
  expect_error(lingam_parce(dat$data, independence = "foo"), "should be one of")
  expect_error(lingam_parce(dat$data, ind_corr = -1), "ind_corr")
  expect_error(lingam_parce(dat$data, reg_method = "bad_method"), "should be one of")

  pk_bad <- matrix(-1L, 4, 4)
  expect_error(
    lingam_parce(dat$data, prior_knowledge = pk_bad),
    "shape of prior knowledge"
  )
})

test_that("lingam_parce detects the latent confounder as an unresolved block", {
  dat <- generate_parce_sample(n = 1000, seed = 42)
  res <- lingam_parce(dat$data, reg_method = "ols")

  block <- res$causal_order[[1]]
  expect_gt(length(block), 1L)
  expect_setequal(block, dat$confounded_pair)

  x2 <- dat$confounded_pair[1]
  x3 <- dat$confounded_pair[2]
  expect_true(is.na(res$adjacency_matrix[x2, x3]))
  expect_true(is.na(res$adjacency_matrix[x3, x2]))

  # a downstream true edge (x3 -> x0) is recovered with the correct sign
  x0 <- 1L
  expect_gt(res$adjacency_matrix[x0, x3], 0)

  # x0 -> x5 (true coefficient 0.5) is recovered with the correct sign
  x5 <- 6L
  expect_gt(res$adjacency_matrix[x5, x0], 0)
})

test_that("alpha = 0 disables rejection and returns a full causal order", {
  dat <- generate_parce_sample(n = 500, seed = 42)
  res <- lingam_parce(dat$data, alpha = 0, reg_method = "ols")

  expect_true(all(vapply(res$causal_order, length, integer(1)) == 1L))
  expect_false(anyNA(res$adjacency_matrix))
})

test_that("lingam_parce with unconfounded data returns no block and a plausible order", {
  # seed = 1 gives a full order at this n; the independence test is a
  # statistical test and can occasionally reject a genuine direct edge
  # (e.g. seed = 42 spuriously blocks x0/x4, which are directly connected).
  dat <- generate_lingam_sample_6(n = 1000, seed = 1)
  res <- lingam_parce(dat$data, reg_method = "ols")

  expect_true(all(vapply(res$causal_order, length, integer(1)) == 1L))
  expect_false(anyNA(res$adjacency_matrix))

  order_vec <- unlist(res$causal_order)
  # x3 is the true root; it should be placed before its descendants
  x3 <- which(names(dat$data) == "x3")
  x0 <- which(names(dat$data) == "x0")
  expect_lt(which(order_vec == x3), which(order_vec == x0))
})

test_that("independence = 'fcorr' returns a valid ParceLingamResult", {
  dat <- generate_parce_sample(n = 500, seed = 42)
  res <- lingam_parce(dat$data, independence = "fcorr", reg_method = "ols")

  expect_s3_class(res, "ParceLingamResult")
  expect_equal(res$independence, "fcorr")
  expect_setequal(unlist(res$causal_order), 1:6)
})

test_that("estimate_total_effect_parce warns and returns NA for a confounded variable", {
  dat <- generate_parce_sample(n = 1000, seed = 42)
  res <- lingam_parce(dat$data, reg_method = "ols")

  x2 <- dat$confounded_pair[1]
  expect_warning(
    eff <- estimate_total_effect_parce(dat$data, res, from_index = x2, to_index = 1),
    "unresolved"
  )
  expect_true(is.na(eff))

  # a well-identified pair (x0 -> x5) returns a numeric estimate
  x0 <- 1L
  x5 <- 6L
  eff2 <- estimate_total_effect_parce(dat$data, res, from_index = x0, to_index = x5)
  expect_true(is.numeric(eff2) && !is.na(eff2))
})

test_that("estimate_total_effect_parce works with method = 'ols'", {
  dat <- generate_parce_sample(n = 1000, seed = 42)
  res <- lingam_parce(dat$data, reg_method = "ols")

  eff <- estimate_total_effect_parce(dat$data, res,
                                     from_index = 1L, to_index = 6L,
                                     method = "ols")
  expect_true(is.numeric(eff) && !is.na(eff))
})

test_that("get_error_independence_p_values_parce returns NA for confounded pairs and valid p-values elsewhere", {
  dat <- generate_parce_sample(n = 1000, seed = 42)
  res <- lingam_parce(dat$data, reg_method = "ols")

  pvals <- get_error_independence_p_values_parce(dat$data, res)
  expect_equal(dim(pvals), c(6L, 6L))

  x2 <- dat$confounded_pair[1]
  x3 <- dat$confounded_pair[2]
  expect_true(is.na(pvals[x2, x3]))

  finite_vals <- pvals[!is.na(pvals)]
  expect_true(all(finite_vals >= 0 & finite_vals <= 1))
})

test_that("print.ParceLingamResult runs without error and shows blocks", {
  dat <- generate_parce_sample(n = 500, seed = 42)
  res <- lingam_parce(dat$data, reg_method = "ols")

  expect_output(print(res), "Bottom-Up ParceLiNGAM")
  expect_output(print(res), "\\(")
})

test_that("lingam_parce with prior_knowledge runs without error", {
  dat <- generate_lingam_sample_6(n = 500, seed = 1)
  pk <- make_prior_knowledge(6, exogenous_variables = 4, labels = names(dat$data))
  res <- lingam_parce(dat$data, prior_knowledge = pk, reg_method = "ols")

  expect_s3_class(res, "ParceLingamResult")
})

test_that("lingam_parce_bootstrap returns a usable BootstrapResult", {
  dat <- generate_parce_sample(n = 500, seed = 42)

  bs <- lingam_parce_bootstrap(dat$data,
    n_sampling = 8L, reg_method = "ols", seed = 42, verbose = FALSE
  )

  expect_s3_class(bs, "BootstrapResult")
  expect_null(bs$causal_orders)
  expect_false(anyNA(bs$adjacency_matrices))
  expect_false(anyNA(bs$total_effects))

  p <- get_probabilities(bs)
  expect_equal(dim(p), c(6L, 6L))
  expect_false(anyNA(p))

  dc <- get_causal_direction_counts(bs)
  expect_s3_class(dc, "data.frame")

  dag <- get_directed_acyclic_graph_counts(bs)
  expect_named(dag, c("dag", "count"))

  te <- get_total_causal_effects(bs)
  expect_s3_class(te, "data.frame")
})

test_that("lingam_parce_bootstrap is reproducible with the same seed", {
  dat <- generate_parce_sample(n = 400, seed = 1)

  bs1 <- lingam_parce_bootstrap(dat$data, n_sampling = 5L, reg_method = "ols",
                                 seed = 99L, verbose = FALSE)
  bs2 <- lingam_parce_bootstrap(dat$data, n_sampling = 5L, reg_method = "ols",
                                 seed = 99L, verbose = FALSE)

  expect_equal(bs1$adjacency_matrices, bs2$adjacency_matrices)
})

test_that("parallel lingam_parce_bootstrap is reproducible and NA-free", {
  skip_on_cran()
  skip_if_not(parallel::detectCores() >= 2L, "requires >= 2 cores")

  dat <- generate_parce_sample(n = 400, seed = 1)

  bs1 <- lingam_parce_bootstrap(dat$data, n_sampling = 6L, reg_method = "ols",
                                 seed = 7L, verbose = FALSE,
                                 parallel = TRUE, n_cores = 2L)
  bs2 <- lingam_parce_bootstrap(dat$data, n_sampling = 6L, reg_method = "ols",
                                 seed = 7L, verbose = FALSE,
                                 parallel = TRUE, n_cores = 2L)

  expect_equal(bs1$adjacency_matrices, bs2$adjacency_matrices)
  expect_equal(bs1$total_effects, bs2$total_effects)
  expect_false(anyNA(bs1$adjacency_matrices))
})
