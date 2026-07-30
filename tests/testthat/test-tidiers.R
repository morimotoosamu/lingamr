test_that("tidy.LingamResult returns long data.frame of edges", {
  res <- fit_direct_500()
  td  <- tidy(res)

  expect_s3_class(td, "data.frame")
  expect_named(td, c("from", "to", "estimate"))
  expect_type(td$from, "character")
  expect_type(td$to, "character")
  expect_type(td$estimate, "double")
  # number of rows == number of non-zero edges
  expect_equal(nrow(td), sum(abs(res$adjacency_matrix) > 0))
})

test_that("tidy.LingamResult threshold filters edges", {
  res <- fit_direct_500()

  td_all  <- tidy(res, threshold = 0)
  td_high <- tidy(res, threshold = 100)  # exclude all edges

  expect_gte(nrow(td_all), nrow(td_high))
  expect_equal(nrow(td_high), 0L)
})

test_that("tidy.LingamResult from/to follow j -> i convention", {
  # pass the true structure directly to confirm the convention: B["x0","x3"] = 3 is x3 -> x0
  dat  <- sample6_100()
  fake <- fake_lingam_result(dat$true_adjacency)
  td <- tidy(fake)

  edge <- td[td$from == "x3" & td$to == "x0", ]
  expect_equal(nrow(edge), 1L)
  expect_equal(edge$estimate, 3.0)
})

test_that("tidy.LingamResult returns 0-row data.frame for empty graph", {
  fake <- fake_lingam_result(
    matrix(0, 3, 3, dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
  )
  td <- tidy(fake)

  expect_equal(nrow(td), 0L)
  expect_named(td, c("from", "to", "estimate"))
})

test_that("glance.LingamResult returns one-row summary", {
  res <- fit_direct_300()
  g   <- glance(res)

  expect_s3_class(g, "data.frame")
  expect_equal(nrow(g), 1L)
  expect_true(all(c("n_variables", "n_edges", "causal_order") %in% names(g)))
  expect_equal(g$n_variables, 6L)
})

test_that("tidy.BootstrapResult returns direction counts", {
  bs  <- bs_direct_300_15()
  td  <- tidy(bs)

  expect_s3_class(td, "data.frame")
  expect_true(all(c("from", "to", "count", "proportion") %in% names(td)))
})

test_that("tidy/glance for ParceLingamResult keep NA entries", {
  fake <- fake_parce_result()

  td <- tidy(fake)
  expect_named(td, c("from", "to", "estimate"))
  expect_equal(sum(is.na(td$estimate)), 2L)
  expect_equal(td$estimate[td$from == "a" & td$to == "b"], 1.5)

  g <- glance(fake)
  expect_equal(nrow(g), 1L)
  expect_equal(g$n_edges, 1L)
  expect_equal(g$n_na_entries, 2L)
  expect_match(g$causal_order, "(a, c)", fixed = TRUE)
})

test_that("tidy/glance for RCDResult count confounded pairs", {
  B <- matrix(0, 3, 3, dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
  B["c", "b"] <- -0.8
  B["a", "b"] <- NA
  B["b", "a"] <- NA
  fake <- fake_rcd_result(
    ancestors_list = list(integer(0), integer(0), 2L),
    adjacency_matrix = B
  )

  td <- tidy(fake)
  expect_named(td, c("from", "to", "estimate"))
  expect_equal(sum(is.na(td$estimate)), 2L)

  g <- glance(fake)
  expect_equal(g$n_edges, 1L)
  expect_equal(g$n_confounded_pairs, 1L)
  expect_false("causal_order" %in% names(g))
})

test_that("tidy/glance for LiMResult work", {
  fake <- fake_lim_result()

  td <- tidy(fake)
  expect_named(td, c("from", "to", "estimate"))
  expect_equal(nrow(td), 1L)

  g <- glance(fake)
  expect_equal(g$n_edges, 1L)
  expect_equal(g$n_discrete, 1L)
  expect_equal(g$causal_order, "x1 -> x2 -> x3")
})

test_that("tidy/glance for MultiGroupLingamResult stack groups", {
  mg  <- generate_multi_group_sample(n = c(300, 300), seed = 42)
  res <- lingam_multi_group(mg$data_list, reg_method = "ols")

  td <- tidy(res)
  expect_named(td, c("group", "from", "to", "estimate"))
  expect_setequal(unique(td$group), names(res$adjacency_matrices))

  # each group's block equals tidy() of that group's LingamResult
  g1_name <- names(res$adjacency_matrices)[1]
  block <- td[td$group == g1_name, c("from", "to", "estimate")]
  expect_equal(block, tidy(get_group_result(res, 1)), ignore_attr = TRUE)

  gl <- glance(res)
  expect_equal(nrow(gl), 1L)
  expect_equal(gl$n_groups, length(res$adjacency_matrices))
  expect_true(all(c("n_variables", "causal_order") %in% names(gl)))
})

test_that("tidy.MultiGroupBootstrapResult adds a group column", {
  mg <- generate_multi_group_sample(n = c(300, 300), seed = 42)
  bs <- lingam_multi_group_bootstrap(mg$data_list,
    n_sampling = 10L, reg_method = "ols", seed = 42, verbose = FALSE
  )

  td <- tidy(bs)
  expect_true(all(c("group", "from", "to", "count", "proportion") %in% names(td)))
  expect_setequal(unique(td$group), names(bs))
})

test_that("tidy.ImputationBootstrapResult collapses imputations", {
  skip_if_not_installed("mice")

  d <- make_missing_sample6(n = 200)
  bs <- bootstrap_with_imputation(d$X,
    n_sampling = 5L, n_repeats = 2L, seed = 42, verbose = FALSE
  )

  td <- tidy(bs)
  expect_s3_class(td, "data.frame")
  expect_true(all(c("from", "to", "count", "proportion") %in% names(td)))
})

test_that("tidy/glance for ResitResult work", {
  fake <- fake_resit_result()

  td <- tidy(fake)
  expect_named(td, c("from", "to", "estimate"))
  expect_equal(nrow(td), 2L)
  expect_true(all(td$estimate == 1))

  g <- glance(fake)
  expect_equal(
    names(g),
    c("n_variables", "n_edges", "regressor", "causal_order")
  )
  expect_equal(g$n_variables, 3L)
  expect_equal(g$n_edges, 2L)
  expect_equal(g$regressor, "gam")
  expect_equal(g$causal_order, "x0 -> x1 -> x2")
})

test_that("tidy/glance for CAMUVResult keep NA pairs and count them", {
  fake <- fake_camuv_result()

  td <- tidy(fake)
  expect_named(td, c("from", "to", "estimate"))
  # 1 identified edge (a -> b) + the a/c NA pair in both directions
  expect_equal(sum(!is.na(td$estimate)), 1L)
  expect_equal(sum(is.na(td$estimate)), 2L)

  g <- glance(fake)
  expect_equal(
    names(g),
    c("n_variables", "n_edges", "n_confounded_pairs", "regressor")
  )
  expect_equal(g$n_variables, 3L)
  expect_equal(g$n_edges, 1L)
  expect_equal(g$n_confounded_pairs, 1L)
  expect_equal(g$regressor, "gam")
  expect_false("causal_order" %in% names(g))
})
