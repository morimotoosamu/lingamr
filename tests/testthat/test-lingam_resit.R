# RESIT のテスト。
#
# 検定ベースのアルゴリズムなので Python との数値一致ではなく構造回復で判定する
# （リポジトリ規約）。generate_resit_sample() の DGP は seed 1..50 の事前スイープで
# 全 seed が「4辺完全回復・余分な辺ゼロ・順序正解」を返すことを確認済み（n = 300）。
# regressor = "gam" を使うテストは mgcv（Suggests）に依存するため
# skip_if_not_installed("mgcv") を先頭に置く。

# ── 構造 ──────────────────────────────────────────────────────────────────────

test_that("lingam_resit returns a well-formed ResitResult", {
  skip_if_not_installed("mgcv")
  res <- fit_resit_300()

  expect_s3_class(res, "ResitResult")
  expect_named(res, c("adjacency_matrix", "causal_order", "regressor"))

  B <- res$adjacency_matrix
  expect_equal(dim(B), c(4L, 4L))
  expect_equal(dimnames(B), list(paste0("x", 0:3), paste0("x", 0:3)))
  expect_true(all(B %in% c(0, 1)))

  expect_equal(sort(res$causal_order), 1:4)
  expect_equal(res$regressor, "gam")
})

# ── バリデーション ────────────────────────────────────────────────────────────

test_that("lingam_resit validates its inputs", {
  skip_if_not_installed("mgcv")
  X <- resit_sample_300()$data

  expect_error(
    lingam_resit(data.frame(a = letters[1:10], b = letters[1:10])),
    "numeric"
  )
  Xna <- X
  Xna[1, 1] <- NA
  expect_error(lingam_resit(Xna), "missing values")
  expect_error(lingam_resit(X[, 1, drop = FALSE]), "at least 2 variables")
  expect_error(lingam_resit(X[1:5, ]), "at least 6 observations")
  expect_error(lingam_resit(X, alpha = -1), "alpha must be")
  expect_error(lingam_resit(X, alpha = c(0.01, 0.05)), "alpha must be")

  Xconst <- X
  Xconst$x0 <- 1
  expect_error(lingam_resit(Xconst), "constant")

  bad_pk <- matrix(2, 4, 4)
  expect_error(lingam_resit(X, prior_knowledge = bad_pk), "prior_knowledge")
  expect_error(
    lingam_resit(X, prior_knowledge = matrix(-1, 3, 3)),
    "shape of prior knowledge"
  )
})

test_that("lingam_resit validates the regressor argument", {
  X <- resit_sample_300()$data

  expect_error(lingam_resit(X, regressor = "forest"), "'arg'")
  expect_error(lingam_resit(X, regressor = 1L), "regressor must be")
  expect_error(lingam_resit(X, regressor = NULL), "regressor must be")

  # a regressor returning the wrong length is caught with a clear message
  bad_reg <- function(X, y) rep(0, 3)
  expect_error(lingam_resit(X, regressor = bad_reg), "fitted values")
  # ... as is one returning NA
  na_reg <- function(X, y) rep(NA_real_, length(y))
  expect_error(lingam_resit(X, regressor = na_reg), "fitted values")
})

# ── アルゴリズム挙動（回復） ──────────────────────────────────────────────────

test_that("lingam_resit recovers the nonlinear structure with the gam regressor", {
  skip_if_not_installed("mgcv")
  dat <- resit_sample_300()
  res <- fit_resit_300()
  B <- res$adjacency_matrix

  # all 4 true edges, no extra edges (seed-swept DGP, see header comment)
  expect_equal(B["x1", "x0"], 1)
  expect_equal(B["x2", "x0"], 1)
  expect_equal(B["x2", "x1"], 1)
  expect_equal(B["x3", "x2"], 1)
  expect_equal(sum(B), 4)

  # causal order consistent with every true edge (from before to)
  ord <- res$causal_order
  true_edges <- rbind(c(2L, 1L), c(3L, 1L), c(3L, 2L), c(4L, 3L)) # (to, from)
  for (r in seq_len(nrow(true_edges))) {
    expect_lt(
      which(ord == true_edges[r, 2]),
      which(ord == true_edges[r, 1])
    )
  }

  # DAG: adjacency reordered by the causal order is strictly lower-triangular
  reordered <- B[ord, ord]
  expect_true(all(reordered[upper.tri(reordered, diag = TRUE)] == 0))
})

test_that("lingam_resit is deterministic for the same input", {
  skip_if_not_installed("mgcv")
  X <- generate_resit_sample(n = 100, seed = 3)$data

  res1 <- lingam_resit(X)
  res2 <- lingam_resit(X)
  expect_identical(res1, res2)
})

# ── ユーザー関数 regressor ────────────────────────────────────────────────────

test_that("a user-supplied linear regressor runs end to end", {
  dat <- generate_lingam_sample_6(n = 300, seed = 1)
  lin_reg <- function(X, y) stats::lm.fit(cbind(1, X), y)$fitted.values

  res <- lingam_resit(dat$data, regressor = lin_reg)

  expect_s3_class(res, "ResitResult")
  expect_equal(res$regressor, "user function")
  B <- res$adjacency_matrix
  expect_true(all(B %in% c(0, 1)))
  # the output is always a DAG consistent with the estimated order
  ord <- res$causal_order
  reordered <- B[ord, ord]
  expect_true(all(reordered[upper.tri(reordered, diag = TRUE)] == 0))
})

test_that("a user function calling mgcv matches regressor = \"gam\" exactly", {
  skip_if_not_installed("mgcv")
  X <- generate_resit_sample(n = 200, seed = 1)$data

  gam_user <- function(X, y) {
    d <- ncol(X)
    df <- as.data.frame(X)
    names(df) <- paste0("V", seq_len(d))
    df$.y <- y
    f <- stats::as.formula(paste(
      ".y ~", paste(sprintf("s(V%d, k = 10)", seq_len(d)), collapse = " + ")
    ))
    as.vector(stats::fitted(mgcv::gam(f, data = df)))
  }

  res_str <- lingam_resit(X)
  res_fn <- lingam_resit(X, regressor = gam_user)

  expect_identical(res_str$adjacency_matrix, res_fn$adjacency_matrix)
  expect_identical(res_str$causal_order, res_fn$causal_order)
  expect_equal(res_fn$regressor, "user function")
})

# ── prior knowledge ───────────────────────────────────────────────────────────

test_that("prior knowledge can forbid an edge", {
  skip_if_not_installed("mgcv")
  X <- generate_resit_sample(n = 200, seed = 1)$data

  # forbid the true edge x2 -> x3 (pk[to, from] = 0); the parent pruning by
  # prior knowledge guarantees B["x3", "x2"] = 0 regardless of the data
  pk <- matrix(-1, 4, 4)
  pk[4, 3] <- 0
  res <- lingam_resit(X, prior_knowledge = pk)

  expect_equal(res$adjacency_matrix["x3", "x2"], 0)
})

test_that("prior knowledge can force an order against the data", {
  skip_if_not_installed("mgcv")
  X <- generate_resit_sample(n = 200, seed = 1)$data

  # assert a path x3 -> x0 (pk[to = x0, from = x3] = 1): x3 can never be
  # picked as a sink while x0 is unresolved, so x3 must precede x0
  pk <- matrix(-1, 4, 4)
  pk[1, 4] <- 1
  res <- lingam_resit(X, prior_knowledge = pk)

  ord <- res$causal_order
  expect_lt(which(ord == 4), which(ord == 1))
})

test_that("an all-unknown prior knowledge matrix matches no prior knowledge", {
  skip_if_not_installed("mgcv")
  X <- generate_resit_sample(n = 200, seed = 1)$data

  res_none <- lingam_resit(X)
  res_na <- lingam_resit(X, prior_knowledge = matrix(-1, 4, 4))

  expect_identical(res_none$adjacency_matrix, res_na$adjacency_matrix)
  expect_identical(res_none$causal_order, res_na$causal_order)
})

# ── print ─────────────────────────────────────────────────────────────────────

test_that("print.ResitResult shows the order and the 0/1 note", {
  res <- fake_resit_result()

  expect_output(print(res), "RESIT Result")
  expect_output(print(res), "x0 -> x1 -> x2")
  expect_output(print(res), "0/1 edge indicators, not coefficients")
  expect_output(print(res), "Regressor : gam")

  out <- capture.output(v <- withVisible(print(res)))
  expect_false(v$visible)
  expect_identical(v$value, res)
})

# ── bootstrap ─────────────────────────────────────────────────────────────────

test_that("lingam_resit_bootstrap returns a BootstrapResult without total effects", {
  skip_if_not_installed("mgcv")
  X <- generate_resit_sample(n = 150, seed = 2)$data

  bs <- lingam_resit_bootstrap(X, n_sampling = 3L, seed = 42, verbose = FALSE)

  expect_s3_class(bs, "BootstrapResult")
  expect_equal(dim(bs$adjacency_matrices), c(3L, 4L, 4L))
  expect_true(all(bs$adjacency_matrices %in% c(0, 1)))
  expect_null(bs$total_effects)
  expect_error(get_total_causal_effects(bs), "no total effects")

  probs <- get_probabilities(bs)
  expect_true(all(probs >= 0 & probs <= 1))

  # causal_orders is populated, so order-stability queries work
  expect_equal(dim(bs$causal_orders), c(3L, 4L))
  stab <- get_causal_order_stability(bs)
  expect_s3_class(stab, "causal_order_stability")
  expect_equal(stab$n_sampling, 3L)
})

test_that("lingam_resit_bootstrap is reproducible with the same seed", {
  skip_if_not_installed("mgcv")
  X <- generate_resit_sample(n = 150, seed = 2)$data

  bs1 <- lingam_resit_bootstrap(X, n_sampling = 2L, seed = 7, verbose = FALSE)
  bs2 <- lingam_resit_bootstrap(X, n_sampling = 2L, seed = 7, verbose = FALSE)

  expect_equal(bs1$adjacency_matrices, bs2$adjacency_matrices)
  expect_equal(bs1$causal_orders, bs2$causal_orders)
})

test_that("lingam_resit_bootstrap validates its inputs", {
  X <- generate_resit_sample(n = 150, seed = 2)$data

  expect_error(lingam_resit_bootstrap(X, n_sampling = 0), "n_sampling")
  expect_error(
    lingam_resit_bootstrap(X[1:5, ], n_sampling = 2L),
    "at least 6 observations"
  )
  expect_error(
    lingam_resit_bootstrap(X, n_sampling = 2L, alpha = -1),
    "alpha must be"
  )
  expect_error(
    lingam_resit_bootstrap(X, n_sampling = 2L, regressor = "forest"),
    "'arg'"
  )
})

test_that("lingam_resit_bootstrap works in parallel", {
  skip_if_not_installed("mgcv")
  skip_on_cran()
  X <- generate_resit_sample(n = 150, seed = 2)$data

  bs <- lingam_resit_bootstrap(X,
    n_sampling = 2L, seed = 42, verbose = FALSE,
    parallel = TRUE, n_cores = 2L
  )

  expect_s3_class(bs, "BootstrapResult")
  expect_equal(dim(bs$adjacency_matrices), c(2L, 4L, 4L))
  expect_null(bs$total_effects)
})
