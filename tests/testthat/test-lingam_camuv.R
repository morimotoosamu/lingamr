# CAM-UV のテスト。
#
# 検定ベースのアルゴリズムなので Python との数値一致ではなく構造回復で判定する
# （リポジトリ規約。上流は pygam の LinearGAM、本移植は mgcv::gam なので
# 数値一致はそもそも定義できない）。generate_camuv_sample() の DGP は
# seed 1..10 の事前スイープで全 seed が「3辺完全回復・NAペア2組・余分な辺ゼロ」
# を返すことを確認済み（n = 500）。
# regressor = "gam" を使うテストは mgcv（Suggests）に依存するため
# skip_if_not_installed("mgcv") を先頭に置く。

# ── 構造 ──────────────────────────────────────────────────────────────────────

test_that("lingam_camuv returns a well-formed CAMUVResult", {
  skip_if_not_installed("mgcv")
  dat <- camuv_sample_500()
  res <- fit_camuv_500()

  expect_s3_class(res, "CAMUVResult")
  expect_named(
    res,
    c("adjacency_matrix", "parents_list", "confounded_pairs", "regressor")
  )

  B <- res$adjacency_matrix
  expect_equal(dim(B), c(6L, 6L))
  expect_equal(dimnames(B), list(names(dat$data), names(dat$data)))
  expect_true(all(B %in% c(0, 1) | is.na(B)))

  expect_length(res$parents_list, 6L)
  expect_equal(names(res$parents_list), names(dat$data))

  expect_true(is.matrix(res$confounded_pairs))
  expect_equal(colnames(res$confounded_pairs), c("var1", "var2"))
  expect_equal(res$regressor, "gam")

  # parents_list / confounded_pairs は adjacency_matrix と整合する
  for (i in seq_len(6L)) {
    expect_equal(which(!is.na(B[i, ]) & B[i, ] == 1), res$parents_list[[i]],
      ignore_attr = TRUE
    )
  }
  na_idx <- which(is.na(B), arr.ind = TRUE)
  na_pairs <- unique(t(apply(na_idx, 1, sort)))
  expect_equal(
    na_pairs[order(na_pairs[, 1]), , drop = FALSE],
    unname(res$confounded_pairs[order(res$confounded_pairs[, 1]), ,
      drop = FALSE
    ]),
    ignore_attr = TRUE
  )
})

# ── バリデーション ────────────────────────────────────────────────────────────

test_that("lingam_camuv validates its inputs", {
  X <- camuv_sample_500()$data

  expect_error(
    lingam_camuv(data.frame(a = letters[1:10], b = letters[1:10])),
    "numeric"
  )
  Xna <- X
  Xna[1, 1] <- NA
  expect_error(lingam_camuv(Xna), "missing values")
  expect_error(lingam_camuv(X[, 1, drop = FALSE]), "at least 2 variables")
  expect_error(lingam_camuv(X[1:5, ]), "at least 6 observations")
  expect_error(lingam_camuv(X, alpha = -1), "alpha must be")
  expect_error(lingam_camuv(X, alpha = c(0.01, 0.05)), "alpha must be")
  expect_error(
    lingam_camuv(X, num_explanatory_vals = 0L),
    "num_explanatory_vals"
  )
  expect_error(lingam_camuv(X, independence = "lingam"), "should be one of")
  expect_error(lingam_camuv(X, ind_corr = -1), "ind_corr must be")
  expect_error(
    lingam_camuv(X, independence = "fcorr", num_explanatory_vals = 3L),
    "fcorr"
  )

  Xconst <- X
  Xconst$x0 <- 1
  expect_error(lingam_camuv(Xconst), "constant")
})

test_that("lingam_camuv validates prior_knowledge pairs", {
  X <- camuv_sample_500()$data

  expect_error(lingam_camuv(X, prior_knowledge = "x0"), "prior_knowledge")
  expect_error(lingam_camuv(X, prior_knowledge = list()), "prior_knowledge")
  expect_error(
    lingam_camuv(X, prior_knowledge = list(c(1, 2, 3))),
    "prior_knowledge"
  )
  expect_error(
    lingam_camuv(X, prior_knowledge = matrix(1:9, 3, 3)),
    "prior_knowledge"
  )
  expect_error(
    lingam_camuv(X, prior_knowledge = list(c(0, 2))),
    "between 1 and 6"
  )
  expect_error(
    lingam_camuv(X, prior_knowledge = list(c(1, 7))),
    "between 1 and 6"
  )
  expect_error(
    lingam_camuv(X, prior_knowledge = list(c(2, 2))),
    "distinct"
  )
})

test_that("lingam_camuv validates the regressor argument", {
  X <- camuv_sample_500()$data

  expect_error(lingam_camuv(X, regressor = "forest"), "'arg'")
  expect_error(lingam_camuv(X, regressor = 1L), "regressor must be")

  # a regressor returning the wrong length is caught with a clear message
  bad_reg <- function(X, y) rep(0, 3)
  expect_error(lingam_camuv(X, regressor = bad_reg), "fitted values")
})

# ── アルゴリズム挙動（回復） ──────────────────────────────────────────────────

test_that("lingam_camuv recovers direct edges and unobserved-path pairs", {
  skip_if_not_installed("mgcv")
  dat <- camuv_sample_500()
  res <- fit_camuv_500()
  B <- res$adjacency_matrix

  # 3 true direct edges: x0 -> x1, x0 -> x3, x2 -> x4
  expect_equal(B["x1", "x0"], 1)
  expect_equal(B["x3", "x0"], 1)
  expect_equal(B["x4", "x2"], 1)
  expect_equal(sum(B > 0, na.rm = TRUE), 3)

  # UCP pair {x2, x5} and UBP pair {x3, x4} are NA in both directions
  expect_true(is.na(B["x2", "x5"]) && is.na(B["x5", "x2"]))
  expect_true(is.na(B["x3", "x4"]) && is.na(B["x4", "x3"]))
  expect_equal(sum(is.na(B)), 4)

  # ... and reported as confounded pairs, matching the generator's oracle
  expect_equal(
    unname(res$confounded_pairs[order(res$confounded_pairs[, 1]), ]),
    unname(dat$confounded_pairs[order(dat$confounded_pairs[, 1]), ])
  )
})

test_that("prior_knowledge forbids the specified causal direction", {
  skip_if_not_installed("mgcv")
  skip_on_cran()
  dat <- camuv_sample_500()

  # forbid x0 (col 1) from being a cause of x1 (col 2)
  res <- lingam_camuv(dat$data, prior_knowledge = list(c(1L, 2L)))
  B <- res$adjacency_matrix

  expect_false(!is.na(B["x1", "x0"]) && B["x1", "x0"] == 1)
  expect_false(2L %in% res$parents_list[["x0"]] &&
                 1L %in% res$parents_list[["x1"]])
})

test_that("independence = 'fcorr' returns a valid CAMUVResult", {
  skip_if_not_installed("mgcv")
  skip_on_cran()
  dat <- generate_camuv_sample(n = 200, seed = 1)

  res <- lingam_camuv(dat$data, independence = "fcorr")

  expect_s3_class(res, "CAMUVResult")
  expect_equal(dim(res$adjacency_matrix), c(6L, 6L))
  expect_true(all(res$adjacency_matrix %in% c(0, 1) |
                    is.na(res$adjacency_matrix)))
})

test_that("lingam_camuv is deterministic for the same input", {
  skip_if_not_installed("mgcv")
  skip_on_cran()
  X <- generate_camuv_sample(n = 150, seed = 3)$data

  res1 <- lingam_camuv(X)
  res2 <- lingam_camuv(X)
  expect_identical(res1, res2)
})

# ── ユーザー関数 regressor ────────────────────────────────────────────────────

test_that("a user-supplied linear regressor runs end to end", {
  dat <- generate_camuv_sample(n = 150, seed = 1)
  lin_reg <- function(X, y) stats::lm.fit(cbind(1, X), y)$fitted.values

  res <- lingam_camuv(dat$data, regressor = lin_reg)

  expect_s3_class(res, "CAMUVResult")
  expect_equal(res$regressor, "user function")
})

# ── print ─────────────────────────────────────────────────────────────────────

test_that("print.CAMUVResult prints parent sets and confounded pairs", {
  fake <- fake_camuv_result()

  out <- capture.output(print(fake))
  expect_true(any(grepl("CAM-UV Result", out)))
  expect_true(any(grepl("Regressor : gam", out)))
  expect_true(any(grepl("P(b) = {a}", out, fixed = TRUE)))
  expect_true(any(grepl("a -- c", out)))
  expect_true(any(grepl("row = to, col = from", out)))

  # empty confounded pairs print the "none" line instead
  fake0 <- fake_camuv_result(
    adjacency_matrix = matrix(0, 3, 3,
      dimnames = list(c("a", "b", "c"), c("a", "b", "c"))
    ),
    parents_list = list(a = integer(0), b = integer(0), c = integer(0)),
    confounded_pairs = matrix(integer(0), 0, 2,
      dimnames = list(NULL, c("var1", "var2"))
    )
  )
  out0 <- capture.output(print(fake0))
  expect_true(any(grepl("No pairs with an unobserved", out0)))

  # print() returns its argument invisibly (capture.output absorbs the text)
  ret <- NULL
  invisible(capture.output(ret <- withVisible(print(fake))))
  expect_false(ret$visible)
  expect_identical(ret$value, fake)
})
