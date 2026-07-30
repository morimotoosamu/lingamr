# 複数のテストファイルで共有するフィクスチャ集。
#
# 方針:
# - generate_*() は内部で set.seed() するため呼び出し時点の RNG 状態に依存せず、
#   フィット関数（lingam_direct / lingam_rcd / lingam_parce / lingam_multi_group /
#   lingam_var / lingam_high_dim）と seed 固定の bootstrap も決定的。
#   そこで高コストなデータ生成・フィットを memo() でメモ化し、ファイル横断で共有する。
# - lingam_lim だけは ambient RNG を消費する（optim 初期値の runif）ため、
#   test-lingam_lim.R はこのフィクスチャを一切使わない（除外対象）。
#   その他の除外: test-lingam_var_snapshot.R（golden-value）、test-hsic.R、
#   test-generate_samples.R、test-prior_knowledge.R、再現性テスト（同 seed 2回比較）。
# - フィクスチャの引数（n / seed / reg_method / n_sampling など）を変更すると
#   golden 的な期待値を持つテストが壊れるので変更しないこと。
#   条件が異なるケースは各テストファイル側で直接呼び出す。
# - 返り値はすべて通常の list / S3 オブジェクト（copy-on-modify）なので、
#   テスト内でローカルに書き換えても他テストへ波及しない。環境や R6 は返さない。

# 初回呼び出し時のみ fn() を評価し、以後はキャッシュを返すメモ化ラッパ
memo <- function(fn) {
  val <- NULL
  computed <- FALSE
  function() {
    if (!computed) {
      val <<- fn()
      computed <<- TRUE
    }
    val
  }
}

# --- データフィクスチャ -------------------------------------------------------

sample6_100      <- memo(function() generate_lingam_sample_6(n = 100, seed = 1))
sample6_200      <- memo(function() generate_lingam_sample_6(n = 200, seed = 1))
sample6_300      <- memo(function() generate_lingam_sample_6(n = 300, seed = 1))
sample6_500_s42  <- memo(function() generate_lingam_sample_6(n = 500, seed = 42))
sample6_800_s42  <- memo(function() generate_lingam_sample_6(n = 800, seed = 42))
sample6_2000_s42 <- memo(function() generate_lingam_sample_6(n = 2000, seed = 42))

rcd_sample_300    <- memo(function() generate_rcd_sample(n = 300, seed = 42))
parce_sample_500  <- memo(function() generate_parce_sample(n = 500, seed = 42))
parce_sample_1000 <- memo(function() generate_parce_sample(n = 1000, seed = 42))
mg_sample_300     <- memo(function() generate_multi_group_sample(n = c(300, 300), seed = 1))

vars_1000_s42 <- memo(function() generate_varlingam_sample(n = 1000, seed = 42))
vars_1500_s42 <- memo(function() generate_varlingam_sample(n = 1500, seed = 42))
vars_2000_s42 <- memo(function() generate_varlingam_sample(n = 2000, seed = 42))

# x5 列に NA を混ぜた sample6（多重代入テスト用）。自前で set.seed するので決定的。
make_missing_sample6 <- function(n = 300, seed = 1, na_frac = 0.1) {
  sample6 <- generate_lingam_sample_6(n = n, seed = seed)
  X <- sample6$data
  set.seed(seed)
  na_idx <- sample.int(n, size = round(na_frac * n))
  X$x5[na_idx] <- NA
  list(X = X, na_idx = na_idx)
}

# --- フィットフィクスチャ（既定条件のみキャッシュ） ---------------------------

fit_direct_200  <- memo(function() lingam_direct(sample6_200()$data, reg_method = "ols"))
fit_direct_300  <- memo(function() lingam_direct(sample6_300()$data, reg_method = "ols"))
fit_direct_500  <- memo(function() lingam_direct(sample6_500_s42()$data, reg_method = "ols"))
fit_direct_2000 <- memo(function() lingam_direct(sample6_2000_s42()$data, reg_method = "ols"))

fit_rcd_300    <- memo(function() lingam_rcd(rcd_sample_300()$data))
fit_parce_500  <- memo(function() lingam_parce(parce_sample_500()$data, reg_method = "ols"))
fit_parce_1000 <- memo(function() lingam_parce(parce_sample_1000()$data, reg_method = "ols"))
fit_mg_300     <- memo(function() lingam_multi_group(mg_sample_300()$data_list, reg_method = "ols"))

# VAR-LiNGAM の共通引数ラッパ（非キャッシュ。小さい n の一回きりフィット用）
fit_var_default <- function(data) {
  lingam_var(data, lags = 1, reg_method = "ols", criterion = NULL, prune = FALSE)
}

fit_var_1000 <- memo(function() fit_var_default(vars_1000_s42()$data))
fit_var_1500 <- memo(function() fit_var_default(vars_1500_s42()$data))
fit_var_2000 <- memo(function() fit_var_default(vars_2000_s42()$data))

# --- bootstrap フィクスチャ（seed 固定なので内部 set.seed により決定的） ------

bs_direct_200_10 <- memo(function() {
  lingam_direct_bootstrap(sample6_200()$data, n_sampling = 10L, reg_method = "ols", seed = 1L)
})
bs_direct_300_15 <- memo(function() {
  lingam_direct_bootstrap(sample6_300()$data, n_sampling = 15L, reg_method = "ols", seed = 42L)
})
bs_direct_800_20 <- memo(function() {
  lingam_direct_bootstrap(sample6_800_s42()$data, n_sampling = 20L, reg_method = "ols", seed = 1L)
})

# --- fake 結果オブジェクト（tidiers / autoplot の手組みオブジェクトを共通化） --

# 3 変数の隣接行列: b <- a (1.5)、a/c ペアは潜在交絡を表す NA
fake_na_adjacency_3 <- function() {
  B <- matrix(0, 3, 3, dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
  B["b", "a"] <- 1.5
  B["a", "c"] <- NA
  B["c", "a"] <- NA
  B
}

fake_lingam_result <- function(adjacency_matrix,
                               causal_order = seq_len(ncol(adjacency_matrix))) {
  structure(
    list(adjacency_matrix = adjacency_matrix, causal_order = causal_order),
    class = "LingamResult"
  )
}

fake_parce_result <- function(adjacency_matrix = fake_na_adjacency_3(),
                              causal_order = list(2L, c(1L, 3L)),
                              p_values = matrix(0, 3, 3),
                              independence = "hsic") {
  structure(
    list(adjacency_matrix = adjacency_matrix, causal_order = causal_order,
         p_values = p_values, independence = independence),
    class = "ParceLingamResult"
  )
}

fake_rcd_result <- function(ancestors_list,
                            adjacency_matrix = fake_na_adjacency_3()) {
  structure(
    list(adjacency_matrix = adjacency_matrix, ancestors_list = ancestors_list),
    class = "RCDResult"
  )
}

fake_lim_result <- function() {
  B <- matrix(0, 3, 3, dimnames = list(paste0("x", 1:3), paste0("x", 1:3)))
  B["x2", "x1"] <- 1.2
  structure(
    list(adjacency_matrix = B, causal_order = 1:3,
         is_continuous = c(TRUE, FALSE, TRUE)),
    class = "LiMResult"
  )
}
