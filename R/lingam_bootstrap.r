# =============================================================================
# Bootstrap for Direct LiNGAM - R Implementation
# Based on the Python implementation from the LiNGAM Project
# https://sites.google.com/view/sshimizu06/lingam
# https://github.com/cdt15/lingam
# =============================================================================

# =============================================================================
# ユーティリティ関数
# =============================================================================

#' DAG中の全パスを探索（from_index → to_index）
#' @param adjacency_matrix 隣接行列
#' @param from_index 始点インデックス (1-based)
#' @param to_index 終点インデックス (1-based)
#' @param min_causal_effect 最小因果効果の閾値
#' @return list(paths, effects)
#' @keywords internal
find_all_paths <- function(adjacency_matrix, from_index, to_index, min_causal_effect = 0.0) {
  B <- adjacency_matrix
  B[is.na(B)] <- 0

  # 有効なエッジのみ残す
  B[abs(B) <= min_causal_effect] <- 0

  p <- ncol(B)
  paths <- list()
  effects <- c()

  # DFS で全パスを探索
  # B[i, j] != 0 は j → i を意味する
  dfs <- function(current, target, visited, path, effect) {
    if (current == target && length(path) > 1) {
      paths[[length(paths) + 1]] <<- path
      effects[length(effects) + 1] <<- effect
      return()
    }
    for (next_node in 1:p) {
      if (next_node %in% visited) next
      if (B[next_node, current] == 0) next
      dfs(
        next_node, target, c(visited, next_node),
        c(path, next_node), effect * B[next_node, current]
      )
    }
  }

  dfs(from_index, to_index, c(from_index), c(from_index), 1.0)

  return(list(paths = paths, effects = effects))
}


#' 総合因果効果を計算
#' @param adjacency_matrix 隣接行列
#' @param from_index 原因変数のインデックス (1-based)
#' @param to_index 結果変数のインデックス (1-based)
#' @return 総合因果効果
#' @keywords internal
calculate_total_effect <- function(adjacency_matrix, from_index, to_index) {
  result <- find_all_paths(adjacency_matrix, from_index, to_index, min_causal_effect = 0.0)
  if (length(result$effects) == 0) {
    return(0.0)
  }
  return(sum(result$effects))
}


# =============================================================================
# ブートストラップ実行関数
# =============================================================================

#' Direct LiNGAM のブートストラップ
#'
#' @param X 数値行列 (n_samples x n_features)
#' @param n_sampling ブートストラップの反復回数
#' @param prior_knowledge 事前知識行列 (NULL可)
#' @param apply_prior_knowledge_softly 事前知識のソフト適用 (logical)
#' @param measure 独立性の評価尺度 ("pwling" or "kernel")
#' @param reg_method 回帰手法 ("ols", "lasso", "adaptive_lasso")
#' @param lambda ラムダ選択 ("lambda.min", "lambda.1se", "AIC", "BIC","oracle")
#' @param seed 乱数シード (NULL可)
#' @param verbose 進捗を表示するか (logical)
#' @return BootstrapResult (list)
#' @export
#' @examples
#' LiNGAM_sample_1000 <- generate_lingam_sample_6()
#'
#' # Fast example with OLS
#' bs <- lingam_direct_bootstrap(LiNGAM_sample_1000$data,
#'   n_sampling = 10L,
#'   reg_method = "ols",
#'   seed = 42
#' )
#' get_probabilities(bs)
#'
#' \donttest{
#' # With LASSO (requires glmnet)
#' bs_lasso <- lingam_direct_bootstrap(LiNGAM_sample_1000$data,
#'   n_sampling = 30L,
#'   seed = 42
#' )
#' }
lingam_direct_bootstrap <- function(X,
                             n_sampling,
                             prior_knowledge = NULL,
                             apply_prior_knowledge_softly = FALSE,
                             measure = "pwling",
                             reg_method = "adaptive_lasso",
                             lambda = "BIC",
                             seed = NULL,
                             verbose = TRUE) {
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be a numeric matrix.")
  n_sampling <- as.integer(n_sampling)
  if (n_sampling <= 0) stop("n_sampling must be > 0.")
  if (!is.null(seed)) set.seed(seed)
  n_samples <- nrow(X)
  n_features <- ncol(X)
  adjacency_matrices <- array(0, dim = c(n_sampling, n_features, n_features))
  total_effects <- array(0, dim = c(n_sampling, n_features, n_features))
  resampled_indices <- vector("list", n_sampling)
  if (verbose) {
    message(sprintf("Bootstrap: %d iterations, method=%s", n_sampling, reg_method))
    t_start <- proc.time()
  }
  for (i in seq_len(n_sampling)) {
    if (verbose && (i %% 10 == 0 || i == 1)) {
      message(sprintf("  iteration %d / %d", i, n_sampling))
    }
    idx <- sample(n_samples, replace = TRUE)
    resampled_X <- X[idx, , drop = FALSE]
    resampled_indices[[i]] <- idx
    result <- lingam_direct(
      resampled_X,
      prior_knowledge = prior_knowledge,
      apply_prior_knowledge_softly = apply_prior_knowledge_softly,
      measure = measure,
      reg_method = reg_method,
      lambda = lambda
    )
    adjacency_matrices[i, , ] <- result$adjacency_matrix
    total_effects[i, , ] <- estimate_all_total_effects(
      resampled_X, result,
      method = reg_method,
      lambda = lambda
    )
  }
  if (verbose) {
    elapsed <- (proc.time() - t_start)["elapsed"]
    message(sprintf("Completed in %.1f seconds.", elapsed))
  }
  create_bootstrap_result(adjacency_matrices, total_effects, resampled_indices)
}


# =============================================================================
# BootstrapResult オブジェクト
# =============================================================================

#' BootstrapResult を作成
#' @param adjacency_matrices array (n_sampling x n_features x n_features)
#' @param total_effects array (n_sampling x n_features x n_features)
#' @param resampled_indices list of index vectors
#' @return BootstrapResult (list with class attribute)
#' @keywords internal
create_bootstrap_result <- function(adjacency_matrices, total_effects, resampled_indices = NULL) {
  obj <- list(
    adjacency_matrices = adjacency_matrices,
    total_effects      = total_effects,
    resampled_indices  = resampled_indices
  )
  class(obj) <- "BootstrapResult"
  return(obj)
}


#' BootstrapResult の内容を表示
#'
#' @param x BootstrapResult オブジェクト
#' @param ... 追加の引数 (S3メソッド互換用)
#' @method print BootstrapResult
#' @export
#' @examples
#' LiNGAM_sample_1000 <- generate_lingam_sample_6()
#'
#' bs_model <- LiNGAM_sample_1000$data |>
#'   lingam_direct_bootstrap(n_sampling = 30L, seed = 42)
#'
#' bs_model |>
#'   print()
print.BootstrapResult <- function(x, ...) {
  n_sampling <- dim(x$adjacency_matrices)[1]
  n_features <- dim(x$adjacency_matrices)[2]
  cat(sprintf("BootstrapResult: %d samplings, %d features\n", n_sampling, n_features))
}


# =============================================================================
# BootstrapResult のメソッド群
# =============================================================================

#' 因果方向のカウント・割合・因果効果を取得
#'
#' @param result BootstrapResult オブジェクト
#' @param n_directions 上位何件を返すか (NULL = 全て)
#' @param min_causal_effect 因果効果の最小閾値 (NULL = 0)
#' @param split_by_causal_effect_sign 因果効果の符号で分割するか
#' @param labels 変数名ベクトル (NULL可。指定するとfrom_name, to_name列を追加)
#' @return A data frame containing the following columns:
#' * `from`, `to`: 1-based indices of the causal (from) and effect (to) variables.
#' * `count`: Number of bootstrap samples in which this specific causal direction was identified.
#' * `proportion`: The frequency of the direction's occurrence (count / n_sampling), representing its bootstrap probability.
#' * `mean_effect`: The average value of the estimated causal effects across samples where this direction was identified.
#' * `median_effect`: The median value of the estimated causal effects, providing a robust estimate of the effect size.
#' * `sd_effect`: The standard deviation of the causal effect estimates, indicating the stability of the effect size.
#' * `ci_lower`, `ci_upper`: The lower (2.5%) and upper (97.5%) bounds of the bootstrap confidence interval for the causal effect.
#' * `sign` (optional): The sign of the causal effect (1 for positive, -1 for negative), included if `split_by_causal_effect_sign = TRUE`.
#' * `from_name`, `to_name` (optional): Character labels for the variables, included if `labels` were provided.
#' @importFrom stats sd median quantile
#' @export
#' @examples
#' LiNGAM_sample_1000 <- generate_lingam_sample_6()
#'
#' bs_model <- LiNGAM_sample_1000$data |>
#'   lingam_direct_bootstrap(n_sampling = 30L, seed = 42)
#'
#' bs_model |>
#'   get_causal_direction_counts(labels = names(LiNGAM_sample_1000$data))
get_causal_direction_counts <- function(result,
                                        n_directions = NULL,
                                        min_causal_effect = NULL,
                                        split_by_causal_effect_sign = FALSE,
                                        labels = NULL) {
  stopifnot(inherits(result, "BootstrapResult"))

  if (is.null(min_causal_effect)) min_causal_effect <- 0.0
  if (min_causal_effect < 0) stop("min_causal_effect must be >= 0.")

  am_array <- result$adjacency_matrices
  am_array[is.na(am_array)] <- 0
  n_sampling <- dim(am_array)[1]
  n_features <- dim(am_array)[2]

  # 全ブートストラップから方向と効果量を収集
  directions_list <- list()
  for (s in seq_len(n_sampling)) {
    am <- am_array[s, , ]
    idx <- which(abs(am) > min_causal_effect, arr.ind = TRUE)
    if (nrow(idx) == 0) next

    effects <- sapply(seq_len(nrow(idx)), function(k) am[idx[k, 1], idx[k, 2]])

    if (split_by_causal_effect_sign) {
      directions_list[[length(directions_list) + 1]] <- data.frame(
        to     = idx[, 1],
        from   = idx[, 2],
        sign   = as.integer(sign(effects)),
        effect = effects
      )
    } else {
      directions_list[[length(directions_list) + 1]] <- data.frame(
        to     = idx[, 1],
        from   = idx[, 2],
        effect = effects
      )
    }
  }

  # 空の場合
  if (length(directions_list) == 0) {
    empty <- data.frame(
      from = integer(0), to = integer(0),
      count = integer(0), proportion = numeric(0),
      mean_effect = numeric(0), median_effect = numeric(0),
      sd_effect = numeric(0), ci_lower = numeric(0), ci_upper = numeric(0)
    )
    if (split_by_causal_effect_sign) empty$sign <- integer(0)
    if (!is.null(labels)) {
      empty$from_name <- character(0)
      empty$to_name <- character(0)
    }
    return(empty)
  }

  directions <- do.call(rbind, directions_list)

  # グループキーの構築
  if (split_by_causal_effect_sign) {
    group_key <- paste(directions$from, directions$to, directions$sign, sep = "_")
  } else {
    group_key <- paste(directions$from, directions$to, sep = "_")
  }

  # グループごとに集計
  unique_keys <- unique(group_key)
  results_list <- lapply(unique_keys, function(key) {
    mask <- group_key == key
    subset_df <- directions[mask, ]
    effects <- subset_df$effect

    row <- data.frame(
      from = subset_df$from[1],
      to = subset_df$to[1],
      count = length(effects),
      proportion = length(effects) / n_sampling,
      mean_effect = base::mean(effects),
      median_effect = stats::median(effects),
      sd_effect = if (length(effects) > 1) stats::sd(effects) else 0,
      ci_lower = stats::quantile(effects, 0.025, names = FALSE),
      ci_upper = stats::quantile(effects, 0.975, names = FALSE),
      stringsAsFactors = FALSE
    )

    if (split_by_causal_effect_sign) {
      row$sign <- subset_df$sign[1]
    }

    row
  })

  agg <- do.call(rbind, results_list)

  # 降順ソート
  agg <- agg[order(-agg$count), ]
  rownames(agg) <- NULL

  # 変数名の付与
  if (!is.null(labels)) {
    agg$from_name <- labels[agg$from]
    agg$to_name <- labels[agg$to]
  }

  # 上位 n_directions 件
  if (!is.null(n_directions)) {
    n_directions <- min(n_directions, nrow(agg))
    agg <- agg[seq_len(n_directions), ]
  }

  return(agg)
}


#' DAG カウントを取得
#'
#' @param result BootstrapResult オブジェクト
#' @param n_dags 上位何件を返すか (NULL = 全て)
#' @param min_causal_effect 因果効果の最小閾値 (NULL = 0)
#' @param split_by_causal_effect_sign 因果効果の符号で分割するか
#' @return list(dag = list of data.frames, count = integer vector)
#' @export
#' @examples
#' LiNGAM_sample_1000 <- generate_lingam_sample_6()
#'
#' bs_model <- LiNGAM_sample_1000$data |>
#'   lingam_direct_bootstrap(n_sampling = 30L, seed = 42)
#'
#' bs_model |>
#'   get_directed_acyclic_graph_counts()
get_directed_acyclic_graph_counts <- function(result,
                                              n_dags = NULL,
                                              min_causal_effect = NULL,
                                              split_by_causal_effect_sign = FALSE) {
  stopifnot(inherits(result, "BootstrapResult"))

  if (is.null(min_causal_effect)) min_causal_effect <- 0.0
  if (min_causal_effect < 0) stop("min_causal_effect must be >= 0.")

  am_array <- result$adjacency_matrices
  am_array[is.na(am_array)] <- 0
  n_sampling <- dim(am_array)[1]

  # 各 DAG を文字列キーに変換
  dag_keys <- character(n_sampling)
  dag_list <- vector("list", n_sampling)

  if (split_by_causal_effect_sign) {
    sign_array <- sign(am_array)
    sign_array[abs(am_array) <= min_causal_effect] <- 0
    dag_keys <- apply(sign_array, MARGIN = 1, FUN = paste, collapse = ",")
    dag_list <- lapply(seq_len(n_sampling), function(s) sign_array[s, , ])
  } else {
    bin_array <- abs(am_array) > min_causal_effect
    mode(bin_array) <- "integer" # 文字列化のために integer に変換
    dag_keys <- apply(bin_array, MARGIN = 1, FUN = paste, collapse = ",")
    dag_list <- lapply(seq_len(n_sampling), function(s) bin_array[s, , ])
  }

  # カウント
  tbl <- sort(table(dag_keys), decreasing = TRUE)
  if (!is.null(n_dags)) {
    n_dags <- min(n_dags, length(tbl))
    tbl <- tbl[seq_len(n_dags)]
  }

  # 結果の構築
  dags_result <- list()
  counts_result <- as.integer(tbl)

  for (i in seq_along(tbl)) {
    key <- names(tbl)[i]
    match_idx <- which(dag_keys == key)[1]
    mat <- dag_list[[match_idx]]

    # ★修正: 両方とも which(mat != 0) で統一
    idx <- which(mat != 0, arr.ind = TRUE)

    if (nrow(idx) == 0) {
      if (split_by_causal_effect_sign) {
        dags_result[[i]] <- data.frame(
          from = integer(0), to = integer(0),
          sign = integer(0)
        )
      } else {
        dags_result[[i]] <- data.frame(from = integer(0), to = integer(0))
      }
    } else if (split_by_causal_effect_sign) {
      dags_result[[i]] <- data.frame(
        from = idx[, 2],
        to   = idx[, 1],
        sign = sapply(seq_len(nrow(idx)), function(k) mat[idx[k, 1], idx[k, 2]])
      )
    } else {
      dags_result[[i]] <- data.frame(
        from = idx[, 2],
        to   = idx[, 1]
      )
    }
  }

  return(list(dag = dags_result, count = counts_result))
}


#' ブートストラップ確率を取得
#'
#' @param result BootstrapResult オブジェクト
#' @param min_causal_effect 因果効果の最小閾値 (NULL = 0)
#' @return 確率行列 (n_features x n_features)
#' @export
#' @examples
#' LiNGAM_sample_1000 <- generate_lingam_sample_6()
#'
#' bs_model <- LiNGAM_sample_1000$data |>
#'   lingam_direct_bootstrap(n_sampling = 30L, seed = 42)
#'
#' bs_model |>
#'   get_probabilities()
get_probabilities <- function(result, min_causal_effect = NULL) {
  stopifnot(inherits(result, "BootstrapResult"))

  if (is.null(min_causal_effect)) min_causal_effect <- 0.0
  if (min_causal_effect < 0) stop("min_causal_effect must be >= 0.")

  am_array <- result$adjacency_matrices
  am_array[is.na(am_array)] <- 0

  bp <- apply(abs(am_array) > min_causal_effect, MARGIN = c(2, 3), FUN = mean)

  return(bp)
}


#' 総合因果効果リストを取得
#'
#' @param result BootstrapResult オブジェクト
#' @param min_causal_effect 因果効果の最小閾値 (NULL = 0)
#' @return data.frame (from, to, effect, probability)
#' @importFrom stats median
#' @export
#' @examples
#' LiNGAM_sample_1000 <- generate_lingam_sample_6()
#'
#' bs_model <- LiNGAM_sample_1000$data |>
#'   lingam_direct_bootstrap(n_sampling = 30L, seed = 42)
#'
#' bs_model |>
#'   get_total_causal_effects()
get_total_causal_effects <- function(result, min_causal_effect = NULL) {
  stopifnot(inherits(result, "BootstrapResult"))

  if (is.null(min_causal_effect)) min_causal_effect <- 0.0
  if (min_causal_effect < 0) stop("min_causal_effect must be >= 0.")

  te_array <- result$total_effects
  te_array[is.na(te_array)] <- 0

  probs <- apply(abs(te_array) > min_causal_effect, MARGIN = c(2, 3), FUN = mean)

  # 確率 > 0 の因果方向
  idx <- which(abs(probs) > 0, arr.ind = TRUE)
  if (nrow(idx) == 0) {
    return(data.frame(
      from = integer(0), to = integer(0),
      effect = numeric(0), probability = numeric(0)
    ))
  }

  from_vec <- idx[, 2]
  to_vec <- idx[, 1]

  prob_vec <- probs[idx]

  median_mat <- apply(te_array, MARGIN = c(2, 3), FUN = function(x) {
    nonzero <- x[abs(x) > 0]
    if (length(nonzero) == 0) {
      return(0)
    }
    stats::median(nonzero)
  })

  # 中央値も idx を使って一発で抽出
  effect_vec <- median_mat[idx]

  # 確率の降順でソート
  ord <- order(-prob_vec)
  data.frame(
    from        = from_vec[ord],
    to          = to_vec[ord],
    effect      = effect_vec[ord],
    probability = prob_vec[ord]
  )
}


#' 指定した2変数間の全パスとブートストラップ確率を取得
#'
#' @param result BootstrapResult オブジェクト
#' @param from_index 始点インデックス (1-based)
#' @param to_index 終点インデックス (1-based)
#' @param min_causal_effect 因果効果の最小閾値 (NULL = 0)
#' @return data.frame (path, effect, probability)
#' @importFrom stats median
#' @export
#' @examples
#' LiNGAM_sample_1000 <- generate_lingam_sample_6()
#'
#' bs_model <- LiNGAM_sample_1000$data |>
#'   lingam_direct_bootstrap(n_sampling = 30L, seed = 42)
#' bs_model |>
#'   get_paths(1, 6)
get_paths <- function(result, from_index, to_index, min_causal_effect = NULL) {
  stopifnot(inherits(result, "BootstrapResult"))

  if (is.null(min_causal_effect)) min_causal_effect <- 0.0
  if (min_causal_effect < 0) stop("min_causal_effect must be >= 0.")

  am_array <- result$adjacency_matrices
  n_sampling <- dim(am_array)[1]

  # 全パスの収集
  # 事前にリストに蓄積し、最後に unlist
  paths_collector <- vector("list", n_sampling)
  effects_collector <- vector("list", n_sampling)

  for (s in seq_len(n_sampling)) {
    am <- am_array[s, , ]
    res <- find_all_paths(am, from_index, to_index, min_causal_effect)
    if (length(res$paths) > 0) {
      paths_collector[[s]] <- vapply(res$paths, paste, "", collapse = "_")
      effects_collector[[s]] <- res$effects
    }
  }
  paths_str_list <- unlist(paths_collector)
  effects_list <- unlist(effects_collector)

  if (length(paths_str_list) == 0) {
    return(data.frame(path = character(0), effect = numeric(0), probability = numeric(0)))
  }

  # カウント
  tbl <- table(paths_str_list)
  tbl <- sort(tbl, decreasing = TRUE)
  probs <- as.numeric(tbl) / n_sampling

  # 各パスの中央値効果
  path_strs <- names(tbl)
  effects_median <- sapply(path_strs, function(ps) {
    stats::median(effects_list[paths_str_list == ps])
  })

  # パス文字列をリストに変換
  path_list <- lapply(path_strs, function(ps) {
    as.integer(strsplit(ps, "_")[[1]])
  })

  data.frame(
    path        = I(path_list),
    effect      = as.numeric(effects_median),
    probability = probs,
    row.names   = NULL
  )
}


# =============================================================================
# 結果の表示・可視化ヘルパー
# =============================================================================

#' ブートストラップ確率を DiagrammeR で描画
#'
#' @param result BootstrapResult オブジェクト
#' @param labels 変数名ベクトル (NULL可)
#' @param min_causal_effect 表示する最小因果効果
#' @param min_probability 表示する最小確率
#' @param rankdir レイアウト方向
#' @param shape ノード形状
#' @return grViz オブジェクト
#' @export
#' @examples
#' LiNGAM_sample_1000 <- generate_lingam_sample_6()
#'
#' bs_model <- LiNGAM_sample_1000$data |>
#'   lingam_direct_bootstrap(n_sampling = 30L, seed = 42)
#' bs_model |>
#'   plot_bootstrap_probabilities()
plot_bootstrap_probabilities <- function(result,
                                         labels = NULL,
                                         min_causal_effect = NULL,
                                         min_probability = 0.5,
                                         rankdir = "TB",
                                         shape = "circle") {
  stopifnot(inherits(result, "BootstrapResult"))

  if (!requireNamespace("DiagrammeR", quietly = TRUE)) {
    stop("Package 'DiagrammeR' is required. Please install it.", call. = FALSE)
  }

  bp <- get_probabilities(result, min_causal_effect)
  p <- ncol(bp)
  if (is.null(labels)) labels <- paste0("x", seq_len(p) - 1)

  # エッジの生成
  edge_lines <- c()
  for (i in 1:p) {
    for (j in 1:p) {
      if (bp[i, j] >= min_probability) {
        edge_lines <- c(
          edge_lines,
          sprintf(
            "  %s -> %s [label = ' %.2f', penwidth = %.1f]",
            labels[j], labels[i], bp[i, j], bp[i, j] * 3 + 0.5
          )
        )
      }
    }
  }

  if (length(edge_lines) == 0) {
    message("No edges exceed the specified threshold. Please lower 'min_probability'.")
    return(invisible(NULL))
  }

  dot <- paste0(
    "digraph bootstrap_result {\n",
    sprintf("  graph [rankdir = %s, fontsize = 14,\n", rankdir),
    "         label = 'Bootstrap Probabilities',\n",
    "         labelloc = t, fontname = 'Helvetica-Bold']\n",
    sprintf("  node [shape = %s, style = filled, fillcolor = lightyellow,\n", shape),
    "        fontname = Helvetica, fontsize = 14, width = 0.6]\n",
    "  edge [fontname = Helvetica, fontsize = 10, fontcolor = blue, color = gray40]\n\n",
    paste(edge_lines, collapse = "\n"), "\n",
    "}\n"
  )

  DiagrammeR::grViz(dot)
}


#' ブートストラップ結果から因果効果の代表値の隣接行列を作成
#'
#' @param result BootstrapResult オブジェクト
#' @param stat 代表値 ("mean" or "median")
#' @param min_causal_effect 因果効果の最小閾値（これ以下はゼロ扱い）(NULL = 0)
#' @param min_probability この確率未満のエッジはゼロにする (NULL = 0)
#' @param labels 変数名ベクトル (NULL可)
#' @return 隣接行列 (n_features x n_features)
#' @export
#' @examples
#' LiNGAM_sample_1000 <- generate_lingam_sample_6()
#'
#' bs_model <- LiNGAM_sample_1000$data |>
#'   lingam_direct_bootstrap(n_sampling = 30L, seed = 42)
#' bs_model |>
#'   get_adjacency_matrix_summary()
get_adjacency_matrix_summary <- function(result,
                                         stat = "median",
                                         min_causal_effect = NULL,
                                         min_probability = NULL,
                                         labels = NULL) {
  stopifnot(inherits(result, "BootstrapResult"))
  if (!(stat %in% c("mean", "median"))) {
    stop("'stat' must be either 'mean' or 'median'.")
  }

  if (is.null(min_causal_effect)) min_causal_effect <- 0.0
  if (is.null(min_probability)) min_probability <- 0.0

  am_array <- result$adjacency_matrices
  am_array[is.na(am_array)] <- 0
  n_sampling <- dim(am_array)[1]
  n_features <- dim(am_array)[2]

  B <- matrix(0, nrow = n_features, ncol = n_features)

  for (i in 1:n_features) {
    for (j in 1:n_features) {
      if (i == j) next

      # 全ブートストラップでの (i,j) 要素を取得
      vals <- am_array[, i, j]

      # 閾値を超える値のみ抽出
      significant <- vals[abs(vals) > min_causal_effect]

      # 確率の計算
      prob <- length(significant) / n_sampling

      if (prob < min_probability || length(significant) == 0) {
        B[i, j] <- 0
      } else {
        B[i, j] <- if (stat == "mean") mean(significant) else median(significant)
      }
    }
  }

  # 変数名の付与
  if (!is.null(labels)) {
    rownames(B) <- labels
    colnames(B) <- labels
  }

  return(B)
}
