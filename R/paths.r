# =============================================================================
# Direct LiNGAM - DAG のパス列挙と経路効果
# Based on the Python implementation from the LiNGAM Project
# https://sites.google.com/view/sshimizu06/lingam
# https://github.com/cdt15/lingam
#
# License: MIT + file LICENSE
#
# Original work:
#   Copyright (c) 2019 T.Ikeuchi, G.Haraoka, M.Ide, W.Kurebayashi, S.Shimizu
#
# Portions of this work:
#   Copyright (c) 2026 O.Morimoto
# =============================================================================


#' DAG 中の全パスを深さ優先探索で列挙する
#'
#' `B[i, j]` が j → i を表す隣接行列を受け取り、`from_index` から
#' `to_index` に至る全パスとそれぞれの経路効果（係数の積）を返す。
#'
#' @param adjacency_matrix 隣接行列 (n x n)。`B[i,j]` は j → i の係数。
#' @param from_index 始点インデックス (1-based)
#' @param to_index 終点インデックス (1-based)
#' @param min_causal_effect このしきい値以下の係数は存在しないエッジとみなす
#' @return list(paths, effects)
#' @keywords internal
find_all_paths <- function(adjacency_matrix, from_index, to_index, min_causal_effect = 0.0) {
  B <- adjacency_matrix
  B[is.na(B)] <- 0
  B[abs(B) <= min_causal_effect] <- 0

  p <- ncol(B)
  paths <- list()
  effects <- c()

  dfs <- function(current, target, visited, path, effect) {
    if (current == target && length(path) > 1) {
      paths[[length(paths) + 1]] <<- path
      effects[length(effects) + 1] <<- effect
      return()
    }
    for (next_node in seq_len(p)) {
      if (next_node %in% visited) next
      if (B[next_node, current] == 0) next
      dfs(
        next_node, target, c(visited, next_node),
        c(path, next_node), effect * B[next_node, current]
      )
    }
  }

  dfs(from_index, to_index, c(from_index), c(from_index), 1.0)
  list(paths = paths, effects = effects)
}


#' 隣接行列から2変数間の総合因果効果を計算する
#'
#' `find_all_paths()` で列挙した全経路効果の総和を返す。
#' パスが存在しない場合は 0 を返す。
#'
#' @param adjacency_matrix 隣接行列 (n x n)。`B[i,j]` は j → i の係数。
#' @param from_index 原因変数のインデックス (1-based)
#' @param to_index 結果変数のインデックス (1-based)
#' @return 総合因果効果（スカラー）
#' @keywords internal
calculate_total_effect <- function(adjacency_matrix, from_index, to_index) {
  result <- find_all_paths(adjacency_matrix, from_index, to_index, min_causal_effect = 0.0)
  if (length(result$effects) == 0) 0.0 else sum(result$effects)
}
