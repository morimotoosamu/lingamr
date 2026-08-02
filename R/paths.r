# =============================================================================
# Direct LiNGAM - DAG path enumeration and path effects
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


#' Enumerate all paths in a DAG via depth-first search
#'
#' Takes an adjacency matrix where `B[i, j]` represents j -> i, and returns
#' all paths from `from_index` to `to_index` together with each path effect
#' (the product of the coefficients).
#'
#' @param adjacency_matrix Adjacency matrix (n x n). `B[i,j]` is the coefficient of j -> i.
#' @param from_index Start index (1-based)
#' @param to_index End index (1-based)
#' @param min_causal_effect Coefficients at or below this threshold are treated as nonexistent edges
#' @return list(paths, effects)
#' @keywords internal
find_all_paths <- function(adjacency_matrix, from_index, to_index, min_causal_effect = 0.0) {
  B <- adjacency_matrix
  B[is.na(B)] <- 0
  B[abs(B) <= min_causal_effect] <- 0

  p <- ncol(B)
  acc <- new.env(parent = emptyenv())
  acc$paths <- list()
  acc$effects <- list()

  dfs <- function(current, target, visited, path, effect) {
    if (current == target && length(path) > 1) {
      n <- length(acc$paths) + 1L
      acc$paths[[n]] <- path
      acc$effects[[n]] <- effect
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
  list(paths = acc$paths, effects = unlist(acc$effects, use.names = FALSE))
}


#' Compute the total causal effect between two variables from an adjacency matrix
#'
#' Returns the sum of all path effects enumerated by `find_all_paths()`.
#' Returns 0 if no path exists.
#'
#' @param adjacency_matrix Adjacency matrix (n x n). `B[i,j]` is the coefficient of j -> i.
#' @param from_index Index of the cause variable (1-based)
#' @param to_index Index of the effect variable (1-based)
#' @return Total causal effect (scalar)
#' @keywords internal
calculate_total_effect <- function(adjacency_matrix, from_index, to_index) {
  result <- find_all_paths(adjacency_matrix, from_index, to_index, min_causal_effect = 0.0)
  if (length(result$effects) == 0) 0.0 else sum(result$effects)
}


#' Compute all pairwise total causal effects at once
#'
#' Equivalent to calling [calculate_total_effect()] for every (from, to)
#' pair, but computed as the truncated Neumann series `B + B^2 + ... +
#' B^(p-1)`. In a DAG every walk is a simple path (a repeated node would
#' imply a cycle), so the series contains exactly the same path-product
#' terms the DFS enumerates, and entries with no connecting path are exactly
#' zero. Only valid for acyclic `B` -- with cycles the series would count
#' non-simple walks that [find_all_paths()] excludes.
#'
#' @param adjacency_matrix Adjacency matrix (p x p). `B[i,j]` is the
#'   coefficient of j -> i; NAs are treated as absent edges.
#' @return p x p matrix `TE` with `TE[to, from]` = total effect of from on to
#' @keywords internal
calculate_total_effects_all <- function(adjacency_matrix) {
  B <- adjacency_matrix
  B[is.na(B)] <- 0
  p <- ncol(B)
  TE <- matrix(0, p, p)
  P <- B
  for (len in seq_len(max(p - 1L, 1L))) {
    TE <- TE + P
    if (len == p - 1L || !any(P != 0)) break
    P <- P %*% B
  }
  TE
}
