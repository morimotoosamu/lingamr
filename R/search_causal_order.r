# =============================================================================
# Direct LiNGAM - 因果順序の探索（pwling / kernel）と事前知識の処理
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


#' 事前知識から部分順序を抽出
#' @param pk 事前知識行列 (NaN = 不明)
#' @return matrix (n x 2), 各行は (from, to) の部分順序
#' @keywords internal
extract_partial_orders <- function(pk) {
  # パスがあるペア (pk == 1)
  path_idx <- which(pk == 1, arr.ind = TRUE)
  # パスがないペア (pk == 0)
  no_path_idx <- which(pk == 0, arr.ind = TRUE)

  # --- パスありペアの矛盾チェック ---
  if (nrow(path_idx) > 0) {
    check_pairs <- rbind(path_idx, path_idx[, 2:1, drop = FALSE])
    dup <- duplicated(check_pairs) | duplicated(check_pairs, fromLast = TRUE)
    if (any(dup)) {
      bad <- unique(check_pairs[dup, , drop = FALSE])
      stop(paste(
        "The prior knowledge contains inconsistencies at indices:",
        paste(apply(bad, 1, function(r) paste0("(", r[1], ",", r[2], ")")),
          collapse = ", "
        )
      ))
    }
  }

  # --- パスなしペアの重複除去 ---
  if (nrow(no_path_idx) > 0) {
    check_pairs2 <- rbind(no_path_idx, no_path_idx[, 2:1, drop = FALSE])
    # 双方向に 0 が入っているペアを見つけて除外
    pair_key <- paste(check_pairs2[, 1], check_pairs2[, 2], sep = ",")
    tbl <- table(pair_key)
    dup_keys <- names(tbl[tbl > 1])
    # no_path_idx のうち、双方向に存在するものを除外
    no_path_key <- paste(no_path_idx[, 1], no_path_idx[, 2], sep = ",")
    keep <- !(no_path_key %in% dup_keys)
    no_path_idx <- no_path_idx[keep, , drop = FALSE]
  }

  # path_pairs と no_path_pairs[:, [1,0]] を結合
  combined <- matrix(nrow = 0, ncol = 2)
  if (nrow(path_idx) > 0) {
    combined <- rbind(combined, path_idx)
  }
  if (nrow(no_path_idx) > 0) {
    combined <- rbind(combined, no_path_idx[, 2:1, drop = FALSE])
  }

  if (nrow(combined) == 0) {
    return(matrix(nrow = 0, ncol = 2))
  }

  combined <- unique(combined)
  # [to, from] -> [from, to]
  result <- combined[, 2:1, drop = FALSE]
  colnames(result) <- NULL
  return(result)
}


#' 残差 (xi を xj に回帰したときの残差)
#' 残差ベクトルの計算
#' @param xi 対象変数ベクトル
#' @param xj 説明変数ベクトル
#' @param standardized データが標準化済みか (default: FALSE)
#' @return 回帰後の残差ベクトル
#' @keywords internal
residual_vec <- function(xi, xj, standardized = FALSE) {
  if (standardized) {
    # 高速版: mean = 0 を仮定
    beta <- sum(xi * xj) / sum(xj * xj)
  } else {
    # 汎用版: 中心化を含む
    xi_c <- xi - mean(xi)
    xj_c <- xj - mean(xj)
    beta <- sum(xi_c * xj_c) / sum(xj_c * xj_c)
  }
  xi - beta * xj
}


#' エントロピーの最大エントロピー近似
#' @param u 入力ベクトル
#' @return 近似エントロピー値
#' @keywords internal
entropy_approx <- function(u) {
  n <- length(u)
  k1 <- 79.047
  k2 <- 7.4129
  gamma <- 0.37457
  (1 + log(2 * pi)) / 2 -
    k1 * (sum(log(cosh(u))) / n - gamma)^2 -
    k2 * (sum(u * exp(-u^2 / 2)) / n)^2
}


#' 相互情報量の差
#' @param xi_std 標準化された xi
#' @param xj_std 標準化された xj
#' @param ri_j xi を xj で回帰した残差
#' @param rj_i xj を xi で回帰した残差
#' @return 相互情報量の差
#' @keywords internal
diff_mutual_info <- function(xi_std, xj_std, ri_j, rj_i) {
  sd_ri_j <- sd_pop(ri_j)
  sd_rj_i <- sd_pop(rj_i)
  (entropy_approx(xj_std) + entropy_approx(ri_j / sd_ri_j)) -
    (entropy_approx(xi_std) + entropy_approx(rj_i / sd_rj_i))
}


#' 候補変数の探索
#' @param U 現在の未確定変数の集合
#' @param Aknw 事前知識行列
#' @param apply_prior_knowledge_softly ソフト適用の有無
#' @param partial_orders 抽出された部分順序
#' @return list(Uc, Vj)
#' @keywords internal
search_candidate <- function(U, Aknw, apply_prior_knowledge_softly, partial_orders) {
  # 事前知識なし

  if (is.null(Aknw)) {
    return(list(Uc = U, Vj = integer(0)))
  }

  # --- ハード適用 ---
  if (!apply_prior_knowledge_softly) {
    if (!is.null(partial_orders) && nrow(partial_orders) > 0) {
      Uc <- setdiff(U, partial_orders[, 2])
      if (length(Uc) == 0) Uc <- U
      return(list(Uc = Uc, Vj = integer(0)))
    } else {
      return(list(Uc = U, Vj = integer(0)))
    }
  }

  # --- ソフト適用 ---
  # 外生変数の探索
  Uc <- integer(0)
  for (j in U) {
    index <- setdiff(U, j)
    # NA (不明) を含む行は sum が NA → isTRUE() で FALSE になり候補から外れる
    # (Python 版の NaN.sum() == 0 が False になる挙動と同じ)
    if (isTRUE(sum(Aknw[j, index]) == 0)) {
      Uc <- c(Uc, j)
    }
  }

  # 内生変数の探索 → 候補の絞り込み
  if (length(Uc) == 0) {
    U_end <- integer(0)
    for (j in U) {
      index <- setdiff(U, j)
      s <- sum(Aknw[j, index], na.rm = TRUE)
      if (!is.na(s) && s > 0) {
        U_end <- c(U_end, j)
      }
    }
    # シンク特徴量
    for (i in U) {
      index <- setdiff(U, i)
      if (isTRUE(sum(Aknw[index, i]) == 0)) {
        U_end <- c(U_end, i)
      }
    }
    Uc <- setdiff(U, unique(U_end))
    if (length(Uc) == 0) Uc <- U
  }

  # V^(j) の構築
  Vj <- integer(0)
  for (i in U) {
    if (i %in% Uc) next
    if (isTRUE(sum(Aknw[i, Uc]) == 0)) {
      Vj <- c(Vj, i)
    }
  }

  return(list(Uc = Uc, Vj = Vj))
}


#' pwling による因果順序の探索
#' @param X データ行列
#' @param U 全変数インデックス
#' @param Uc 候補変数のインデックス
#' @param Vj 事前知識に基づく変数集合
#' @return 選ばれた変数のインデックス
#' @keywords internal
search_causal_order_pwling <- function(X, U, Uc, Vj) {
  if (length(Uc) == 1) return(Uc[1])
  n <- nrow(X)
  p <- ncol(X)

  # --- 一括で標準化し、各列のエントロピーを事前計算（ペアに依存しないため）---
  X_std <- matrix(0, nrow = n, ncol = p)
  H <- numeric(p)
  for (k in U) {
    xk <- X[, k]
    xk <- xk - sum(xk) / n
    X_std[, k] <- xk / sqrt(sum(xk * xk) / n)
    H[k] <- entropy_approx(X_std[, k])
  }

  # 標準化済みなので相関行列は crossprod / n（BLAS で一括計算）。
  # 回帰係数 beta と残差の母標準偏差 sqrt(1 - r^2) はここから解析的に得られる。
  R <- crossprod(X_std[, U, drop = FALSE]) / n
  pos <- integer(p)
  pos[U] <- seq_along(U)

  in_Uc <- logical(p)
  in_Uc[Uc] <- TRUE
  in_Vj <- logical(p)
  in_Vj[Vj] <- TRUE

  M_acc <- numeric(p)
  for (i in Uc) {
    xi_std <- X_std[, i]
    for (j in U) {
      if (i == j) next
      # diff_mutual_info は反対称（dm_ji = -dm_ij）なので、
      # 両方が候補のペアは i < j の側で1回だけ計算して両者に加算する
      if (in_Uc[j] && j < i) next
      xj_std <- X_std[, j]
      r_ij <- R[pos[i], pos[j]]
      sd_r <- sqrt(max(0, 1 - r_ij^2))

      H_ri_j <- if (in_Vj[i] && in_Uc[j]) {
        H[i]
      } else {
        entropy_approx((xi_std - r_ij * xj_std) / sd_r)
      }
      H_rj_i <- if (in_Vj[j] && in_Uc[i]) {
        H[j]
      } else {
        entropy_approx((xj_std - r_ij * xi_std) / sd_r)
      }

      dm <- (H[j] + H_ri_j) - (H[i] + H_rj_i)
      M_acc[i] <- M_acc[i] + min(0, dm)^2
      if (in_Uc[j]) M_acc[j] <- M_acc[j] + min(0, -dm)^2
    }
  }
  return(Uc[which.max(-M_acc[Uc])])
}


#' カーネル法の相互情報量：変数1側の前計算
#'
#' `kernel_mi_core()` で使う行列 `E1 = tmp1^-1 K1`（`tmp1 = K1 + n*kappa/2 * I`）
#' を計算する。候補変数ごとに1回だけ呼べばよく、ペアごとの再計算を避けられる。
#'
#' @param x 変数1のベクトル
#' @param kappa 正則化パラメータ
#' @param sigma ガウスカーネルの幅
#' @return 行列 E1 (n x n)
#' @keywords internal
kernel_mi_prepare <- function(x, kappa, sigma) {
  n <- length(x)
  K <- exp(-1 / (2 * sigma^2) * outer(x, x, "-")^2)
  c0 <- n * kappa / 2
  tmp <- K
  diag(tmp) <- diag(tmp) + c0
  # tmp と K は可換なので E = tmp^-1 K = I - c0 * tmp^-1（対称）
  E <- -c0 * chol2inv(chol(tmp))
  diag(E) <- diag(E) + 1
  E
}


#' カーネル法の相互情報量：本体
#'
#' 求める量は 2n x 2n 行列の logdet の差だが、ブロック構造と Schur 補行列により
#' n x n の Cholesky 分解だけで等価に計算できる：
#' `MI = -1/2 * (logdet(tmp2^2 - K2 K1 tmp1^-2 K1 K2) - logdet(tmp2^2))`
#'
#' @param E1 `kernel_mi_prepare()` で前計算した変数1側の行列
#' @param x2 変数2のベクトル
#' @param kappa 正則化パラメータ
#' @param sigma ガウスカーネルの幅
#' @return 相互情報量
#' @keywords internal
kernel_mi_core <- function(E1, x2, kappa, sigma) {
  n <- length(x2)
  K2 <- exp(-1 / (2 * sigma^2) * outer(x2, x2, "-")^2)
  tmp2 <- K2
  diag(tmp2) <- diag(tmp2) + n * kappa / 2
  W <- E1 %*% K2                       # = tmp1^-1 K1 K2
  S <- crossprod(tmp2) - crossprod(W)  # Schur 補行列（tmp2 は対称なので crossprod = tmp2^2）
  logdet_S <- 2 * sum(log(diag(chol(S))))
  logdet_tmp2_sq <- 4 * sum(log(diag(chol(tmp2))))
  (-1 / 2) * (logdet_S - logdet_tmp2_sq)
}


#' カーネル法による相互情報量
#' @param x1 変数1
#' @param x2 変数2
#' @param param パラメータベクトル (kappa, sigma)
#' @return 相互情報量
#' @keywords internal
mutual_information_kernel <- function(x1, x2, param) {
  kappa <- param[1]
  sigma <- param[2]
  E1 <- kernel_mi_prepare(x1, kappa, sigma)
  kernel_mi_core(E1, x2, kappa, sigma)
}


#' カーネル法による因果順序の探索
#' @param X データ行列
#' @param U 全変数
#' @param Uc 候補変数
#' @param Vj 事前知識集合
#' @return 選ばれた変数のインデックス
#' @keywords internal
search_causal_order_kernel <- function(X, U, Uc, Vj) {
  if (length(Uc) == 1) {
    return(Uc[1])
  }

  n <- nrow(X)
  if (n > 1000) {
    param <- c(2e-3, 0.5)
  } else {
    param <- c(2e-2, 1.0)
  }

  kappa <- param[1]
  sigma <- param[2]
  Tkernels <- numeric(length(Uc))

  for (idx in seq_along(Uc)) {
    j <- Uc[idx]
    # 候補 j 側のカーネル行列・逆行列は内側ループで不変なので1回だけ計算する
    E1 <- kernel_mi_prepare(X[, j], kappa, sigma)
    Tkernel <- 0
    for (i in U) {
      if (i == j) next
      ri_j <- if (j %in% Vj && i %in% Uc) {
        X[, i]
      } else {
        residual_vec(X[, i], X[, j])
      }
      Tkernel <- Tkernel + kernel_mi_core(E1, ri_j, kappa, sigma)
    }
    Tkernels[idx] <- Tkernel
  }

  return(Uc[which.min(Tkernels)])
}
