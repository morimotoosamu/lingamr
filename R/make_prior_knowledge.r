#' 事前知識行列を作成
#'
#' @param n_variables 変数の数
#' @param exogenous_variables 外生変数 (1-based index または変数名, NULL可)
#'   指定した変数は他のどの変数からも影響を受けないとする
#' @param sink_variables シンク変数 (1-based index または変数名, NULL可)
#'   指定した変数は他のどの変数にも影響を与えないとする
#' @param paths 有向パスが存在する変数ペア (NULL可)
#'   list(c(from, to), ...) の形式。インデックスまたは変数名で指定
#' @param no_paths 有向パスが存在しない変数ペア (NULL可)
#'   list(c(from, to), ...) の形式。インデックスまたは変数名で指定
#' @param labels 変数名ベクトル (NULL可)
#'   変数名で指定する場合は必須。data.frame の colnames() 等を渡す
#'
#' @return 事前知識行列 (n_variables x n_variables)
#'   -1: 不明, 0: パスなし, 1: パスあり
#'
#' @examples
#' # インデックスで指定
#' pk <- make_prior_knowledge(6, exogenous_variables = c(4))
#'
#' # 変数名で指定
#' pk <- make_prior_knowledge(6,
#'   exogenous_variables = "x3",
#'   sink_variables = c("x1", "x4"),
#'   paths = list(c("x3", "x0"), c("x3", "x2")),
#'   no_paths = list(c("x5", "x2")),
#'   labels = c("x0", "x1", "x2", "x3", "x4", "x5")
#' )
make_prior_knowledge <- function(n_variables,
                                 exogenous_variables = NULL,
                                 sink_variables = NULL,
                                 paths = NULL,
                                 no_paths = NULL,
                                 labels = NULL) {
  # --- 基本バリデーション ---
  if (!is.numeric(n_variables) || n_variables < 2) {
    stop("n_variables must be an integer >= 2.")
  }
  n_variables <- as.integer(n_variables)

  # --- labels のバリデーション ---
  if (!is.null(labels)) {
    if (length(labels) != n_variables) {
      stop(sprintf(
        "Length of 'labels' (%d) must equal n_variables (%d).",
        length(labels), n_variables
      ))
    }
    if (anyDuplicated(labels) > 0) {
      stop("'labels' must not contain duplicate names.")
    }
  }

  # --- 変数名 → インデックス変換ヘルパー ---
  resolve_variable <- function(x, arg_name) {
    if (is.character(x)) {
      if (is.null(labels)) {
        stop(sprintf(
          "'%s' contains variable name(s), but 'labels' is not specified.\n  Provide labels or use integer indices.",
          arg_name
        ))
      }
      pos <- match(x, labels)
      if (any(is.na(pos))) {
        bad <- x[is.na(pos)]
        stop(sprintf(
          "Variable name(s) %s in '%s' not found in labels.\n  Available: %s",
          paste0("'", bad, "'", collapse = ", "),
          arg_name,
          paste(labels, collapse = ", ")
        ))
      }
      return(pos)
    } else if (is.numeric(x)) {
      x <- as.integer(x)
      if (any(x < 1) || any(x > n_variables)) {
        stop(sprintf("'%s' indices must be between 1 and %d.", arg_name, n_variables))
      }
      return(x)
    } else {
      stop(sprintf("'%s' must be integer indices or character variable names.", arg_name))
    }
  }

  resolve_pairs <- function(pairs, arg_name) {
    if (is.null(pairs)) {
      return(NULL)
    }
    if (!is.list(pairs)) stop(sprintf("'%s' must be a list of length-2 vectors.", arg_name))

    resolved <- lapply(seq_along(pairs), function(i) {
      p <- pairs[[i]]
      if (length(p) != 2) {
        stop(sprintf("Each element of '%s' must be a length-2 vector (from, to).", arg_name))
      }
      # ペア内で型が混在していても各要素を個別に解決
      from_idx <- resolve_variable(p[1], sprintf("%s[[%d]][1]", arg_name, i))
      to_idx <- resolve_variable(p[2], sprintf("%s[[%d]][2]", arg_name, i))
      c(from_idx, to_idx)
    })
    return(resolved)
  }

  # --- 引数の解決 ---
  if (!is.null(exogenous_variables)) {
    exogenous_variables <- resolve_variable(exogenous_variables, "exogenous_variables")
  }
  if (!is.null(sink_variables)) {
    sink_variables <- resolve_variable(sink_variables, "sink_variables")
  }
  paths <- resolve_pairs(paths, "paths")
  no_paths <- resolve_pairs(no_paths, "no_paths")

  # --- 事前知識行列の初期化（全て -1 = 不明）---
  prior_knowledge <- matrix(-1L, nrow = n_variables, ncol = n_variables)

  if (!is.null(labels)) {
    rownames(prior_knowledge) <- labels
    colnames(prior_knowledge) <- labels
  }

  # --- no_paths の設定: (from, to) → prior_knowledge[to, from] = 0 ---
  if (!is.null(no_paths)) {
    for (pair in no_paths) {
      prior_knowledge[pair[2], pair[1]] <- 0L
    }
  }

  # --- paths の設定: (from, to) → prior_knowledge[to, from] = 1 ---
  if (!is.null(paths)) {
    for (pair in paths) {
      prior_knowledge[pair[2], pair[1]] <- 1L
    }
  }

  # --- sink_variables: 他の変数に影響を与えない → 列を全て 0 ---
  if (!is.null(sink_variables)) {
    for (var in sink_variables) {
      prior_knowledge[, var] <- 0L
    }
  }

  # --- exogenous_variables: 他の変数から影響を受けない → 行を全て 0 ---
  if (!is.null(exogenous_variables)) {
    for (var in exogenous_variables) {
      prior_knowledge[var, ] <- 0L
    }
  }

  # --- 対角要素は -1 ---
  diag(prior_knowledge) <- -1L

  return(prior_knowledge)
}
