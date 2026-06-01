#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

# ggplot2::aes() 内で参照する data frame の列名を、R CMD check の
# 「no visible binding for global variable」NOTE 回避のため宣言する。
utils::globalVariables("residual")
