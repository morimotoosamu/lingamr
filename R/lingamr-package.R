#' @details
#' tutorial: `vignette("lingamr", package = "lingamr")`
#' @references
#' Shimizu, S., et al. (2011). DirectLiNGAM: A direct method for learning a
#' linear non-Gaussian structural equation model. *Journal of Machine Learning
#' Research*, 12, 1225-1248.
#'
#' Hyvärinen, A., Zhang, K., Shimizu, S., & Hoyer, P. O. (2010). Estimation of a
#' structural vector autoregression model using non-Gaussianity. *Journal of
#' Machine Learning Research*, 11, 1709-1731.
#'
#' Moneta, A., Entner, D., Hoyer, P. O., & Coad, A. (2013). Causal inference by
#' independent component analysis: Theory and applications. *Oxford Bulletin of
#' Economics and Statistics*, 75(5), 705-730. VARLiNGAM R code:
#' <https://sites.google.com/site/dorisentner/publications/VARLiNGAM>
#'
#' Python implementation (DirectLiNGAM, VAR-LiNGAM): cdt15/lingam,
#' <https://github.com/cdt15/lingam>
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

# ggplot2::aes() 内で参照する data frame の列名を、R CMD check の
# 「no visible binding for global variable」NOTE 回避のため宣言する。
utils::globalVariables(c("residual", "x", "y", "xend", "yend", "label", "name", "variable"))
