# R/imports.r
# Single collection point for all @importFrom stats declarations used across
# the package. Do not add per-file @importFrom stats tags elsewhere; add the
# function to one of the lines below instead. Each @importFrom tag must stay
# on a single line (roxygen2 rejects wrapped ones).

#' @importFrom stats coef cor.test cov deviance dist lm lm.fit logLik
#' @importFrom stats median quantile sd var
#' @importFrom stats binomial glm optim pchisq plogis
#' @importFrom stats rexp rlnorm rnorm rt runif
NULL
