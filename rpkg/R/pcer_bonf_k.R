#' Partial conditional error for a Bonferroni local test at size k (k = |I| is the size of the intersection hypothesis)
#' @description Partial conditional error for a Bonferroni local test at size k (k = |I| is the size of the intersection hypothesis)
#'
#' @param z1 A numeric vector giving first stage z-values computed from first stage p-values.
#' @param v A numeric vector giving the proportions of pre-planned measurements collected up to the interim analysis.
#' @param alpha significance level
#' @export
#' @details eWHORM simulations
#' @author Sonja Zehetmayer
#' 

pcer_bonf_k <- function(z1, v, alpha, k) {
  ck <- qnorm(1 - alpha / k)  # one-sided critical value for local level alpha/k
  1 - pnorm((ck - sqrt(v) * z1) / sqrt(1 - v))
}

