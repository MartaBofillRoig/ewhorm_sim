#' help function to calculate quantiles (e.g., 1-alpha lower confidence interval or median) of concordance with inverse normal method to combine concordance of stage 1 and stage 2.
#' @description help function to calculate 1-alpha confidence interval of concordance with pooled variances (modified Hanley-McNeil approach) according to Newcombe, 2006.
#'
#' @param thetahat1 estimated empirical concordance of stage 1
#' @param thetahat2 estimated empirical concordance of stage 2
#' @param theta confidence limit - to be found by uniroot function
#' @param N1 sample size of group 1
#' @param N2 sample size of group 2
#' @param value quantile - significance level or median
#' @keywords internal
#' @returns variance of concordance
#' @export
#' @details eWHORM simulations
#' @author Sonja Zehetmayer
#' 


CIinvn_help<-function(theta,thetahat1,thetahat2,value,N1,N2) 
{
  N<-N1+N2
  (sqrt(N1/N)*(thetahat1-theta)/sqrt(var2(theta,N1/3,N1/3))+sqrt((N2/N))*(thetahat2-theta)/sqrt(var2(theta,N2/3,N2/3)))-qnorm(1-value)
  #(1-pnorm(sqrt(N1/N)*(thetahat1-theta)/sqrt(var2(theta,m1,m1))+sqrt((N2/N))*(thetahat2-theta)/sqrt(var2(theta,m2,m2))))-alpha#qnorm(1-alpha)
}

