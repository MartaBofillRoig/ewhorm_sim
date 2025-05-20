#' help function to calculate 1-alpha confidence interval of concordance with weighted variances according to sample sizes.
#' @description help function to calculate 1-alpha confidence interval of concordance with pooled variances (modified Hanley-McNeil approach) according to Newcombe, 2006.
#'
#' @param thetahat estimated empirical concordance
#' @param theta confidence limit - to be found by uniroot function
#' @param N1 sample size of stage 1
#' @param N2 sample size of stage 2
#' @param alpha significance level
#' @keywords internal
#' @returns variance of concordance
#' @export
#' @details eWHORM simulations
#' @author Sonja Zehetmayer
#' 
CI2_help<-function(theta,thetahat,alpha,N1,N2) #pooled variance
{
  N<-N1+N2
  var<-sqrt(N1/N)*var2(theta,N1/3,N1/3)+sqrt(N2/N)*var2(theta,N2/3,N2/3)
  sqrt(var)*qnorm(alpha)-theta+thetahat
}
