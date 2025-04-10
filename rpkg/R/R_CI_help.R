#' help function to calculate 1-alpha confidence interval of concordance with pooled variances (modified Hanley-McNeil approach) according to Newcombe, 2006.
#' @description help function to calculate 1-alpha confidence interval of concordance with pooled variances (modified Hanley-McNeil approach) according to Newcombe, 2006.
#'
#' @param thetahat estimated empirical concordance
#' @param theta confidence limit - to be found by uniroot function
#' @param N1 sample size of stage 1
#' @param N2 sample size of stage 2
#' @param alpha significance level
#' @keywords internal
#' @export
#' @details eWHORM simulations
#' @author Sonja Zehetmayer
#' 

CI_help<-function(theta,thetahat,N1,N2,alpha) #pooled variance
{
  sqrt(var2(theta,N1/3,N2/3))*qnorm(alpha)-theta+thetahat
}
