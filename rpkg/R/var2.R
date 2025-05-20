#' modified Hanley–McNeil approach for variance  estimation  (Newcombe, 2006)
#' @description #modified Hanley–McNeil approach for variance  estimation  (Newcombe, 2006)
#'
#' @param theta concordance
#' @param m sample size of group 1
#' @param n sample size of group 2
#' @keywords internal
#' @returns variance of concordance
#' @export
#' @details eWHORM simulations
#' @author Sonja Zehetmayer
#' 

var2<-function(theta,m,n)  
{
  mstar<-1/2*(m+n)-1
  nstar<-1/2*(m+n)-1
  theta*(1-theta)*(1+nstar*(1-theta)/(2-theta)+mstar*theta/(1+theta))/(m*n)
}
