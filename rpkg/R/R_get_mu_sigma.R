#' Function to transform mean and standard deviation from "original" baseline values to parameters of the normal distribution of the log values.
#' @description Function to transform mean and standard deviation from "original" baseline values to parameters of the normal distribution of the log values.
#'
#' @param mu_raw_0 num mean value of original data at baseline
#' @param sd_raw_0 num standard deviation of original data at baseline
#' @param reductrate_6 vector of reduction rates after 6 months for each dose (between 0 and 1)
#' @param reductrate_12 vector of reduction rates after 12 months for each dose (between 0 and 1)
#' @param rho correlation between baseline and 6 months follow-up observations and between 6 months and 12 months follow-up observations
#' @keywords internal
#' @returns list of mean values and variance covariance matrix
#' @export
#' @details eWHORM simulations
#' @importFrom gtools combinations
#' @author Marta Bofill Roig
#'



get_mu_sigma = function(mu_raw_0, sd_raw_0 , reductrate_6 , reductrate_12, rho )
{
  
  
  # define parameters
  mu_raw_6 <- (1-reductrate_6)*mu_raw_0 #mean value six months after baseline
  mu_log_0 <- log(mu_raw_0^2/sqrt(mu_raw_0^2 + sd_raw_0^2))  # calculate the mean for the log transformation for the baseline
  
  sd_log_0 <- sqrt(log(1+(sd_raw_0^2/mu_raw_0^2) )) # sd(log(x0))  Baseline
  sd_raw_6 <- mu_raw_6*sqrt(exp(sd_log_0^2) - 1)   # calculating sd for the log after six month
  
  
  #12 month 
  mu_raw_12 <- (1-reductrate_12)*mu_raw_0 #mean value 12 months after baseline
  sd_raw_12 <- mu_raw_12*sqrt(exp(sd_log_0^2) - 1)   # calculating sd for the log after 12 month
  
  
  
  #Obtain the vector mu for each time point
  mu_log_0<-rep(mu_log_0,4)
  mu_log_6<-log(mu_raw_6^2/sqrt(mu_raw_6^2 + sd_raw_6^2))
  mu_log_12<-log(mu_raw_12^2/sqrt(mu_raw_12^2 + sd_raw_12^2))
  
  
  # Variance-covariance matrix
  # Elements of the matrix
  var0 <- sd_log_0^2 # cov(log(X0)
  var6 <- sd_log_0^2 # cov(log(X06m) 
  var12<- sd_log_0^2 # cov(log(X12m)
  
  #cov_X0_X0 <- var0
  cov_X0_X6m <- rho*sd_log_0^2 # cov(log(X0), log(X6m)) 
  cov_X0_X12m <- rho^2*sd_log_0^2 #cov(log(X0), log(X12m)) 
  cov_X6_X0 <- rho*sd_log_0^2 # cov(log(X6), log(X0m)) 
  #cov_X6_X6 <- var6
  cov_X6_X12m <- rho*sd_log_0^2 #cov(log(X6), log(X12m)) 
  cov_X12m_X0 <- rho^2*sd_log_0^2 #cov(log(X12m), log(X0)) 
  cov_X12m_X6 <- rho*sd_log_0^2 #cov(log(X12m), log(X6m)) 
  #cov_X12m_X12m <- var12  #cov(log(X12m), log(X12m)) 
  
  
  #variance covariance matrix
  sg <- matrix(c(var0, cov_X0_X6m, cov_X0_X12m, cov_X6_X0, var6, cov_X6_X12m,cov_X12m_X0, cov_X12m_X6, var12), ncol=3)
  
  
  return(list(mu_raw_0, c(mu_log_0), c(mu_log_6), c(mu_log_12), sg))
}
