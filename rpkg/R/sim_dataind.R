#' Simulate data from a multi-arm trial with shared control
#' @description Function to simulate trial data (1-stage, multiple arms)
#'
#' @param n_arms number of arms (including control)
#' @param N total sample size
#' @param mu_6m 6-month mean response per arm (vector of length `n_arm`)
#' @param mu_12m 12-month mean response per arm (vector of length `n_arm`)
#' @param sigma covariance matrix between 6- and 12-month responses assumed equal across arms (matrix of dim 2x2)
#' @param rmonth mean number of patients recruited per month (recruitment times assumed exponential distributed / number of patients follows a Poisson distribution)
#' @returns simulated data consisting of the responses at 6 and 12 months, treatment arm, and recruitment time for each subject.
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats model.matrix
#' @export
#' @details eWHORM simulations
#' @author Marta Bofill Roig

# n_arms=4; N=100; mu_0m =c(0,0,0,0); mu_6m =c(1,2,3,4); mu_12m=c(1,2,3,4); sg=matrix(c(1,0,0,0,1,0,0,0,1), ncol=3); rmonth=1

sim_dataind <- function(n_arms, N, mu_0m, mu_6m, mu_12m, sg, rmonth){

  treatments <- factor(c(sample(rep(1:n_arms, floor(N/n_arms))), sample(1:n_arms, N-floor(N/n_arms)*n_arms, replace=T)),
                      # sample(1:n_arms, N, replace = TRUE),
                       levels = 1:n_arms,
                       labels = c("Placebo", "Low", "Medium", "High")[1:n_arms])
  X <- model.matrix(~ treatments - 1)
  y <- X %*% matrix(c(mu_0m, mu_6m, mu_12m), nrow=n_arms, byrow = F) + rmvnorm(n=N, mean = c(0,0,0), sigma = sg )

  # Treatment indicator for dataframe
  max_col_indices <- apply(X, 1, get_max_col_index)
  treat = unname(max_col_indices)

  treat = factor(treat, levels = 1:n_arms,
               labels = c("Placebo", "Low", "Medium", "High")[1:n_arms])

  # Recruitment
  # time = sample(1:ceiling(N/rmonth), N, replace = T)
  time = rexp(N, rate=rmonth)

  # Output
  data = data.frame(y_0m=y[,1],y_6m=y[,2], y_12m=y[,3], treat=treat, recruit_time = time)

  return(data)
}

#sim_dataind (n_arms = 4, N=60,  mu_0m=c(0,0,0,0), mu_6m=c(1,1,1,1), mu_12m=c(1,2,3,4), sg=matrix(c(1,0,0,0,1,0,0,0,1),3), rmonth=1)
