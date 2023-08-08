#' Simulate data from a multi-arm trial with shared control
#' @description Function to simulate trial data (1-stage, multiple arms)
#'
#' @param n_arms number of arms (including control)
#' @param N total sample size
#' @param mu_6m 6-month mean response per arm (vector of length `n_arm`)
#' @param mu_12m 12-month mean response per arm (vector of length `n_arm`)
#' @param sigma covariance matrix between 6- and 12-month responses assumed equal across arms (matrix of dim 2x2)
#' @param rmonth recruitment per month (recruitment speed assumed constant over time)
#' @returns simulated data consisting of the responses at 6 and 12 months, treatment arm, and recruitment time for each subject.
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats model.matrix
#' @export
#' @details eWHORM simulations
#' @author Marta Bofill Roig

sim_data <- function(n_arms, N, mu_6m, mu_12m, sigma, rmonth){

  treatments <- factor(c(sample(rep(1:n_arms, floor(N/n_arms))), sample(1:n_arms, N-floor(N/n_arms)*n_arms, replace=T)),
                      # sample(1:n_arms, N, replace = TRUE),
                       levels = 1:n_arms,
                       labels = c("Placebo", "Low", "Medium", "High")[1:n_arms])
  X <- model.matrix(~ treatments - 1)
  y <- X %*% matrix(c(mu_6m, mu_12m), nrow=n_arms, byrow = F) + rmvnorm(n=N, mean = c(0,0), sigma = sigma )

  # Treatment indicator for dataframe
  max_col_indices <- apply(X, 1, get_max_col_index)
  treat = unname(max_col_indices)

  treat = factor(treat, levels = 1:n_arms,
               labels = c("Placebo", "Low", "Medium", "High")[1:n_arms])

  # Recruitment
  time = sample(1:ceiling(N/rmonth), N, replace = T)

  # Output
  data = data.frame(y_6m=y[,1], y_12m=y[,2], treat=treat, recruit_time = time)

  return(data)
}
