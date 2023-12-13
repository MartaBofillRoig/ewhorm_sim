

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
  time = sample(1:ceiling(N/rmonth), N, replace = T)
  
  # Output
  data = data.frame(y_0m=y[,1],y_6m=y[,2], y_12m=y[,3], treat=treat, recruit_time = time)
  
  return(data)
}


#Set values for testing: No rational just for them
##################################################


rmonth <- 12
n_arms <- 4
rho <- 0.5
sd0 <- 575
sdx6m <- 300
sdx12m <- 250

# Mean values
#############

mu0 <- c(650, 650, 650, 650)
mu6 <- c((1 - r0)*650, (1 - r1)*650, (1 - r2)*650, (1 - r3)*650)
mu12 <- c((1 - r0)*650, (1 - r1)*650, (1 - r2)*650, (1 - r3)*650)


# Variance-covariance matrix
############################

#1. Elements of the matrix

varb <- sd0*sd0 # cov(X0, X0) <- rho*sd0*sd0
var6m <- sdx6m*sdx6m # cov(X06m, X06m) <- rho*sdx6m*sdx6m
var12m <- sdx12m*sdx12m # cov(X12m, X12m) <- rho*sdx12m*sdx12m



cov_X0_X6m <- rho*sd0*sdx6m # cov(X0, X6m) <- rho*sdx0*sdx6m
cov_X0_X12m <- rho*sd0*sdx12m #cov(X0, X12m) <- rho*sdx0*sdx12m
cov_X6_X12m <- rho*sdx6m*sdx12m #cov(X6, X12m) <- rho*sdx6m*sdx12m


#2. Matrix


sg <- matrix(c(varb,cov_X0_X6m, cov_X0_X12m, cov_X0_X6m,var6m, cov_X6_X12m,
               cov_X0_X12m, cov_X6_X12m, var12m), ncol=3)


#Now call the function to simulate the data
###########################################


rb = sim_dataind(n_arms = 4, N = 150, mu_0m = mu0, mu_6m = mu6, mu_12m = mu12, sg, rmonth)



# Check on the values you get: mean sd per group
################################################

rb %>% group_by(treat) %>% 
       summarize(baseline_mean_sd = paste0(round(mean(y_0m)), " (", round(sd(y_0m)),")"),
                 sixmont_mean_sd = paste0(round(mean(y_6m)), " (", round(sd(y_6m)),")"),
                 twlvmont_mean_sd = paste0(round(mean(y_12m)), " (", round(sd(y_12m)),")"))










