
#Libraries
##########

library(dplyr) # alternatively tidyverse
library(tidyr) # alternatively tidyverse
library(future)
library(furrr)
library(gMCP)
library(mvtnorm)


#FUNCTION 1: Do the calculation to obtain the mu and the sigma
##############################################################

get_mu_sigma = function(baseline_mean = 650, baseline_sd = 575, r0 = 0.15, r1 = 0.6, r2 = 0.8, 
                        r3 = 0.9,  rho = 0.5, N = 150, n1 = 90){

n2 = N - n1

mn = baseline_mean
sd = baseline_sd 


########################################
# PART 1: defining parameters
########################################

# Mean and standard deviation at baseline
#########################################

#mn <- 650; sd <- 575

#Reduction rate
###############

#r0 <- 0.15; r1 <- 0.6; r2 <- 0.8; r3 <- 0.9;
reduct_rate <- as.numeric(c(r0, r1, r2, r3))

#Correlation between the difference 
###################################


#The mean after six months from the common baseline mn
######################################################

dtmn <- as.data.frame(reduct_rate)

mean_6m <- data.frame(apply(dtmn, 2, function(r){ (1 - r)*mn}))
names(mean_6m) <- "mean_6month"
sd0 <- sd #we assumed the same sd for all the group at the baseline only

# calculate the mean for the log transformation for the baseline
################################################################

moy0 <- log(mn^2/sqrt(mn^2 + sd0^2)) # mean(log(x0))  Baseline


# calculating sd for the log after six month
############################################

sde0 <- sqrt(log(1+(sd0^2/mn^2) )) # sd(log(x0))  Baseline
sd_6m <- data.frame(apply(mean_6m, 2, function(m){ m*sqrt(exp(sde0) - 1)}))
names(sd_6m) <- "sd_6month"


#Twelve month : THIS IS THE SAME CODE TO SIX MONTH SO I CAN SKIP THE CODE 
###############AT THIS POINT UNLESS WE HAVE NEW VALUES FOR THE REDUCTION
###############RATE AND SD
###############

#Assume also that the decrease remained the same
################################################

#Start with mean of original scale
##################################

mean_12m <- mean_6m
names(mean_12m) <- "mean_12month"

#Standard deviation of the log(X)
#################################

sd_12m <- sd_6m


##############################################
# PART 2 : Setting the all for the simulation
##############################################

#Let's start with the mean for the log: Combine mean and sd in the same data.fr

mean_sd_6m <- cbind(mn = mean_6m, sd = sd_6m)

#Please write a function to calculate directly the mean log
###########################################################

mslg <- function(mean_6month, sd_6month){
  
  moy <- log(mean_6month^2/sqrt(mean_6month^2 + sd_6month^2))
  
  return(moy)
} #good!

#Now apply it
#############

mulog <- apply(mean_sd_6m, 1, function(m){ 
  do.call(mslg, as.list(m))
})

#Obtain the vector mu for each time point
#########################################

mu0 <-  rep(moy0, 4) # since at the baseline they are the same
mu6 <- mulog # for each arm
mu12 <- mulog



# Variance-covariance matrix
############################

#1. Elements of the matrix


varb <- 1*sde0^2 # cov(log(X0)
var6m <- 1*sde0^2 # cov(log(X06m) 
var12m <- 1*sde0^2 # cov(log(X12m)

cov_X0_X0 <- varb
cov_X0_X6m <- rho*sde0^2 # cov(log(X0), log(X6m)) 
cov_X0_X12m <- rho^2*sde0^2 #cov(log(X0), log(X12m)) 
cov_X6_X0 <- rho*sde0^2 # cov(log(X6), log(X0m)) 
cov_X6_X6 <- var6m
cov_X6_X12m <- rho*sde0^2 #cov(log(X6), log(X12m)) 
cov_X12m_X0 <- rho^2*sde0^2 #cov(log(X12m), log(X0)) 
cov_X12m_X6 <- rho*sde0^2 #cov(log(X12m), log(X6m)) 
cov_X12m_X12m <- var12m  #cov(log(X12m), log(X12m)) 


#2. Matrix: note the correlation matrix but the variance covariance matrix


sg <- matrix(c(varb, cov_X0_X6m, cov_X0_X12m, cov_X6_X0, var6m, cov_X6_X12m,
               cov_X12m_X0, cov_X12m_X6, cov_X12m_X12m), ncol=3)



return(list(mn, mu0, mu6, mu12, sg))
}


#FUNCTION 2: Using mu and sigma obtained above then run the function below
##########################################################################

do_pce_baseline = function(mu0, mu6, mu12, sg_m, n1, n2, alpha, alpha1, rmonth = 11){
  
  
  #Set the number of trials to run and other parameters for future plan
  n_trials <- 10000
  n_cores <- availableCores()-1 
  plan(multisession, workers = n_cores)
  
  # Run the simulations in parallel using future_map
  
  
  # First case alpha = 0.07
  
  alp = alpha
  bdr = alpha1
  
  results_list <- future_map(1:n_trials, function(i) sim_trial_pceind(n_arms = 4, 
                                N1 = n1 , N2 = n2, mu_0m = mu0,   mu_6m = mu6, 
                                mu_12m = mu12,  sg = sg_m, 
                                rmonth, alpha1 = bdr, alpha = alp,
                                v = c(1/2, 1/2, 0), sim_out=T), 
                             .options=furrr_options(seed = TRUE))
  
  
  
  # Summary results
  res_stage2 <- matrix(unlist(lapply(results_list, function(element) element$stage2_arms)),
                       ncol = 3, byrow = T) 
  
  # percentage of cases in which doses 1, 2 and 3 were present in the second stage
  armsel1 <- c(sum(res_stage2[,1]), sum(res_stage2[,2]), sum(res_stage2[,3]))/n_trials
  
  
  acc <- matrix(unlist(lapply(results_list, function(element) element$simdec_output)),
                ncol = 3, byrow = T) 
  
  armacc1 <- c(sum(acc[,1], na.rm = T), sum(acc[,2], na.rm = T), sum(acc[,3]))/n_trials
  
  
  
  return(list(armsel1, armacc1))
}


#FUNCTION 3: function 1 AND function 2
########################################


simul_res = function(baseline_mean, baseline_sd,  r0, r1, r2, r3,  rho, N, n1, alpha, alpha1){

#Maybe this is not optimal but
  
sz = as.numeric(N)
n1 = as.numeric(n1)
r0 = as.numeric(r0)


alp = alpha
bdr = alpha1

#Compute mu et sigma
####################

val = get_mu_sigma(baseline_mean,  baseline_sd,  r0, r1, r2, r3,  rho, N , n1) # rho is the correlation coefficient


# extract the matrix and the mu for each time
#############################################

mtr1 = matrix(unlist(val[5]), nrow=3,byrow = T) #matrix I
m0 = as.numeric(unlist(val[2])) #mu0
m6 = as.numeric(unlist(val[3])) #mu6m
m12 = as.numeric(unlist(val[4])) #mu12m

#do pce and also the simulation with futur_map
##############################################

aa = do_pce_baseline(mu0 = m0, mu6 = m6, mu12 = m12, sg_m = mtr1,
                     n1 = n1, n2 = N - n1, alpha = alp, alpha1 = bdr)


return(aa)

}


# FUNCTION 4: Do the simulation per scenario given a table x
###########################################################


sim_it = function(x){
  
  #Do the simulation for each scenario in your table
  ##################################################
  
  
  result_list <- apply(x, 1, function(row) {
    do.call(simul_res, as.list(row))
  })
  
  b = length(unique(x$rho))
  
  # Now extract element in each list (scenarios) wich corresponds to proportion of selection
  ##########################################################################################
  
  
  res1 = data.frame(t(sapply(lapply(result_list, function(x) x[[1]]), #extract the first element in each list: proportion of sel
                             function(x) unlist(x)))) %>% # each element beacomes a column in a data.frame
    mutate(Hypothesis = paste0(rep(unique(x[, "alpha1"]), b), "- rho = ", x$rho) )#, # 
  
  scen = c("H1", "H2", "H3") # Hypothesis: Low, medium, High
  names(res1)[1:3] = scen
  
  
  res2 = data.frame(t(sapply(lapply(result_list, function(x) x[[2]]), #extract the second element in each list: proportion of sel
                             function(x) unlist(x)))) %>% 
    mutate(Hypothesis = paste0(rep(unique(x[, "alpha1"]), b), "- rho = ", x$rho) )#
  
  
  names(res2)[1:3] = scen
  
  all_scenario <- list(res1, res2)
  
  return(all_scenario)
  
}


