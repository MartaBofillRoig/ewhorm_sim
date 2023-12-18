


# Importing packages
####################

library(dplyr)
library(mvtnorm)
library(future)
library(furrr)
library(gMCP)



########################################
# PART 1: defining parameters
########################################



# Mean and standard deviation at baseline
#########################################

mn <- 650; sd <- 575

#Reduction rate
###############

r0 <- 0.15; r1 <- 0.6; r2 <- 0.8; r3 <- 0.9;

#Correlation between the difference, this will change to be between the outcome for each time points
####################################################################################################

rho <- 0.5

#The mean after six months from the common baseline mn
######################################################

mn0 <- (1 - r0)*mn
mn1 <- (1 - r1)*mn
mn2 <- (1 - r2)*mn
mn3 <- (1 - r3)*mn

sd0 <- sd #we assumed the same sd for all the group at the baseline only

moy0 <- log(mn^2/sqrt(mn^2 + sd0^2)) # mean(log(x0))  Baseline
sde0 <- sqrt(log(1+(sd0^2/mn^2) )) # sd(log(x0))  Baseline


#As we assumed that the estimated sigma for the log is the same for every other, 
#we can deduce the sdi for each of the other group, to estimate a good mean for the log

sd06m <- mn0*sqrt(exp(sde0) - 1) # sd((x6m)) 
sd16m <- mn1*sqrt(exp(sde0) - 1) # sd((x6m)) LOW D. for twelve months
sd26m <- mn2*sqrt(exp(sde0) - 1) # sd((x6m)) Medium D. for twelve months
sd36m <- mn3*sqrt(exp(sde0) - 1) # sd((x6m)) HIGH for twelve months


#Twelve month
#############

#Assume also that the decrease remained the same
################################################

#Start with mean of original scale
##################################

mn0 <- (1 - r0)*mn
mn1 <- (1 - r1)*mn
mn2 <- (1 - r2)*mn
mn3 <- (1 - r3)*mn


sd0 <- sd #we assumed the same sd for all the group at the baseline only

moy0 <- log(mn^2/sqrt(mn^2 + sd0^2)) # mean(log(x0))  BASELINE
sde0 <- sqrt(log(1+(sd0^2/mn^2) )) # sd(log(x0))  Baseline

#Standard deviation of the log(X)
#################################

sd012m <- mn0*sqrt(exp(sde0) - 1) # sd((x12m)) 
sd112m <- mn1*sqrt(exp(sde0) - 1) # sd((x12m)) LOW D. for twelve months
sd212m <- mn2*sqrt(exp(sde0) - 1) # sd((x12m)) Medium D. for twelve months
sd312m <- mn3*sqrt(exp(sde0) - 1) # sd((x12m)) HIGH for twelve months


#Now let's calculate the mean for the log
#########################################

moy1 <- log(mn0^2/sqrt(mn0^2 + sd012m^2)) # mean(log(x0)) PLACEBO for twelve months
moy2 <- log(mn1^2/sqrt(mn1^2+sd112m^2)) # mean(log(x6m)) LOW D. for twelve months
moy3 <-  log(mn2^2/sqrt(mn2^2+sd212m^2))# mean(log(x6m)) Medium D. for twelve months
moy4 <-  log(mn3^2/sqrt(mn3^2+sd312m^2))# mean(log(x6m)) HIGH for twelve months

mu0 <-  c(moy0, moy0, moy0, moy0)
mu6 <- c(moy1, moy2, moy3, moy4)
mu12 <- c(moy1, moy2, moy3, moy4)




#Variance of the diff of the logs
#################################

rho <- 0.5; vardiff <- sde0^2  #sde1^2 + sde2^2 - 2* rho*sde1*sde2


# Variance-covariance matrix
############################

#1. Elements of the matrix

varb <- 1*sde0^2 # cov(log(X0), log(X0)) <- rho*sd0*sd0
var6m <- 1*sde0^2 # cov(log(X06m), log(X06m)) <- rho*sdx6m*sdx6m
var12m <- 1*sde0^2 # cov(log(X12m), log(X12m)) <- rho*sdx12m*sdx12m


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

########################################
# PART 2: TESTING THE TWO NEW FUNCTIONS
########################################

## A: testing sim_dataind

# simulate the data
###################

rb <- sim_dataind(n_arms = 4, N = 150, mu_0m = mu0, mu_6m = mu6, mu_12m = mu12, sg, rmonth)


# Check on the values you get: mean sd per group
################################################

rb %>% group_by(treat) %>% 
  summarize(baseline_mean_sd = paste0(round(mean(y_0m)), " (", round(sd(y_0m), 2),")"),
            sixmont_mean_sd = paste0(round(mean(y_6m)), " (", round(sd(y_6m), 2),")"),
            twlvmont_mean_sd = paste0(round(mean(y_12m)), " (", round(sd(y_12m), 2),")"))



# Let's do the simulation using the new function for pce 
########################################################

sim_trial_pceind(n_arms = 4, N1 , N2, mu_0m = c(moy0, moy0, moy0, moy0), 
                             mu_6m = c(moy1, moy2, moy3, moy4), mu_12m = c(moy1, moy2, moy3, moy4),
                             sg, rmonth, alpha1 = 0.1, alpha = 0.025,v = c(1/2,1/2,0), sim_out=T)



set.seed(20)

#Set the number of trials to run and other parameters for future plan
n_trials <- 10000
n_cores <- availableCores()-1 
plan(multisession, workers = n_cores)

# Run the simulations in parallel using future_map

# First case alpha = 0.1

results_list <- future_map(1:n_trials, function(i) sim_trial_pceind(n_arms = 4, N1 , N2, mu_0m = c(moy0, 
                                                                    moy0, moy0, moy0), 
                                                                    mu_6m = c(moy1, moy2, moy3, moy4), 
                                                                    mu_12m = c(moy1, moy2, moy3, moy4),
                                                                    sg, rmonth, alpha1 = 0.1, alpha = 0.025,
                                                                    v = c(1/2,1/2,0), sim_out=T), .options=furrr_options(seed = TRUE))



# Summary results
res_stage2 <- matrix(unlist(lapply(results_list, function(element) element$stage2_arms)),
                    ncol = 3, byrow = T) 

# percentage of cases in which doses 1, 2 and 3 were present in the second stage
armsel1 <- c(sum(res_stage2[,1]), sum(res_stage2[,2]), sum(res_stage2[,3]))/n_trials


acc <- matrix(unlist(lapply(results_list, function(element) element$simdec_output)),
             ncol = 3, byrow = T) 

armacc1 <- c(sum(acc[,1], na.rm = T), sum(acc[,2], na.rm = T), sum(acc[,3]))/n_trials




# Old function based on the difference only
###########################################


# Then mean(log(xi) - log(x0))

mu1 <-  moy0 - moy1
mu2 <-  moy0 - moy2
mu3 <-  moy0 - moy3
mu4 <-  moy0 - moy4

# The mu for the simulation: this is about the effect in log
############################################################


mu <- c(mu1, mu2, mu3, mu4)


#Run it here
############

results_list2 <- future_map(1:n_trials, function(i) sim_trial_pce(n_arms=4, N1=n1, N2=n2,
                                                                 mu_6m=mu, mu_12m=mu, sigma=sg_m, 
                                                                 rmonth=2, alpha1=0.1, alpha=0.025, 
                                                                 sim_out=T), .options=furrr_options(seed = TRUE))



# Summary results

res_stage2_old <- matrix(unlist(lapply(results_list2, 
                            function(element) element$stage2_arms)),
                                       ncol = 3, byrow = T) 

# percentage of cases in which doses 1, 2 and 3 were present in the second stage
armsel2 <- c(sum(res_stage2_old[,1]), sum(res_stage2_old[,2]), 
             sum(res_stage2_old[,3]))/n_trials


# percentage of cases in which doses 1, 2 and 3 were rejected at the final analsis (second stage)

acc <- matrix(unlist(lapply(results_list2, 
                           function(element) element$simdec_output)),
                                      ncol = 3, byrow = T) 

armacc2 <- c(sum(acc[,1], na.rm = T), sum(acc[,2], na.rm = T), sum(acc[,3], na.rm = T))/n_trials



# Comparison of the output
##########################

#Proportion of selected

armsel1

armsel2

#Proportion of rejected

armacc1
armacc2



