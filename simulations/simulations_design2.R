##########################
# eWHORM simulations
# August 2023
# Marta Bofill Roig
########################## 

rm(list = ls()) 
# Remove Package
remove.packages("ewhorm")

# server
# setwd("~/GitHub/ewhorm_sim/simulations")
# install.packages('~/GitHub/ewhorm_sim/ewhorm_0.1.tar.gz',repos=NULL)

# local
setwd("C:/Users/mbofi/Dropbox/CeMSIIS/GitHub/ewhorm_sim/simulations")
# source("C:/Users/mbofi/Dropbox/CeMSIIS/GitHub/ewhorm_sim/R/get_hyp_mat.R")
# source("C:/Users/mbofi/Dropbox/CeMSIIS/GitHub/ewhorm_sim/R/get_max_col_index.R")
# source("C:/Users/mbofi/Dropbox/CeMSIIS/GitHub/ewhorm_sim/R/sim_trial.R")
# source("C:/Users/mbofi/Dropbox/CeMSIIS/GitHub/ewhorm_sim/R/sim_data.R")
install.packages('C:/Users/mbofi/Dropbox/CeMSIIS/GitHub/ewhorm_sim/ewhorm_0.1.tar.gz',repos=NULL)
library(ewhorm)

# packges needed for this script
library(future) 
library(purrr)
library(furrr) 
# underlying dependencies
require(mvtnorm)#sim_data function 
require(multcomp)#aux functions
require(gtools)#aux functions


##########################################################
##########################################################


mu = c(0,0,0,0)
sg_m=matrix(c(1,.9,.9,1),nrow=2,byrow = T)

db <- ewhorm::sim_data(n_arms = 4,
                       N = 30 * 4,
                       mu_6m = mu,
                       mu_12m= mu,
                       sigma=diag(1,2),
                       rmonth =12)

summary(db)


##########################################################
##########################################################
set.seed(32)
# evaluate trial duration with respect to the rmonth, also assumptions regarding the break between stages
v=c(0,0,0)

sim_trial_pce(n_arms=4, N1=30*4, N2=30*2, mu_6m=mu, mu_12m=mu, sigma=sg_m, rmonth=2, alpha1=0.1, alpha=0.05) 
sim_trial_pce(n_arms=4, N1=30*4, N2=30*2, mu_6m=mu, mu_12m=mu, sigma=sg_m, rmonth=2, 
              alpha1=0.1, alpha=0.05, v=v)

sim_trial_pce(n_arms=4, N1=30*4, N2=30*2, mu_6m=mu, mu_12m=mu, sigma=sg_m, rmonth=2, 
              alpha1=0.1, alpha=0.05, v=v, sim_out=T)


##########################################################
##########################################################

# Run simulations
# Set a specific seed for reproducibility
set.seed(5267)  
# Set the number of trials to run and other parameters
n_trials <- 100000
# n_trials <- 100
# n_cores <- detectCores()-1  # Adjust the number of cores based on your machine's capabilities
n_cores <- availableCores()-1
# Set up the "multicore" future plan
# plan(multicore, workers = n_cores)
plan(multisession, workers = n_cores)

# Run the simulations in parallel using future_map
results_list <- future_map(1:n_trials, function(i) sim_trial_pce(n_arms=4, N1=30*4, N2=30*2, mu_6m=mu, mu_12m=mu, sigma=sg_m, rmonth=2, alpha1=0.1, alpha=0.05, sim_out=T), .options=furrr_options(seed = TRUE))


# Using simplified output
res_stage2 = matrix(unlist(lapply(results_list, function(element) element$stage2_arms)),
                      ncol = 3, byrow = T) 
c(sum(res_stage2[,1]), sum(res_stage2[,1]), sum(res_stage2[,3]))/n_trials
# 

res_decision = matrix(unlist(lapply(results_list, function(element) element$simdec_output)),
                      ncol = 3, byrow = T)  
res_decision4 = ifelse(res_decision[,1]+res_decision[,2]+res_decision[,3] >= 1 
                          & is.na(res_decision[,1]+res_decision[,2]+res_decision[,3])==F,
                          1, 0) 
res_decision <- cbind(res_decision, res_decision4) 

sum(res_decision[,4])/n_trials
