##########################
# eWHORM simulations
# August 2023
# Marta Bofill Roig
########################## 

install.packages('C:/Users/mbofi/Dropbox/CeMSIIS/GitHub/ewhorm_sim/ewhorm_0.1.tar.gz',repos=NULL)
library(ewhorm)

# packges needed for this script
library(future) 
library(furrr) 

# Set a specific seed for reproducibility
set.seed(522)

##########################################################
##########################################################
# Settings under the global null

# Means in each group
mu = c(0,0,0,0)
# Matrix var-cov between 6- and 12- months
sg_m=matrix(c(1,.7,.7,1),nrow=2,byrow = T)

# Simulation trial (1 stage only)
db <- ewhorm::sim_data(n_arms = 4,
                       N = 30 * 4,
                       mu_6m = mu,
                       mu_12m= mu,
                       sigma = diag(1,2),
                       rmonth = 12)

summary(db)

##########################################################
##########################################################

# Simulation and analysis 
sim_trial_pce(n_arms=4, N1=30*4, N2=30*2, 
              mu_6m=mu, mu_12m=mu, sigma=sg_m, 
              rmonth=12, alpha1=0.1, alpha=0.05) 

# Simplified output for simulations
sim_trial_pce(n_arms=4, N1=30*4, N2=30*2, 
              mu_6m=mu, mu_12m=mu, sigma=sg_m, 
              rmonth=12, alpha1=0.1, alpha=0.05,  
              sim_out=T) 
# $stage2_arms: arms included in stage 2
# $simdec_output: decision for each treatment arm (1: rejected; NA in arm 3 if the highest dose was not included in stage2)

##########################################################
##########################################################
# Run simulations

# Set the number of trials to run and other parameters for future plan
n_trials <- 100000
n_cores <- availableCores()-1 
plan(multisession, workers = n_cores)

v0=c(0,0,0)

# Run the simulations in parallel using future_map
results_list <- future_map(1:n_trials, function(i) sim_trial_pce(n_arms=4, N1=30*4, N2=30*2, mu_6m=mu, mu_12m=mu, sigma=sg_m, rmonth=2, alpha1=0.000001, alpha=0.05, v=v0, sim_out=T), .options=furrr_options(seed = TRUE))


# Summary results
res_stage2 = matrix(unlist(lapply(results_list, function(element) element$stage2_arms)),
                    ncol = 3, byrow = T) 

# percentage of cases in which doses 1, 2 and 3 were present in the second stage
c(sum(res_stage2[,1]), sum(res_stage2[,1]), sum(res_stage2[,3]))/n_trials

# Type 1 error rates
res_decision = matrix(unlist(lapply(results_list, function(element) element$simdec_output)),
                      ncol = 3, byrow = T)  
colSums(res_decision, na.rm = T)/n_trials

# Global null
res_decision4 = ifelse(res_decision[,1]+res_decision[,2]+res_decision[,3] >= 1 
                       & is.na(res_decision[,1]+res_decision[,2]+res_decision[,3])==F,
                       1, 0) 
res_decision <- cbind(res_decision, res_decision4) 
sum(res_decision[,4])/n_trials 

