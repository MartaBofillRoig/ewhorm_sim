##########################
# eWHORM simulations
# July 2023
# Marta Bofill Roig
########################## 

rm(list = ls()) 
# Remove Package
# remove.packages("ewhorm") 

# server
setwd("~/GitHub/ewhorm_sim/simulations")
install.packages('~/GitHub/ewhorm_sim/ewhorm_0.1.tar.gz',repos=NULL)

# local
# setwd("C:/Users/mbofi/Dropbox/CeMSIIS/GitHub/ewhorm_sim/simulations")
# source("C:/Users/mbofi/Dropbox/CeMSIIS/GitHub/ewhorm_sim/R/get_hyp_mat.R")
# source("C:/Users/mbofi/Dropbox/CeMSIIS/GitHub/ewhorm_sim/R/get_max_col_index.R")
# source("C:/Users/mbofi/Dropbox/CeMSIIS/GitHub/ewhorm_sim/R/sim_trial.R")
# source("C:/Users/mbofi/Dropbox/CeMSIIS/GitHub/ewhorm_sim/R/sim_data.R")
# install.packages('C:/Users/mbofi/Dropbox/CeMSIIS/GitHub/ewhorm_sim/ewhorm_0.1.tar.gz',repos=NULL)
library(ewhorm)

# packges needed for this script
library(future) 
library(purrr)
library(furrr) 
# underlying dependencies
require(mvtnorm)#sim_data function 
require(multcomp)#aux functions
require(gtools)#aux functions

# Example
# n_arms=4
# N1=120
# N2=60
# mu = c(0,0,0,0)
# mu_6m = mu
# mu_12m= mu
# sg_m=matrix(c(1,.9,.9,1),nrow=2,byrow = T)
# rmonth=12
# alpha1=0.5
# alpha=0.05
# p_safety=c(0.9,0.8,0.7)
# safety=T
# promising=T

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
# evaluate trial duration with respect to the rmonth, also assumptions regarding the break between stages

sim_trial(n_arms=4, N1=30*4, N2=30*2, mu_6m=mu, mu_12m=mu, sigma=diag(0.5,2), alpha1=0.5, alpha=0.05,rmonth = 12)

##########################################################
##########################################################
# Run simulations
# Set a specific seed for reproducibility
set.seed(55)  
# Set the number of trials to run and other parameters
n_trials <- 100000
# 00 
# n_cores <- detectCores()-1  # Adjust the number of cores based on your machine's capabilities
n_cores <- availableCores()-1
# Set up the "multicore" future plan
# plan(multicore, workers = n_cores)
plan(multisession, workers = n_cores)

# Run the simulations in parallel using future_map
results_list <- future_map(1:n_trials, function(i) sim_trial(n_arms=4, N1=30*4, N2=30*2, mu_6m=mu, mu_12m=mu, sigma=sg_m, alpha1=.5, alpha=0.05, rmonth = 12), .options=furrr_options(seed = TRUE))

# Extract the two sets of results from the list
combined_pvalue_values <- sapply(results_list, function(x) x$combined_pvalue)
selected_dose_values <- sapply(results_list, function(x) x$selected_dose)
safety_values <- sapply(results_list, function(x) x$safety)
pvalue_stage1_values <- sapply(results_list, function(x) x$pvalue_stage1)
pvalue_stage2_values <- sapply(results_list, function(x) x$pvalue_stage2)

recruit_times1_values <- sapply(results_list, function(x) x$recruit_time1)
recruit_times2_values <- sapply(results_list, function(x) x$recruit_time2)

summary(recruit_times1_values)
summary(recruit_times2_values)
# Calculate the means

summary(combined_pvalue_values)
summary(pvalue_stage1_values)
summary(pvalue_stage2_values)

mean_result1 <- sum((combined_pvalue_values<0.05))/length(combined_pvalue_values)
summary_result2 <- table(as.factor(selected_dose_values))
mean_safety <- mean(safety_values)


# Print the means
cat("Type 1 error:", mean_result1, "\n")
cat("Selected dose:", summary_result2, "\n") 
cat("Safety selected dose:", mean_safety, "\n")

##########################################################
##########################################################
# Run simulations
# Set a specific seed for reproducibility
set.seed(4235)  
# Set the number of trials to run and other parameters
n_trials <- 100000
# 00 
# n_cores <- detectCores()-1  # Adjust the number of cores based on your machine's capabilities
n_cores <- availableCores()-1
# Set up the "multicore" future plan
# plan(multicore, workers = n_cores)
plan(multisession, workers = n_cores)


sg_m=matrix(c(1,1,1,1),nrow=2,byrow = T)
mu1 = c(0,0,0,1.3)
# Run the simulations in parallel using future_map
results_list <- future_map(1:n_trials, function(i) sim_trial(n_arms=4, N1=30*4, N2=30*2, mu_6m=mu1, mu_12m=mu1, sigma=sg_m, alpha1=.1, alpha=0.05,rmonth = 12, promising = T), .options=furrr_options(seed = TRUE))

# Extract the two sets of results from the list
combined_pvalue_values <- sapply(results_list, function(x) x$combined_pvalue)
selected_dose_values <- sapply(results_list, function(x) x$selected_dose) 
pvalue_stage1_values <- sapply(results_list, function(x) x$pvalue_stage1)
pvalue_stage2_values <- sapply(results_list, function(x) x$pvalue_stage2) 

mean_result1 <- sum((combined_pvalue_values<0.05))/length(combined_pvalue_values)
summary_result2 <- table(as.factor(selected_dose_values)) 


# Print the means
cat("Type 1 error:", mean_result1, "\n")
cat("Selected dose:", summary_result2, "\n")  


##########################################################
##########################################################

set <- list(
  scenario1 = c(n_arms=4, N1=30*4, N2=30*2, mu_6m=mu, mu_12m=mu, sd_y=0.1, alpha1=0.5, alpha=0.05),
  scenario2 = c(n_arms=4, N1=30*4, N2=30*2, mu_6m=mu, mu_12m=mu, sd_y=0.1, alpha1=0.5, alpha=0.05)
)

# sim_trial(c(n_arms=4, N1=30*4, N2=30*2, mu_6m=mu, mu_12m=mu, sd_y=0.1, alpha1=0.5, alpha=0.05))
results <- future_map(set, ~ future_map2(replicate(5, .x), .y, sim_trial))


