##########################
# eWHORM simulations
# July 2023
# Marta Bofill Roig
##########################

rm(list = ls())

setwd("C:/Users/mbofi/Dropbox/CeMSIIS/GitHub/ewhorm_sim/2023-07-simulations")
source("C:/Users/mbofi/Dropbox/CeMSIIS/GitHub/ewhorm_sim/R/aux-functions.R")
source("C:/Users/mbofi/Dropbox/CeMSIIS/GitHub/ewhorm_sim/R/sim_trial.R")
source("C:/Users/mbofi/Dropbox/CeMSIIS/GitHub/ewhorm_sim/R/sim_data.R")

# packges needed for this script
library(future) 
library(purrr)
library(furrr)
# library(parallel)

# Settings 
n_arms=4; N1=30*4; N2=30*2; sd_y=0.1; alpha1=0.5
# mu=c(0,1,2,5); 
mu=c(0,0,0,0);
mu_6m=mu; mu_12m=mu

# underlying dependencies
require(mvtnorm)#sim_data function


##########################################################
##########################################################
# evaluate trial duration with respect to the rmonth, also assumptions regarding the break between stages
mu=c(0,0,0,0); sigma=matrix(c(0.1,0,0,0.1), nrow = 2, byrow = T)
y=sim_data(n_arms=4, N=30*4, mu_6m=mu, mu_12m=mu+c(0,1,1,2), sigma=sigma, rmonth=10) 
y
summary(y$recruit_time)

##########################################################
##########################################################
# Run simulations
# Set a specific seed for reproducibility
set.seed(421)  
# Set the number of trials to run and other parameters
n_trials <- 100000 
# n_cores <- detectCores()-1  # Adjust the number of cores based on your machine's capabilities
n_cores <- availableCores()-1
# Set up the "multicore" future plan
# plan(multicore, workers = n_cores)
plan(multisession, workers = n_cores)

# Run the simulations in parallel using future_map
results_list <- future_map(1:n_trials, function(i) sim_trial(n_arms=4, N1=30*4, N2=30*2, mu_6m=mu, mu_12m=mu, sd_y=0.1, alpha1=0.5, alpha=0.05), .options=furrr_options(seed = TRUE))

# Extract the two sets of results from the list
result1_values <- sapply(results_list, function(x) x$result1)
result2_values <- sapply(results_list, function(x) x$result2)

# Calculate the means
mean_result1 <- mean(result1_values)
summary_result2 <- table(as.factor(result2_values))

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


