
#Importing libraries
####################

library(dplyr)
library(tidyr)
library(future)
library(furrr)
library(gMCP)
library(mvtnorm)

#Call the function and run from where they are
##############################################

source("Simul_from_table/all_functions.R")

# Now let's apply it to all the scenarios
#########################################

# The correlation

rh <- 0.5

rh2 <- 0.7

# Define manually the dataset
#############################


#case 1:
########

cs1 <- data.frame(baseline_mean = 650, baseline_sd = 575, r0 = rep(0.15, 3), 
                 r1 = rep(0.65, 3), r2 = rep(0.75, 3), r3 = rep(0.90, 3), 
                 rho = rep(rh, 3),
                 N = rep(150, 3), n1 = rep(90, 3), alpha = rep(0.025, 3), 
                 alpha1 = c(0.1, 0.2, 0.5))

cs0 <- data.frame(baseline_mean = 650, baseline_sd = 575, r0 = rep(0.15, 3), 
                 r1 = rep(0.65, 3), r2 = rep(0.75, 3), r3 = rep(0.90, 3), 
                 rho = rep(rh2, 3),
                 N = rep(150, 3), n1 = rep(90, 3), alpha = rep(0.025, 3), 
                 alpha1 = c(0.1, 0.2, 0.5))

cs <- rbind(cs1, cs0) %>% 
  mutate_at(vars(baseline_mean:alpha1), list(~as.numeric(.)))

# Do the simulation

s10 <- sim_it(cs)



#Case 2: 
########

cs2 <- data.frame(baseline_mean = 650, baseline_sd = 575, r0 = rep(0.15, 3), 
                 r1 = rep(0.7, 3), r2 = rep(0.65, 3), r3 = rep(0.90, 3), 
                 rho = rep(rh, 3),
                 N = rep(150, 3), n1 = rep(90, 3), alpha = rep(0.025, 3),
                 alpha1 = c(0.1, 0.2, 0.5) )

cs20 <- data.frame(baseline_mean = 650, baseline_sd = 575, r0 = rep(0.15, 3), 
                  r1 = rep(0.7, 3), r2 = rep(0.65, 3), r3 = rep(0.90, 3), 
                  rho = rep(rh2, 3),
                  N = rep(150, 3), n1 = rep(90, 3), alpha = rep(0.025, 3),
                  alpha1 = c(0.1, 0.2, 0.5) )

cs02 <- rbind(cs2, cs20)


# Do the simulation
###################

s20 <- sim_it(cs02)



# Functions exist to handle this output: need a discussion on the shape


# Please add the other parameters like: test and scenario when low dose is promising
