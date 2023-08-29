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
sim_trial_pce(n_arms=4, N1=30*4, N2=30*2, mu_6m=mu, mu_12m=mu, sigma=sg_m, rmonth=2, alpha1=0.1, alpha=0.05)

sim_trial_pce(n_arms=4, N1=30*4, N2=30*2, mu_6m=mu, mu_12m=mu, sigma=diag(0.5,2), alpha1=0.5, alpha=0.05,rmonth = 12)
