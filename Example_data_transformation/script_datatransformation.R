
# Importing packages
####################

library(dplyr)
library(mvtnorm)


########################################
# PART 1: defining parameters
########################################

# Mean and standard deviation at baseline
#########################################

mn <- 650; sd <- 575

#Reduction rate
###############

r0 <- 0.15; r1 <- 0.6; r2 <- 0.8; r3 <- 0.9;
reduct_rate <- c(r0, r1, r2, r3)

#Correlation between the difference 
###################################

rho <- 0.5

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

rho <- 0.5

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

