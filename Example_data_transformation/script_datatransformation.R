\subsubsection{Code}

\begin{lstlisting}[language=R]
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


#Correlation between the difference 
###################################

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


#Twelve month : THIS IS THE SAME CODE TO SIX MONTH SO I CAN SKIP THE CODE AT THIS POINT UNLESS WE HAVE NEW VALUES FOR THE REDUCTION RATE AND SD
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


##############################################
# PART 2 : Setting the all for the simulation
##############################################

#Let's start with the mean for the log
######################################

moy1 <- log(mn0^2/sqrt(mn0^2 + sd012m^2)) # mean(log(x0)) PLACEBO for twelve months
moy2 <- log(mn1^2/sqrt(mn1^2+sd112m^2)) # mean(log(x6m)) LOW D. for twelve months
moy3 <-  log(mn2^2/sqrt(mn2^2+sd212m^2))# mean(log(x6m)) Medium D. for twelve months
moy4 <-  log(mn3^2/sqrt(mn3^2+sd312m^2))# mean(log(x6m)) HIGH for twelve months

mu0 <-  c(moy0, moy0, moy0, moy0)
mu6 <- c(moy1, moy2, moy3, moy4)
mu12 <- c(moy1, moy2, moy3, moy4)



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
cov_X12m_X0 <- rho*sde0^2 #cov(log(X12m), log(X0)) 
cov_X12m_X6 <- rho*sde0^2 #cov(log(X12m), log(X6m)) 
cov_X12m_X12m <- var12m  #cov(log(X12m), log(X12m)) 


#2. Matrix: note the correlation matrix but the variance covariance matrix


sg <- matrix(c(varb, cov_X0_X6m, cov_X0_X12m, cov_X6_X0, var6m, cov_X6_X12m,
               cov_X12m_X0, cov_X12m_X6, cov_X12m_X12m), ncol=3)

\end{lstlisting}