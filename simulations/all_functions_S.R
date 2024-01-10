
#Libraries
##########

library(dplyr) # alternatively tidyverse
library(tidyr) # alternatively tidyverse
library(future)
library(furrr)
library(gMCP)
library(mvtnorm)


#r0 <- 0.15; r1 <- 0.6; r2 <- 0.8; r3 <- 0.9;
#mu_raw_0 <- 650; sd_raw_0 <- 575;N = 150; n1 = 90
#get_mu_sigma(mu_raw_0 = 650, sd_raw_0 = 575, r0 = 0.15, r1 = 0.6, r2 = 0.8, r3 = 0.9,  rho = 0.5)

#FUNCTION 1: Do the calculation to obtain the mu and the sigma
##############################################################

get_mu_sigma = function(mu_raw_0 = 650, sd_raw_0 = 575, r0 = 0.15, r1 = 0.6, r2 = 0.8, 
                        r3 = 0.9,  rho = 0.5)
{

########################################
# PART 1: defining parameters
########################################

reduct_rate <- as.numeric(c(r0, r1, r2, r3)) #Reduction rate

mu_raw_6 <- (1-reduct_rate)*mu_raw_0 ##The mean after six months from the common baseline
mu_log_0 <- log(mu_raw_0^2/sqrt(mu_raw_0^2 + sd_raw_0^2))  # calculate the mean for the log transformation for the baseline

sd_log_0 <- sqrt(log(1+(sd_raw_0^2/mu_raw_0^2) )) # sd(log(x0))  Baseline
sd_raw_6 <- mu_raw_6*sqrt(exp(sd_log_0) - 1)   # calculating sd for the log after six month


#Twelve month : THIS IS THE SAME CODE As SIX MONTHs
mu_raw_12<-mu_raw_6
sd_raw_12<-sd_raw_6

##############################################
# PART 2 : Setting the all for the simulation
##############################################
#Obtain the vector mu for each time point
#########################################

mu_log_0<-rep(mu_log_0,4)
mu_log_6<-log(mu_raw_6^2/sqrt(mu_raw_6^2 + sd_raw_6^2))
mu_log_12<-log(mu_raw_12^2/sqrt(mu_raw_12^2 + sd_raw_12^2))


# Variance-covariance matrix
############################

#1. Elements of the matrix


var0 <- sd_log_0^2 # cov(log(X0)
var6 <- sd_log_0^2 # cov(log(X06m) 
var12<- sd_log_0^2 # cov(log(X12m)

#cov_X0_X0 <- var0
cov_X0_X6m <- rho*sd_log_0^2 # cov(log(X0), log(X6m)) 
cov_X0_X12m <- rho^2*sd_log_0^2 #cov(log(X0), log(X12m)) 
cov_X6_X0 <- rho*sd_log_0^2 # cov(log(X6), log(X0m)) 
#cov_X6_X6 <- var6
cov_X6_X12m <- rho*sd_log_0^2 #cov(log(X6), log(X12m)) 
cov_X12m_X0 <- rho^2*sd_log_0^2 #cov(log(X12m), log(X0)) 
cov_X12m_X6 <- rho*sd_log_0^2 #cov(log(X12m), log(X6m)) 
#cov_X12m_X12m <- var12  #cov(log(X12m), log(X12m)) 


#2. Matrix: not the correlation matrix but the variance covariance matrix
sg <- matrix(c(var0, cov_X0_X6m, cov_X0_X12m, cov_X6_X0, var6, cov_X6_X12m,cov_X12m_X0, cov_X12m_X6, var12), ncol=3)


return(list(mu_raw_0, c(mu_log_0), c(mu_log_6), c(mu_log_12), sg))
}


#FUNCTION 2: Using mu and sigma obtained above then run the function below
##########################################################################

# n_trials=1000; n_arms=4; N1=900;N2=60 mu_0m =c(10,10,10,10); mu_6m =c(10,10,10,10); mu_12m=c(10,10,10,10); sg=matrix(c(1,0,0,0,1,0,0,0,1), ncol=3); 
#rmonth=1;alpha1=.1; alpha=0.025;
#do_pce_baseline(n_trials=10000,n_arms = 4,N1 = 90 , N2 = 60, mu_0m = mu_0m,   mu_6m = mu_6m, mu_12m = mu_12m,  sg = sg, rmonth=1, alpha1 = .1, alpha = .025, v = c(1/2, 1/2, 0), sim_out=T,sel_scen=0, side=T,test="t")
  
do_pce_baseline = function(n_trials,n_arms = 4,N1 = N1 , N2 = N1, mu_0m = mu_0,   mu_6m = mu_6m, 
                           mu_12m = mu12m,  sg = sg, rmonth, alpha1 , alpha ,
                           v = c(1/2, 1/2, 0), sim_out=sim_out,sel_scen=sel_scen, side=side,test=test)
{

  #Set the number of trials to run and other parameters for future plan
  #n_trials <- 10000
  n_cores <- availableCores()-1 
  plan(multisession, workers = n_cores)
  
  # Run the simulations in parallel using future_map
  
  # First case alpha = 0.07
  
  #alp = alpha
  #bdr = alpha1
  
  results_list <- future_map(1:n_trials, function(i) sim_trial_pceind_test(n_arms = n_arms,N1 = N1 , N2 = N2, mu_0m = mu_0m,   mu_6m = mu_6m, 
                                                                          mu_12m = mu_12m,  sg = sg, rmonth=rmonth, alpha1=alpha1 , alpha =alpha,
                                                                          v = c(1/2, 1/2, 0), sim_out=sim_out,sel_scen=sel_scen, side=side,test=test), 
                                                                          .options=furrr_options(seed = TRUE))
  
   # Summary results
  res_stage2 <- matrix(unlist(lapply(results_list, function(element) element$stage2_arms)),ncol = 3, byrow = T) 
  
  # percentage of cases in which doses 1, 2 and 3 were present in the second stage
  armsel1 <- c(sum(res_stage2[,1]), sum(res_stage2[,2]), sum(res_stage2[,3]))/n_trials
  
  
  acc <- matrix(unlist(lapply(results_list, function(element) element$simdec_output)),ncol = 3, byrow = T) 
  
  armacc1 <- c(sum(acc[,1], na.rm = T)/sum(res_stage2[,1]), sum(acc[,2], na.rm = T)/sum(res_stage2[,2]), sum(acc[,3], na.rm = T)/sum(res_stage2[,3]))
  
  return(list(armsel1, armacc1))
}


#FUNCTION 3: function 1 AND function 2
########################################

#simul_res (mu_raw_0 = 650, sd_raw_0 = 575, r0 = 0.15, r1 = 0.6, r2 = 0.8, r3 = 0.9,  rho = 0.5, n_trials=1000,n_arms = 4,N1 = 90 , N2 = 60,
                  #   rmonth=11, alpha1=.1 , alpha=.025, v = c(1/2, 1/2, 0), sim_out=T,sel_scen=0, side=T,test="t")
  

simul_res = function(mu_raw_0 = 650, sd_raw_0 = 575, r0 = 0.15, r1 = 0.6, r2 = 0.8, r3 = 0.9,  rho = 0.5,
                     n_trials=n_trials,n_arms = 4,N1 = N1 , N2 = N1, rmonth=rmonth, alpha1=alpha1 , alpha=alpha ,
                     #v = c(1/2, 1/2, 0), 
                     sim_out=sim_out,sel_scen=sel_scen, side=side,test=test)
{
  
  ##Maybe this is not optimal but
#  
#  sz = as.numeric(N)
#  n1 = as.numeric(n1)
#  r0 = as.numeric(r0)
 
  #Compute mu et sigma
  ####################
  v = c(1/2, 1/2, 0)
  
  val = get_mu_sigma(mu_raw_0,sd_raw_0,  r0, r1, r2, r3,  rho) # rho is the correlation coefficient
 
  
  # extract the matrix and the mu for each time
  #############################################
  
  mtr1 = matrix(unlist(val[5]), nrow=3,byrow = T) #matrix I
  m0 = as.numeric(unlist(val[2])) #mu0
  m6 = as.numeric(unlist(val[3])) #mu6m
  m12 = as.numeric(unlist(val[4])) #mu12m
  
  #do pce and also the simulation with futur_map
  ##############################################
  
  aa = do_pce_baseline(n_trials=n_trials,n_arms = 4,N1 = N1 , N2 = N1, mu_0m = m0,   mu_6m = m6, 
                       mu_12m = mu12,  sg = mtr1, rmonth=rmonth, alpha1=alpha1 , alpha=alpha ,
                       v = c(1/2, 1/2, 0), sim_out=sim_out,sel_scen=sel_scen, side=side,test=test)
  
  return(aa)
  
}

oo<-mapply(simul_res, mu_raw_0 = 650, sd_raw_0 = 575, r0 = 0.15, r1 = 0.6, r2 = 0.8, r3 = 0.9,  rho = c(0.5,.7), n_trials=1000,n_arms = 4,
       N1 = 90 , N2 = 60,  rmonth=11, alpha1=c(.1,.2,.5) , alpha=.025, sim_out=T,sel_scen=0, side=T,test="t")
           
dfr<-data.frame(650, 575, 0.15, 0.6, 0.8, 0.9,  0.5, 1000,4,
90 , 60, 11, c(.1,.5), .025, T,0, T,"t")

oo1<-mapply(simul_res, 650, 575, 0, 0, 0, 0,  c(0.5,.7,0.5,0.7), 10000,4,
           90 , 60, 11, c(.1,.1,.5,.5), .025, T,0, T,"t")


par(mfrow=c(2,2))
plot(unlist(oo1[,1]),ylim=c(0,1))
plot(unlist(oo1[,2]),ylim=c(0,1))
plot(unlist(oo1[,3]),ylim=c(0,1))
plot(unlist(oo1[,4]),ylim=c(0,1))


oo2<-mapply(simul_res, dfr)


# FUNCTION 4: Do the simulation per scenario given a table x
###########################################################


sim_it = function(x)
  {
  
  #Do the simulation for each scenario in your table
  ##################################################
  
  
  result_list <- apply(x, 1, function(row) {
    do.call(simul_res, as.list(row))
  })
#}

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


rh <- 0.5

rh2 <- 0.7

# Define manually the dataset
#############################


#case 1:
########

cs1 <- data.frame(mu_raw_0 = 650, sd_raw_0 = 575, r0 = rep(0.15, 3), 
                  r1 = rep(0.65, 3), r2 = rep(0.75, 3), r3 = rep(0.90, 3), 
                  rho = rep(rh, 3),
                  N = rep(150, 3), n1 = rep(90, 3), alpha = rep(0.025, 3), 
                  alpha1 = c(0.1, 0.2, 0.5))

cs0 <- data.frame(mu_raw_0 = 650, sd_raw_0 = 575, r0 = rep(0.15, 3), 
                  r1 = rep(0.65, 3), r2 = rep(0.75, 3), r3 = rep(0.90, 3), 
                  rho = rep(rh2, 3),n_trials=10000,n_arms=4,
                  N1 = rep(90, 3), N2 = rep(90, 3), rmonth=11, 
                  alpha1 = c(0.1, 0.2, 0.5),alpha = rep(0.025, 3), 
                  sim_out=T,sel_scen=0, side=T,test="t")

cs <- rbind(cs1, cs0) %>% 
  mutate_at(vars(baseline_mean:alpha1), list(~as.numeric(.)))

# Do the simulation

s10 <- sim_it(cs0)



mu_raw_0 = 650, sd_raw_0 = 575, r0 = 0.15, r1 = 0.6, r2 = 0.8, r3 = 0.9,  rho = 0.5,
n_trials=n_trials,n_arms = 4,N1 = N1 , N2 = N1, rmonth=rmonth, alpha1=alpha1 , alpha=alpha ,
v = c(1/2, 1/2, 0), sim_out=sim_out,sel_scen=sel_scen, side=side,test=test