
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
#get_mu_sigma(650, 575, 0, 0, 0, 0,.5)
#FUNCTION 1: Do the calculation to obtain the mu and the sigma
##############################################################

#get_mu_sigma

#FUNCTION 2: Using mu and sigma obtained above then run the function below
##########################################################################

# n_trials=1000; n_arms=4; N1=900;N2=60 mu_0m =c(10,10,10,10); mu_6m =c(10,10,10,10); mu_12m=c(10,10,10,10); sg=matrix(c(1,0,0,0,1,0,0,0,1), ncol=3); 
#rmonth=1;alpha1=.1; alpha=0.025;
#do_pce_baseline(n_trials=10000,n_arms = 4,N1 = 90 , N2 = 60, mu_0m = mu_0m,   mu_6m = mu_6m, mu_12m = mu_12m,  sg = sg, rmonth=1, alpha1 = .1, alpha = .025, sim_out=T,sel_scen=0, side=T,test="t",dropout=.1,rr=rep(0,4))
#do_pce_baseline(n_trials=10000,n_arms = 4, N1=60 , N=90, mu_0m=rep(6.18,4), mu_6m=rep(6.18,4), mu_12m=rep(6.18,4), sg=matrix(c(.57,.28,.14,.28,.57,.28,.14,.28,.57),3), 
#                                        rmonth=1, alpha1 = 0.1, alpha = 0.025,sim_out=T,sel_scen=0, side=T,test="t",dropout=.1,rr=rep(0,4))

#do_pce_baseline(n_trials=n_trials,n_arms = n_arms, N1=N1 , N2=N2, mu_0m=mu_0m, mu_6m=mu_6m, mu_12m=mu_12m, sg=sg, 
#                rmonth=rmonth, alpha1 = alpha1, alpha = alpha,sim_out=sim_out,sel_scen=sel_scen, side=side,test=test,dropout=.1,rr=rep(0,4))

do_pce_baseline = function(n_trials,n_arms = 4,N1, N2, mu_0m, mu_6m, mu_12m,  sg, rmonth, alpha1 , alpha , sim_out,sel_scen, side,test,dropout,rr,bound)
{

  #Set the number of trials to run and other parameters for future plan
  #n_trials <- 10000
  n_cores <- availableCores()-1 
  plan(multisession, workers = n_cores)
  
  # Run the simulations in parallel using future_map
  results_list <- future_map(1:n_trials, function(i)#sim_trial_pceind_test (n_arms = 4, N1=60 , N2=90, mu_0m=rep(6.18,4), 
                                                    #                        mu_6m=rep(6.18,4), mu_12m=rep(6.18,4), 
                                                    #                        sg=matrix(c(.57,0,0,0,.57,0,0,0,.57),3), 
                                                    #                        rmonth=1, alpha1 = 0.1, alpha = 0.025,sim_out=T,sel_scen=0, side=T,test="t"), 
                                          #.options=furrr_options(seed = TRUE))
    
    
    sim_trial_pceind_test(n_arms = n_arms,N1 = N1 , N2 = N2, mu_0m = mu_0m,   mu_6m = mu_6m, 
                                                                          mu_12m = mu_12m,  sg = sg, rmonth=rmonth, alpha1=alpha1 , alpha =alpha,
                                                                          sim_out=sim_out,sel_scen=sel_scen, side=side,test=test,dropout=dropout,rr=rr,bound=bound), 
                                                                          .options=furrr_options(seed = TRUE))
  
  # percentage of cases in which doses 1, 2 and 3 were present in the second stage
  sel_stage2 <- matrix(unlist(lapply(results_list, function(element) element$stage2_arms)),ncol = 3, byrow = T) 
  armsel <- c(sum(sel_stage2[,1]), sum(sel_stage2[,2]), sum(sel_stage2[,3]))/n_trials
  
  #conditional power
  h_condpow <- matrix(unlist(lapply(results_list, function(element) element$simdec_output)),ncol = 3, byrow = T) 
  disjpow<-sum(apply(h_condpow,1,sum,na.rm=TRUE)>0)/n_trials
  condpow <- c(sum(h_condpow[,1], na.rm = T), sum(h_condpow[,2], na.rm = T), sum(h_condpow[,3], na.rm = T))/(armsel*n_trials)#sum(sel_stage2[,3])
  pow <-     c(sum(h_condpow[,1], na.rm = T), sum(h_condpow[,2], na.rm = T), sum(h_condpow[,3], na.rm = T))/n_trials
  
  
  #power of multiarmed trials 1 and 2
  h_ma1<-matrix(unlist(lapply(results_list, function(element) element$decision_ma1)),ncol = 3, byrow = T)
  h_ma2<-matrix(unlist(lapply(results_list, function(element) element$decision_ma2)),ncol = 3, byrow = T)
  h_ma1a<-matrix(unlist(lapply(results_list, function(element) element$decision_ma1a)),ncol = 3, byrow = T)
  
  pow_ma1<-c(apply(h_ma1,2,sum,na.rm=TRUE)/n_trials,sum(apply(h_ma1,1,sum,na.rm=TRUE)>0)/n_trials)
  pow_ma2<-c(apply(h_ma2,2,sum,na.rm=TRUE)/n_trials,sum(apply(h_ma2,1,sum,na.rm=TRUE)>0)/n_trials)
  pow_ma1a<-c(apply(h_ma1a,2,sum,na.rm=TRUE)/n_trials,sum(apply(h_ma1a,1,sum,na.rm=TRUE)>0)/n_trials)

  #recruittime
  recruit1 <- matrix(unlist(lapply(results_list, function(element) element$recruit_time1)),ncol = 1, byrow = T) 
  recruittime1<-sum(recruit1, na.rm = T)/n_trials
  recruit2 <- matrix(unlist(lapply(results_list, function(element) element$recruit_time2)),ncol = 1, byrow = T) 
  recruittime2<-sum(recruit2, na.rm = T)/n_trials

  
  return(c(armsel, condpow,pow,disjpow,pow_ma1,pow_ma2,pow_ma1a,recruittime1,recruittime2))
}

#FUNCTION 3: function 1 AND function 2
########################################

  

simul_res = function(mu_raw_0, sd_raw_0 , r0_6,r1_6,r2_6,r3_6, r0_12,r1_12,r2_12,r3_12,  rho ,
                     n_trials,n_arms = 4,N1 , N2, rmonth, alpha1 , alpha ,
                     sim_out1,sel_scen, side1,test1,dropout,rr0,rr1,rr2,rr3,bound)
{
  
  ##Maybe this is not optimal but
#  
#  sz = as.numeric(N)
#  n1 = as.numeric(n1)
#  r0 = as.numeric(r0)

  #sim_out<-ifelse(sim_out1==F,0,1)
  #if (sim_out1==1) sim_out=="T"
  if (sim_out1==1) sim_out<-T
  if (sim_out1==0) sim_out<-F
  if (side1==1) side<-T
  if (side1==0) side<-F
  if (test1==0) test<-"l"
  if (test1==1) test<-"m"
  if (test1==2) test<-"t"
  if (test1==3) test<-"w"
  if (test1==4) test<-"w1"
  if (test1==5) test<-"f"
  
   #N2<-N-N1
  
  #Compute mu et sigma
  ####################
  reductrate_6<-c(r0_6,r1_6,r2_6,r3_6) 
  reductrate_12<-c(r0_12,r1_12,r2_12,r3_12) 
  
  rr<-c(rr0,rr1,rr2,rr3) 
  
  #if (geom==0)
  val = get_mu_sigma(mu_raw_0,sd_raw_0,  reductrate_6,reductrate_12,  rho) # rho is the correlation coefficient
                    
  # extract the matrix and the mu for each time
  #############################################
  
  mtr1 = matrix(unlist(val[5]), nrow=3,byrow = T) #matrix I
  m0 = as.numeric(unlist(val[2])) #mu0
  m6 = as.numeric(unlist(val[3])) #mu6m
  m12 = as.numeric(unlist(val[4])) #mu12m
  
  #do pce and also the simulation with futur_map
  ##############################################
  
  aa = #do_pce_baseline(n_trials=n_trials,n_arms = 4,N1 = N1, N2 = N2, mu_0m = m0,   mu_6m = m6, mu_12m = m12,  sg = mtr1, 
      #                  rmonth=rmonth, alpha1=alpha1 , alpha=alpha ,   sim_out=sim_out,sel_scen=sel_scen, side=side,test=test)
        do_pce_baseline(n_trials=n_trials,n_arms = 4,N1=N1 , N2=N2, mu_0m=m0, mu_6m=m6, mu_12m=m12, sg=mtr1,
                    rmonth=rmonth, alpha1=alpha1 , alpha=alpha , sim_out=sim_out,sel_scen=sel_scen, side=side,test=test,dropout=dropout,rr=rr,bound=bound)
  
  #return(list(mtr1,m0,m6,m12,aa))
  return(aa)
  
}

simul_res (mu_raw_0 = 650, sd_raw_0 = 575, r0_6=0, r1_6=0,r2_6=0,r3_6=0,r0_12=0,r1_12=0,r2_12=0,r3_12=0, rho = 0.5, 
           n_trials=10000,n_arms = 4,N1 = 60 , N2 = 150,
                rmonth=1, alpha1=.1 , alpha=.025, sim_out1=1,sel_scen=0, side1=1,test1=2,dropout=.1,rr0=0,rr1=0,rr2=0,rr3=0,bound=0)

r0_6=0; r1_6=0;r2_6=0;r3_6=0;r0_12=0;r1_12=0;r2_12=0;r3_12=0 rho = 0.5, 

simul_res (650, 575, 0, 0,0,0,0,0,0,0, 0.5,10000,4, 60 , 150,1, .1 , .025, 1,0, 1,2,.1,,rr0=0,rr1=0,rr2=0,rr3=0)


do_pce_baseline(n_trials=10000,n_arms = 4, N1=60 , N=90, mu_0m=rep(6.18,4), mu_6m=rep(6.18,4), mu_12m=rep(6.18,4), sg=matrix(c(.57,.28,.14,.28,.57,.28,.14,.28,.57),3), 
                rmonth=1, alpha1 = 0.1, alpha = 0.025,sim_out=T,sel_scen=0, side=T,test="t",dropout=.1,rr=rep(0,4))


do_pce_baseline (n_trials=10000,n_arms = 4,N1 = 90 , N2 = 60, mu_0m = m0,   mu_6m = m6, mu_12m = m12,  sg = sg, rmonth=1, alpha1 = .1, alpha = .025, sim_out=T,sel_scen=0, side=T,test="t",dropout=.1,rep(0,4))






#oo<-mapply(simul_res, mu_raw_0 = 650, sd_raw_0 = 575, r0 = 0.15, r1 = 0.6, r2 = 0.8, r3 = 0.9,  rho = c(0.5,.7), n_trials=1000,n_arms = 4,
#       N1 = 90 , N2 = 60,  rmonth=11, alpha1=c(.1,.2,.5) , alpha=.025, sim_out=T,sel_scen=0, side=T,test="t")
           
#dfr<-data.frame(650, 575, 0.15, 0.6, 0.8, 0.9,  0.5, 1000,4,
#90 , 60, 11, c(.1,.5), .025, T,0, T,"t")

#

#oo1<-mapply(simul_res, 650, 575, 0, 0, 0, 0,  c(0,0.5,1,0,0.5,1), 10000,4, 60 , 90, 1, c(.1,.1,.1,.5,.5,.5), .025, T,0, T,"t")


#oo2<-mapply(simul_res, 650, 575, 0.15, 0.20, 0.25, 0.3,  c(0,0.5,1,0,0.5,1), 10000,4, 60 , 90, 1, c(.1,.1,.1,.5,.5,.5), .025, T,0, T,"t")

#oo3<-mapply(simul_res, 650, 575, 0.15, 0.3, 0.6, 0.9,  c(0,0.5,1,0,0.5,1), 10000,4, 60 , 90, 1, c(.1,.1,.1,.5,.5,.5), .025, T,0, T,"t")



#oo2<-mapply(simul_res, dfr)


# FUNCTION 4: Do the simulation per scenario given a table x
###########################################################


#sim_it = function(x)
#  {
#  
#  #Do the simulation for each scenario in your table
#  ##################################################
#  
#  
#  result_list <- apply(x, 1, function(row) {
#    do.call(simul_res, as.list(row))
#  })
##}##
#
#  b = length(unique(x$rho))
#  
#  # Now extract element in each list (scenarios) wich corresponds to proportion of selection
#  ##########################################################################################
#  
#  
#  res1 = data.frame(t(sapply(lapply(result_list, function(x) x[[1]]), #extract the first element in each list: proportion of sel
#                             function(x) unlist(x)))) %>% # each element beacomes a column in a data.frame
#    mutate(Hypothesis = paste0(rep(unique(x[, "alpha1"]), b), "- rho = ", x$rho) )#, # 
#  
#  scen = c("H1", "H2", "H3") # Hypothesis: Low, medium, High
#  names(res1)[1:3] = scen
#  
#  
#  res2 = data.frame(t(sapply(lapply(result_list, function(x) x[[2]]), #extract the second element in each list: proportion of sel
#                             function(x) unlist(x)))) %>% 
#    mutate(Hypothesis = paste0(rep(unique(x[, "alpha1"]), b), "- rho = ", x$rho) )#
#  
#  
#  names(res2)[1:3] = scen
#  
#  all_scenario <- list(res1, res2)
#  
#  return(all_scenario)
#  
#}


# Define manually the dataset
#############################


#case 1:
########

#cs1 <- data.frame(mu_raw_0 = 650, sd_raw_0 = 575, r0 = rep(0.15, 3), 
#                  r1 = rep(0.65, 3), r2 = rep(0.75, 3), r3 = rep(0.90, 3), 
#                  rho = rep(rh, 3),
#                  N = rep(150, 3), n1 = rep(90, 3), alpha = rep(0.025, 3), 
#                  alpha1 = c(0.1, 0.2, 0.5))#

#cs0 <- data.frame(mu_raw_0 = 650, sd_raw_0 = 575, r0 = rep(0.15, 3), 
#                  r1 = rep(0.65, 3), r2 = rep(0.75, 3), r3 = rep(0.90, 3), 
#                  rho = rep(rh2, 3),n_trials=10000,n_arms=4,
#                  N1 = rep(90, 3), N2 = rep(90, 3), rmonth=11, 
#                  alpha1 = c(0.1, 0.2, 0.5),alpha = rep(0.025, 3), 
#                  sim_out=T,sel_scen=0, side=T,test="t")

#cs <- rbind(cs1, cs0) %>% 
#  mutate_at(vars(baseline_mean:alpha1), list(~as.numeric(.)))

# Do the simulation

#s10 <- sim_it(cs0)



#mu_raw_0 = 650, sd_raw_0 = 575, r0 = 0.15, r1 = 0.6, r2 = 0.8, r3 = 0.9,  rho = 0.5,
#n_trials=n_trials,n_arms = 4,N1 = N1 , N2 = N1, rmonth=rmonth, alpha1=alpha1 , alpha=alpha ,
#v = c(1/2, 1/2, 0), sim_out=sim_out,sel_scen=sel_scen, side=side,test=test