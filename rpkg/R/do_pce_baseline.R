#' Help function to simulate trial
#' @description Help function to simulate trial
#'
#' @param n_trials number of simulation runs
#' @param n_arms number of arms (including control)
#' @param N1 sample size stage 1
#' @param N2 sample size stage 2
#' @param mu_0m Baseline value per arm (vector of length `n_arm`)
#' @param mu_6m 6-month value per arm (vector of length `n_arm`)
#' @param mu_12m 12-month value per arm (vector of length `n_arm`)
#' @param sg covariance matrix between 6- and 12-month mean differences assumed equal across arms (matrix of dim 2x2)
#' @param alpha1 significance level for dose selection (futility boundary)
#' @param alpha significance level for selected dose vs control comparison
#' @param sel_scen choose between two different options in case that in interim analysis low dose is promising, but median dose not: 0: do not continue with low dose or median dose; 1: continue with low and median doses
#' @param side TRUE/FALSE referring to the side for 1-side testing (if TRUE then lower = side)
#' @param test defines type of analysis: "t" calculates a t-test, "l" a linear model with baseline values as covariables, "w" Wilcoxon test of differences, and "w1" Wilcoxon test of follow-up values
#' @param dropout dropoutrate, between 0 and 1
#' @param rr responder rate for each dose (vector of length `n_arm`), which gives the proportion of patients with value 0 at follow-up
#' @param bound lower bound to define total responder in simulation study
#' @keywords internal
#' @returns A vector consisting of summary measures for data simulation with n_trials repetitions: frequency selected for stage 2, conditional power, power, disunctive power, power of MA1, power of MA2, power of ma1a, concordances
#' @export
#' @import dplyr
#' @import tidyr
#' @import future
#' @import furrr
#' @import gMCP
#' @import mvtnorm
#' @details eWHORM simulations
#' @author Marta Bofill Roig, Sonja Zehetmayer


# library(dplyr)
# library(tidyr)
# library(future)
# library(furrr)
# library(gMCP)
# library(mvtnorm)


do_pce_baseline = function(n_trials,n_arms = 4,N1, N2, mu_0m, mu_6m, mu_12m,  sg, alpha1 , alpha , sel_scen, side,test,dropout,rr,bound)
{
  
  #Set the number of trials to run and other parameters for future plan
  n_cores <- availableCores()-1 
  plan(multisession, workers = n_cores)
  
  # Run the simulations in parallel using future_map
  results_list <- future_map(1:n_trials, function(i) 
    sim_trial_pceind_test(n_arms = n_arms,N1 = N1 , N2 = N2, mu_0m = mu_0m,   mu_6m = mu_6m, 
                          mu_12m = mu_12m,  sg = sg, alpha1=alpha1 , alpha =alpha,
                          sel_scen=1, side=side,test=test,dropout=dropout,rr=rr,bound=bound)
    , 
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
  
   
  ##################new BIAS
  #concordance
  #true value:
  true<-sim_trial_pceind_test(n_arms = n_arms,N1 = N1*5000 , N2 = N2*5000, mu_0m = mu_0m,   mu_6m = mu_6m, 
                              mu_12m = mu_12m,  sg = sg, alpha1=alpha1 , alpha =alpha,
                              sel_scen=sel_scen, side=side,test=test,dropout=dropout,rr=rr,bound=bound)
  trueconc<-true$concMA2
  #simul
  concord <- matrix(unlist(lapply(results_list, function(element) element$conc)),ncol = 3, byrow = T)
  sdconcord<-apply(concord,2,sd,na.rm=TRUE)
  meanconcord<-apply(concord,2,mean,na.rm=TRUE)
  biasconcord<-meanconcord-trueconc
  concordMA1 <- matrix(unlist(lapply(results_list, function(element) element$concMA1)),ncol = 3, byrow = T)
  meanconcordMA1<-apply(concordMA1,2,mean,na.rm=TRUE)
  concord11 <- matrix(unlist(lapply(results_list, function(element) element$conc1)),ncol = 3, byrow = T)
  conc1.<-apply(concord11,2,mean,na.rm=TRUE)
  concord21 <- matrix(unlist(lapply(results_list, function(element) element$conc2)),ncol = 3, byrow = T)
  conc2.<-apply(concord21,2,mean,na.rm=TRUE)
  
  concordcond <- matrix(unlist(lapply(results_list, function(element) element$conccond)),ncol = 3, byrow = T)
  meanconcordcond<-apply(concordcond,2,mean,na.rm=TRUE)
  biasconcordcond<-meanconcordcond-trueconc
  
  concordinvn <- matrix(unlist(lapply(results_list, function(element) element$concinvn)),ncol = 3, byrow = T)
  meanconcordinvn<-apply(concordinvn,2,mean,na.rm=TRUE)
  biasconcordinvn<-meanconcordinvn-trueconc
  
  
  #CI
  concordCI <- matrix(unlist(lapply(results_list, function(element) element$concCI)),ncol = 3, byrow = T)
  meanconcordCI<-apply(concordCI,2,mean,na.rm=TRUE)
  coverage<-apply(t(concordCI)<trueconc,1,sum,na.rm=TRUE)/apply(!is.na(concordCI),2,sum)
  
  concordCIcond <- matrix(unlist(lapply(results_list, function(element) element$concCIcond)),ncol = 3, byrow = T)
  meanconcordCIcond<-apply(concordCIcond,2,mean,na.rm=TRUE)
  coveragecond<-apply(t(concordCIcond)<trueconc,1,sum,na.rm=TRUE)/apply(!is.na(concordCIcond),2,sum)
  
  concordCIinvn <- matrix(unlist(lapply(results_list, function(element) element$concCIinvn)),ncol = 3, byrow = T)
  meanconcordCIinvn<-apply(concordCIinvn,2,mean,na.rm=TRUE)
  coverageinvn<-apply(t(concordCIinvn)<trueconc,1,sum,na.rm=TRUE)/apply(!is.na(concordCIinvn),2,sum)
  

  return(c(armsel, condpow,pow,disjpow,pow_ma1,pow_ma2,pow_ma1a,#recruittime1,recruittime2,
           trueconc,meanconcord,sdconcord,biasconcord,meanconcordMA1,meanconcordCI,coverage,meanconcordcond,biasconcordcond,#coverageBH,
           coveragecond,meanconcordCIinvn,
           coverageinvn,meanconcordCIcond,meanconcordinvn,biasconcordinvn))
}
