#' Simulate trial
#' @description Simulate trial
#'
#' @param mu_raw_0 num mean value of original data at baseline
#' @param sd_raw_0 num standard deviation of original data at baseline
#' @param r0_6 reduction rate of control dose at Month 6 with values between 0 and 1
#' @param r1_6 reduction rate of low dose at Month 6 with values between 0 and 1
#' @param r2_6 reduction rate of medium dose at Month 6 with values between 0 and 1
#' @param r3_6 reduction rate of high dose at Month 6 with values between 0 and 1
#' @param r0_12 reduction rate of control dose at Month 12 with values between 0 and 1
#' @param r1_12 reduction rate of low dose at Month 12 with values between 0 and 1
#' @param r2_12 reduction rate of medium dose at Month 12 with values between 0 and 1
#' @param r3_12 reduction rate of high dose at Month 12 with values between 0 and 1
#' @param rho correlation between baseline and 12 months follow-up observations and between 6 months and 12 months follow-up observations
#' @param n_trials number of simulation runs
#' @param n_arms number of arms (including control)
#' @param N1 sample size stage 1
#' @param N2 sample size stage 2
#' @param alpha1 significance level for dose selection (futility boundary)
#' @param alpha significance level for selected dose vs control comparison
#' @param sel_scen choose between two different options in case that in interim analysis low dose is promising, but median dose not: 0: do not continue with low dose or median dose; 1: continue with low and median doses
#' @param side1 TRUE/FALSE referring to the side for 1-side testing (if TRUE then lower = side)
#' @param test1 defines type of analysis: "t" calculates a t-test, "l" a linear model with baseline values as covariables, "w" Wilcoxon test of differences, and "w1" Wilcoxon test of follow-up values
#' @param dropout dropoutrate with values between 0 and 1
#' @param rr0 total responder rate for control dose, which gives the proportion of patients with value 0 at follow-up
#' @param rr1 total responder rate for low dose, which gives the proportion of patients with value 0 at follow-up
#' @param rr2 total responder rate for medium dose, which gives the proportion of patients with value 0 at follow-up
#' @param rr3 total responder rate for high dose, which gives the proportion of patients with value 0 at follow-up
#' @param bound lower bound to define total responder in simulation study
#' @keywords internal
#' @returns A vector consisting of summary measures for data simulation with n_trials repetitions: frequency selected for stage 2, conditional power, power ,disjunctive power, power of MA1, power of MA2, power of ma1a, concordances
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


simul_res = function(mu_raw_0, sd_raw_0 , r0_6,r1_6,r2_6,r3_6, r0_12,r1_12,r2_12,r3_12,  rho ,
                     n_trials,n_arms = 4,N1 , N2, alpha1 , alpha ,
                     sel_scen, side1, test1, dropout, rr0, rr1, rr2, rr3, bound)
{
  
  if (side1==1) side<-T
  if (side1==0) side<-F
  if (test1==0) test<-"l"
  if (test1==2) test<-"t"
  if (test1==3) test<-"w"
  if (test1==4) test<-"w1"
  
  
  #Compute mu et sigma
  ####################
  reductrate_6<-c(r0_6,r1_6,r2_6,r3_6) 
  reductrate_12<-c(r0_12,r1_12,r2_12,r3_12) 
  
  rr<-c(rr0,rr1,rr2,rr3) 
  
  val = get_mu_sigma(mu_raw_0,sd_raw_0,  reductrate_6,reductrate_12,  rho) # rho is the correlation coefficient
  
  # extract the matrix and the mu for each time
  #############################################
  
  mtr1 = matrix(unlist(val[5]), nrow=3,byrow = T) #matrix I
  m0 = as.numeric(unlist(val[2])) #mu0
  m6 = as.numeric(unlist(val[3])) #mu6m
  m12 = as.numeric(unlist(val[4])) #mu12m
  
  #do pce and simulation with future_map
  ##############################################
  
  aa =     do_pce_baseline(n_trials=n_trials,n_arms = 4,N1=N1 , N2=N2, mu_0m=m0, mu_6m=m6, mu_12m=m12, sg=mtr1,
                    alpha1=alpha1 , alpha=alpha , sel_scen=sel_scen, side=side,test=test,dropout=dropout,rr=rr,bound=bound)
  
  return(aa)
}
