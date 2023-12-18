

sim_trial_pceind <- function(n_arms = 4, N1 , N2, mu_0m, mu_6m, mu_12m, sg, rmonth, alpha1 = 0.1, alpha = 0.025,v = c(1/2,1/2,0), sim_out=F)
{
  
  
  side <- T # this because db_stage1$y_6m-db_stage1$y_0m is negative #this to allow to change the side, can we add this as an argument of the function??? #Fabrice
  #na <- n_arms - 1
  
  db_stage1 <- sim_dataind(n_arms = n_arms-1, N = N1, mu_0m = mu_0m[1:n_arms-1], mu_6m = mu_6m[1:n_arms-1], mu_12m = mu_12m[1:n_arms-1], sg = sg, rmonth = rmonth)
  recruit_time1 <- max(db_stage1$recruit_time)
  
  
  
  db_stage1$treat <- relevel(db_stage1$treat, ref = "Placebo")
  
  db_stage1$diff6_0<-db_stage1$y_6m-db_stage1$y_0m
  db_stage1$diff12_0<-db_stage1$y_12m-db_stage1$y_0m
  
  
  #Linear model
  #############
  
  frame.sim <- subset(db_stage1, treat != "High")  #why do I need this?
  LowvsC <- lm(diff6_0 ~ treat, data=subset(frame.sim, treat!="Medium")) # adjust with the baseline please
  MedvsC <- lm(diff6_0  ~ treat , data=subset(frame.sim, treat!="Low"))
  
  
  #Obtain one-sided p-value
  #########################
  
  plow <- pt(coef(summary(LowvsC))[, 3], LowvsC$df, lower = side)[2]
  #t.test(db_stage1$diff6_0[db_stage1$treat!="Medium"]~db_stage1$treat[db_stage1$treat!="Medium"],alternative="less")
  
  pmed <- pt(coef(summary(MedvsC))[, 3], MedvsC$df, lower = side)[2]
  #t.test(db_stage1$diff6_0[db_stage1$treat!="Low"]~db_stage1$treat[db_stage1$treat!="Low"],alternative="less")
  pval.surr <- c(plow, pmed)  #pvalue of surrogate endpoint stage 1
  
  
  
  #######################################
  # decisions based on pvalues from linear model at 6 month
  # 1:Placebo, 2:low, 3:Median, 4:High
  if(sum(pval.surr < alpha1) == 2){  ## both have p<alpha1, both low and high are selected for second stage
    sc <- 2
    sel <- 2:3
  }
  if( (sum(pval.surr < alpha1) == 1 & which.min(pval.surr) == 1) ){  #low has p<alpha1, median has p>alpha1, both are NOT selected for second stage
    sc <- 1
    sel <- 0 
  } 
  if( (sum(pval.surr < alpha1) == 1 & which.min(pval.surr) == 2)){  #low has p>alpha1, median has p<alpha1, median is selected for second stage
    sc <- 3
    sel <- 3
  } 
  if(sum(pval.surr < alpha1) == 0){ # both have p>alpha1, both are NOT selected for second stage
    sc <- 0; sel <- 0 
  }
  
  
  #######################################
  # pvalues ttest 12 months, stage 1
  
  
  pval <- c()
  
  for(j in 1:(n_arms-2)){
    
    sub1 <- subset(db_stage1,(db_stage1$treat==levels(db_stage1$treat)[1])+ (db_stage1$treat==levels(db_stage1$treat)[j+1])==1)
    mod1 <- lm(diff12_0 ~ treat, sub1)
    res1 <- summary(mod1)
    pval[j] <- pt(coef(res1)[2,3], mod1$df, lower.tail = side)
  }
  z <- qnorm(1-pval)
  
  decision_s1 <- c()
  decision_s1[1] = ifelse(pval[1]<alpha1, "Reject", "Accept")
  decision_s1[2] = ifelse(pval[2]<alpha1, "Reject", "Accept")
  
  decision_stage1 = data.frame(decision_s1, row.names = levels(db_stage1$treat)[2:3])
  
  
  #preplanning adaptive conditional error
  #######################################
  
  N=N1+N2
  graph_bh <- BonferroniHolm(3)
  
  # the package assumes that wj are equal for all j
  z1 <- c(z,0)
  preplan <- doInterim(graph=graph_bh,z1=z1,v=v,alpha=alpha)
  # preplan@Aj
  # preplan@BJ
  
  #######################################
  # stage2
  # sc=2 --> Arm A and B continue to stage 2
  # sc=1 --> Arm A and B stop and start Arm C #if Arm B should be drop then drop Arm A too, start the second stage with only the Arm C Arm A or B continue to stage 2
  # sc=0 --> Arm A and B stop and start Arm C
  # sc=3 --> Arm A stop continue with Arm B and start Arm C
  
  
  if(sc==2){  # Arm A and B continue to stage 2  
    
    db_stage2 <- sim_dataind(n_arms=3, N = N2, mu_0m =mu_0m[c(1,2,3)], mu_6m =mu_6m[c(1,2,3)], mu_12m=mu_12m[c(1,2,3)], sg=sg, rmonth=rmonth)
    levels(db_stage2$treat) <- levels(db_stage1$treat)[c(1,2,3)]
    recruit_time2 <- max(db_stage2$recruit_time)
    
    db_stage2$diff6_0<-db_stage2$y_6m-db_stage2$y_0m
    db_stage2$diff12_0<-db_stage2$y_12m-db_stage2$y_0m
    
    
    
    pval2 <- c()
    
    for(j in 1:2){
      
      sub2 <- subset(db_stage2,(db_stage2$treat==levels(db_stage2$treat)[1])+(db_stage2$treat==levels(db_stage2$treat)[j+1])==1)
      mod2 <- lm(y_12m ~ diff12_0, sub2) #are we using this model or should we use individual models?
      res2 <- summary(mod2)
      pval2[j] <- pt(coef(res2)[2,3], mod2$df, lower.tail = side)
    }
    
    
    Avalues <- c(preplan@BJ[7]/2, #H123
                 preplan@BJ[6], #H12
                 preplan@BJ[5], #H13
                 preplan@BJ[3], #H23
                 preplan@BJ[2], #H2
                 preplan@BJ[4]  #H1
    )
    
    decision <- c()
    decision[1] <- ifelse(sum(pval2[1] <= Avalues[c(1,2,3,6)])==4, "Reject", "Accept") 
    decision[2] <- ifelse(sum(pval2[2] <= Avalues[c(1,2,4,5)])==4, "Reject", "Accept")
    
    decision_stage2 <- data.frame(decision, row.names = levels(db_stage2$treat)[2:3])
    
    pval2 <- data.frame(pval2, row.names = levels(db_stage2$treat)[2:3])
    
    stage2_arms <- c(1,1,0)
    simdec_output <- c(ifelse(decision[1]=="Reject", 1, 0),
                       ifelse(decision[2]=="Reject", 1, 0),
                       NA)
    
    decision_intersection = ifelse(sum(pval2 <= Avalues[1]) > 0, "Reject", "Accept")
  }
  
  
  
  if(sc == 1 ){# this means that low dose was the only one to be selected, so we drop both and start 
    #with high dose only in the second stage. Or none of the drug was accepted to cont
    
    #if(sel == 0){
    
    db_stage2 <- sim_dataind(n_arms=2, N = N2, mu_0m = mu_0m[c(1,4)],mu_6m = mu_6m[c(1,4)], mu_12m=mu_12m[c(1,4)], sg=sg, rmonth=rmonth)
    
    levels(db_stage2$treat) <- c(levels(db_stage1$treat)[c(1)],"High") # c(levels(db_stage1$treat)[c(1,sel)],"High")
    recruit_time2 <- max(db_stage2$recruit_time)
    
    db_stage2$diff6_0<-db_stage2$y_6m-db_stage2$y_0m
    db_stage2$diff12_0<-db_stage2$y_12m-db_stage2$y_0m
    
    
    pval2 <- c()
    
    mod2 <- lm(y_12m ~ treat, db_stage2) #are we using this model or should we use individual models?
    res2 <- summary(mod2)
    pval2 <- pt(coef(res2)[2,3], mod2$df, lower.tail = side)
    
    
    
    #}
    
    # if(sel=="Low"){ this should not be the case
    #   Avalues <- c(preplan@BJ[7]/2, #H123
    #                preplan@BJ[6]/2, #H12
    #                preplan@BJ[5], #H13
    #                preplan@BJ[3], #H23
    #                preplan@BJ[2], #H2
    #                preplan@BJ[1]  #H3
    #   )
    # }
    #if(sel == 3){# "Medium"){
    Avalues <- c(preplan@BJ[7]/2, #H123
                 preplan@BJ[6], #H12
                 preplan@BJ[5]/2, #H13
                 preplan@BJ[3], #H23
                 preplan@BJ[2], #H2
                 preplan@BJ[1]  #H3
    )
    
    decision <- c()
    
    dec <- sum(pval2 <= Avalues[c(1,3, 4, 6)])
    decision <- ifelse(dec == 4, "Reject", "Accept")
    
    decision_stage2 <- data.frame(decision, row.names = levels(db_stage2$treat)[2])
    
    pval2 <- data.frame(pval2, row.names = levels(db_stage2$treat)[2])
    
    stage2_arms <- c(0,0,1)
    simdec_output <- c(0, 0,
                       ifelse(decision[1]=="Reject", 1, 0))
  }  
  
  if(sc == 3){
    
    db_stage2 <- sim_dataind(n_arms=3, N=N2, mu_0m=mu_0m[c(1,3,4)], mu_6m=mu_6m[c(1,3,4)], mu_12m=mu_12m[c(1,3,4)], sg=sg, rmonth=rmonth)
    
    levels(db_stage2$treat) <- c(levels(db_stage1$treat)[c(1,3)],"High")
    recruit_time2 <- max(db_stage2$recruit_time)
    
    db_stage2$diff6_0<-db_stage2$y_6m-db_stage2$y_0m
    db_stage2$diff12_0<-db_stage2$y_12m-db_stage2$y_0m
    
    
    pval2 <- c()
    
    for(j in 1:2)
    {
      sub2 <- subset(db_stage2,(db_stage2$treat==levels(db_stage2$treat)[1])+ (db_stage2$treat==levels(db_stage2$treat)[j+1])==1)
      mod2 <- lm(diff12_0 ~ treat, sub2) #are we using this model or should we use individual models?
      res2 <- summary(mod2)
      pval2[j] <- pt(coef(res2)[2,3], mod2$df, lower.tail = side)
    }
    
    
    
    Avalues <- c(preplan@BJ[7]/2, #H123
                 preplan@BJ[6], #H12
                 preplan@BJ[5]/2, #H13
                 preplan@BJ[3], #H23
                 preplan@BJ[2], #H2
                 preplan@BJ[1]  #H3
    )
    
    
    # pval2[1] <= Avalues[c(1,2,4,5)] #p2
    # pval2[2] <= Avalues[c(1,3,4,6)] #p3
    # decision_stage2 = matrix(c(pval2[1] <= Avalues[c(1,2,4,5)],
    #                            pval2[2] <= Avalues[c(1,3,4,6)]),
    #                            byrow = F, ncol = 2)
    
    
    decision <- c()
    
    dec1 <- sum(pval2[1] <= Avalues[c(1,2,4,5)])#Avalues[c(1,2,4,5)])
    dec2 <- sum(pval2[2] <= Avalues[c(1,3,4,6)])
    
    decision[1] <- ifelse(dec1==4, "Reject", "Accept")
    decision[2] <- ifelse(dec2==4, "Reject", "Accept")
    
    #
    decision_stage2 <- data.frame(decision, row.names = levels(db_stage2$treat)[2:3])
    
    pval2 <- data.frame(pval2, row.names = levels(db_stage2$treat)[2:3])
    
    stage2_arms <- c(0,1,1)
    simdec_output <- c(0,
                       ifelse(decision[1]=="Reject", 1, 0),
                       ifelse(decision[2]=="Reject", 1, 0))
    
    decision_intersection <- ifelse(sum(pval2 <= Avalues[1]) > 0, "Reject", "Accept")
    
  } 
  
  
  
  if(sc == 0){ # in that case please you should start dose 3
    
    db_stage2 <- sim_dataind(n_arms=2, N=N2, mu_0m=mu_0m[c(1,4)], mu_6m=mu_6m[c(1,4)], mu_12m=mu_12m[c(1,4)], sg=sg, rmonth=rmonth)
    
    levels(db_stage2$treat) <- c(levels(db_stage1$treat)[c(1)],"High")
    recruit_time2 <- max(db_stage2$recruit_time)
    
    db_stage2$diff6_0<-db_stage2$y_6m-db_stage2$y_0m
    db_stage2$diff12_0<-db_stage2$y_12m-db_stage2$y_0m
    
    mod2 <- lm(diff12_0 ~ treat, db_stage2) #are we using this model or should we use individual models?
    res2 <- summary(mod2)
    pval2 <- pt(coef(res2)[2,3], mod2$df, lower.tail = side)
    
    Avalues <- c(preplan@BJ[7], #H123
                 preplan@BJ[5], #H13
                 preplan@BJ[3], #H23
                 preplan@BJ[1]  #H3
    )
    
    decision_stage2 = matrix(pval2 <= Avalues, ncol = 1) #p3
    
    
    decision <- c()
    decision[1] <- ifelse(sum(pval2 <= Avalues)==4, "Reject", "Accept")
    
    #
    decision_stage2 <- data.frame(decision, 
                                  row.names = levels(db_stage2$treat)[2])
    
    pval2 <- data.frame(pval2, row.names = levels(db_stage2$treat)[2])
    
    stage2_arms <- c(0,0,1)
    simdec_output <- c(0, 0,
                       ifelse(decision[1]=="Reject", 1, 0))
    
    decision_intersection <- ifelse(sum(pval2 <= Avalues[1]) > 0, "Reject", "Accept")
  }
  
  #######################################
  if(sim_out==F){
    return(list(pvalue_stage1=pval,
                pvalue_stage2=pval2,
                sc=sc,
                decision_stage1=decision_stage1,
                decision_stage2=decision_stage2,
                critical_values=preplan,
                selected_dose=sel,
                recruit_time1=recruit_time1,
                recruit_time2=recruit_time2,
                pval.surr=pval.surr))
  }else{
    res_intersection=ifelse(decision_intersection == "Reject", 1,0)
    return(list(stage2_arms=stage2_arms,
                selected_dose=sel,
                simdec_output=simdec_output,
                res_intersection=res_intersection))
  }
  
}

# Importing packages
####################

library(dplyr)
library(mvtnorm)
library(future)
library(furrr)
library(gMCP)



########################################
# PART 1: defining parameters
########################################



# Mean and standard deviation at baseline
#########################################

mn <- 650; sd <- 575

#Reduction rate
###############

r0 <- 0.15; r1 <- 0.6; r2 <- 0.8; r3 <- 0.9;

#Correlation between the difference, this will change to be between the outcome for each time points
####################################################################################################

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


#Twelve month
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


#Now let's calculate the mean for the log
#########################################

moy1 <- log(mn0^2/sqrt(mn0^2 + sd012m^2)) # mean(log(x0)) PLACEBO for twelve months
moy2 <- log(mn1^2/sqrt(mn1^2+sd112m^2)) # mean(log(x6m)) LOW D. for twelve months
moy3 <-  log(mn2^2/sqrt(mn2^2+sd212m^2))# mean(log(x6m)) Medium D. for twelve months
moy4 <-  log(mn3^2/sqrt(mn3^2+sd312m^2))# mean(log(x6m)) HIGH for twelve months

mu0 <-  c(moy0, moy0, moy0, moy0)
mu6 <- c(moy1, moy2, moy3, moy4)
mu12 <- c(moy1, moy2, moy3, moy4)




#Variance of the diff of the logs
#################################

rho <- 0.5; vardiff <- sde0^2  #sde1^2 + sde2^2 - 2* rho*sde1*sde2


# Variance-covariance matrix
############################

#1. Elements of the matrix

varb <- 1*sde0^2 # cov(log(X0), log(X0)) <- rho*sd0*sd0
var6m <- 1*sde0^2 # cov(log(X06m), log(X06m)) <- rho*sdx6m*sdx6m
var12m <- 1*sde0^2 # cov(log(X12m), log(X12m)) <- rho*sdx12m*sdx12m


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

########################################
# PART 2: TESTING THE TWO NEW FUNCTIONS
########################################

## A: testing sim_dataind

# simulate the data
###################

rb <- sim_dataind(n_arms = 4, N = 150, mu_0m = mu0, mu_6m = mu6, mu_12m = mu12, sg, rmonth)


# Check on the values you get: mean sd per group
################################################

rb %>% group_by(treat) %>% 
  summarize(baseline_mean_sd = paste0(round(mean(y_0m)), " (", round(sd(y_0m), 2),")"),
            sixmont_mean_sd = paste0(round(mean(y_6m)), " (", round(sd(y_6m), 2),")"),
            twlvmont_mean_sd = paste0(round(mean(y_12m)), " (", round(sd(y_12m), 2),")"))



# Let's do the simulation using the new function for pce 
########################################################

sim_trial_pceind(n_arms = 4, N1 , N2, mu_0m = c(moy0, moy0, moy0, moy0), 
                             mu_6m = c(moy1, moy2, moy3, moy4), mu_12m = c(moy1, moy2, moy3, moy4),
                             sg, rmonth, alpha1 = 0.1, alpha = 0.025,v = c(1/2,1/2,0), sim_out=T)



set.seed(20)

#Set the number of trials to run and other parameters for future plan
n_trials <- 10000
n_cores <- availableCores()-1 
plan(multisession, workers = n_cores)

# Run the simulations in parallel using future_map

# First case alpha = 0.1

results_list <- future_map(1:n_trials, function(i) sim_trial_pceind(n_arms = 4, N1 , N2, mu_0m = c(moy0, 
                                                                    moy0, moy0, moy0), 
                                                                    mu_6m = c(moy1, moy2, moy3, moy4), 
                                                                    mu_12m = c(moy1, moy2, moy3, moy4),
                                                                    sg, rmonth, alpha1 = 0.1, alpha = 0.025,
                                                                    v = c(1/2,1/2,0), sim_out=T), .options=furrr_options(seed = TRUE))



# Summary results
res_stage2 <- matrix(unlist(lapply(results_list, function(element) element$stage2_arms)),
                    ncol = 3, byrow = T) 

# percentage of cases in which doses 1, 2 and 3 were present in the second stage
armsel1 <- c(sum(res_stage2[,1]), sum(res_stage2[,2]), sum(res_stage2[,3]))/n_trials


acc <- matrix(unlist(lapply(results_list, function(element) element$simdec_output)),
             ncol = 3, byrow = T) 

armacc1 <- c(sum(acc[,1], na.rm = T), sum(acc[,2], na.rm = T), sum(acc[,3]))/n_trials




# Old function based on the difference only
###########################################


# Then mean(log(xi) - log(x0))

mu1 <-  moy0 - moy1
mu2 <-  moy0 - moy2
mu3 <-  moy0 - moy3
mu4 <-  moy0 - moy4

# The mu for the simulation: this is about the effect in log
############################################################


mu <- c(mu1, mu2, mu3, mu4)


#Run it here
############

results_list2 <- future_map(1:n_trials, function(i) sim_trial_pce(n_arms=4, N1=n1, N2=n2,
                                                                 mu_6m=mu, mu_12m=mu, sigma=sg_m, 
                                                                 rmonth=2, alpha1=0.1, alpha=0.025, 
                                                                 sim_out=T), .options=furrr_options(seed = TRUE))



# Summary results

res_stage2_old <- matrix(unlist(lapply(results_list2, 
                            function(element) element$stage2_arms)),
                                       ncol = 3, byrow = T) 

# percentage of cases in which doses 1, 2 and 3 were present in the second stage
armsel2 <- c(sum(res_stage2_old[,1]), sum(res_stage2_old[,2]), 
             sum(res_stage2_old[,3]))/n_trials


# percentage of cases in which doses 1, 2 and 3 were rejected at the final analsis (second stage)

acc <- matrix(unlist(lapply(results_list2, 
                           function(element) element$simdec_output)),
                                      ncol = 3, byrow = T) 

armacc2 <- c(sum(acc[,1], na.rm = T), sum(acc[,2], na.rm = T), sum(acc[,3], na.rm = T))/n_trials



# Comparison of the output
##########################

#Proportion of selected

armsel1

armsel2

#Proportion of rejected

armacc1
armacc2



