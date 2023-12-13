#' Simulate data from a multi-arm multi-stage trial with shared control and two initial doses, where an additional dose could be added after the interim analysis
#' @description Function to simulate trial data (2-stages, with dose selection). The analyses are performed using partial conditional error.
#'
#' @param n_arms number of arms (including control)
#' @param N1 sample size stage 1
#' @param N2 sample size stage 2
#' @param mu_6m Mean difference between baseline and 6-month per arm (vector of length `n_arm`)
#' @param mu_12m Mean difference between baseline and 12-month per arm (vector of length `n_arm`)
#' @param sigma covariance matrix between 6- and 12-month mean differences assumed equal across arms (matrix of dim 2x2)
#' @param alpha1 significance level for dose selection (futility analysis)
#' @param alpha significance level for selected dose vs control comparison
#' @param v vector giving the proportions of pre-planned measurements collected up to the interim analysis
#' @param rmonth recruitment per month (recruitment speed assumed constant over time)
#' @param sim_out Option for simplified output for simulations (if `sim_out=TRUE` simplified version, the value is `FALSE` by default)
#' @returns A list consisting of pvalues at stage 1, pvalues at stage 2, the decision at stages 1 and 2, the selected dose at stage 1, and the time at which the last patient was recruited in stage 1 and 2.
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats runif
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @importFrom stats aov
#' @importFrom stats t.test
#' @importFrom stats coef
#' @importFrom stats lm
#' @importFrom stats pt
#' @importFrom multcomp glht
#' @importFrom multcomp mcp
#' @importFrom gMCP doInterim
#' @importFrom gMCP BonferroniHolm
#' @export
#' @details eWHORM simulations
#' @author Marta Bofill Roig



# Function to simulate trial data (2-stages, with dose selection)
# individual observations are simulated

# n_arms=4; N1=100; mu_0m =c(0,0,0,0); mu_6m =c(1,1,1,1); mu_12m=c(1,2,3,4); sg=matrix(c(1,0,0,0,1,0,0,0,1), ncol=3); rmonth=1
#sim_trial_pceind (n_arms = 4, N1=60 , N=90, mu_0m=c(0,0,0,0), mu_6m=c(1,1,1,1), mu_12m=c(1,2,3,4), sg=matrix(c(1,0,0,0,1,0,0,0,1),3), rmonth=1, alpha1 = 0.1, alpha = 0.025,v = c(1/2,1/2,0), sim_out=T)

sim_trial_pceind <- function(n_arms = 4, N1 , N2, mu_0m, mu_6m, mu_12m, sg, rmonth, alpha1 = 0.1, alpha = 0.025,v = c(1/2,1/2,0), sim_out=F)
  {


  side <- F #this to allow to change the side, can we add this as an argument of the function??? #Fabrice
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

