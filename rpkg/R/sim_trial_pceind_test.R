#' Simulate data from a multi-arm multi-stage trial with shared control and two initial doses, where an additional dose could be added after the interim analysis; individual observations are simulated
#' @description Function to simulate trial data (2-stages, with dose selection). The analyses are performed using partial conditional error rates.
#'
#' @param n_arms number of arms (including control)
#' @param N1 sample size stage 1
#' @param N2 sample size stage 2
#' @param mu_0m Baseline value per arm (vector of length `n_arm`)
#' @param mu_6m 6-month value per arm (vector of length `n_arm`)
#' @param mu_12m 12-month value per arm (vector of length `n_arm`)
#' @param sigma covariance matrix between 6- and 12-month mean differences assumed equal across arms (matrix of dim 2x2)
#' @param alpha1 significance level for dose selection (futility boundary)
#' @param alpha significance level for selected dose vs control comparison
#' @param sel_scen choose between two different options in case that in interim analysis low dose is promising, but median dose not: 0: do not continue with low dose or median dose; 1: continue with low and median doses
#' @param side TRUE/FALSE referring to the side for 1-side testing (if TRUE then lower = side)
#' @param test defines type of analysis: "t" calculates a t-test, "l" a linear model with baseline values as covariables, "w" Wilcoxon test of differences, and "w1" Wilcoxon test of follow-up values
#' @param dropout dropoutrate, between 0 and 1
#' @param rr responder rate for each dose (vector of length `n_arm`), which gives the proportion of patients with value 0 at follow-up
#' @param bound lower bound to define total responder in simulation study
#' @returns A list consisting of pvalues at stage 1, pvalues at stage 2, the decision at stages 1 and 2, the selected dose at stage 1, concordance values, corresponding bias and confidence limits and corresponding values on the multi-armed trials.
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

#library(lme4)
#library(lmerTest)
library(multcomp)
library(gMCP)
# 



sim_trial_pceind_test <- function(n_arms = 4, N1 , N2, mu_0m, mu_6m, mu_12m, sg, alpha1 , alpha = 0.025, sel_scen, side=T,test,dropout,rr,bound)
{
  N1orig<-N1
  N1<-floor(N1*(1-dropout))
  db_stage1 <- sim_dataind(n_arms = n_arms-1, N = N1, mu_0m = mu_0m[1:n_arms-1], mu_6m = mu_6m[1:n_arms-1], mu_12m = mu_12m[1:n_arms-1], sg = sg, rr=rr[1:n_arms-1],bound=bound)
  
  db_stage1$treat <- relevel(db_stage1$treat, ref = "Placebo")
  
  db_stage1$diff6_0<-db_stage1$y_6m-db_stage1$y_0m
  db_stage1$diff12_0<-db_stage1$y_12m-db_stage1$y_0m
  
  
  #interim anlysis

  #Obtain one-sided p-value

  plow <- t.test(db_stage1$y_6m[db_stage1$treat!="Medium"]~db_stage1$treat[db_stage1$treat!="Medium"],alternative="greater")$p.value
  pmed<-t.test(db_stage1$y_6m[db_stage1$treat!="Low"]~db_stage1$treat[db_stage1$treat!="Low"],alternative="greater")$p.value
  
  
  if (test=="w"){
    plow <- wilcox.test(db_stage1$diff6_0[db_stage1$treat!="Medium"]~db_stage1$treat[db_stage1$treat!="Medium"],alternative="greater")$p.value
    pmed<-wilcox.test(db_stage1$diff6_0[db_stage1$treat!="Low"]~db_stage1$treat[db_stage1$treat!="Low"],alternative="greater")$p.value
    }
  
  if (test=="w1"){
    plow <- wilcox.test(db_stage1$y_6m[db_stage1$treat!="Medium"]~db_stage1$treat[db_stage1$treat!="Medium"],alternative="greater")$p.value
    pmed<-wilcox.test(db_stage1$y_6m[db_stage1$treat!="Low"]~db_stage1$treat[db_stage1$treat!="Low"],alternative="greater")$p.value
  }
  
    pval.surr <- c(plow, pmed)  #pvalue of surrogate endpoint stage 1
  
  
  
  #######################################
  # decisions based on pvalues from linear model at 6 month
  # 1:Placebo, 2:low, 3:Median, 4:High
  if(sum(pval.surr < alpha1) == 2){  ## both have p<alpha1, both low and medium are selected for second stage
    sc <- 2
    sel <- 2:3
  }
  if( (sum(pval.surr < alpha1) == 1 & which.min(pval.surr) == 1) ){  #low has p<alpha1, median has p>alpha1, both are %NOT% selected for second stage
    sc <- 1
    sel <- 0
  }
  if( (sum(pval.surr < alpha1) == 1 & which.min(pval.surr) == 2)){  #low has p>alpha1, medium has p<alpha1, median is selected for second stage
    sc <- 3
    sel <- 3
  }
  if(sum(pval.surr < alpha1) == 0){ # both have p>alpha1, both are NOT selected for second stage
    sc <- 0; sel <- 0
  }
  
  
  #######################################
  # pvalues ttest 12 months, stage 1
  
  
  pval1 <- c()  #12months p-value of first stage
  
  if (test=="t"){
    p12low <- t.test(db_stage1$diff12_0[db_stage1$treat!="Medium"]~db_stage1$treat[db_stage1$treat!="Medium"],alternative="greater")$p.value
    p12med<-t.test(db_stage1$diff12_0[db_stage1$treat!="Low"]~db_stage1$treat[db_stage1$treat!="Low"],alternative="greater")$p.value
    pval1<-cbind(p12low,p12med)
  }
  
  if (test=="l"){  #two individual models, due to robustness in case of different variances
    
    for(j in 1:(n_arms-2)){
      
      sub1 <- subset(db_stage1,(db_stage1$treat==levels(db_stage1$treat)[1])+ (db_stage1$treat==levels(db_stage1$treat)[j+1])==1)
      mod1 <- lm(y_12m ~ treat+y_0m, sub1)
      res1 <- summary(mod1)
      pval1[j] <- pt(coef(res1)[2,3], mod1$df, lower.tail = side)
    }
  }
  
  
  if (test=="w"){
    p12low <- wilcox.test(db_stage1$diff12_0[db_stage1$treat!="Medium"]~db_stage1$treat[db_stage1$treat!="Medium"],alternative="greater")$p.value
    p12med<-wilcox.test(db_stage1$diff12_0[db_stage1$treat!="Low"]~db_stage1$treat[db_stage1$treat!="Low"],alternative="greater")$p.value
    pval1<-cbind(p12low,p12med)
    #pval
  }
  
  conc1 <- vector(length=3)  
    
  if (test=="w1"){
    p12low <- wilcox.test(db_stage1$y_12m[db_stage1$treat!="Medium"]~db_stage1$treat[db_stage1$treat!="Medium"],alternative="greater")$p.value
    p12med<-wilcox.test(db_stage1$y_12m[db_stage1$treat!="Low"]~db_stage1$treat[db_stage1$treat!="Low"],alternative="greater")$p.value
    pval1<-cbind(p12low,p12med)
    conc1low <- wilcox.test(db_stage1$y_12m[db_stage1$treat!="Medium"]~db_stage1$treat[db_stage1$treat!="Medium"],alternative="greater")$statistic/(N1/3*N1/3)
    conc1med<-wilcox.test(db_stage1$y_12m[db_stage1$treat!="Low"]~db_stage1$treat[db_stage1$treat!="Low"],alternative="greater")$statistic/(N1/3*N1/3)
    conc1<-cbind(conc1low,conc1med,NA)
  }
  
 
  z <- qnorm(1-pmax(pval1,1e-15))
  
  decision_s1 <- c()
  decision_s1[1] = ifelse(pval1[1]<alpha1, "Reject", "Accept")
  decision_s1[2] = ifelse(pval1[2]<alpha1, "Reject", "Accept")
  
  decision_stage1 = data.frame(decision_s1, row.names = levels(db_stage1$treat)[2:3])
  
  
  #preplanning adaptive conditional error
  #######################################
  N=N1+N2
  graph_bh <- BonferroniHolm(3)
  
  v = c(N1/N,N1/N,0)
  # the package assumes that wj are equal for all j
  z1 <- c(z,0)
  preplan <- doInterim(graph=graph_bh,z1=z1,v=v,alpha=alpha)

  #######################################
  # stage2
  # sc=2 --> Arm A and B continue to stage 2
  # sc=1 --> Arm A and B stop and start Arm C #if Arm B should be drop then drop Arm A too, start the second stage with only the Arm C Arm A or B continue to stage 2
  # sc=0 --> Arm A and B stop and start Arm C
  # sc=3 --> Arm A stop continue with Arm B and start Arm C
  
  conc2<-vector(length=3)
  N2orig<-N2
  
  if(sc==2){  # Arm A and B continue to stage 2
    N2<-floor(N2*(1-dropout)) #floor(N2/(1+dropout))
    db_stage2 <- sim_dataind(n_arms=3, N = N2, mu_0m =mu_0m[c(1,2,3)], mu_6m =mu_6m[c(1,2,3)], mu_12m=mu_12m[c(1,2,3)], sg=sg, rr=rr[c(1,2,3)],bound=bound)
    levels(db_stage2$treat) <- levels(db_stage1$treat)[c(1,2,3)]
    
    db_stage2$diff6_0<-db_stage2$y_6m-db_stage2$y_0m
    db_stage2$diff12_0<-db_stage2$y_12m-db_stage2$y_0m
    

    pval2 <- c()
    
    if (test=="l"){
      
      for(j in 1:2){
        sub2 <- subset(db_stage2,(db_stage2$treat==levels(db_stage2$treat)[1])+(db_stage2$treat==levels(db_stage2$treat)[j+1])==1)
        mod2 <- lm(y_12m~treat+y_0m, sub2) 
        res2 <- summary(mod2)
        pval2[j] <- pt(coef(res2)[2,3], mod2$df, lower.tail = side)
      }
    }
    
    
    if (test=="t"){
      p12low2 <- t.test(db_stage2$diff12_0[db_stage2$treat!="Medium"]~db_stage2$treat[db_stage2$treat!="Medium"],alternative="greater")$p.value
      p12med2<-t.test(db_stage2$diff12_0[db_stage2$treat!="Low"]~db_stage2$treat[db_stage2$treat!="Low"],alternative="greater")$p.value
      pval2<-cbind(p12low2,p12med2)
    }
    
       
    if (test=="w"){
      p12low2 <- wilcox.test(db_stage2$diff12_0[db_stage2$treat!="Medium"]~db_stage2$treat[db_stage2$treat!="Medium"],alternative="greater")$p.value
      p12med2<-wilcox.test(db_stage2$diff12_0[db_stage2$treat!="Low"]~db_stage2$treat[db_stage2$treat!="Low"],alternative="greater")$p.value
      pval2<-cbind(p12low2,p12med2)
    }
    

    if (test=="w1"){
      p12low2 <- wilcox.test(db_stage2$y_12m[db_stage2$treat!="Medium"]~db_stage2$treat[db_stage2$treat!="Medium"],alternative="greater")$p.value
      p12med2<-wilcox.test(db_stage2$y_12m[db_stage2$treat!="Low"]~db_stage2$treat[db_stage2$treat!="Low"],alternative="greater")$p.value
      pval2<-cbind(p12low2,p12med2)
      conc2low <- wilcox.test(db_stage2$y_12m[db_stage2$treat!="Medium"]~db_stage2$treat[db_stage2$treat!="Medium"],alternative="greater")$statistic/(N2/3*N2/3)
      conc2med<-wilcox.test(db_stage2$y_12m[db_stage2$treat!="Low"]~db_stage2$treat[db_stage2$treat!="Low"],alternative="greater")$statistic/(N2/3*N2/3)
      conc2<-cbind(conc2low,conc2med,NA)
    }
    

    #Avalues <- c(preplan@BJ[7], #H123
    #             preplan@BJ[6], #H12
    #             preplan@BJ[5], #H13
    #             preplan@BJ[3], #H23
    #             preplan@BJ[2], #H2
    #             preplan@BJ[4]  #H1
    #)
    
    #decision <- c()
    #decision[1] <- ifelse(sum(pval2[1] <= Avalues[c(1,2,3,6)])==4, "Reject", "Accept")
    #decision[2] <- ifelse(sum(pval2[2] <= Avalues[c(1,2,4,5)])==4, "Reject", "Accept")
    
    #############################
    
    #global Null H123
    dec123<-max(pval2[1]<preplan@BJ[7]*preplan@Aj[7,1]/(preplan@Aj[7,1]+preplan@Aj[7,2]),
                pval2[2]<preplan@BJ[7]*preplan@Aj[7,2]/(preplan@Aj[7,1]+preplan@Aj[7,2]),preplan@BJ[7]>1)
    #H12
    dec12<-max(pval2[1]<preplan@BJ[6]*preplan@Aj[6,1]/(preplan@Aj[6,1]+preplan@Aj[6,2]),
               pval2[2]<preplan@BJ[6]*preplan@Aj[6,2]/(preplan@Aj[6,1]+preplan@Aj[6,2]),preplan@BJ[6]>1)
    #H13
    dec13<-max(pval2[1]<preplan@BJ[5])
    #H23
    dec23<-max(pval2[2]<preplan@BJ[3])
    #H1
    dec1<-min(pval2[1]<preplan@BJ[4],dec123,dec12,dec13)
    #H2
    dec2<-min(pval2[2]<preplan@BJ[2],dec123,dec12,dec23)

    stage2_arms <- c(1,1,0)
    simdec_output <- c(dec1,dec2,NA)
    
    decision_intersection = dec123
  }
  

  if(sc == 1 ){# this means that low dose was the only promising in the interim analysis
    
    if(sel_scen == 0){ #so we drop both low and median doses and start with high dose only in the second stage.

      N2<-floor(N2/(1+dropout))
      db_stage2 <- sim_dataind(n_arms=2, N = N2, mu_0m = mu_0m[c(1,4)],mu_6m = mu_6m[c(1,4)], mu_12m=mu_12m[c(1,4)], sg=sg, rr=rr[c(1,4)],bound=bound)
      levels(db_stage2$treat) <- c(levels(db_stage1$treat)[c(1)],"High") # c(levels(db_stage1$treat)[c(1,sel)],"High")
      db_stage2$diff6_0<-db_stage2$y_6m-db_stage2$y_0m
      db_stage2$diff12_0<-db_stage2$y_12m-db_stage2$y_0m
      
      pval2 <- c()
      
      
      if (test=="l"){
        
        mod2 <- lm(y_12m~treat+y_0m, db_stage2) 
        res2 <- summary(mod2)
        pval2 <- pt(coef(res2)[2,3], mod2$df, lower.tail = side)
      }
      
      
      if (test=="t"){
        pval2 <- t.test(db_stage2$diff12_0~db_stage2$treat,alternative="greater")$p.value
      }
      
      
      if (test=="w"){
        pval2 <- wilcox.test(db_stage2$diff12_0~db_stage2$treat,alternative="greater")$p.value
      }
      
      if (test=="w1"){
        pval2 <- wilcox.test(db_stage2$y_12m~db_stage2$treat,alternative="greater")$p.value
        conc2<-  c(NA,NA,wilcox.test(db_stage2$y_12m~db_stage2$treat,alternative="greater")$statistic/(N2/2*N2/2))
      }      
      
      
      
      dec123<-max(pval2[1]<preplan@BJ[7],preplan@BJ[7]>1)
      #H12
      dec13<-max(pval2[1]<preplan@BJ[5])
      dec23<-max(pval2[1]<preplan@BJ[3])
      #H1
      dec3<-min(pval2[1]<preplan@BJ[1],dec123,dec13,dec23)
      #
      
      stage2_arms <- c(0,0,1)
      simdec_output <- c(0,0,dec3)
      
      decision_intersection = dec123
      
    }
    
    if(sel_scen == 1){ #continue with low and median doses
      
      N2<-floor(N2/(1+dropout))
      db_stage2 <- sim_dataind(n_arms=3, N = N2, mu_0m =mu_0m[c(1,2,3)], mu_6m =mu_6m[c(1,2,3)], mu_12m=mu_12m[c(1,2,3)], sg=sg, rr=rr[c(1,2,3)],bound=bound)
      levels(db_stage2$treat) <- levels(db_stage1$treat)[c(1,2,3)]
      
      db_stage2$diff6_0<-db_stage2$y_6m-db_stage2$y_0m
      db_stage2$diff12_0<-db_stage2$y_12m-db_stage2$y_0m
      
      
      pval2 <- c()
      
      if (test=="l"){
        
        for(j in 1:2){
          
          sub2 <- subset(db_stage2,(db_stage2$treat==levels(db_stage2$treat)[1])+(db_stage2$treat==levels(db_stage2$treat)[j+1])==1)
          mod2 <- lm(y_12m~treat+y_0m, sub2) 
          res2 <- summary(mod2)
          pval2[j] <- pt(coef(res2)[2,3], mod2$df, lower.tail = side)
        }
      }
      
      
      if (test=="t"){
        p12low2 <- t.test(db_stage2$diff12_0[db_stage2$treat!="Medium"]~db_stage2$treat[db_stage2$treat!="Medium"],alternative="greater")$p.value
        p12med2<-t.test(db_stage2$diff12_0[db_stage2$treat!="Low"]~db_stage2$treat[db_stage2$treat!="Low"],alternative="greater")$p.value
        pval2<-cbind(p12low2,p12med2)
      }
      
      if (test=="w"){
        p12low2 <- wilcox.test(db_stage2$diff12_0[db_stage2$treat!="Medium"]~db_stage2$treat[db_stage2$treat!="Medium"],alternative="greater")$p.value
        p12med2<-wilcox.test(db_stage2$diff12_0[db_stage2$treat!="Low"]~db_stage2$treat[db_stage2$treat!="Low"],alternative="greater")$p.value
        pval2<-cbind(p12low2,p12med2)
      }
      
      if (test=="w1"){
        p12low2 <- wilcox.test(db_stage2$y_12m[db_stage2$treat!="Medium"]~db_stage2$treat[db_stage2$treat!="Medium"],alternative="greater")$p.value
        p12med2<-wilcox.test(db_stage2$y_12m[db_stage2$treat!="Low"]~db_stage2$treat[db_stage2$treat!="Low"],alternative="greater")$p.value
        pval2<-cbind(p12low2,p12med2)
        conc2low <- wilcox.test(db_stage2$y_12m[db_stage2$treat!="Medium"]~db_stage2$treat[db_stage2$treat!="Medium"],alternative="greater")$statistic/(N2/3*N2/3)
        conc2med<-wilcox.test(db_stage2$y_12m[db_stage2$treat!="Low"]~db_stage2$treat[db_stage2$treat!="Low"],alternative="greater")$statistic/(N2/3*N2/3)
        conc2<-cbind(conc2low,conc2med,NA)
      }
      
      
      #global Null H123
      dec123<-max(pval2[1]<preplan@BJ[7]*preplan@Aj[7,1]/(preplan@Aj[7,1]+preplan@Aj[7,2]),
                  pval2[2]<preplan@BJ[7]*preplan@Aj[7,2]/(preplan@Aj[7,1]+preplan@Aj[7,2]),preplan@BJ[7]>1)
      #H12
      dec12<-max(pval2[1]<preplan@BJ[6]*preplan@Aj[6,1]/(preplan@Aj[6,1]+preplan@Aj[6,2]),
                 pval2[2]<preplan@BJ[6]*preplan@Aj[6,2]/(preplan@Aj[6,1]+preplan@Aj[6,2]),preplan@BJ[6]>1)
      #H13
      dec13<-max(pval2[1]<preplan@BJ[5])
      #H23
      dec23<-max(pval2[2]<preplan@BJ[3])
      #H1
      dec1<-min(pval2[1]<preplan@BJ[4],dec123,dec12,dec13)
      #H2
      dec2<-min(pval2[2]<preplan@BJ[2],dec123,dec12,dec23)

      stage2_arms <- c(1,1,0)
      simdec_output <- c(dec1,dec2,NA)
      
      decision_intersection = dec123
    }
  }
  
  if(sc == 3){
    
    N2<-floor(N2/(1+dropout))
    db_stage2 <- sim_dataind(n_arms=3, N=N2, mu_0m=mu_0m[c(1,3,4)], mu_6m=mu_6m[c(1,3,4)], mu_12m=mu_12m[c(1,3,4)], sg=sg, rr=rr[c(1,3,4)],bound=bound)
    
    levels(db_stage2$treat) <- c(levels(db_stage1$treat)[c(1,3)],"High")
    
    db_stage2$diff6_0<-db_stage2$y_6m-db_stage2$y_0m
    db_stage2$diff12_0<-db_stage2$y_12m-db_stage2$y_0m
    
    
    pval2 <- c()
    
    
    if (test=="l"){
      
      for(j in 1:2){
        
        sub2 <- subset(db_stage2,(db_stage2$treat==levels(db_stage2$treat)[1])+(db_stage2$treat==levels(db_stage2$treat)[j+1])==1)
        mod2 <- lm(y_12m~treat+y_0m, sub2) 
        res2 <- summary(mod2)
        pval2[j] <- pt(coef(res2)[2,3], mod2$df, lower.tail = side)
      }
    }
    
    
    if (test=="t"){
      p12med2 <- t.test(db_stage2$diff12_0[db_stage2$treat!="High"]~db_stage2$treat[db_stage2$treat!="High"],alternative="greater")$p.value
      p12hi2<-t.test(db_stage2$diff12_0[db_stage2$treat!="Medium"]~db_stage2$treat[db_stage2$treat!="Medium"],alternative="greater")$p.value
      pval2<-cbind(p12med2,p12hi2)
    }
    
    if (test=="w"){
      p12med2 <- wilcox.test(db_stage2$diff12_0[db_stage2$treat!="High"]~db_stage2$treat[db_stage2$treat!="High"],alternative="greater")$p.value
      p12hi2<-wilcox.test(db_stage2$diff12_0[db_stage2$treat!="Medium"]~db_stage2$treat[db_stage2$treat!="Medium"],alternative="greater")$p.value
      pval2<-cbind(p12med2,p12hi2)
    }
    
    if (test=="w1"){
      p12med2 <- wilcox.test(db_stage2$y_12m[db_stage2$treat!="High"]~db_stage2$treat[db_stage2$treat!="High"],alternative="greater")$p.value
      p12hi2<-wilcox.test(db_stage2$y_12m[db_stage2$treat!="Medium"]~db_stage2$treat[db_stage2$treat!="Medium"],alternative="greater")$p.value
      pval2<-cbind(p12med2,p12hi2)
      conc2med<-wilcox.test(db_stage2$y_12m[db_stage2$treat!="High"]~db_stage2$treat[db_stage2$treat!="High"],alternative="greater")$statistic/(N2/3*N2/3)
      conc2hi <- wilcox.test(db_stage2$y_12m[db_stage2$treat!="Medium"]~db_stage2$treat[db_stage2$treat!="Medium"],alternative="greater")$statistic/(N2/3*N2/3)
      conc2<-cbind(NA,conc2med,conc2hi)
    }
    
    
    
    dec123<-max(pval2[1]<preplan@Aj[7,2],#preplan@BJ[7]*(preplan@Aj[7,2]/(preplan@Aj[7,2]+preplan@Aj[7,3]),
                pval2[2]<(preplan@Aj[7,1]+preplan@Aj[7,3]),#preplan@BJ[7]*preplan@Aj[7,3]/(preplan@Aj[7,2]+preplan@Aj[7,3]),
                preplan@BJ[7]>1)
    #H23
    dec23<-max(pval2[1]<preplan@Aj[3,2],#preplan@BJ[3]*preplan@Aj[3,2]/(preplan@Aj[3,2]+preplan@Aj[3,3]),
               pval2[2]<preplan@Aj[3,3])#,BJ[3]*preplan@Aj[3,3]/(preplan@Aj[3,2]+preplan@Aj[3,3]),preplan@BJ[6]>1)
    #H12
    dec12<-max(pval2[1]<preplan@BJ[6])
    #H13
    dec13<-max(pval2[2]<preplan@BJ[5])
    
    #H2
    dec2<-min(pval2[1]<preplan@BJ[2],dec123,dec23,dec12)
    #H3
    dec3<-min(pval2[2]<preplan@BJ[1],dec123,dec23,dec13)
    #
    
    
    stage2_arms <- c(0,1,1)
    simdec_output <- c(0,dec2,dec3)
    
    decision_intersection = dec123
    }

  
  if(sc == 0){ # in that case start dose 3
    
    N2<-floor(N2/(1+dropout))
    db_stage2 <- sim_dataind(n_arms=2, N=N2, mu_0m=mu_0m[c(1,4)], mu_6m=mu_6m[c(1,4)], mu_12m=mu_12m[c(1,4)], sg=sg, rr=rr[c(1,4)],bound=bound)
    
    levels(db_stage2$treat) <- c(levels(db_stage1$treat)[c(1)],"High")
    
    db_stage2$diff6_0<-db_stage2$y_6m-db_stage2$y_0m
    db_stage2$diff12_0<-db_stage2$y_12m-db_stage2$y_0m
    
    
    if (test=="l"){
      mod2 <- lm(y_12m~treat+y_0m, db_stage2) 
      res2 <- summary(mod2)
      pval2 <- pt(coef(res2)[2,3], mod2$df, lower.tail = side)
    }
    

    if (test=="t"){
      pval2 <- t.test(db_stage2$diff12_0~db_stage2$treat,alternative="greater")$p.value
    }
    
    if (test=="w"){
      pval2 <- wilcox.test(db_stage2$diff12_0~db_stage2$treat,alternative="greater")$p.value
    }
    
    if (test=="w1"){
      pval2 <- wilcox.test(db_stage2$y_12m~db_stage2$treat,alternative="greater")$p.value
      conc2<-c(NA,NA,wilcox.test(db_stage2$y_12m~db_stage2$treat,alternative="greater")$statistic/(N2/2*N2/2))
    }
    
    
    dec123<-max(pval2[1]<preplan@BJ[7],preplan@BJ[7]>1)
    #H12
    dec13<-max(pval2[1]<preplan@BJ[5])
    dec23<-max(pval2[1]<preplan@BJ[3])
    #H1
    dec3<-min(pval2[1]<preplan@BJ[1],dec123,dec13,dec23)

    stage2_arms <- c(0,0,1)
    simdec_output <- c(0,0,dec3)
    
    decision_intersection = dec123
  }
  
  
  #Computation of concordance and confidence interval for test ==w1
  conc<-c(NA,NA,NA)       #unconditional concordance
  concCI<-c(NA,NA,NA)     #unconditional CI concordance
  conccond<-c(NA,NA,NA)   #conditional concordance
  concCIcond<-c(NA,NA,NA) #unconditional CI concordance
  concCIinvn<-c(NA,NA,NA) #concordance inverse normal method
  concinvn<-c(NA,NA,NA)   #concordance CI inverse normal method
  
  if (test=="w1")
  {
    if (sc==2) 
    {            
      conc    <-c((conc1[1:2]*N1/3+conc2[1:2]*N2/3)/(N1/3+N2/3),NA) #low and medium
      conccond<-c((conc1[1:2]*N1/3+conc2[1:2]*N2/3)/(N1/3+N2/3),NA) #low and medium
      concCI<-c(uniroot(CI2_help, lower = 0, upper = .99, tol = .Machine$double.eps^0.25, maxiter = 10000,thetahat=conc[1],alpha=alpha,N1=N1,N2=N2)$root, 
                uniroot(CI2_help, lower = 0, upper = .99, tol = .Machine$double.eps^0.25, maxiter = 10000,thetahat=conc[2],alpha=alpha,N1=N1,N2=N2)$root, 
                NA)
      concCIcond<-concCI
      concCIinvn<-c(uniroot(CIinvn_help, lower = 0, upper = .99, tol = .Machine$double.eps^0.25, maxiter = 10000,
                            thetahat1=conc1[1],thetahat2=conc2[1],value=alpha,N1=N1,N2=N2orig)$root,
                    uniroot(CIinvn_help, lower = 0, upper = .99, tol = .Machine$double.eps^0.25, maxiter = 10000,
                            thetahat1=conc1[2],thetahat2=conc2[2],value=alpha,N1=N1,N2=N2orig)$root,
                NA)
      #Median estimate:
      concinvn<-c(uniroot(CIinvn_help, lower = 0, upper = .99, tol = .Machine$double.eps^0.25, maxiter = 10000,
                          thetahat1=conc1[1],thetahat2=conc2[1],value=0.5,N1=N1,N2=N2orig)$root,
                    uniroot(CIinvn_help, lower = 0, upper = .99, tol = .Machine$double.eps^0.25, maxiter = 10000,
                            thetahat1=conc1[2],thetahat2=conc2[2],value=0.5,N1=N1,N2=N2orig)$root,
                    NA)
      
    }
  
    if (sc==1) 
    { 
      conc<-c((conc1[1:2]*N1/3+conc2[1:2]*N2/3)/(N1/3+N2/3),NA) #low and medium
      conccond<-c((conc1[1:2]*N1/3+conc2[1:2]*N2/3)/(N1/3+N2/3),NA) #low and medium
      concCI<-c(uniroot(CI2_help, lower = 0, upper = .99, tol = .Machine$double.eps^0.25, maxiter = 10000,thetahat=conc[1],alpha=alpha,N1=N1,N2=N2)$root,
                uniroot(CI2_help, lower = 0, upper = .99, tol = .Machine$double.eps^0.25, maxiter = 10000,thetahat=conc[2],alpha=alpha,N1=N1,N2=N2)$root,
                NA)
      concCIcond<-concCI
      concCIinvn<-c(uniroot(CIinvn_help, lower = 0, upper = .99, tol = .Machine$double.eps^0.25, maxiter = 10000,
                            thetahat1=conc1[1],thetahat2=conc2[1],value=alpha,N1=N1,N2=N2orig)$root,
                    uniroot(CIinvn_help, lower = 0, upper = .99, tol = .Machine$double.eps^0.25, maxiter = 10000,
                            thetahat1=conc1[2],thetahat2=conc2[2],value=alpha,N1=N1,N2=N2orig)$root,
                       NA)
      #Median estimate:
      concinvn<-c(uniroot(CIinvn_help, lower = 0, upper = .99, tol = .Machine$double.eps^0.25, maxiter = 10000,
                          thetahat1=conc1[1],thetahat2=conc2[1],value=0.5,N1=N1,N2=N2orig)$root,
                  uniroot(CIinvn_help, lower = 0, upper = .99, tol = .Machine$double.eps^0.25, maxiter = 10000,
                          thetahat1=conc1[2],thetahat2=conc2[2],value=0.5,N1=N1,N2=N2orig)$root,
                  NA)

    }
  
    if (sc==3) 
    {            
      conc<-c(conc1[1],(conc1[2]*N1/3+conc2[2]*N2/3)/(N1/3+N2/3),conc2[3]) #medium and high
      conccond<-c(NA,conc[2:3]) #medium and high
      concCI<-c(uniroot(CI_help, lower = 0, upper = .99, tol = .Machine$double.eps^0.25, maxiter = 10000,thetahat=conc[1],N1=N1,N2=N2,alpha=alpha)$root,
                uniroot(CI2_help, lower = 0, upper = .99, tol = .Machine$double.eps^0.25, maxiter = 10000,thetahat=conc[2],alpha=alpha,N1=N1,N2=N2)$root,
                uniroot(CI_help, lower = 0, upper = .99, tol = .Machine$double.eps^0.25, maxiter = 10000,thetahat=conc[3],N1=N2,N2=N2,alpha=alpha)$root)
      concCIcond<-c(NA,
                    uniroot(CI2_help, lower = 0, upper = .99, tol = .Machine$double.eps^0.25, maxiter = 10000,thetahat=conc[2],alpha=alpha,N1=N1,N2=N2)$root,
                    uniroot(CI_help, lower = 0, upper = .99, tol = .Machine$double.eps^0.25, maxiter = 10000,thetahat=conc[3],N1=N2,N2=N2,alpha=alpha)$root)
      concCIinvn<-c(NA,
                    uniroot(CIinvn_help, lower = 0, upper = .99, tol = .Machine$double.eps^0.25, maxiter = 10000,
                            thetahat1=conc1[2],thetahat2=conc2[2],value=alpha,N1=N1,N2=N2orig)$root,
                    uniroot(CI_help, lower = 0, upper = .99, tol = .Machine$double.eps^0.25, maxiter = 10000,thetahat=conc[3],N1=N2,N2=N2,alpha=alpha)$root)
      #Median estimate:
      concinvn<-c(NA,
                  uniroot(CIinvn_help, lower = 0, upper = .99, tol = .Machine$double.eps^0.25, maxiter = 10000,
                          thetahat1=conc1[2],thetahat2=conc2[2],value=0.5,N1=N1,N2=N2orig)$root,
                  conc2[3])#uniroot(CI_help, lower = 0, upper = .99, tol = .Machine$double.eps^0.25, maxiter = 10000,thetahat=conc[3],N1=N2,N2=N2,alpha=alpha)$root)
      
    }
  
    if (sc==0) 
    {            
      conc<-c(conc1[1:2],conc2[3]) #high
      conccond<-c(NA,NA,conc2[3]) #high
      concCI<-c(uniroot(CI_help, lower = .01, upper = .99, tol = .Machine$double.eps^0.25, maxiter = 10000,thetahat=conc[1],N1=N1,N2=N1,alpha=alpha)$root,
               uniroot(CI_help, lower = .01, upper = .99, tol = .Machine$double.eps^0.25, maxiter = 10000,thetahat=conc[2],N1=N1,N2=N1,alpha=alpha)$root,
               uniroot(CI_help, lower = .01, upper = .99, tol = .Machine$double.eps^0.25, maxiter = 10000,thetahat=conc[3],N1=N2,N2=N2,alpha=alpha)$root)
      concCIcond<-c(NA,
                NA,
                uniroot(CI_help, lower = .01, upper = .99, tol = .Machine$double.eps^0.25, maxiter = 10000,thetahat=conc[3],N1=N2,N2=N2,alpha=alpha)$root)  
      concCIinvn<-concCIcond
      concinvn<-conccond
    }
  }
  
  
  ################################################
  # Multi-arm trial 1
  ####
  #Part 1 for low and median dose
  db_stage_ma1a <- sim_dataind(n_arms = n_arms-1, N = floor((N1+N2)/5*3), mu_0m = mu_0m[1:n_arms-1], mu_6m = mu_6m[1:n_arms-1], mu_12m = mu_12m[1:n_arms-1], sg = sg, rr=rr[1:n_arms-1],bound=bound)

  db_stage_ma1a$treat <- relevel(db_stage_ma1a$treat, ref = "Placebo")
  
  db_stage_ma1a$diff6_0<-db_stage_ma1a$y_6m-db_stage_ma1a$y_0m
  db_stage_ma1a$diff12_0<-db_stage_ma1a$y_12m-db_stage_ma1a$y_0m
  
  
  #Part 2 for high dose
  db_stage_ma1b <- sim_dataind(n_arms = 2, N = floor((N1+N2)/5*2), mu_0m = mu_0m[c(1,n_arms)], mu_6m = mu_6m[c(1,n_arms)], mu_12m = mu_12m[c(1,n_arms)], sg = sg, rr=rr[c(1,n_arms)],bound=bound)
  
  db_stage_ma1b$treat <- relevel(db_stage_ma1b$treat, ref = "Placebo")
  
  db_stage_ma1b$diff6_0<-db_stage_ma1b$y_6m-db_stage_ma1b$y_0m
  db_stage_ma1b$diff12_0<-db_stage_ma1b$y_12m-db_stage_ma1b$y_0m
  
  
  
  if (test=="l")
  {
  mod_ma1a <- aov(y_12m ~ treat+y_0m, db_stage_ma1a)
  model_dunnett_ma1a = summary(glht(model = mod_ma1a, linfct=mcp(treat="Dunnett"), alternative = "less"))
  pval_dunnett_ma1a = model_dunnett_ma1a$test$pvalues
  
  #part2
  mod_ma1b <- aov(y_12m ~ treat+y_0m, db_stage_ma1b)
  model_dunnett_ma1b = summary(glht(model = mod_ma1b, linfct=mcp(treat="Dunnett"), alternative = "less"))
  pval_dunnett_ma1b = model_dunnett_ma1b$test$pvalues
  
  decision_ma1<-c((pval_dunnett_ma1a<=(alpha/3*2))*1,(pval_dunnett_ma1b<=(alpha/3))*1)
  
  decision_ma1a<-c((pval_dunnett_ma1a<=(alpha/3*2))*1,(1-max(decision_ma1[1:2]))*((pval_dunnett_ma1b<=(alpha/3))*1))
  }
  
  if (test=="t"){
    pma1loa<-t.test(db_stage_ma1a$diff12_0[db_stage_ma1a$treat!="Medium"]~db_stage_ma1a$treat[db_stage_ma1a$treat!="Medium"],alternative="greater")$p.value
    pma1mea<-t.test(db_stage_ma1a$diff12_0[db_stage_ma1a$treat!="Low"]~db_stage_ma1a$treat[db_stage_ma1a$treat!="Low"],alternative="greater")$p.value
    pvalma1a<-cbind(pma1loa,pma1mea)
    
    pma1hi<-t.test(db_stage_ma1b$diff12_0~db_stage_ma1b$treat,alternative="greater")$p.value
    decision_ma1<-c((p.adjust(pvalma1a,"holm")<=(alpha/3*2))*1,(p.adjust(pma1hi,"holm")<=(alpha/3))*1)
    
    decision_ma1a<-c((p.adjust(pvalma1a,"holm")<=(alpha/3*2))*1,(1-max(decision_ma1[1:2]))*((p.adjust(pma1hi,"holm")<=(alpha/3))*1))
    
  }
  
  if (test=="w"){
    pma1loa<-wilcox.test(db_stage_ma1a$diff12_0[db_stage_ma1a$treat!="Medium"]~db_stage_ma1a$treat[db_stage_ma1a$treat!="Medium"],alternative="greater")$p.value
    pma1mea<-wilcox.test(db_stage_ma1a$diff12_0[db_stage_ma1a$treat!="Low"]~db_stage_ma1a$treat[db_stage_ma1a$treat!="Low"],alternative="greater")$p.value
    pvalma1a<-cbind(pma1loa,pma1mea)
    
    pma1hi<-wilcox.test(db_stage_ma1b$diff12_0~db_stage_ma1b$treat,alternative="greater")$p.value
    
    decision_ma1<-c((p.adjust(pvalma1a,"holm")<=(alpha/3*2))*1,(p.adjust(pma1hi,"holm")<=(alpha/3))*1)
    decision_ma1a<-c((p.adjust(pvalma1a,"holm")<=(alpha/3*2))*1,(1-max(decision_ma1[1:2]))*((p.adjust(pma1hi,"holm")<=(alpha/3))*1))
  }
  
  
  if (test=="w1"){
    pma1loa<-wilcox.test(db_stage_ma1a$y_12m[db_stage_ma1a$treat!="Medium"]~db_stage_ma1a$treat[db_stage_ma1a$treat!="Medium"],alternative="greater")$p.value
    pma1mea<-wilcox.test(db_stage_ma1a$y_12m[db_stage_ma1a$treat!="Low"]~db_stage_ma1a$treat[db_stage_ma1a$treat!="Low"],alternative="greater")$p.value
    pvalma1a<-cbind(pma1loa,pma1mea)
    
    pma1hi<-wilcox.test(db_stage_ma1b$y_12m~db_stage_ma1b$treat,alternative="greater")$p.value
    
    decision_ma1<-c((p.adjust(pvalma1a,"holm")<=(alpha/3*2))*1,(p.adjust(pma1hi,"holm")<=(alpha/3))*1)
    decision_ma1a<-c((p.adjust(pvalma1a,"holm")<=(alpha/3*2))*1,(1-max(decision_ma1[1:2]))*((p.adjust(pma1hi,"holm")<=(alpha/3))*1))
  }
  
  concMA1lo<-wilcox.test(db_stage_ma1a$y_12m[db_stage_ma1a$treat!="Medium"]~db_stage_ma1a$treat[db_stage_ma1a$treat!="Medium"],alternative="greater")$statistic/(floor((N1+N2)/15*3)*floor((N1+N2)/15*3))
  concMA1me<-wilcox.test(db_stage_ma1a$y_12m[db_stage_ma1a$treat!="Low"]~db_stage_ma1a$treat[db_stage_ma1a$treat!="Low"],alternative="greater")$statistic/(floor((N1+N2)/15*3)*floor((N1+N2)/15*3))
  concMA1hi<-wilcox.test(db_stage_ma1b$y_12m~db_stage_ma1b$treat,alternative="greater")$statistic/(floor((N1+N2)/10*2)*floor((N1+N2)/10*2))
  concMA1<-vector(length=3)
  concMA1<-c(concMA1lo,concMA1me,concMA1hi)
  
  
  
  
  ################################################
  # Multi-arm trial 2
  db_stage_ma2 <- sim_dataind(n_arms = n_arms, N = N1+N2, mu_0m = mu_0m, mu_6m = mu_6m, mu_12m = mu_12m, sg = sg, rr=rr,bound=bound)
  
  db_stage_ma2$treat <- relevel(db_stage_ma2$treat, ref = "Placebo")
  
  db_stage_ma2$diff6_0<-db_stage_ma2$y_6m-db_stage_ma2$y_0m
  db_stage_ma2$diff12_0<-db_stage_ma2$y_12m-db_stage_ma2$y_0m
  
  if (test=="l"){
    mod_ma2 <- aov(y_12m ~ treat+y_0m, db_stage_ma2)
    model_dunnett_ma2 = summary(glht(model = mod_ma2, linfct=mcp(treat="Dunnett"), alternative = "less"))
    pval_dunnett_ma2 = model_dunnett_ma2$test$pvalues
  
    decision_ma2<-(pval_dunnett_ma2<=alpha)*1
  }
  
  if (test=="t"){
    pma1loa<-t.test(db_stage_ma2$diff12_0[db_stage_ma2$treat %in% c("Low","Placebo")]~db_stage_ma2$treat[db_stage_ma2$treat %in% c("Low","Placebo")],alternative="greater")$p.value
    pma1mea<-t.test(db_stage_ma2$diff12_0[db_stage_ma2$treat%in% c("Medium","Placebo")]~db_stage_ma2$treat[db_stage_ma2$treat %in% c("Medium","Placebo")],alternative="greater")$p.value
    pma1hia<-t.test(db_stage_ma2$diff12_0[db_stage_ma2$treat %in% c("High","Placebo")]~db_stage_ma2$treat[db_stage_ma2$treat %in% c("High","Placebo")],alternative="greater")$p.value
                                                                                                                                                                                                    
    pvalma1a<-cbind(pma1loa,pma1mea)
    pvalma1a
    decision_ma2<-(p.adjust(c(pma1loa,pma1mea,pma1hia),"holm")<alpha)*1
  }
  
  if (test=="w"){
    pma1loa<-wilcox.test(db_stage_ma2$diff12_0[db_stage_ma2$treat %in% c("Low","Placebo")]~db_stage_ma2$treat[db_stage_ma2$treat %in% c("Low","Placebo")],alternative="greater")$p.value
    pma1mea<-wilcox.test(db_stage_ma2$diff12_0[db_stage_ma2$treat%in% c("Medium","Placebo")]~db_stage_ma2$treat[db_stage_ma2$treat %in% c("Medium","Placebo")],alternative="greater")$p.value
    pma1hia<-wilcox.test(db_stage_ma2$diff12_0[db_stage_ma2$treat %in% c("High","Placebo")]~db_stage_ma2$treat[db_stage_ma2$treat %in% c("High","Placebo")],alternative="greater")$p.value
    decision_ma2<-(p.adjust(c(pma1loa,pma1mea,pma1hia),"holm")<alpha)*1
  }
  
  if (test=="w1"){
    pma1loa<-wilcox.test(db_stage_ma2$y_12m[db_stage_ma2$treat %in% c("Low","Placebo")]~db_stage_ma2$treat[db_stage_ma2$treat %in% c("Low","Placebo")],alternative="greater")$p.value
    pma1mea<-wilcox.test(db_stage_ma2$y_12m[db_stage_ma2$treat%in% c("Medium","Placebo")]~db_stage_ma2$treat[db_stage_ma2$treat %in% c("Medium","Placebo")],alternative="greater")$p.value
    pma1hia<-wilcox.test(db_stage_ma2$y_12m[db_stage_ma2$treat %in% c("High","Placebo")]~db_stage_ma2$treat[db_stage_ma2$treat %in% c("High","Placebo")],alternative="greater")$p.value
    decision_ma2<-(p.adjust(c(pma1loa,pma1mea,pma1hia),"holm")<alpha)*1
  }
  

  concMAlo<-wilcox.test(db_stage_ma2$y_12m[db_stage_ma2$treat %in% c("Low","Placebo")]~db_stage_ma2$treat[db_stage_ma2$treat %in% c("Low","Placebo")],alternative="greater")$statistic/(floor((N1+N2)/4)*floor((N1+N2)/4))
  concMAme<-wilcox.test(db_stage_ma2$y_12m[db_stage_ma2$treat%in% c("Medium","Placebo")]~db_stage_ma2$treat[db_stage_ma2$treat %in% c("Medium","Placebo")],alternative="greater")$statistic/(floor((N1+N2)/4)*floor((N1+N2)/4))
  concMAhi<-wilcox.test(db_stage_ma2$y_12m[db_stage_ma2$treat %in% c("High","Placebo")]~db_stage_ma2$treat[db_stage_ma2$treat %in% c("High","Placebo")],alternative="greater")$statistic/(floor((N1+N2)/4)*floor((N1+N2)/4))
  concMA2<-vector(length=3)
  concMA2<-c(concMAlo,concMAme,concMAhi)
  
  #######################################
  
    res_intersection=ifelse(decision_intersection == "Reject", 1,0)
    return(list(stage2_arms=stage2_arms,
                selected_dose=sel,
                simdec_output=simdec_output,
                res_intersection=res_intersection,
                decision_ma1=decision_ma1,
                decision_ma2=decision_ma2,
                decision_ma1a=decision_ma1a,
                conc1=conc1,
                conc2=conc2,
                conc=conc,
                concMA1=concMA1,
                concMA2=concMA2,
                concCI=concCI,
                conccond=conccond,
                concCIcond=concCIcond,
                concCIinvn=concCIinvn,
                concinvn=concinvn
                ))
  
}


