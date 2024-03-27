#' Simulate data from a multi-arm multi-stage trial with shared control and two initial doses, where an additional dose could be added after the interim analysis
#' @description Function to simulate trial data (2-stages, with dose selection). The analyses are performed using partial conditional error.
#'
#' @param n_arms number of arms (including control)
#' @param N1 sample size stage 1
#' @param N2 sample size stage 2
#' @param mu_0m Baseline value per arm (vector of length `n_arm`)
#' @param mu_6m 6-month value per arm (vector of length `n_arm`)
#' @param mu_12m 12-month value per arm (vector of length `n_arm`)
#' @param sigma covariance matrix between 6- and 12-month mean differences assumed equal across arms (matrix of dim 2x2)
#' @param alpha1 significance level for dose selection (futility analysis)
#' @param alpha significance level for selected dose vs control comparison
#' @param v vector giving the proportions of pre-planned measurements collected up to the interim analysis
#' @param rmonth recruitment per month (recruitment speed assumed constant over time)
#' @param sel_scen choose between two different options in case that in interim analysis low dose is promising, but median dose not: 0: do not continue with low dose or median dose; 1: continue with low and median doses
#' @param sim_out Option for simplified output for simulations (if `sim_out=TRUE` simplified version, the value is `FALSE` by default)
#' @param side TRUE/FALSE referring to the side for 1-side testing (if TRUE then lower = side)
#' @param analysis defines type of analysis: "t" calculates a t-test, "l" a linear model with baseline values as covariables, and "m" defines a mixture model, "w" Wilcoxon test of differences
#' @param dropout dropoutrate, between 0 and 1
#' @param rr responder rate for each dose, which gives the proportion of patients with value 0 at follow-up
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

library(lme4)
library(lmerTest)
library(multcomp)
# Function to simulate trial data (2-stages, with dose selection)
# individual observations are simulated

#n_arms = 4; N1=120; N2=80; mu_0m=c(6.9,6.9,6.9,6.9); mu_6m=c(6.9,5.9,5.9,5.9); mu_12m=c(6.9,5.9,5.9,5.9); sg=matrix(c(1,.6,.4,0.6,1,.6,0.4,0.6,1),3); rmonth=1; alpha1 = 0.1; alpha = 0.025; sim_out=T;test="w";dropout=0.1;rr=c(0,.2,.4,.6)
#sim_trial_pceind_test (n_arms, N1 , N2, mu_0m, mu_6m, mu_12m, sg, rmonth, alpha1, alpha, sim_out,sel_scen,side,test,dropout,rr)
#sim_trial_pceind_test(n_arms = 4,N1 = 90 , N2 = 60, mu_0m = mu_0m,   mu_6m = mu_6m, mu_12m = mu_12m,  sg = sg, rmonth=1, alpha1 = .1, alpha = .025, sim_out=T,sel_scen=0, side=T,test="l",dropout=.1,rr=rep(0,4))
sim_trial_pceind_test <- function(n_arms = 4, N1 , N2, mu_0m, mu_6m, mu_12m, sg, rmonth, alpha1 , alpha = 0.025, sim_out=F, sel_scen, side=T,test,dropout,rr,bound)
{
  N1<-floor(N1/(1+dropout))
  db_stage1 <- sim_dataind(n_arms = n_arms-1, N = N1, mu_0m = mu_0m[1:n_arms-1], mu_6m = mu_6m[1:n_arms-1], mu_12m = mu_12m[1:n_arms-1], sg = sg, rmonth = rmonth,rr=rr[1:n_arms-1],bound=bound)
  recruit_time1 <- max(db_stage1$recruit_time)
  
  
  
  db_stage1$treat <- relevel(db_stage1$treat, ref = "Placebo")
  
  db_stage1$diff6_0<-db_stage1$y_6m-db_stage1$y_0m
  db_stage1$diff12_0<-db_stage1$y_12m-db_stage1$y_0m
  
  
  #interim anlysis
  #############
  
  #frame.sim <- subset(db_stage1, treat != "High")  #why do I need this?
  #LowvsC <- lm(diff6_0 ~ treat, data=subset(frame.sim, treat!="Medium")) # adjust with the baseline please
  #MedvsC <- lm(diff6_0  ~ treat , data=subset(frame.sim, treat!="Low"))
  
  
  #Obtain one-sided p-value
  #########################
  
  #plow <- pt(coef(summary(LowvsC))[, 3], LowvsC$df, lower = side)[2]
  plow <- t.test(db_stage1$diff6_0[db_stage1$treat!="Medium"]~db_stage1$treat[db_stage1$treat!="Medium"],alternative="greater")$p.value
  
  #pmed <- pt(coef(summary(MedvsC))[, 3], MedvsC$df, lower = side)[2]
  pmed<-t.test(db_stage1$diff6_0[db_stage1$treat!="Low"]~db_stage1$treat[db_stage1$treat!="Low"],alternative="greater")$p.value
  
  
  if (test=="w"){
    plow <- wilcox.test(db_stage1$diff6_0[db_stage1$treat!="Medium"]~db_stage1$treat[db_stage1$treat!="Medium"],alternative="greater")$p.value
    pmed<-wilcox.test(db_stage1$diff6_0[db_stage1$treat!="Low"]~db_stage1$treat[db_stage1$treat!="Low"],alternative="greater")$p.value
    }
  
  if (test=="w1"){
    plow <- wilcox.test(db_stage1$y_6m[db_stage1$treat!="Medium"]~db_stage1$treat[db_stage1$treat!="Medium"],alternative="greater")$p.value
    pmed<-wilcox.test(db_stage1$y_6m[db_stage1$treat!="Low"]~db_stage1$treat[db_stage1$treat!="Low"],alternative="greater")$p.value
  }
  
  if (test=="f"){  #we select according to wilcoxon also for fisher option
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
  if( (sum(pval.surr < alpha1) == 1 & which.min(pval.surr) == 1) ){  #low has p<alpha1, median has p>alpha1, both are NOT selected for second stage
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
    #pval
  }
  
  if (test=="l"){  #two individual models, due to robustness in case of different variances
    
    for(j in 1:(n_arms-2)){
      
      sub1 <- subset(db_stage1,(db_stage1$treat==levels(db_stage1$treat)[1])+ (db_stage1$treat==levels(db_stage1$treat)[j+1])==1)
      mod1 <- lm(diff12_0 ~ treat+y_0m, sub1)
      res1 <- summary(mod1)
      pval1[j] <- pt(coef(res1)[2,3], mod1$df, lower.tail = side)
    }
  }
  
  if (test=="m"){
    db_stage1$patID<-c(1:N1)
    db_stage1long<-reshape(data=db_stage1, direction = "long",v.names="diffmf", idvar = "patID", times=c(6,12),varying=c("diff6_0","diff12_0"),timevar="month")
    db_stage1long$month<-factor(db_stage1long$month,c(6,12))
    
    md0 <- lmer(diffmf ~ y_0m + treat+ month+(1 | patID), data = db_stage1long)
    pval1<-summary(md0)$coefficients[3:4,5]
  }
  
  if (test=="w"){
    p12low <- wilcox.test(db_stage1$diff12_0[db_stage1$treat!="Medium"]~db_stage1$treat[db_stage1$treat!="Medium"],alternative="greater")$p.value
    p12med<-wilcox.test(db_stage1$diff12_0[db_stage1$treat!="Low"]~db_stage1$treat[db_stage1$treat!="Low"],alternative="greater")$p.value
    pval1<-cbind(p12low,p12med)
    #pval
  }
  
  if (test=="w1"){
    p12low <- wilcox.test(db_stage1$y_12m[db_stage1$treat!="Medium"]~db_stage1$treat[db_stage1$treat!="Medium"],alternative="greater")$p.value
    p12med<-wilcox.test(db_stage1$y_12m[db_stage1$treat!="Low"]~db_stage1$treat[db_stage1$treat!="Low"],alternative="greater")$p.value
    pval1<-cbind(p12low,p12med)
    #pval
  }
  
  if (test=="f"){
    low1<-table(db_stage1$y_12m[db_stage1$treat!="Medium"]<=bound,db_stage1$treat[db_stage1$treat!="Medium"])[,c(1,2)]
    med1<-table(db_stage1$y_12m[db_stage1$treat!="Low"]<=bound,db_stage1$treat[db_stage1$treat!="Low"])[,c(1,3)]
    pvallow1<-try(fisher.test(low1,alternative="greater")$p.value,silent=TRUE)
    if (is.numeric(pvallow1)==FALSE) {pvallow1<-1}
    pvalmed1<-try(fisher.test(med1,alternative="greater")$p.value,silent=TRUE)
    if (is.numeric(pvalmed1)==FALSE) {pvalmed1<-1}
    pval1<-c(pvallow1,pvalmed1)
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
  # preplan@Aj
  # preplan@BJ
  
  #######################################
  # stage2
  # sc=2 --> Arm A and B continue to stage 2
  # sc=1 --> Arm A and B stop and start Arm C #if Arm B should be drop then drop Arm A too, start the second stage with only the Arm C Arm A or B continue to stage 2
  # sc=0 --> Arm A and B stop and start Arm C
  # sc=3 --> Arm A stop continue with Arm B and start Arm C
  
  
  if(sc==2){  # Arm A and B continue to stage 2
    N2<-floor(N2/(1+dropout))
    db_stage2 <- sim_dataind(n_arms=3, N = N2, mu_0m =mu_0m[c(1,2,3)], mu_6m =mu_6m[c(1,2,3)], mu_12m=mu_12m[c(1,2,3)], sg=sg, rmonth=rmonth,rr=rr[c(1,2,3)],bound=bound)
    levels(db_stage2$treat) <- levels(db_stage1$treat)[c(1,2,3)]
    recruit_time2 <- max(db_stage2$recruit_time)
    
    db_stage2$diff6_0<-db_stage2$y_6m-db_stage2$y_0m
    db_stage2$diff12_0<-db_stage2$y_12m-db_stage2$y_0m
    
    
    
    pval2 <- c()
    
    if (test=="l"){
      
      for(j in 1:2){
        
        sub2 <- subset(db_stage2,(db_stage2$treat==levels(db_stage2$treat)[1])+(db_stage2$treat==levels(db_stage2$treat)[j+1])==1)
        mod2 <- lm(diff12_0~treat+y_0m, sub2) 
        res2 <- summary(mod2)
        pval2[j] <- pt(coef(res2)[2,3], mod2$df, lower.tail = side)
      }
    }
    
    
    if (test=="t"){
      p12low2 <- t.test(db_stage2$diff12_0[db_stage2$treat!="Medium"]~db_stage2$treat[db_stage2$treat!="Medium"],alternative="greater")$p.value
      p12med2<-t.test(db_stage2$diff12_0[db_stage2$treat!="Low"]~db_stage2$treat[db_stage2$treat!="Low"],alternative="greater")$p.value
      pval2<-cbind(p12low2,p12med2)
    }
    
    
    if (test=="m"){
      db_stage2$patID<-c(1:N2)
      db_stage2long<-reshape(data=db_stage2, direction = "long",v.names="diffmf", idvar = "patID", times=c(6,12),varying=c("diff6_0","diff12_0"),timevar="month")
      db_stage2long$month<-factor(db_stage2long$month,c(6,12))
      
      md0 <- lmer(diffmf ~ y_0m + treat+ month+(1 | patID), data = db_stage2long)
      pval2<-summary(md0)$coefficients[3:4,5]
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
    }
    
    if (test=="f"){
      low2<-table(db_stage2$y_12m[db_stage2$treat!="Medium"]<=bound,db_stage2$treat[db_stage2$treat!="Medium"])[,c(1,2)]
      med2<-table(db_stage2$y_12m[db_stage2$treat!="Low"]<=bound,db_stage2$treat[db_stage2$treat!="Low"])[,c(1,3)]
      pvallow2<-try(fisher.test(low2,alternative="greater")$p.value,silent=TRUE)
      if (is.numeric(pvallow2)==FALSE) {pvallow2<-1}
      pvalmed2<-try(fisher.test(med2,alternative="greater")$p.value,silent=TRUE)
      if (is.numeric(pvalmed2)==FALSE) {pvalmed2<-1}
      pval2<-c(pvallow2,pvalmed2)
      #pval2<-c(fisher.test(low2,alternative="greater")$p.value,fisher.test(med2,alternative="greater")$p.value)
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
    #
    
    
    
    #decision_stage2 <- data.frame(c(dec1,dec2), row.names = levels(db_stage2$treat)[2:3])
    
    #pval2df <- data.frame(c(pval2), row.names = levels(db_stage2$treat)[2:3])
    
    stage2_arms <- c(1,1,0)
    simdec_output <- c(dec1,dec2,NA)#c(ifelse(decision[1]=="Reject", 1, 0),                       ifelse(decision[2]=="Reject", 1, 0),                       NA)
    
    decision_intersection = dec123#ifelse(sum(pval2 <= Avalues[1]) > 0, "Reject", "Accept")
  }
  
  
  
  if(sc == 1 ){# this means that low dose was the only promising in the interim analysis
    
    if(sel_scen == 0){ #so we drop both low and median doses and start with high dose only in the second stage.
      N2<-floor(N2/(1+dropout))
      
      db_stage2 <- sim_dataind(n_arms=2, N = N2, mu_0m = mu_0m[c(1,4)],mu_6m = mu_6m[c(1,4)], mu_12m=mu_12m[c(1,4)], sg=sg, rmonth=rmonth,rr=rr[c(1,4)],bound=bound)
      
      levels(db_stage2$treat) <- c(levels(db_stage1$treat)[c(1)],"High") # c(levels(db_stage1$treat)[c(1,sel)],"High")
      recruit_time2 <- max(db_stage2$recruit_time)
      
      db_stage2$diff6_0<-db_stage2$y_6m-db_stage2$y_0m
      db_stage2$diff12_0<-db_stage2$y_12m-db_stage2$y_0m
      
      
      pval2 <- c()
      
      #mod2 <- lm(diff12_0 ~ treat, db_stage2) #are we using this model or should we use individual models?
      #res2 <- summary(mod2)
      #pval2 <- pt(coef(res2)[2,3], mod2$df, lower.tail = side)
      
      if (test=="l"){
        
        mod2 <- lm(diff12_0~treat+y_0m, db_stage2) 
        res2 <- summary(mod2)
        pval2 <- pt(coef(res2)[2,3], mod2$df, lower.tail = side)
      }
      
      
      
      if (test=="t"){
        pval2 <- t.test(db_stage2$diff12_0~db_stage2$treat,alternative="greater")$p.value
      }
      
      
      if (test=="m"){
        db_stage2$patID<-c(1:N2)
        db_stage2long<-reshape(data=db_stage2, direction = "long",v.names="diffmf", idvar = "patID", times=c(6,12),varying=c("diff6_0","diff12_0"),timevar="month")
        db_stage2long$month<-factor(db_stage2long$month,c(6,12))
        
        md0 <- lmer(diffmf ~ y_0m + treat+ month+(1 | patID), data = db_stage2long)
        pval2<-summary(md0)$coefficients[3,5]
      }
      
      if (test=="w"){
        pval2 <- wilcox.test(db_stage2$diff12_0~db_stage2$treat,alternative="greater")$p.value
      }
      
      if (test=="w1"){
        pval2 <- wilcox.test(db_stage2$y_12m~db_stage2$treat,alternative="greater")$p.value
      }      
      
      if (test=="f"){
        high2<-table(db_stage2$y_12m<=bound,db_stage2$treat)
        pval2<-try(fisher.test(high2,alternative="greater")$p.value,silent=TRUE)
        if (is.numeric(pval2)==FALSE) {pval2<-1}
        }
      
      
      
      dec123<-max(pval2[1]<preplan@BJ[7],preplan@BJ[7]>1)
      #H12
      dec13<-max(pval2[1]<preplan@BJ[5])
      dec23<-max(pval2[1]<preplan@BJ[3])
      #H1
      dec3<-min(pval2[1]<preplan@BJ[1],dec123,dec13,dec23)
      #
      
      
      
      #decision_stage2 <- data.frame(c(dec1,dec2), row.names = levels(db_stage2$treat)[2:3])
      
      #pval2df <- data.frame(c(pval2), row.names = levels(db_stage2$treat)[2:3])
      
      stage2_arms <- c(0,0,1)
      simdec_output <- c(0,0,dec3)#c(ifelse(decision[1]=="Reject", 1, 0),                       ifelse(decision[2]=="Reject", 1, 0),                       NA)
      
      decision_intersection = dec123#ifelse(sum(pval2 <= Avalues[1]) > 0, "Reject", "Accept")
      
    }
    
    if(sel_scen == 1){ #continue with low and median doses
      
      N2<-floor(N2/(1+dropout))
      db_stage2 <- sim_dataind(n_arms=3, N = N2, mu_0m =mu_0m[c(1,2,3)], mu_6m =mu_6m[c(1,2,3)], mu_12m=mu_12m[c(1,2,3)], sg=sg, rmonth=rmonth,rr=rr[c(1,2,3)],bound=bound)
      levels(db_stage2$treat) <- levels(db_stage1$treat)[c(1,2,3)]
      recruit_time2 <- max(db_stage2$recruit_time)
      
      db_stage2$diff6_0<-db_stage2$y_6m-db_stage2$y_0m
      db_stage2$diff12_0<-db_stage2$y_12m-db_stage2$y_0m
      
      
      
      pval2 <- c()
      
      
      if (test=="l"){
        
        for(j in 1:2){
          
          sub2 <- subset(db_stage2,(db_stage2$treat==levels(db_stage2$treat)[1])+(db_stage2$treat==levels(db_stage2$treat)[j+1])==1)
          mod2 <- lm(diff12_0~treat+y_0m, sub2) 
          res2 <- summary(mod2)
          pval2[j] <- pt(coef(res2)[2,3], mod2$df, lower.tail = side)
        }
      }
      
      
      if (test=="t"){
        p12low2 <- t.test(db_stage2$diff12_0[db_stage2$treat!="Medium"]~db_stage2$treat[db_stage2$treat!="Medium"],alternative="greater")$p.value
        p12med2<-t.test(db_stage2$diff12_0[db_stage2$treat!="Low"]~db_stage2$treat[db_stage2$treat!="Low"],alternative="greater")$p.value
        pval2<-cbind(p12low2,p12med2)
      }
      
      if (test=="m"){#zweiseitig!!!!
        db_stage2$patID<-c(1:N2)
        db_stage2long<-reshape(data=db_stage2, direction = "long",v.names="diffmf", idvar = "patID", times=c(6,12),varying=c("diff6_0","diff12_0"),timevar="month")
        db_stage2long$month<-factor(db_stage2long$month,c(6,12))
        
        md0 <- lmer(diffmf ~ y_0m + treat+ month+(1 | patID), data = db_stage2long)
        pval2<-summary(md0)$coefficients[3:4,5]
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
      }
      
      if (test=="f"){
        low2<-table(db_stage2$y_12m[db_stage2$treat!="Medium"]<=bound,db_stage2$treat[db_stage2$treat!="Medium"])[,c(1,2)]
        med2<-table(db_stage2$y_12m[db_stage2$treat!="Low"]<=bound,db_stage2$treat[db_stage2$treat!="Low"])[,c(1,3)]
        pvallow2<-try(fisher.test(low2,alternative="greater")$p.value,silent=TRUE)
        if (is.numeric(pvallow2)==FALSE) {pvallow2<-1}
        pvalmed2<-try(fisher.test(med2,alternative="greater")$p.value,silent=TRUE)
        if (is.numeric(pvalmed2)==FALSE) {pvalmed2<-1}
        pval2<-c(pvallow2,pvalmed2)
        #pval2<-c(fisher.test(low2,alternative="greater")$p.value,fisher.test(med2,alternative="greater")$p.value)
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
      #
      
      
      
      #decision_stage2 <- data.frame(c(dec1,dec2), row.names = levels(db_stage2$treat)[2:3])
      
      #pval2df <- data.frame(c(pval2), row.names = levels(db_stage2$treat)[2:3])
      
      stage2_arms <- c(1,1,0)
      simdec_output <- c(dec1,dec2,NA)#c(ifelse(decision[1]=="Reject", 1, 0),                       ifelse(decision[2]=="Reject", 1, 0),                       NA)
      
      decision_intersection = dec123#ifelse(sum(pval2 <= Avalues[1]) > 0, "Reject", "Accept")
      
    }
  }
  
  if(sc == 3){
    
    N2<-floor(N2/(1+dropout))
    db_stage2 <- sim_dataind(n_arms=3, N=N2, mu_0m=mu_0m[c(1,3,4)], mu_6m=mu_6m[c(1,3,4)], mu_12m=mu_12m[c(1,3,4)], sg=sg, rmonth=rmonth,rr=rr[c(1,3,4)],bound=bound)
    
    levels(db_stage2$treat) <- c(levels(db_stage1$treat)[c(1,3)],"High")
    recruit_time2 <- max(db_stage2$recruit_time)
    
    db_stage2$diff6_0<-db_stage2$y_6m-db_stage2$y_0m
    db_stage2$diff12_0<-db_stage2$y_12m-db_stage2$y_0m
    
    
    pval2 <- c()
    
    #for(j in 1:2)
    #  {
    #sub2 <- subset(db_stage2,(db_stage2$treat==levels(db_stage2$treat)[1])+ (db_stage2$treat==levels(db_stage2$treat)[j+1])==1)
    #mod2 <- lm(diff12_0 ~ treat, sub2) #are we using this model or should we use individual models?
    #res2 <- summary(mod2)
    #pval2[j] <- pt(coef(res2)[2,3], mod2$df, lower.tail = side)
    #  }
    
    
    #
    if (test=="l"){
      
      for(j in 1:2){
        
        sub2 <- subset(db_stage2,(db_stage2$treat==levels(db_stage2$treat)[1])+(db_stage2$treat==levels(db_stage2$treat)[j+1])==1)
        mod2 <- lm(diff12_0~treat+y_0m, sub2) 
        res2 <- summary(mod2)
        pval2[j] <- pt(coef(res2)[2,3], mod2$df, lower.tail = side)
      }
    }
    
    
    if (test=="t"){
      p12med2 <- t.test(db_stage2$diff12_0[db_stage2$treat!="High"]~db_stage2$treat[db_stage2$treat!="High"],alternative="greater")$p.value
      p12hi2<-t.test(db_stage2$diff12_0[db_stage2$treat!="Medium"]~db_stage2$treat[db_stage2$treat!="Medium"],alternative="greater")$p.value
      pval2<-cbind(p12med2,p12hi2)
    }
    
    if (test=="m"){
      db_stage2$patID<-c(1:N2)
      db_stage2long<-reshape(data=db_stage2, direction = "long",v.names="diffmf", idvar = "patID", times=c(6,12),varying=c("diff6_0","diff12_0"),timevar="month")
      db_stage2long$month<-factor(db_stage2long$month,c(6,12))
      
      md0 <- lmer(diffmf ~ y_0m + treat+ month+(1 | patID), data = db_stage2long)
      pval2<-summary(md0)$coefficients[3:4,5]
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
    }
    
    if (test=="f"){
      med2<-table(db_stage2$y_12m[db_stage2$treat!="High"]<=bound,db_stage2$treat[db_stage2$treat!="High"])[,c(1,2)]
      high2<-table(db_stage2$y_12m[db_stage2$treat!="Medium"]<=bound,db_stage2$treat[db_stage2$treat!="Medium"])[,c(1,3)]
      pvalmed2<-try(fisher.test(med2,alternative="greater")$p.value,silent=TRUE)
      if (is.numeric(pvalmed2)==FALSE) {pvalmed2<-1}
      pvalhi2<-try(fisher.test(high2,alternative="greater")$p.value,silent=TRUE)
      if (is.numeric(pvalhi2)==FALSE) {pvalhi2<-1}
      pval2<-c(pvalmed2,pvalhi2)#pval2<-c(fisher.test(med2,alternative="greater")$p.value,fisher.test(high2,alternative="greater")$p.value)
      
      
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
    
    #H1
    dec2<-min(pval2[1]<preplan@BJ[2],dec123,dec23,dec12)
    dec3<-min(pval2[2]<preplan@BJ[1],dec123,dec23,dec13)
    #
    
    
    
    #decision_stage2 <- data.frame(c(dec1,dec2), row.names = levels(db_stage2$treat)[2:3])
    
    #pval2df <- data.frame(c(pval2), row.names = levels(db_stage2$treat)[2:3])
    
    stage2_arms <- c(0,1,1)
    simdec_output <- c(0,dec2,dec3)#c(ifelse(decision[1]=="Reject", 1, 0),                       ifelse(decision[2]=="Reject", 1, 0),                       NA)
    
    decision_intersection = dec123#ifelse(sum(pval2 <= Avalues[1]) > 0, "Reject", "Accept")
    
    
    }
  
  
  
  if(sc == 0){ # in that case start dose 3
    
    N2<-floor(N2/(1+dropout))
    db_stage2 <- sim_dataind(n_arms=2, N=N2, mu_0m=mu_0m[c(1,4)], mu_6m=mu_6m[c(1,4)], mu_12m=mu_12m[c(1,4)], sg=sg, rmonth=rmonth,rr=rr[c(1,4)],bound=bound)
    
    levels(db_stage2$treat) <- c(levels(db_stage1$treat)[c(1)],"High")
    recruit_time2 <- max(db_stage2$recruit_time)
    
    db_stage2$diff6_0<-db_stage2$y_6m-db_stage2$y_0m
    db_stage2$diff12_0<-db_stage2$y_12m-db_stage2$y_0m
    
    #mod2 <- lm(diff12_0 ~ treat, db_stage2) #are we using this model or should we use individual models?
    #res2 <- summary(mod2)
    #pval2 <- pt(coef(res2)[2,3], mod2$df, lower.tail = side)
    
    #
    if (test=="l"){
      
      mod2 <- lm(diff12_0~treat+y_0m, db_stage2) 
      res2 <- summary(mod2)
      pval2 <- pt(coef(res2)[2,3], mod2$df, lower.tail = side)
    }
    
    
    
    if (test=="t"){
      pval2 <- t.test(db_stage2$diff12_0~db_stage2$treat,alternative="greater")$p.value
    }
    
    if (test=="m"){
      db_stage2$patID<-c(1:N2)
      db_stage2long<-reshape(data=db_stage2, direction = "long",v.names="diffmf", idvar = "patID", times=c(6,12),varying=c("diff6_0","diff12_0"),timevar="month")
      db_stage2long$month<-factor(db_stage2long$month,c(6,12))
      
      md0 <- lmer(diffmf ~ y_0m + treat+ month+(1 | patID), data = db_stage2long)
      pval2<-summary(md0)$coefficients[3,5]
    }
    
    if (test=="w"){
      pval2 <- wilcox.test(db_stage2$diff12_0~db_stage2$treat,alternative="greater")$p.value
    }
    
    if (test=="w1"){
      pval2 <- wilcox.test(db_stage2$y_12m~db_stage2$treat,alternative="greater")$p.value
    }
    
    if (test=="f"){
      high2<-table(db_stage2$y_12m<=bound,db_stage2$treat)
      pval2<-try(fisher.test(high2,alternative="greater")$p.value,silent=TRUE)
      if (is.numeric(pval2)==FALSE) {pval2<-1}
    }
    
    
    dec123<-max(pval2[1]<preplan@BJ[7],preplan@BJ[7]>1)
    #H12
    dec13<-max(pval2[1]<preplan@BJ[5])
    dec23<-max(pval2[1]<preplan@BJ[3])
    #H1
    dec3<-min(pval2[1]<preplan@BJ[1],dec123,dec13,dec23)
    #
    
    
    
    #decision_stage2 <- data.frame(c(dec1,dec2), row.names = levels(db_stage2$treat)[2:3])
    
    #pval2df <- data.frame(c(pval2), row.names = levels(db_stage2$treat)[2:3])
    
    stage2_arms <- c(0,0,1)
    simdec_output <- c(0,0,dec3)#c(ifelse(decision[1]=="Reject", 1, 0),                       ifelse(decision[2]=="Reject", 1, 0),                       NA)
    
    decision_intersection = dec123#ifelse(sum(pval2 <= Avalues[1]) > 0, "Reject", "Accept")
  }
  
  
  ################################################
  # Multi-arm trial 1
  ####
  #Part 1 for low and median dose
  db_stage_ma1a <- sim_dataind(n_arms = n_arms-1, N = floor((N/5*3)/(1+dropout)), mu_0m = mu_0m[1:n_arms-1], mu_6m = mu_6m[1:n_arms-1], mu_12m = mu_12m[1:n_arms-1], sg = sg, rmonth = rmonth,rr=rr[1:n_arms-1],bound=bound)
  #recruit_time_ma1a <- max(db_stage_ma1$recruit_timea)
  
  db_stage_ma1a$treat <- relevel(db_stage_ma1a$treat, ref = "Placebo")
  
  db_stage_ma1a$diff6_0<-db_stage_ma1a$y_6m-db_stage_ma1a$y_0m
  db_stage_ma1a$diff12_0<-db_stage_ma1a$y_12m-db_stage_ma1a$y_0m
  
  
  #Part 2 for high dose
  db_stage_ma1b <- sim_dataind(n_arms = 2, N = floor((N/5*2)/(1+dropout)), mu_0m = mu_0m[c(1,n_arms)], mu_6m = mu_6m[c(1,n_arms)], mu_12m = mu_12m[c(1,n_arms)], sg = sg, rmonth = rmonth,rr=rr[c(1,n_arms)],bound=bound)
  recruit_time_ma1b <- max(db_stage_ma1b$recruit_time)
  
  db_stage_ma1b$treat <- relevel(db_stage_ma1b$treat, ref = "Placebo")
  
  db_stage_ma1b$diff6_0<-db_stage_ma1b$y_6m-db_stage_ma1b$y_0m
  db_stage_ma1b$diff12_0<-db_stage_ma1b$y_12m-db_stage_ma1b$y_0m
  
  
  
  if (test=="l")
  {
  mod_ma1a <- aov(diff12_0 ~ treat+y_0m, db_stage_ma1a)
  model_dunnett_ma1a = summary(glht(model = mod_ma1a, linfct=mcp(treat="Dunnett"), alternative = "less"))
  pval_dunnett_ma1a = model_dunnett_ma1a$test$pvalues
  
  #part2
  mod_ma1b <- aov(diff12_0 ~ treat+y_0m, db_stage_ma1b)
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
    #decision_ma1<-c((pvalma1a<=(alpha/3*2))*1,(pma1hi<=(alpha/3))*1)
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
  
  
  if (test=="f"){
    tab<-table(db_stage_ma1a$y_12m<=bound,db_stage_ma1a$treat)
    pvalma1a1<-try(fisher.test(tab[,c(1,2)],alternative="greater")$p.value,silent=TRUE)
    pvalma1a2<-try(fisher.test(tab[,c(1,3)],alternative="greater")$p.value,silent=TRUE)
    if (is.numeric(pvalma1a1)==FALSE) {pvalma1a1<-1}
    if (is.numeric(pvalma1a2)==FALSE) {pvalma1a2<-1}
    
    
    pvalma1a<-c(pvalma1a1,pvalma1a2)
    
    tabH<-table(db_stage_ma1b$y_12m<=bound,db_stage_ma1b$treat)
    pma1hi<-try(fisher.test(tabH[,c(1,2)],alternative="greater")$p.value,silent=TRUE)
    if (is.numeric(pma1hi)==FALSE) {pma1hi<-1}
    
    
    decision_ma1<-c((p.adjust(pvalma1a,"holm")<=(alpha/3*2))*1,(p.adjust(pma1hi,"holm")<=(alpha/3))*1)
    decision_ma1a<-c((p.adjust(pvalma1a,"holm")<=(alpha/3*2))*1,(1-max(decision_ma1[1:2]))*((p.adjust(pma1hi,"holm")<=(alpha/3))*1))
  }
  
  
  
  ################################################
  # Multi-arm trial 2
  db_stage_ma2 <- sim_dataind(n_arms = n_arms, N = floor(N/(1+dropout)), mu_0m = mu_0m, mu_6m = mu_6m, mu_12m = mu_12m, sg = sg, rmonth = rmonth,rr=rr,bound=bound)
  recruit_time_ma2 <- max(db_stage_ma2$recruit_time)
  
  db_stage_ma2$treat <- relevel(db_stage_ma2$treat, ref = "Placebo")
  
  db_stage_ma2$diff6_0<-db_stage_ma2$y_6m-db_stage_ma2$y_0m
  db_stage_ma2$diff12_0<-db_stage_ma2$y_12m-db_stage_ma2$y_0m
  
  if (test=="l")
  {
  mod_ma2 <- aov(diff12_0 ~ treat+y_0m, db_stage_ma2)
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
  
  if (test=="f"){
    tab2<-table(db_stage_ma2$y_12m<=bound,db_stage_ma2$treat)
    
    pma1loa<-try(fisher.test(tab2[,c(1,2)],alternative="greater")$p.value,silent=TRUE)
    pma1mea<-try(fisher.test(tab2[,c(1,3)],alternative="greater")$p.value,silent=TRUE)
    pma1hia<-try(fisher.test(tab2[,c(1,4)],alternative="greater")$p.value,silent=TRUE)
    
    if (is.numeric(pma1loa)==FALSE) {pma1loa<-1}
    if (is.numeric(pma1mea)==FALSE) {pma1mea<-1}
    if (is.numeric(pma1hia)==FALSE) {pma1hia<-1}
    
    decision_ma2<-(p.adjust(c(pma1loa,pma1mea,pma1hia),"holm")<alpha)*1
  }
  
  
  
  #######################################
  if(sim_out==F){
    res_intersection=ifelse(decision_intersection == "Reject", 1,0)
    return(list(pval.surr=pval.surr,
                pvalue_stage1=pval1,
                pvalue_stage2=pval2,
                sc=sc,
                decision_stage1=decision_stage1,
                decision_stage2=decision_stage2,
                critical_values=preplan,
                recruit_time1=recruit_time1,
                recruit_time2=recruit_time2,
                stage2_arms=stage2_arms,
                selected_dose=sel,
                simdec_output=simdec_output,
                res_intersection=res_intersection
    ))
  }else{
    res_intersection=ifelse(decision_intersection == "Reject", 1,0)
    return(list(stage2_arms=stage2_arms,
                selected_dose=sel,
                simdec_output=simdec_output,
                res_intersection=res_intersection,
                decision_ma1=decision_ma1,
                decision_ma2=decision_ma2,
                decision_ma1a=decision_ma1a,
                recruit_time1=recruit_time1,
                recruit_time2=recruit_time2))
  }
  #par(mfrow=c(2,3))
  #boxplot(db_stage1$y_0m~db_stage1$treat)
  #boxplot(db_stage1$diff12_0~db_stage1$treat)
  #boxplot(db_stage1$y_12m~db_stage1$treat)
  #boxplot(db_stage2$y_0m~db_stage2$treat)
  #boxplot(db_stage2$diff12_0~db_stage2$treat)
  #boxplot(db_stage2$y_12m~db_stage2$treat)
}

