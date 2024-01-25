r0_no=0
r1_no=0
r2_no=0
r3_no=0

r0_lo=0
r1_lo=0
r2_lo=0
r3_lo=0.6

r0_me=0
r1_me=0
r2_me=0.5
r3_me=0.6


r0_hi=0
r1_hi=0.6
r2_hi=0.6
r3_hi=0.6





sim_trial_Trichuris <- function(n_arms = 4, N , mu,sigma,reductrate, rmonth, alpha = 0.025, side=T,test="t",dropout=0,rr=c(0,0,0,0))
{
  N1<-floor(N1/(1+dropout))
  val<-get_mu_sigma(mu,sigma, reductrate_6=reductrate,reductrate_12 = reductrate, rho=.5) # rho is the correlation coefficient
  db_stage1 <- sim_dataind(n_arms = n_arms, N = N1, mu_0m = val[[2]][1:n_arms], mu_6m = val[[3]][1:n_arms], mu_12m = val[[4]][1:n_arms], sg = val[[5]], rmonth = rmonth,rr=rr[1:n_arms])
  recruit_time1 <- max(db_stage1$recruit_time)
  
  
  
  db_stage1$treat <- relevel(db_stage1$treat, ref = "Placebo")
  
  db_stage1$diff6_0<-db_stage1$y_6m-db_stage1$y_0m
  db_stage1$diff12_0<-db_stage1$y_12m-db_stage1$y_0m
  
  
  
  pval1 <- c()  #12months p-value of first stage
  
  if (test=="t"){
    p12low <- t.test(db_stage1$diff6_0[db_stage1$treat!="Medium"]~db_stage1$treat[db_stage1$treat!="Medium"],alternative="greater")$p.value
    p12med<-t.test(db_stage1$diff6_0[db_stage1$treat!="Low"]~db_stage1$treat[db_stage1$treat!="Low"],alternative="greater")$p.value
    pval1<-cbind(p12low,p12med)
    #pval
  }
  
  if (test=="l"){  #two individual models, due to robustness in case of different variances
    
    for(j in 1:(n_arms-2)){
      
      sub1 <- subset(db_stage1,(db_stage1$treat==levels(db_stage1$treat)[1])+ (db_stage1$treat==levels(db_stage1$treat)[j+1])==1)
      mod1 <- lm(diff6_0 ~ treat+y_0m, sub1)
      res1 <- summary(mod1)
      pval1[j] <- pt(coef(res1)[2,3], mod1$df, lower.tail = side)
    }
  }
  
  #if (test=="m"){
  #  db_stage1$patID<-c(1:N1)
  #  db_stage1long<-reshape(data=db_stage1, direction = "long",v.names="diffmf", idvar = "patID", times=c(6,12),varying=c("diff6_0","diff12_0"),timevar="month")
  #  db_stage1long$month<-factor(db_stage1long$month,c(6,12))
  #  
  #  md0 <- lmer(diffmf ~ y_0m + treat+ month+(1 | patID), data = db_stage1long)
  #  pval1<-summary(md0)$coefficients[3:4,5]
  #}
  
  if (test=="w"){
    p12low <- wilcox.test(db_stage1$diff6_0[db_stage1$treat!="Medium"]~db_stage1$treat[db_stage1$treat!="Medium"],alternative="greater")$p.value
    p12med<-wilcox.test(db_stage1$diff6_0[db_stage1$treat!="Low"]~db_stage1$treat[db_stage1$treat!="Low"],alternative="greater")$p.value
    pval1<-cbind(p12low,p12med)
    #pval
  }
  
  if (test=="w1"){
    p12low <- wilcox.test(db_stage1$y_6m[db_stage1$treat!="Medium"]~db_stage1$treat[db_stage1$treat!="Medium"],alternative="greater")$p.value
    p12med<-wilcox.test(db_stage1$y_6m[db_stage1$treat!="Low"]~db_stage1$treat[db_stage1$treat!="Low"],alternative="greater")$p.value
    pval1<-cbind(p12low,p12med)
    #pval
  }
  
  p.adjust(pval1,method="holm")
}

sim_trial_Trichuris(n_arms = 3, N=150 , mu=50,sigma=50,reductrate=c(r0_hi,r1_hi,r2_hi,r3_hi), rmonth=3, alpha = 0.025, side=T,test="t",dropout=0,rr=c(0,0,0,0))
n_arms = 3; N=150; mu=50;sigma=1;reductrate=c(0,0,1,1); rmonth=3; alpha = 0.025; side=T;test="t";dropout=0;rr=c(0,0,0,0);
