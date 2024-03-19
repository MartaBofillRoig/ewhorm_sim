library(dplyr) # alternatively tidyverse
library(tidyr) # alternatively tidyverse
library(future)
library(furrr)
library(mvtnorm)



#################################################
#
#   3 arms
#
#



sim_trial_Trichuris <- function(n_arms , N , mu,sigma,reductrate, rmonth, alpha , side,test,dropout,rr,rho)
{
  N<-floor(N/(1+dropout))
  val<-get_mu_sigma(mu,sigma, reductrate_6=reductrate,reductrate_12 = reductrate, rho=rho) # rho is the correlation coefficient
  db_stage1 <- sim_dataind(n_arms = n_arms, N = N, mu_0m = val[[2]][1:n_arms], mu_6m = val[[3]][1:n_arms], mu_12m = val[[4]][1:n_arms], sg = val[[5]], rmonth = rmonth,rr=rr[1:n_arms])
  recruit_time1 <- max(db_stage1$recruit_time)
  
  
  
  db_stage1$treat <- relevel(db_stage1$treat, ref = "Placebo")
  
  db_stage1$diff6_0<-db_stage1$y_6m-db_stage1$y_0m
  db_stage1$diff12_0<-db_stage1$y_12m-db_stage1$y_0m
  
  
  
  pval1 <- c()  #12months p-value of first stage
  
  if (test=="t"){
    p12low <- t.test(db_stage1$diff6_0[db_stage1$treat %in% c("Low","Placebo")]~db_stage1$treat[db_stage1$treat %in% c("Low","Placebo")],alternative="greater")$p.value
    p12med<-t.test(db_stage1$diff6_0[db_stage1$treat %in% c("Medium","Placebo")]~db_stage1$treat[db_stage1$treat %in% c("Medium","Placebo")],alternative="greater")$p.value
    pval1<-cbind(p12low,p12med)
    #pval
  }
  
  if (test=="l"){  #two individual models, due to robustness in case of different variances
    
    for(j in 1:2)#(n_arms-2))
    {
      
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
    p12low <- wilcox.test(db_stage1$diff6_0[db_stage1$treat %in% c("Low","Placebo")]~db_stage1$treat[db_stage1$treat %in% c("Low","Placebo")],alternative="greater")$p.value
    p12med<-wilcox.test(db_stage1$diff6_0[db_stage1$treat %in% c("Medium","Placebo")]~db_stage1$treat[db_stage1$treat %in% c("Medium","Placebo")],alternative="greater")$p.value
    pval1<-cbind(p12low,p12med)
  }
  
  if (test=="w1"){
    p12low <- wilcox.test(db_stage1$y_6m[db_stage1$treat %in% c("Low","Placebo")]~db_stage1$treat[db_stage1$treat %in% c("Low","Placebo")],alternative="greater")$p.value
    p12med<-wilcox.test(db_stage1$y_6m[db_stage1$treat %in% c("Medium","Placebo")]~db_stage1$treat[db_stage1$treat %in% c("Medium","Placebo")],alternative="greater")$p.value
    pval1<-cbind(p12low,p12med)
    #pval
  }
  
  p.adjust(pval1,method="holm")
}








do_Trichuris = function(n_trials,n_arms , N , mu,sigma,reductrate0,reductrate1,reductrate2, rmonth, alpha , side1,test1,dropout,rr0,rr1,rr2,rho)
{
  
  #Set the number of trials to run and other parameters for future plan
  #n_trials <- 10000
  
  reductrate<-c(reductrate0,reductrate1,reductrate2)
  rr<-c(rr0,rr1,rr2)
  
  if (side1==1) side<-T
  if (side1==0) side<-F
  if (test1==0) test<-"l"
  if (test1==1) test<-"m"
  if (test1==2) test<-"t"
  if (test1==3) test<-"w"
  if (test1==4) test<-"w1"
  
  n_cores <- availableCores()-1 
  plan(multisession, workers = n_cores)
  
  
  # Run the simulations in parallel using future_map
  results_list <- future_map(1:n_trials, function(i)
    sim_trial_Trichuris (n_arms=n_arms , N=N , mu=mu,sigma=sigma,reductrate=reductrate, rmonth=rmonth, 
                         alpha=alpha , side=side,test=test,dropout=dropout,rr=rr,rho),.options=furrr_options(seed = TRUE))
  
  pow<-c(apply(matrix(unlist(results_list),ncol=2,byrow=T)<=alpha,2,sum)/n_trials,sum((apply(matrix(unlist(results_list),ncol=2,byrow=T)<=alpha,1,sum)>0))/n_trials)
  
  
  return(pow)
}

#do_Trichuris(n_trials=1000,n_arms = 3, N=150 , mu=50,sigma=50,reductrate0=0,reductrate1=0,reductrate2=0, rmonth=3, alpha = 0.025, side1=T,test1="w1",dropout=0,rr0=0,rr1=0,rr2=0)
#do_Trichuris(n_trials=1000,n_arms = 3, N=150 , mu=50,sigma=50,reductrate=c(r0_hi,r1_hi,r2_hi,r3_hi), rmonth=3, alpha = 0.025, side=T,test="w1",dropout=0,rr=c(0,0,0,0))
#do_Trichuris(n_trials=1000,n_arms = 3, N=150 , mu=50,sigma=50,reductrate=c(r0_hi,r1_hi,r2_hi,r3_hi), rmonth=3, alpha = 0.025, side=T,test="w1",dropout=0,rr=c(0,0,0,0))
#do_Trichuris(n_trials=1000,n_arms = 3, N=150 , mu=50,sigma=50,reductrate=c(r0_hi,r1_hi,r2_hi,r3_hi), rmonth=3, alpha = 0.025, side=T,test="w1",dropout=0,rr=c(0,0,0,0))





mu<-1000
sigma<-3000#2000
rmonth=1





###################################################################################################################################################
#
#         Results 3 arms (2 doses against Placebo) (13.02.2024)
#
###################################################################################################################################################



##################################
#
# Simulation 2
#
##################################





###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
#
#         Results 3 arms (2 doses against Placebo) (19.02.2024)   ------    N=120
#
###################################################################################################################################################



##################################
#
# Simulation 2
#
##################################

#Trich




rho<-.5

n_arms = 3
N = 120

alpha=.025
side1=1
dropout=.1
#rr0=0
#rr1=.50
#rr2=.5
n_trials<-10000#0
rmonth=3
#no_effect

#simulation.global.null<-c(mu_raw_0, sd_raw_0, r0_6,r1_6,r2_6,r3_6, r0_12,r1_12,r2_12,r3_12,  rho ,   #data.matrix(expand.grid
#            n_trials,n_arms,N1 , N, rmonth, alpha1 , alpha,  sim_out1,sel_scen, side1, test1)#)

test1=0
rr0=0
rr1=0
rr2=0


reduct0<-0
reduct1<-0
reduct2<-c(0,.1,.2,.3,.4,.5,.6,.7,.8)

Trichuris_no_effect<-mapply(do_Trichuris,n_trials,n_arms, N , mu,sigma,reduct0,reduct1,reduct2, rmonth, alpha, side1,test1,dropout,rr0,rr1,rr2,rho)

reduct0<-0
reduct1<-.2
reduct2<-c(0,.1,.2,.3,.4,.5,.6,.7,.8)
Trichuris_low_effect<-mapply(do_Trichuris,n_trials,n_arms, N , mu,sigma,reduct0,reduct1,reduct2, rmonth, alpha, side1,test1,dropout,rr0,rr1,rr2,rho)

reduct0<-0
reduct1<-.4
reduct2<-c(0,.1,.2,.3,.4,.5,.6,.7,.8)
Trichuris_med_effect<-mapply(do_Trichuris,n_trials,n_arms, N , mu,sigma,reduct0,reduct1,reduct2, rmonth, alpha, side1,test1,dropout,rr0,rr1,rr2,rho)

reduct0<-0
reduct1<-.6
reduct2<-c(0,.1,.2,.3,.4,.5,.6,.7,.8)
Trichuris_high_effect<-mapply(do_Trichuris,n_trials,n_arms, N , mu,sigma,reduct0,reduct1,reduct2, rmonth, alpha, side1,test1,dropout,rr0,rr1,rr2,rho)

##############################################################################################

#
#
###
#
#
#       Total responder rho=0.5
test1=0


#reduct0<-0
#reduct1<-0
#reduct2<-0
n_trials<-10000#0
#Trichuris_no_effect<-mapply(do_Trichuris,n_trials,n_arms, N , mu,sigma,reduct0,reduct1,reduct2,reduct3, rmonth, alpha, side1,test1,dropout,rr0,rr1,rr2,rr3)

reduct0<-0
reduct1<-0
reduct2<-c(0,.1,.2,.3,.4,.5,.6,.7,.8)
rr0<-0.1
rr1<-0.1
rr2<-c(0.1,.1,.2-0.1,.3-0.1,.4-0.1,.5-0.1,.6-0.1,.7-0.1,.8-0.1)

Trichuris_no_effect0<-mapply(do_Trichuris,n_trials,n_arms, N , mu,sigma,reduct0,reduct1,reduct2, rmonth, alpha, side1,test1,dropout,rr0,rr1,rr2,rho)

reduct0<-0
reduct1<-.2
reduct2<-c(0,.1,.2,.3,.4,.5,.6,.7,.8)
rr0<-0.1
rr1<-0.1
rr2<-c(0.1,.1,.2-0.1,.3-0.1,.4-0.1,.5-0.1,.6-0.1,.7-0.1,.8-0.1)


Trichuris_low_effect0<-mapply(do_Trichuris,n_trials,n_arms, N , mu,sigma,reduct0,reduct1,reduct2, rmonth, alpha, side1,test1,dropout,rr0,rr1,rr2,rho)

reduct0<-0
reduct1<-.4
reduct2<-c(0,.1,.2,.3,.4,.5,.6,.7,.8)
rr0<-0.1
rr1<-0.3
rr2<-c(0.1,.1,.2-0.1,.3-0.1,.4-0.1,.5-0.1,.6-0.1,.7-0.1,.8-0.1)


Trichuris_med_effect0<-mapply(do_Trichuris,n_trials,n_arms, N , mu,sigma,reduct0,reduct1,reduct2, rmonth, alpha, side1,test1,dropout,rr0,rr1,rr2,rho)

reduct0<-0
reduct1<-.6
reduct2<-c(0,.1,.2,.3,.4,.5,.6,.7,.8)

rr0<-0.1
rr1<-0.5
rr2<-c(0.1,.1,.2-0.1,.3-0.1,.4-0.1,.5-0.1,.6-0.1,.7-0.1,.8-0.1)

Trichuris_high_effect0<-mapply(do_Trichuris,n_trials,n_arms, N , mu,sigma,reduct0,reduct1,reduct2, rmonth, alpha, side1,test1,dropout,rr0,rr1,rr2,rho)

######################################



######################################################


test1=4


#reduct0<-0
#reduct1<-0
#reduct2<-0
n_trials<-10000#0
#Trichuris_no_effect<-mapply(do_Trichuris,n_trials,n_arms, N , mu,sigma,reduct0,reduct1,reduct2,reduct3, rmonth, alpha, side1,test1,dropout,rr0,rr1,rr2,rr3)

reduct0<-0
reduct1<-0
reduct2<-c(0,.1,.2,.3,.4,.5,.6,.7,.8)
rr0<-0.1
rr1<-0.1
rr2<-c(0.1,.1,.2-0.1,.3-0.1,.4-0.1,.5-0.1,.6-0.1,.7-0.1,.8-0.1)

Trichuris_no_effect4<-mapply(do_Trichuris,n_trials,n_arms, N , mu,sigma,reduct0,reduct1,reduct2, rmonth, alpha, side1,test1,dropout,rr0,rr1,rr2,rho)

reduct0<-0
reduct1<-.2
reduct2<-c(0,.1,.2,.3,.4,.5,.6,.7,.8)
rr0<-0.1
rr1<-0.1
rr2<-c(0.1,.1,.2-0.1,.3-0.1,.4-0.1,.5-0.1,.6-0.1,.7-0.1,.8-0.1)


Trichuris_low_effect4<-mapply(do_Trichuris,n_trials,n_arms, N , mu,sigma,reduct0,reduct1,reduct2, rmonth, alpha, side1,test1,dropout,rr0,rr1,rr2,rho)

reduct0<-0
reduct1<-.4
reduct2<-c(0,.1,.2,.3,.4,.5,.6,.7,.8)
rr0<-0.1
rr1<-0.3
rr2<-c(0.1,.1,.2-0.1,.3-0.1,.4-0.1,.5-0.1,.6-0.1,.7-0.1,.8-0.1)


Trichuris_med_effect4<-mapply(do_Trichuris,n_trials,n_arms, N , mu,sigma,reduct0,reduct1,reduct2, rmonth, alpha, side1,test1,dropout,rr0,rr1,rr2,rho)

reduct0<-0
reduct1<-.6
reduct2<-c(0,.1,.2,.3,.4,.5,.6,.7,.8)

rr0<-0.1
rr1<-0.5
rr2<-c(0.1,.1,.2-0.1,.3-0.1,.4-0.1,.5-0.1,.6-0.1,.7-0.1,.8-0.1)

Trichuris_high_effect4<-mapply(do_Trichuris,n_trials,n_arms, N , mu,sigma,reduct0,reduct1,reduct2, rmonth, alpha, side1,test1,dropout,rr0,rr1,rr2,rho)









##############################################################################################

#
#
###
#
#
#       Total responder rho=0.5
test1=0
rho<-0.25
mu<-500
sigma<-4000#2000


#reduct0<-0
#reduct1<-0
#reduct2<-0
n_trials<-10000#0
#Trichuris_no_effect<-mapply(do_Trichuris,n_trials,n_arms, N , mu,sigma,reduct0,reduct1,reduct2,reduct3, rmonth, alpha, side1,test1,dropout,rr0,rr1,rr2,rr3)

reduct0<-0
reduct1<-0
reduct2<-c(0,.1,.2,.3,.4,.5,.6,.7,.8)
rr0<-0.1
rr1<-0.1
rr2<-c(0.1,.1,.2-0.1,.3-0.1,.4-0.1,.5-0.1,.6-0.1,.7-0.1,.8-0.1)

Trichuris_no_effect0r<-mapply(do_Trichuris,n_trials,n_arms, N , mu,sigma,reduct0,reduct1,reduct2, rmonth, alpha, side1,test1,dropout,rr0,rr1,rr2,rho)

reduct0<-0
reduct1<-.2
reduct2<-c(0,.1,.2,.3,.4,.5,.6,.7,.8)
rr0<-0.1
rr1<-0.1
rr2<-c(0.1,.1,.2-0.1,.3-0.1,.4-0.1,.5-0.1,.6-0.1,.7-0.1,.8-0.1)


Trichuris_low_effect0r<-mapply(do_Trichuris,n_trials,n_arms, N , mu,sigma,reduct0,reduct1,reduct2, rmonth, alpha, side1,test1,dropout,rr0,rr1,rr2,rho)

reduct0<-0
reduct1<-.4
reduct2<-c(0,.1,.2,.3,.4,.5,.6,.7,.8)
rr0<-0.1
rr1<-0.3
rr2<-c(0.1,.1,.2-0.1,.3-0.1,.4-0.1,.5-0.1,.6-0.1,.7-0.1,.8-0.1)


Trichuris_med_effect0r<-mapply(do_Trichuris,n_trials,n_arms, N , mu,sigma,reduct0,reduct1,reduct2, rmonth, alpha, side1,test1,dropout,rr0,rr1,rr2,rho)

reduct0<-0
reduct1<-.6
reduct2<-c(0,.1,.2,.3,.4,.5,.6,.7,.8)

rr0<-0.1
rr1<-0.5
rr2<-c(0.1,.1,.2-0.1,.3-0.1,.4-0.1,.5-0.1,.6-0.1,.7-0.1,.8-0.1)

Trichuris_high_effect0r<-mapply(do_Trichuris,n_trials,n_arms, N , mu,sigma,reduct0,reduct1,reduct2, rmonth, alpha, side1,test1,dropout,rr0,rr1,rr2,rho)

######################################



######################################################


test1=4


#reduct0<-0
#reduct1<-0
#reduct2<-0
n_trials<-10000#0
#Trichuris_no_effect<-mapply(do_Trichuris,n_trials,n_arms, N , mu,sigma,reduct0,reduct1,reduct2,reduct3, rmonth, alpha, side1,test1,dropout,rr0,rr1,rr2,rr3)

reduct0<-0
reduct1<-0
reduct2<-c(0,.1,.2,.3,.4,.5,.6,.7,.8)
rr0<-0.1
rr1<-0.1
rr2<-c(0.1,.1,.2-0.1,.3-0.1,.4-0.1,.5-0.1,.6-0.1,.7-0.1,.8-0.1)

Trichuris_no_effect4r<-mapply(do_Trichuris,n_trials,n_arms, N , mu,sigma,reduct0,reduct1,reduct2, rmonth, alpha, side1,test1,dropout,rr0,rr1,rr2,rho)

reduct0<-0
reduct1<-.2
reduct2<-c(0,.1,.2,.3,.4,.5,.6,.7,.8)
rr0<-0.1
rr1<-0.1
rr2<-c(0.1,.1,.2-0.1,.3-0.1,.4-0.1,.5-0.1,.6-0.1,.7-0.1,.8-0.1)


Trichuris_low_effect4r<-mapply(do_Trichuris,n_trials,n_arms, N , mu,sigma,reduct0,reduct1,reduct2, rmonth, alpha, side1,test1,dropout,rr0,rr1,rr2,rho)

reduct0<-0
reduct1<-.4
reduct2<-c(0,.1,.2,.3,.4,.5,.6,.7,.8)
rr0<-0.1
rr1<-0.3
rr2<-c(0.1,.1,.2-0.1,.3-0.1,.4-0.1,.5-0.1,.6-0.1,.7-0.1,.8-0.1)


Trichuris_med_effect4r<-mapply(do_Trichuris,n_trials,n_arms, N , mu,sigma,reduct0,reduct1,reduct2, rmonth, alpha, side1,test1,dropout,rr0,rr1,rr2,rho)

reduct0<-0
reduct1<-.6
reduct2<-c(0,.1,.2,.3,.4,.5,.6,.7,.8)

rr0<-0.1
rr1<-0.5
rr2<-c(0.1,.1,.2-0.1,.3-0.1,.4-0.1,.5-0.1,.6-0.1,.7-0.1,.8-0.1)

Trichuris_high_effect4r<-mapply(do_Trichuris,n_trials,n_arms, N , mu,sigma,reduct0,reduct1,reduct2, rmonth, alpha, side1,test1,dropout,rr0,rr1,rr2,rho)






#Power
pdf(file ="Trichuris_PowV2_120_x.pdf", width = 9, height = 9, pointsize = 12)
par (mfrow=c(3,3), mar=c(1, 1, 1, 1) + 0.1, oma=c(3, 3, 1, 1))
plot(1:9,Trichuris_no_effect[1,],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col="black")
legend(x=0.3,y=0.75,legend=c("Reduction rate low dose", "0","0.2","0.4","0.6"),lwd=2,col=c("white","black",2,3,5),lty=1,bty="n",cex=1.3)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
mtext(side=3,"Power low Dose",cex=1.3)#,line=2.2)
lines(1:9,Trichuris_low_effect[1,],type="b",lwd=2,col=2)
lines(1:9,Trichuris_med_effect[1,],type="b",lwd=2,col=3,lty=1)
lines(1:9,Trichuris_high_effect[1,],type="b",lwd=2,col=5,lty=1)
#mtext(side=1,c("Reduction rate medium dose"),cex=1,line=2)#,line=2.2)
#axis(1,c(1:9),seq(0,0.8,0.1),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)

plot(1:9,Trichuris_no_effect[2,],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col="black")
mtext(side=3,"Power medium Dose",cex=1.3)#,line=2.2)
lines(1:9,Trichuris_low_effect[2,],type="b",lwd=2,col=2)
lines(1:9,Trichuris_med_effect[2,],type="b",lwd=2,col=3,lty=1)
lines(1:9,Trichuris_high_effect[2,],type="b",lwd=2,col=5,lty=1)
#mtext(side=1,c("Reduction rate medium dose"),cex=1,line=2)#,line=2.2)
#axis(1,c(1:9),seq(0,0.8,0.1),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)

plot(1:9,Trichuris_no_effect[3,],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col="black")
mtext(side=3,"Disjunctive Power",cex=1.3)#,line=2.2)
lines(1:9,Trichuris_low_effect[3,],type="b",lwd=2,col=2)
lines(1:9,Trichuris_med_effect[3,],type="b",lwd=2,col=3,lty=1)
lines(1:9,Trichuris_high_effect[3,],type="b",lwd=2,col=5,lty=1)
mtext(side=4,"Simulation scenario 1",cex=1,line=1)#,line=2.2)
#mtext(side=1,c("Reduction rate medium dose"),cex=1,line=2)#,line=2.2)
#axis(1,c(1:9),seq(0,0.8,0.1),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)




plot(1:9,Trichuris_no_effect0[1,],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col="black")
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
#mtext(side=3,"Power low Dose",cex=1.3)#,line=2.2)
lines(1:9,Trichuris_low_effect0[1,],type="b",lwd=2,col=2)
lines(1:9,Trichuris_med_effect0[1,],type="b",lwd=2,col=3,lty=1)
lines(1:9,Trichuris_high_effect0[1,],type="b",lwd=2,col=5,lty=1)
#mtext(side=1,c("Reduction rate medium dose"),cex=1,line=2)#,line=2.2)
#axis(1,c(1:9),seq(0,0.8,0.1),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
lines(1:9,Trichuris_no_effect4[1,],type="b",lwd=2,col=1,lty=3)
lines(1:9,Trichuris_low_effect4[1,],type="b",lwd=2,col=2,lty=3)
lines(1:9,Trichuris_med_effect4[1,],type="b",lwd=2,col=3,lty=3)
lines(1:9,Trichuris_high_effect4[1,],type="b",lwd=2,col=5,lty=3)



plot(1:9,Trichuris_no_effect0[2,],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col="black",yaxt="n")
#mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
#mtext(side=3,"Power low Dose",cex=1.3)#,line=2.2)
lines(1:9,Trichuris_low_effect0[2,],type="b",lwd=2,col=2)
lines(1:9,Trichuris_med_effect0[2,],type="b",lwd=2,col=3,lty=1)
lines(1:9,Trichuris_high_effect0[2,],type="b",lwd=2,col=5,lty=1)
#mtext(side=1,c("Reduction rate medium dose"),cex=1,line=2)#,line=2.2)
#axis(1,c(1:9),seq(0,0.8,0.1),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
lines(1:9,Trichuris_no_effect4[2,],type="b",lwd=2,col=1,lty=3)
lines(1:9,Trichuris_low_effect4[2,],type="b",lwd=2,col=2,lty=3)
lines(1:9,Trichuris_med_effect4[2,],type="b",lwd=2,col=3,lty=3)
lines(1:9,Trichuris_high_effect4[2,],type="b",lwd=2,col=5,lty=3)


plot(1:9,Trichuris_no_effect0[3,],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col="black",yaxt="n")
#mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
#mtext(side=3,"Power low Dose",cex=1.3)#,line=2.2)
lines(1:9,Trichuris_low_effect0[3,],type="b",lwd=2,col=2)
lines(1:9,Trichuris_med_effect0[3,],type="b",lwd=2,col=3,lty=1)
lines(1:9,Trichuris_high_effect0[3,],type="b",lwd=2,col=5,lty=1)
#mtext(side=1,c("Reduction rate medium dose"),cex=1,line=2)#,line=2.2)
#axis(1,c(1:9),seq(0,0.8,0.1),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
lines(1:9,Trichuris_no_effect4[3,],type="b",lwd=2,col=1,lty=3)
lines(1:9,Trichuris_low_effect4[3,],type="b",lwd=2,col=2,lty=3)
lines(1:9,Trichuris_med_effect4[3,],type="b",lwd=2,col=3,lty=3)
lines(1:9,Trichuris_high_effect4[3,],type="b",lwd=2,col=5,lty=3)
mtext(side=4,"Simulation scenario 2",cex=1,line=1)#,line=2.2)






plot(1:9,Trichuris_no_effect0r[1,],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col="black")
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
#mtext(side=3,"Power low Dose",cex=1.3)#,line=2.2)
lines(1:9,Trichuris_low_effect0r[1,],type="b",lwd=2,col=2)
lines(1:9,Trichuris_med_effect0r[1,],type="b",lwd=2,col=3,lty=1)
lines(1:9,Trichuris_high_effect0r[1,],type="b",lwd=2,col=5,lty=1)
mtext(side=1,c("Reduction rate medium dose"),cex=1,line=2)#,line=2.2)
axis(1,c(1:9),seq(0,0.8,0.1),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
#lines(1:9,Trichuris_no_effect2[1,],type="b",lwd=2,col=1)
#lines(1:9,Trichuris_low_effect2[1,],type="b",lwd=2,col=2)
#lines(1:9,Trichuris_med_effect2[1,],type="b",lwd=2,col=3)
#lines(1:9,Trichuris_high_effect2[1,],type="b",lwd=2,col=5)
#lines(1:9,Trichuris_no_effect3[1,],type="b",lwd=2,col=1)
#lines(1:9,Trichuris_low_effect3[1,],type="b",lwd=2,col=2)
#lines(1:9,Trichuris_med_effect3[1,],type="b",lwd=2,col=3)
#lines(1:9,Trichuris_high_effect3[1,],type="b",lwd=2,col=5)
lines(1:9,Trichuris_no_effect4r[1,],type="b",lwd=2,col=1,lty=3)
lines(1:9,Trichuris_low_effect4r[1,],type="b",lwd=2,col=2,lty=3)
lines(1:9,Trichuris_med_effect4r[1,],type="b",lwd=2,col=3,lty=3)
lines(1:9,Trichuris_high_effect4r[1,],type="b",lwd=2,col=5,lty=3)



plot(1:9,Trichuris_no_effect0r[2,],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col="black",yaxt="n")
#mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
#mtext(side=3,"Power low Dose",cex=1.3)#,line=2.2)
lines(1:9,Trichuris_low_effect0r[2,],type="b",lwd=2,col=2)
lines(1:9,Trichuris_med_effect0r[2,],type="b",lwd=2,col=3,lty=1)
lines(1:9,Trichuris_high_effect0r[2,],type="b",lwd=2,col=5,lty=1)
mtext(side=1,c("Reduction rate medium dose"),cex=1,line=2)#,line=2.2)
axis(1,c(1:9),seq(0,0.8,0.1),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
lines(1:9,Trichuris_no_effect4r[2,],type="b",lwd=2,col=1,lty=3)
lines(1:9,Trichuris_low_effect4r[2,],type="b",lwd=2,col=2,lty=3)
lines(1:9,Trichuris_med_effect4r[2,],type="b",lwd=2,col=3,lty=3)
lines(1:9,Trichuris_high_effect4r[2,],type="b",lwd=2,col=5,lty=3)


plot(1:9,Trichuris_no_effect0r[3,],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col="black",yaxt="n")
#mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
#mtext(side=3,"Power low Dose",cex=1.3)#,line=2.2)
lines(1:9,Trichuris_low_effect0r[3,],type="b",lwd=2,col=2)
lines(1:9,Trichuris_med_effect0r[3,],type="b",lwd=2,col=3,lty=1)
lines(1:9,Trichuris_high_effect0r[3,],type="b",lwd=2,col=5,lty=1)
mtext(side=1,c("Reduction rate medium dose"),cex=1,line=2)#,line=2.2)
axis(1,c(1:9),seq(0,0.8,0.1),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
lines(1:9,Trichuris_no_effect4r[3,],type="b",lwd=2,col=1,lty=3)
lines(1:9,Trichuris_low_effect4r[3,],type="b",lwd=2,col=2,lty=3)
lines(1:9,Trichuris_med_effect4r[3,],type="b",lwd=2,col=3,lty=3)
lines(1:9,Trichuris_high_effect4r[3,],type="b",lwd=2,col=5,lty=3)
mtext(side=4,"Simulation scenario 3",cex=1,line=1)#,line=2.2)

dev.off()










