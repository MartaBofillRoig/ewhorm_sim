##########################################################################
#Example 1

#all doses effective:
mu_sigma<-get_mu_sigma(mu_raw_0=1838,sd_raw_0=2565 , reductrate_6=c(0,.40,.40,.40) , reductrate_12=c(0,0.5,.50,.50), rho=.5 )
mu_sigma

set.seed(24120403)
e1<-sim_trial_pceind_test_example(n_arms = 4,N1 = 120 , N2 = 80, mu_0m = mu_sigma[[2]],   mu_6m = mu_sigma[[3]], 
                      mu_12m = mu_sigma[[4]],  sg = mu_sigma[[5]], #rmonth=1, 
                      alpha1 = .3, alpha = .025, #sim_out=T,
                      sel_scen=1, side=T, test="w1",dropout=.1,rr=c(0.1,.3,.30,.30),bound=0)
#n_arms<-4
db_stage1 <- e1$db_stage1#sim_dataind(n_arms = n_arms-1, N = 120, mu_0m = mu_sigma[[2]][1:n_arms-1], mu_6m = mu_sigma[[3]][1:n_arms-1], 
                         #mu_12m = mu_sigma[[4]][1:n_arms-1], sg = mu_sigma[[5]], rr=c(0.1,.3,.30,.30)[1:n_arms-1],bound=bound)

db_stage2 <- e1$db_stage2



#e1$pval.surr
#[1] 0.002273498 0.004612684

#$pval1
#p12low       p12med
#[1,] 0.00932757 0.0007665152

#$pval2
#p12low2    p12med2
#[1,] 0.0005318758 0.03677223

#pval.surr.1<-e1$pval.surr # c(0.002273498, 0.004612684)
#pval1<-e1$pval1 #c(0.00932757, 0.0007665152)
#pval1.1<-pval1
#pval2<-e1$pval2 #c(0.0005318758, 0.03677223)
#pval2.1<-pval2

round(c(e1$pval.surr,e1$pval1,e1$pval2,e1$concinvn,e1$concCIinvn),3)
e1$simdec_output #test decision 1=reject



alpha<-0.025
z <- qnorm(1-pmax(pval1,1e-15))
N=floor(120*(1-.1))+80#N1+N2
graph_bh <- BonferroniHolm(3)

v = c(N1/N,N1/N,0)
z1 <- c(z,0)

preplan <- doInterim(graph=graph_bh,z1=z1,v=v,alpha=alpha)

#H123
dec123<-max(pval2[1]<preplan@BJ[7]*preplan@Aj[7,1]/(preplan@Aj[7,1]+preplan@Aj[7,2]),
            pval2[2]<preplan@BJ[7]*preplan@Aj[7,2]/(preplan@Aj[7,1]+preplan@Aj[7,2]),preplan@BJ[7]>1)

preplan@BJ[7]*preplan@Aj[7,1]/(preplan@Aj[7,1]+preplan@Aj[7,2])
preplan@BJ[7]*preplan@Aj[7,2]/(preplan@Aj[7,1]+preplan@Aj[7,2])
preplan@BJ[7]

#H12
dec12<-max(pval2[1]<preplan@BJ[6]*preplan@Aj[6,1]/(preplan@Aj[6,1]+preplan@Aj[6,2]),
           pval2[2]<preplan@BJ[6]*preplan@Aj[6,2]/(preplan@Aj[6,1]+preplan@Aj[6,2]),preplan@BJ[6]>1)

preplan@BJ[6]*preplan@Aj[6,1]/(preplan@Aj[6,1]+preplan@Aj[6,2])
preplan@BJ[6]*preplan@Aj[6,2]/(preplan@Aj[6,1]+preplan@Aj[6,2])
preplan@BJ[6]


#H13
dec13<-max(pval2[1]<preplan@BJ[5])
preplan@BJ[5]

#H23
dec23<-max(pval2[2]<preplan@BJ[3])
preplan@BJ[3]

#H1
dec1<-min(pval2[1]<preplan@BJ[4],dec123,dec12,dec13)
preplan@BJ[4]

#H2
dec2<-min(pval2[2]<preplan@BJ[2],dec123,dec12,dec23)
preplan@BJ[2]

stage2_arms <- c(1,1,0)
simdec_output <- c(dec1,dec2,NA)#c(ifelse(decision[1]=="Reject", 1, 0),                       ifelse(decision[2]=="Reject", 1, 0),                       NA)

decision_intersection = dec123#ifelse(sum(pval2 <= Avalues[1]) > 0, "Reject", "Accept")
#




###################################################################
#Example 2

#Trend:
mu_sigma<-get_mu_sigma(mu_raw_0=1838,sd_raw_0=2565 , reductrate_6=c(0,0,.30,.50) , reductrate_12=c(0,0,.40,.60), rho=.5 )
mu_sigma

#set.seed(040412)
set.seed(04041278)
#set.seed(24120403)
e2<-sim_trial_pceind_test_example(n_arms = 4,N1 = 120 , N2 = 80, mu_0m = mu_sigma[[2]],   mu_6m = mu_sigma[[3]], 
                      mu_12m = mu_sigma[[4]],  sg = mu_sigma[[5]], alpha1 = .3, alpha = .025, 
                      sel_scen=1, side=T,test="w1",dropout=.1,rr=c(0.1,.1,.20,.40),bound=0)


round(c(e2$pval.surr,e2$pval1,e2$pval2,e2$concinvn,e2$concCIinvn),3)
e2$simdec_output #test decision 1=reject


#$pval.surr
#[1] 0.46633980 0.07062546

#$pval1
#p12low      p12med
#[1,] 0.9336625 0.02150906

#$pval2
#p12med2      p12hi2
#[1,] 0.005924966 0.005861396

#pval.surr.2<-e2$pval.surr #c(0.46633980, 0.07062546)
#pval1<-e2$pval1#c(0.9336625, 0.02150906)
#pval1.2<-pval1
#pval2<-e2$pval2  #c(0.005924966, 0.005861396)
#pval2.2<-pval2

alpha<-0.025
z <- qnorm(1-pmax(pval1,1e-15))
N=floor(120*(1-.1))+80#N1+N2
graph_bh <- BonferroniHolm(3)

v = c(N1/N,N1/N,0)
z1 <- c(z,0)
preplan <- doInterim(graph=graph_bh,z1=z1,v=v,alpha=alpha)

#H123
dec123<-max(pval2[1]<preplan@Aj[7,2],#preplan@BJ[7]*(preplan@Aj[7,2]/(preplan@Aj[7,2]+preplan@Aj[7,3]),
            pval2[2]<(preplan@Aj[7,1]+preplan@Aj[7,3]),#preplan@BJ[7]*preplan@Aj[7,3]/(preplan@Aj[7,2]+preplan@Aj[7,3]),
            preplan@BJ[7]>1)

preplan@Aj[7,2]
(preplan@Aj[7,1]+preplan@Aj[7,3])
preplan@BJ[7]

#H23
dec23<-max(pval2[1]<preplan@Aj[3,2],#preplan@BJ[3]*preplan@Aj[3,2]/(preplan@Aj[3,2]+preplan@Aj[3,3]),
           pval2[2]<preplan@Aj[3,3])#,BJ[3]*preplan@Aj[3,3]/(preplan@Aj[3,2]+preplan@Aj[3,3]),preplan@BJ[6]>1)

preplan@Aj[3,2]
preplan@Aj[3,3]

#H12
dec12<-max(pval2[1]<preplan@BJ[6])
preplan@BJ[6]

#H13
dec13<-max(pval2[2]<preplan@BJ[5])
preplan@BJ[5]

#H2
dec2<-min(pval2[1]<preplan@BJ[2],dec123,dec23,dec12)
preplan@BJ[2]

#H3
dec3<-min(pval2[2]<preplan@BJ[1],dec123,dec23,dec13)
preplan@BJ[1]
#




#########################################
#
#
#Create Table 5:

rbind(
  round(c(e1$pval.surr,e1$pval1,e1$pval2,e1$concinvn,e1$concCIinvn),3),
  round(c(e2$pval.surr,e2$pval1,e2$pval2,e2$concinvn,e2$concCIinvn),3)
)
#  c(t(pval.surr.1),t(pval1.1),t(pval2.1),0,e1$conccond,e1$concCIcond),
#c(t(pval.surr.2),t(pval1.2),0,t(pval2.2),e2$conccond,e2$concCIcond)
#),3)




#Histogram of stage 1 data;
pdf("plot_Histogram_exampleR.pdf")
par(mfrow=c(3,3))
hist(db_stage1[db_stage1$treat=="Placebo",1],xlim=c(0,10),ylim=c(0,15),breaks=10,xlab=c("log counts"),main=c("Stage 1: 0 months, Placebo"))
hist(db_stage1[db_stage1$treat=="Low",1],xlim=c(0,10),ylim=c(0,15),breaks=10,xlab=c("log counts"),main=c("Stage 1: 0 months, Low dose"))
hist(db_stage1[db_stage1$treat=="Medium",1],xlim=c(0,10),ylim=c(0,15),breaks=10,xlab=c("log counts"),main=c("Stage 1: 0 months, Med dose"))
hist(db_stage1[db_stage1$treat=="Placebo",2],xlim=c(0,10),ylim=c(0,15),breaks=20,xlab=c("log counts"),main=c("Stage 1: 6 months, Placebo"))
hist(db_stage1[db_stage1$treat=="Low",2],xlim=c(0,10),ylim=c(0,15),breaks=20,xlab=c("log counts"),main=c("Stage 1: 6 months, Low dose"))
hist(db_stage1[db_stage1$treat=="Medium",2],xlim=c(0,10),ylim=c(0,15),breaks=20,xlab=c("log counts"),main=c("Stage 1: 6 months, Med dose"))
hist(db_stage1[db_stage1$treat=="Placebo",3],xlim=c(0,10),ylim=c(0,15),breaks=20,xlab=c("log counts"),main=c("Stage 1: 12 months, Placebo"))
hist(db_stage1[db_stage1$treat=="Low",3],xlim=c(0,10),ylim=c(0,15),breaks=20,xlab=c("log counts"),main=c("Stage 1: 12 months, Low dose"))
hist(db_stage1[db_stage1$treat=="Medium",3],xlim=c(0,10),ylim=c(0,15),breaks=20,xlab=c("log counts"),main=c("Stage 1: 12months, Med dose"))
#dev.off()


par(mfrow=c(3,3))
hist(db_stage1[db_stage2$treat=="Placebo",1],xlim=c(0,10),ylim=c(0,15),breaks=10,xlab=c("log counts"),main=c("Stage 2: 0 months, Placebo"))
hist(db_stage1[db_stage2$treat=="Low",1],xlim=c(0,10),ylim=c(0,15),breaks=10,xlab=c("log counts"),main=c("Stage 2: 0 months, Low dose"))
hist(db_stage1[db_stage2$treat=="Medium",1],xlim=c(0,10),ylim=c(0,15),breaks=10,xlab=c("log counts"),main=c("Stage 2: 0 months, Med dose"))
hist(db_stage1[db_stage2$treat=="Placebo",2],xlim=c(0,10),ylim=c(0,15),breaks=20,xlab=c("log counts"),main=c("Stage 2: 6 months, Placebo"))
hist(db_stage1[db_stage2$treat=="Low",2],xlim=c(0,10),ylim=c(0,15),breaks=20,xlab=c("log counts"),main=c("Stage 2: 6 months, Low dose"))
hist(db_stage1[db_stage2$treat=="Medium",2],xlim=c(0,10),ylim=c(0,15),breaks=20,xlab=c("log counts"),main=c("Stage 2: 6 months, Med dose"))
hist(db_stage1[db_stage2$treat=="Placebo",3],xlim=c(0,10),ylim=c(0,15),breaks=20,xlab=c("log counts"),main=c("Stage 2: 12 months, Placebo"))
hist(db_stage1[db_stage2$treat=="Low",3],xlim=c(0,10),ylim=c(0,15),breaks=20,xlab=c("log counts"),main=c("Stage 2: 12 months, Low dose"))
hist(db_stage1[db_stage2$treat=="Medium",3],xlim=c(0,10),ylim=c(0,15),breaks=20,xlab=c("log counts"),main=c("Stage 2: 12months, Med dose"))


par(mfrow=c(3,3))
hist(e2$db_stage1[e2$db_stage1$treat=="Placebo",1],xlim=c(0,10),ylim=c(0,15),breaks=10,xlab=c("log counts"),main=c("Stage 1: 0 months, Placebo"))
hist(e2$db_stage1[e2$db_stage1$treat=="Low",1],xlim=c(0,10),ylim=c(0,15),breaks=10,xlab=c("log counts"),main=c("Stage 1: 0 months, Low dose"))
hist(e2$db_stage1[e2$db_stage1$treat=="Medium",1],xlim=c(0,10),ylim=c(0,15),breaks=10,xlab=c("log counts"),main=c("Stage 1: 0 months, Med dose"))
hist(e2$db_stage1[e2$db_stage1$treat=="Placebo",2],xlim=c(0,10),ylim=c(0,15),breaks=20,xlab=c("log counts"),main=c("Stage 1: 6 months, Placebo"))
hist(e2$db_stage1[e2$db_stage1$treat=="Low",2],xlim=c(0,10),ylim=c(0,15),breaks=20,xlab=c("log counts"),main=c("Stage 1: 6 months, Low dose"))
hist(e2$db_stage1[e2$db_stage1$treat=="Medium",2],xlim=c(0,10),ylim=c(0,15),breaks=20,xlab=c("log counts"),main=c("Stage 1: 6 months, Med dose"))
hist(e2$db_stage1[e2$db_stage1$treat=="Placebo",3],xlim=c(0,10),ylim=c(0,15),breaks=20,xlab=c("log counts"),main=c("Stage 1: 12 months, Placebo"))
hist(e2$db_stage1[e2$db_stage1$treat=="Low",3],xlim=c(0,10),ylim=c(0,15),breaks=20,xlab=c("log counts"),main=c("Stage 1: 12 months, Low dose"))
hist(e2$db_stage1[e2$db_stage1$treat=="Medium",3],xlim=c(0,10),ylim=c(0,15),breaks=20,xlab=c("log counts"),main=c("Stage 1: 12months, Med dose"))
#dev.off()


par(mfrow=c(3,3))
hist(e2$db_stage1[e2$db_stage2$treat=="Placebo",1],xlim=c(0,10),ylim=c(0,15),breaks=10,xlab=c("log counts"),main=c("Stage 2: 0 months, Placebo"))
hist(e2$db_stage1[e2$db_stage2$treat=="Medium",1],xlim=c(0,10),ylim=c(0,15),breaks=10,xlab=c("log counts"),main=c("Stage 2: 0 months, Med dose"))
hist(e2$db_stage1[e2$db_stage2$treat=="High",1],xlim=c(0,10),ylim=c(0,15),breaks=10,xlab=c("log counts"),main=c("Stage 2: 0 months, High dose"))
hist(e2$db_stage1[e2$db_stage2$treat=="Placebo",2],xlim=c(0,10),ylim=c(0,15),breaks=20,xlab=c("log counts"),main=c("Stage 2: 6 months, Placebo"))
hist(e2$db_stage1[e2$db_stage2$treat=="Medium",2],xlim=c(0,10),ylim=c(0,15),breaks=20,xlab=c("log counts"),main=c("Stage 2: 6 months, Med dose"))
hist(e2$db_stage1[e2$db_stage2$treat=="High",2],xlim=c(0,10),ylim=c(0,15),breaks=20,xlab=c("log counts"),main=c("Stage 2: 6 months, High dose"))
hist(e2$db_stage1[e2$db_stage2$treat=="Placebo",3],xlim=c(0,10),ylim=c(0,15),breaks=20,xlab=c("log counts"),main=c("Stage 2: 12 months, Placebo"))
hist(e2$db_stage1[e2$db_stage2$treat=="Medium",3],xlim=c(0,10),ylim=c(0,15),breaks=20,xlab=c("log counts"),main=c("Stage 2: 12months, Med dose"))
hist(e2$db_stage1[e2$db_stage2$treat=="High",3],xlim=c(0,10),ylim=c(0,15),breaks=20,xlab=c("log counts"),main=c("Stage 2: 12months, High dose"))
dev.off()


