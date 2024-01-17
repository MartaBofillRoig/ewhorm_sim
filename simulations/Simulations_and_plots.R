
#oncho

mu<-20
sigma<-28

mu_raw_0 = mu
sd_raw_0 = sigma 
rho = 0.5
n_trials=10000
n_arms = 4
N1<-c(60,60,60,90,90,90)
N = 150
rmonth=1
alpha1<-c(.1,.2,.3,.1,.2,.3)
alpha=.025
sim_out1=1
sel_scen=0
side1=1
test1=2

#no_effect
r0_6=0
r1_6=0
r2_6=0
r3_6=0
r0_12=0
r1_12=0
r2_12=0
r3_12=0

#simulation.global.null<-c(mu_raw_0, sd_raw_0, r0_6,r1_6,r2_6,r3_6, r0_12,r1_12,r2_12,r3_12,  rho ,   #data.matrix(expand.grid
#            n_trials,n_arms,N1 , N, rmonth, alpha1 , alpha,  sim_out1,sel_scen, side1, test1)#)

oncho_no_effect<-mapply(simul_res,mu_raw_0, sd_raw_0 , r0_6,r1_6,r2_6,r3_6, r0_12,r1_12,r2_12,r3_12,  rho ,
       n_trials,4,N1 , N, rmonth, alpha1 , alpha,
       sim_out1,sel_scen, side1,test1)



#low_effect
r0_6=0
r1_6=0.1
r2_6=0.2
r3_6=0.5
r0_12=0
r1_12=0.1
r2_12=0.2
r3_12=0.5

#simulation.global.null<-c(mu_raw_0, sd_raw_0, r0_6,r1_6,r2_6,r3_6, r0_12,r1_12,r2_12,r3_12,  rho ,   #data.matrix(expand.grid
#            n_trials,n_arms,N1 , N, rmonth, alpha1 , alpha,  sim_out1,sel_scen, side1, test1)#)

oncho_low_effect<-mapply(simul_res,mu_raw_0, sd_raw_0 , r0_6,r1_6,r2_6,r3_6, r0_12,r1_12,r2_12,r3_12,  rho ,
                        n_trials,4,N1 , N, rmonth, alpha1 , alpha,
                        sim_out1,sel_scen, side1,test1)


#med_effect
r0_6=0
r1_6=0.1
r2_6=0.5
r3_6=0.7
r0_12=0
r1_12=0.1
r2_12=0.5
r3_12=0.7

#simulation.global.null<-c(mu_raw_0, sd_raw_0, r0_6,r1_6,r2_6,r3_6, r0_12,r1_12,r2_12,r3_12,  rho ,   #data.matrix(expand.grid
#            n_trials,n_arms,N1 , N, rmonth, alpha1 , alpha,  sim_out1,sel_scen, side1, test1)#)

oncho_med_effect<-mapply(simul_res,mu_raw_0, sd_raw_0 , r0_6,r1_6,r2_6,r3_6, r0_12,r1_12,r2_12,r3_12,  rho ,
                         n_trials,4,N1 , N, rmonth, alpha1 , alpha,
                         sim_out1,sel_scen, side1,test1)



#high_effect
r0_6=0
r1_6=0.5
r2_6=0.6
r3_6=0.7
r0_12=0
r1_12=0.5
r2_12=0.6
r3_12=0.7

#simulation.global.null<-c(mu_raw_0, sd_raw_0, r0_6,r1_6,r2_6,r3_6, r0_12,r1_12,r2_12,r3_12,  rho ,   #data.matrix(expand.grid
#            n_trials,n_arms,N1 , N, rmonth, alpha1 , alpha,  sim_out1,sel_scen, side1, test1)#)

oncho_high_effect<-mapply(simul_res,mu_raw_0, sd_raw_0 , r0_6,r1_6,r2_6,r3_6, r0_12,r1_12,r2_12,r3_12,  rho ,
                         n_trials,4,N1 , N, rmonth, alpha1 , alpha,
                         sim_out1,sel_scen, side1,test1)




#datwide<-reshape(dat, idvar = "SubjNr_code", timevar = "FeNO.1_marker", direction = "wide")

#Hyp1
par (mfrow=c(4,4), mar=c(1, 1, 1, 1), oma=c(3, 3, 1, 1))
matplot(1:3,oncho_no_effect[c(1,4,7),],type="p",pch=1:6,xaxt="n",ylim=c(0,1))
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
legend("topright",legend=c(expression(N[1]==60),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3),expression(N[1]==90),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3)),cex=.8,col=c("white",1:3,"white",4:6),lwd=2.5,ncol=2,bty="n",pch=c(1,1:3,1,4:6))
mtext(side=3,"Low Dose",cex=1.3)#,line=2.2)
matplot(1:3,oncho_no_effect[c(2,5,8),],type="p",pch=1:6,xaxt="n",yaxt="n",ylim=c(0,1))
mtext(side=3,"Intermediate Dose",cex=1.3)#,line=2.2)
matplot(1:3,oncho_no_effect[c(3,6,9),],type="p",pch=1:6,xaxt="n",yaxt="n",ylim=c(0,1))
mtext(side=3,"High Dose",cex=1.3)#,line=2.2)
matplot(1:3,oncho_no_effect[c(10,13,17),],type="p",pch=1:6,xaxt="n",yaxt="n",ylim=c(0,1))
mtext(side=3,"Disjunctive Power",cex=1.3)#,line=2.2)
mtext(side=4,"No effective",cex=1,line=1)#,line=2.2)


matplot(1:3,oncho_low_effect[c(1,4,7),],type="p",pch=1:6,xaxt="n",ylim=c(0,1))
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
matplot(1:3,oncho_low_effect[c(2,5,8),],type="p",pch=1:6,xaxt="n",yaxt="n",ylim=c(0,1))
matplot(1:3,oncho_low_effect[c(3,6,9),],type="p",pch=1:6,xaxt="n",yaxt="n",ylim=c(0,1))
matplot(1:3,oncho_low_effect[c(10,13,17),],type="p",pch=1:6,xaxt="n",yaxt="n",ylim=c(0,1))
mtext(side=4,"low effective",cex=1,line=1)#,line=2.2)

matplot(1:3,oncho_med_effect[c(1,4,7),],type="p",pch=1:6,xaxt="n",ylim=c(0,1))
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
matplot(1:3,oncho_med_effect[c(2,5,8),],type="p",pch=1:6,xaxt="n",yaxt="n",ylim=c(0,1))
matplot(1:3,oncho_med_effect[c(3,6,9),],type="p",pch=1:6,xaxt="n",yaxt="n",ylim=c(0,1))
matplot(1:3,oncho_med_effect[c(10,13,17),],type="p",pch=1:6,xaxt="n",yaxt="n",ylim=c(0,1))
mtext(side=4,"Median effective",cex=1,line=1)#,line=2.2)

matplot(1:3,oncho_high_effect[c(1,4,7),],type="p",pch=1:6,xaxt="n",ylim=c(0,1))
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
axis(1,1:3,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1)
matplot(1:3,oncho_high_effect[c(2,5,8),],type="p",pch=1:6,xaxt="n",yaxt="n",ylim=c(0,1))
axis(1,1:3,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1)
matplot(1:3,oncho_high_effect[c(3,6,9),],type="p",pch=1:6,xaxt="n",yaxt="n",ylim=c(0,1))
axis(1,1:3,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1)
matplot(1:3,oncho_high_effect[c(10,13,17),],type="p",pch=1:6,xaxt="n",yaxt="n",ylim=c(0,1))
mtext(side=4,"High effective",cex=1,line=1)
axis(1,1:3,c("Adaptive","MA1","MA2"),padj=-0.3,cex.axis=1)


simul_res(mu_raw_0, sd_raw_0 , r0_6,r1_6,r2_6,r3_6, r0_12,r1_12,r2_12,r3_12,  rho ,   #data.matrix(expand.grid
          n_trials,4,N1 , N, rmonth, alpha1 , alpha,
          sim_out1,sel_scen, side1,test1)

apply(simulation.global.null,1,simul_res)
simul_res(as.numeric(simulation.global.null[1,]))


apply(simulation.global.null,2,simul_res)
mapply(simul_res,simulation.global.null)
sapply(simulation.global.null,simul_res)
lapply(simulation.global.null,simul_res)


simul_res(simulation.global.null[1,])


simul_res(mu_raw_0 = mu, sd_raw_0 = sigma, r0_6 = 0,r1_6 = 0,r2_6 = 0,r3_6 = 0,r0_12 = 0, r1_12 = 0,r2_12 = 0,r3_12 = 0, 
          rho = 0.5, n_trials=10000,n_arms = 4,N1 = 60 , 
          N = 90,rmonth=1, alpha1=.1 , alpha=.025, sim_out=T,sel_scen=0, side=T,test="l")

simul_res(2.0e+01, 2.8e+01, 0.0e+00, 0.0e+00, 0.0e+00, 0.0e+00, 0.0e+00, 0.0e+00, 0.0e+00, 0.0e+00, 5.0e-01, 1.0e+04, 4.0e+00, 6.0e+01, 1.5e+02, 1.0e+00 ,1.0e-01, 2.5e-02,
          TRUE,0,TRUE,"t")


sims<-function(mu_raw_0 = mu, sd_raw_0 = sigma, reductrate_6 = c(0,0,0,0), reductrate_12 = c(0,0,0,0), 
               rho = 0.5, n_trials=10000,n_arms = 4,N1 = c(60,75,90) , 
               N2 = 90,rmonth=1, alpha1=c(.1,.2,.3) , alpha=.025, sim_out=T,sel_scen=0, side=T,test="l")
{
 for (i in 1:length(N1)) 
   for (j in 1:length(alpha1))
     simul_res(mu_raw_0 = mu, sd_raw_0 = sigma, reductrate_6 = c(0,0,0,0), reductrate_12 = c(0,0,0,0), 
               rho = 0.5, n_trials=10000,n_arms = 4,N1 = 60 , 
               N2 = 90,rmonth=1, alpha1=.1 , alpha=.025, sim_out=T,sel_scen=0, side=T,test="l")
}

results[i,]<-simul_res(mu_raw_0 = mu, sd_raw_0 = sigma, reductrate_6 = c(0,0,0,0), reductrate_12 = c(0,0,0,0), 
                      rho = 0.5, n_trials=10000,n_arms = 4,N1 = 60 , 
                      N2 = 90,rmonth=1, alpha1=.1 , alpha=.025, sim_out=T,sel_scen=0, side=T,test="l")


oo.oncho1<-mapply(simul_res, mu, sigma, c(0, 0, 0, 0),c(0,0,0,0),  .5, 10000,4, 60 , 90, 1, c(.1,.2,.5), .025, T,0, T,"l")
oo.oncho1<-simul_res( mu, sigma, c(0, 0, 0, 0),c(0,0,0,0),  .5, 10000,4, 60 , 90, 1, .1, .025, T,0, T,"l")
oo.oncho1<-simul_res (mu_raw_0 = mu, sd_raw_0 = sigma, reductrate_6 = c(0,0,0,0), reductrate_12 = c(0,0,0,0), rho = 0.5, n_trials=10000,n_arms = 4,N1 = 60 , 
           N2 = 90,rmonth=1, alpha1=.1 , alpha=.025, sim_out=T,sel_scen=0, side=T,test="l")



plot(unlist(oo.oncho1)[c(1:6,10)],ylim=c(0,1))
lines(4:6,unlist(oo.oncho1)[7:9],type="p",col="red")
lines(c(4,5,7),unlist(oo.oncho1)[7:9],type="p",col="red")
lines(4:6,unlist(oo.oncho1)[7:9],type="p",col="red")


title(expression(alpha[1]==0.1))
lines(unlist(oo.oncho1[,2]),ylim=c(0,1),type="p",col="red")
lines(unlist(oo.oncho1[,3]),ylim=c(0,1),type="p",col="blue")
legend("topleft",legend=c(expression(rho==0),expression(rho==.5),expression(rho==0.7)),cex=.8,col=c(1,"red","blue"),lwd=2.5,ncol=1,bty="n")
abline(v=3.5)




oo.oncho2<-mapply(simul_res, mu, sigma, 0.15, 0.20, 0.25, 0.3,  c(0,0.5,.7,0,0.5,.7), 10000,4, 60 , 90, 1, c(.1,.1,.1,.5,.5,.5), .025, T,0, T,"t")

oo.oncho3<-mapply(simul_res, mu, sigma, 0.15, 0.3, 0.6, 0.9,  c(0,0.5,.7,0,0.5,.7), 10000,4, 60 , 90, 1, c(.1,.1,.1,.5,.5,.5), .025, T,0, T,"t")

oo.oncho4<-mapply(simul_res, mu, sigma, 0.15, 0.3, 0.3, 0.9,  c(0,0.5,.7,0,0.5,.7), 10000,4, 60 , 90, 1, c(.1,.1,.1,.5,.5,.5), .025, T,0, T,"t")

oo.oncho5<-mapply(simul_res, mu, sigma, 0.15, 0.3, 0.3, 0.9,  c(0,0.5,.7,0,0.5,.7), 10000,4, 60 , 90, 1, c(.1,.1,.1,.5,.5,.5), .025, T,0, T,"l")

#oo.oncho6<-mapply(simul_res, mu, sigma, 0.15, 0.3, 0.3, 0.9,  c(0,0.5,.7,0,0.5,.7), 10000,4, 60 , 90, 1, c(.1,.1,.1,.5,.5,.5), .025, T,1, T,"t")

pdf(file ="Oncho.pdf", width = 7, height = 10, pointsize = 12, paper = "special")

par(mfrow=c(5,2), mar=c(1, 1, 1, 1), oma=c(3, 3, 1, 1))
plot(unlist(oo.oncho1[,1]),ylim=c(0,1))
title(expression(alpha[1]==0.1))
lines(unlist(oo.oncho1[,2]),ylim=c(0,1),type="p",col="red")
lines(unlist(oo.oncho1[,3]),ylim=c(0,1),type="p",col="blue")
legend("topleft",legend=c(expression(rho==0),expression(rho==.5),expression(rho==0.7)),cex=.8,col=c(1,"red","blue"),lwd=2.5,ncol=1,bty="n")
abline(v=3.5)

plot(unlist(oo.oncho1[,4]),ylim=c(0,1))
title(expression(alpha[1]==0.5))
lines(unlist(oo.oncho1[,5]),ylim=c(0,1),type="p",col="red")
lines(unlist(oo.oncho1[,6]),ylim=c(0,1),type="p",col="blue")
abline(v=3.5)

plot(unlist(oo.oncho2[,1]),ylim=c(0,1))
lines(unlist(oo.oncho2[,2]),ylim=c(0,1),type="p",col="red")
lines(unlist(oo.oncho2[,3]),ylim=c(0,1),type="p",col="blue")
abline(v=3.5)
plot(unlist(oo.oncho2[,4]),ylim=c(0,1))
lines(unlist(oo.oncho2[,5]),ylim=c(0,1),type="p",col="red")
lines(unlist(oo.oncho2[,6]),ylim=c(0,1),type="p",col="blue")
abline(v=3.5)

plot(unlist(oo.oncho3[,1]),ylim=c(0,1))
lines(unlist(oo.oncho3[,2]),ylim=c(0,1),type="p",col="red")
lines(unlist(oo.oncho3[,3]),ylim=c(0,1),type="p",col="blue")
abline(v=3.5)
plot(unlist(oo.oncho3[,4]),ylim=c(0,1))
lines(unlist(oo.oncho3[,5]),ylim=c(0,1),type="p",col="red")
lines(unlist(oo.oncho3[,6]),ylim=c(0,1),type="p",col="blue")
abline(v=3.5)

plot(unlist(oo.oncho4[,1]),ylim=c(0,1))
lines(unlist(oo.oncho4[,2]),ylim=c(0,1),type="p",col="red")
lines(unlist(oo.oncho4[,3]),ylim=c(0,1),type="p",col="blue")
abline(v=3.5)
plot(unlist(oo.oncho4[,4]),ylim=c(0,1))
lines(unlist(oo.oncho4[,5]),ylim=c(0,1),type="p",col="red")
lines(unlist(oo.oncho4[,6]),ylim=c(0,1),type="p",col="blue")
abline(v=3.5)



plot(unlist(oo.oncho5[,1]),ylim=c(0,1))
lines(unlist(oo.oncho5[,2]),ylim=c(0,1),type="p",col="red")
lines(unlist(oo.oncho5[,3]),ylim=c(0,1),type="p",col="blue")
abline(v=3.5)
plot(unlist(oo.oncho5[,4]),ylim=c(0,1))
lines(unlist(oo.oncho5[,5]),ylim=c(0,1),type="p",col="red")
lines(unlist(oo.oncho5[,6]),ylim=c(0,1),type="p",col="blue")
abline(v=3.5)

dev.off()


#mans

mu<-510
sigma<-470

oo.mans1<-mapply(simul_res, mu, sigma, 0, 0, 0, 0,  c(0,0.5,.7,0,0.5,.7), 10000,4, 60 , 90, 1, c(.1,.1,.1,.5,.5,.5), .025, T,0, T,"t")


oo.mans2<-mapply(simul_res, mu, sigma, 0.15, 0.20, 0.25, 0.3,  c(0,0.5,.7,0,0.5,.7), 10000,4, 60 , 90, 1, c(.1,.1,.1,.5,.5,.5), .025, T,0, T,"t")

oo.mans3<-mapply(simul_res, mu, sigma, 0.15, 0.3, 0.6, 0.9,  c(0,0.5,.7,0,0.5,.7), 10000,4, 60 , 90, 1, c(.1,.1,.1,.5,.5,.5), .025, T,0, T,"t")

oo.mans4<-mapply(simul_res, mu, sigma, 0.15, 0.3, 0.3, 0.9,  c(0,0.5,.7,0,0.5,.7), 10000,4, 60 , 90, 1, c(.1,.1,.1,.5,.5,.5), .025, T,0, T,"t")

oo.mans5<-mapply(simul_res, mu, sigma, 0.15, 0.3, 0.3, 0.9,  c(0,0.5,.7,0,0.5,.7), 10000,4, 60 , 90, 1, c(.1,.1,.1,.5,.5,.5), .025, T,0, T,"l")

pdf(file ="Mans.pdf", width = 7, height = 10, pointsize = 12, paper = "special")

par(mfrow=c(5,2), mar=c(1, 1, 1, 1), oma=c(3, 3, 1, 1))
plot(unlist(oo.mans1[,1]),ylim=c(0,1))
title(expression(alpha[1]==0.1))
lines(unlist(oo.mans1[,2]),ylim=c(0,1),type="p",col="red")
lines(unlist(oo.mans1[,3]),ylim=c(0,1),type="p",col="blue")
legend("topleft",legend=c(expression(rho==0),expression(rho==.5),expression(rho==0.7)),cex=.8,col=c(1,"red","blue"),lwd=2.5,ncol=1,bty="n")
abline(v=3.5)

plot(unlist(oo.mans1[,4]),ylim=c(0,1))
title(expression(alpha[1]==0.5))
lines(unlist(oo.mans1[,5]),ylim=c(0,1),type="p",col="red")
lines(unlist(oo.mans1[,6]),ylim=c(0,1),type="p",col="blue")
abline(v=3.5)

plot(unlist(oo.mans2[,1]),ylim=c(0,1))
lines(unlist(oo.mans2[,2]),ylim=c(0,1),type="p",col="red")
lines(unlist(oo.mans2[,3]),ylim=c(0,1),type="p",col="blue")
abline(v=3.5)
plot(unlist(oo.mans2[,4]),ylim=c(0,1))
lines(unlist(oo.mans2[,5]),ylim=c(0,1),type="p",col="red")
lines(unlist(oo.mans2[,6]),ylim=c(0,1),type="p",col="blue")
abline(v=3.5)

plot(unlist(oo.mans3[,1]),ylim=c(0,1))
lines(unlist(oo.mans3[,2]),ylim=c(0,1),type="p",col="red")
lines(unlist(oo.mans3[,3]),ylim=c(0,1),type="p",col="blue")
abline(v=3.5)
plot(unlist(oo.mans3[,4]),ylim=c(0,1))
lines(unlist(oo.mans3[,5]),ylim=c(0,1),type="p",col="red")
lines(unlist(oo.mans3[,6]),ylim=c(0,1),type="p",col="blue")
abline(v=3.5)

plot(unlist(oo.mans4[,1]),ylim=c(0,1))
lines(unlist(oo.mans4[,2]),ylim=c(0,1),type="p",col="red")
lines(unlist(oo.mans4[,3]),ylim=c(0,1),type="p",col="blue")
abline(v=3.5)
plot(unlist(oo.mans4[,4]),ylim=c(0,1))
lines(unlist(oo.mans4[,5]),ylim=c(0,1),type="p",col="red")
lines(unlist(oo.mans4[,6]),ylim=c(0,1),type="p",col="blue")
abline(v=3.5)

plot(unlist(oo.mans5[,1]),ylim=c(0,1))
lines(unlist(oo.mans5[,2]),ylim=c(0,1),type="p",col="red")
lines(unlist(oo.mans5[,3]),ylim=c(0,1),type="p",col="blue")
abline(v=3.5)
plot(unlist(oo.mans5[,4]),ylim=c(0,1))
lines(unlist(oo.mans5[,5]),ylim=c(0,1),type="p",col="red")
lines(unlist(oo.mans5[,6]),ylim=c(0,1),type="p",col="blue")
abline(v=3.5)

dev.off()


#loa

mu<-650
sigma<-575

oo.loa1<-mapply(simul_res, mu, sigma, 0, 0, 0, 0,  c(0,0.5,.7,0,0.5,.7), 10000,4, 60 , 90, 1, c(.1,.1,.1,.5,.5,.5), .025, T,0, T,"t")


oo.loa2<-mapply(simul_res, mu, sigma, 0.15, 0.20, 0.25, 0.3,  c(0,0.5,.7,0,0.5,.7), 10000,4, 60 , 90, 1, c(.1,.1,.1,.5,.5,.5), .025, T,0, T,"t")

oo.loa3<-mapply(simul_res, mu, sigma, 0.15, 0.3, 0.6, 0.9,  c(0,0.5,.7,0,0.5,.7), 10000,4, 60 , 90, 1, c(.1,.1,.1,.5,.5,.5), .025, T,0, T,"t")

oo.loa4<-mapply(simul_res, mu, sigma, 0.15, 0.3, 0.3, 0.9,  c(0,0.5,.7,0,0.5,.7), 10000,4, 60 , 90, 1, c(.1,.1,.1,.5,.5,.5), .025, T,0, T,"t")

oo.loa5<-mapply(simul_res, mu, sigma, 0.15, 0.3, 0.3, 0.9,  c(0,0.5,.7,0,0.5,.7), 10000,4, 60 , 90, 1, c(.1,.1,.1,.5,.5,.5), .025, T,0, T,"l")

oo.loa6<-mapply(simul_res, mu, sigma, 0.15, 0.3, 0.3, 0.9,  c(0,0.5,.7,0,0.5,.7), 10000,4, 60 , 90, 1, c(.1,.1,.1,.5,.5,.5), .025, T,1, T,"t")



pdf(file ="Loa.pdf", width = 7, height = 12, pointsize = 12, paper = "special")

par(mfrow=c(6,2), mar=c(1, 1, 1, 1), oma=c(3, 3, 1, 1))
plot(unlist(oo.loa1[,1]),ylim=c(0,1))
abline(v=3.5)
title(expression(alpha[1]==0.1))
lines(unlist(oo.loa1[,2]),ylim=c(0,1),type="p",col="red")
lines(unlist(oo.loa1[,3]),ylim=c(0,1),type="p",col="blue")
legend("topleft",legend=c(expression(rho==0),expression(rho==.5),expression(rho==0.7)),cex=.8,col=c(1,"red","blue"),lwd=2.5,ncol=1,bty="n")
abline(v=3.5)

plot(unlist(oo.loa1[,4]),ylim=c(0,1))
title(expression(alpha[1]==0.5))
lines(unlist(oo.loa1[,5]),ylim=c(0,1),type="p",col="red")
lines(unlist(oo.loa1[,6]),ylim=c(0,1),type="p",col="blue")
abline(v=3.5)

plot(unlist(oo.loa2[,1]),ylim=c(0,1))
lines(unlist(oo.loa2[,2]),ylim=c(0,1),type="p",col="red")
lines(unlist(oo.loa2[,3]),ylim=c(0,1),type="p",col="blue")
abline(v=3.5)
plot(unlist(oo.loa2[,4]),ylim=c(0,1))
lines(unlist(oo.loa2[,5]),ylim=c(0,1),type="p",col="red")
lines(unlist(oo.loa2[,6]),ylim=c(0,1),type="p",col="blue")
abline(v=3.5)

plot(unlist(oo.loa3[,1]),ylim=c(0,1))
lines(unlist(oo.loa3[,2]),ylim=c(0,1),type="p",col="red")
lines(unlist(oo.loa3[,3]),ylim=c(0,1),type="p",col="blue")
abline(v=3.5)
plot(unlist(oo.loa3[,4]),ylim=c(0,1))
lines(unlist(oo.loa3[,5]),ylim=c(0,1),type="p",col="red")
lines(unlist(oo.loa3[,6]),ylim=c(0,1),type="p",col="blue")
abline(v=3.5)

plot(unlist(oo.loa4[,1]),ylim=c(0,1))
lines(unlist(oo.loa4[,2]),ylim=c(0,1),type="p",col="red")
lines(unlist(oo.loa4[,3]),ylim=c(0,1),type="p",col="blue")
abline(v=3.5)
plot(unlist(oo.loa4[,4]),ylim=c(0,1))
lines(unlist(oo.loa4[,5]),ylim=c(0,1),type="p",col="red")
lines(unlist(oo.loa4[,6]),ylim=c(0,1),type="p",col="blue")
abline(v=3.5)


plot(unlist(oo.loa5[,1]),ylim=c(0,1))
lines(unlist(oo.loa5[,2]),ylim=c(0,1),type="p",col="red")
lines(unlist(oo.loa5[,3]),ylim=c(0,1),type="p",col="blue")
abline(v=3.5)
plot(unlist(oo.loa5[,4]),ylim=c(0,1))
lines(unlist(oo.loa5[,5]),ylim=c(0,1),type="p",col="red")
lines(unlist(oo.loa5[,6]),ylim=c(0,1),type="p",col="blue")
abline(v=3.5)

plot(unlist(oo.loa6[,1]),ylim=c(0,1))
lines(unlist(oo.loa6[,2]),ylim=c(0,1),type="p",col="red")
lines(unlist(oo.loa6[,3]),ylim=c(0,1),type="p",col="blue")
abline(v=3.5)
plot(unlist(oo.loa6[,4]),ylim=c(0,1))
lines(unlist(oo.loa6[,5]),ylim=c(0,1),type="p",col="red")
lines(unlist(oo.loa6[,6]),ylim=c(0,1),type="p",col="blue")
abline(v=3.5)

dev.off()