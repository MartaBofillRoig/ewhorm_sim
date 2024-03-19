#save.image(file='eWHORM.RData')
#load('yoursession.RData')


###################################################################
#oncho

mu<-19
sigma<-30
rmonth=1


##################################
#
# Simulation 1
#


r0_6no=0
r1_6no=0
r2_6no=0
r3_6no=0
r0_12no=0
r1_12no=0
r2_12no=0
r3_12no=0

r0_6lo=0
r1_6lo=0
r2_6lo=0
r3_6lo=0.5
r0_12lo=0
r1_12lo=0
r2_12lo=0
r3_12lo=0.6

r0_6me=0
r1_6me=0
r2_6me=0.3
r3_6me=0.5
r0_12me=0
r1_12me=0
r2_12me=0.5
r3_12me=0.6

r0_6hi=0
r1_6hi=0.5
r2_6hi=0.5
r3_6hi=0.5
r0_12hi=0
r1_12hi=0.6
r2_12hi=0.6
r3_12hi=0.6

pdf(file ="Scenarios.pdf", width = 10, height = 10, pointsize = 12, paper = "special")
plot(1:3,c(r1_6lo,r2_6lo,r3_6lo),type="b",lty=2,ylim=c(0,1),ylab="Reduction rate",xlab="Doses",lwd=2,xaxt="n",cex.lab=1.5,cex.axis=1.3)
axis(1,c(1:3),c("low dose","medium dose","high dose"),cex.axis=1.3,cex=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
legend("topleft",legend=c("Only high dose effective","Trend","All doses effective","6 months","6 months","6 months","12 months","12 months","12 months"),lwd=2,col=c("white","white","white",1,2,"blue",1,2,"blue"),bty="n",lty=c(1,1,1,2,2,2),ncol=3,cex=1.3)
lines(1:3,c(r1_12lo,r2_12lo,r3_12lo),type="b",lty=1,lwd=2)
lines(1:3,c(r1_12me,r2_12me,r3_12me),type="b",lty=1,lwd=2,col=2)
lines(1:3,c(r1_6me,r2_6me,r3_6me),type="b",lty=2,lwd=2,col=2)            
lines(1:3,c(r1_12hi,r2_12hi,r3_12hi),type="b",lty=1,lwd=2,col="blue")
lines(1:3,c(r1_6hi,r2_6hi,r3_6hi),type="b",lty=2,lwd=2,col="blue")            
dev.off()

pdf(file ="Scenarios_Trichuris.pdf", width = 10, height = 10, pointsize = 12, paper = "special")
plot(1:3,c(r1_12lo,r2_12lo,r3_12lo),type="b",lty=1,ylim=c(0,1),ylab="Reduction rate",xlab="Doses",lwd=2,xaxt="n",cex.lab=1.5,cex.axis=1.3)
axis(1,c(1:3),c("low dose","medium dose","high dose"),cex.axis=1.3,cex=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
legend("topleft",legend=c("Only high dose effective","Trend","All doses effective"),lwd=2,col=c(1,2,"blue"),bty="n",lty=c(1,1,1),ncol=1,cex=1.3)
lines(1:3,c(r1_12me,r2_12me,r3_12me),type="b",lty=1,lwd=2,col=2)
lines(1:3,c(r1_12hi,r2_12hi,r3_12hi),type="b",lty=1,lwd=2,col="blue")
dev.off()






mu_raw_0 = mu
sd_raw_0 = sigma 
rho = 0.5
n_trials=10000
n_arms = 4
N1<-c(80,80,80,120,120,120)
N = 200-N1

alpha1<-c(.1,.2,.3,.1,.2,.3)
alpha=.025
sim_out1=1
sel_scen=0
side1=1
test1=0
dropout=.1
rr0=0
rr1=0
rr2=0
rr3=0

#no_effect

#simulation.global.null<-c(mu_raw_0, sd_raw_0, r0_6,r1_6,r2_6,r3_6, r0_12,r1_12,r2_12,r3_12,  rho ,   #data.matrix(expand.grid
#            n_trials,n_arms,N1 , N, rmonth, alpha1 , alpha,  sim_out1,sel_scen, side1, test1)#)

oncho_no_effect<-mapply(simul_res,mu_raw_0, sd_raw_0 , r0_6no,r1_6no,r2_6no,r3_6no, r0_12no,r1_12no,r2_12no,r3_12no,  rho ,
       n_trials,4,N1 , N, rmonth, alpha1 , alpha,
       sim_out1,sel_scen, side1,test1,dropout,rr0,rr1,rr2,rr3)



#low_effect
oncho_low_effect<-mapply(simul_res,mu_raw_0, sd_raw_0 , r0_6lo,r1_6lo,r2_6lo,r3_6lo, r0_12lo,r1_12lo,r2_12lo,r3_12lo,  rho ,
                        n_trials,4,N1 , N, rmonth, alpha1 , alpha,
                        sim_out1,sel_scen, side1,test1,dropout,rr0,rr1,rr2,rr3)


#med_effect

oncho_med_effect<-mapply(simul_res,mu_raw_0, sd_raw_0 , r0_6me,r1_6me,r2_6me,r3_6me, r0_12me,r1_12me,r2_12me,r3_12me,  rho ,
                         n_trials,4,N1 , N, rmonth, alpha1 , alpha,
                         sim_out1,sel_scen, side1,test1,dropout,rr0,rr1,rr2,rr3)



#high_effect

oncho_high_effect<-mapply(simul_res,mu_raw_0, sd_raw_0 , r0_6hi,r1_6hi,r2_6hi,r3_6hi, r0_12hi,r1_12hi,r2_12hi,r3_12hi,  rho ,
                         n_trials,4,N1 , N, rmonth, alpha1 , alpha,
                         sim_out1,sel_scen, side1,test1,dropout,rr0,rr1,rr2,rr3)




#Power
pdf(file ="Oncho_Pow.pdf", width = 12, height = 9, pointsize = 12, paper = "special")

par (mfrow=c(3,4), mar=c(1, 1, 1, 1), oma=c(3, 3, 1, 1))
#plot(1:3,oncho_no_effect[7,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1)
#mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
#mtext(side=3,"Low Dose",cex=1.3)#,line=2.2)
#lines(1:3,oncho_no_effect[7,4:6],type="b",lwd=2,col=2)
#lines(1:3,oncho_no_effect[11,1:3],type="b",lwd=2,col=3,lty=2)
#lines(1:3,oncho_no_effect[15,1:3],type="b",lwd=2,col=5,lty=2)

#plot(1:3,oncho_no_effect[8,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
#mtext(side=3,"Medium Dose",cex=1.3)#,line=2.2)
#lines(1:3,oncho_no_effect[8,4:6],type="b",lwd=2,col=2)
#lines(1:3,oncho_no_effect[12,1:3],type="b",lwd=2,col=3,lty=2)
#lines(1:3,oncho_no_effect[16,1:3],type="b",lwd=2,col=5,lty=2)

#plot(1:3,oncho_no_effect[9,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
#mtext(side=3,"High Dose",cex=1.3)#,line=2.2)
#lines(1:3,oncho_no_effect[9,4:6],type="b",lwd=2,col=2)
#lines(1:3,oncho_no_effect[13,1:3],type="b",lwd=2,col=3,lty=2)
#lines(1:3,oncho_no_effect[17,1:3],type="b",lwd=2,col=5,lty=2)

#plot(1:3,oncho_no_effect[10,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
#mtext(side=3,"Disjunctive",cex=1.3)#,line=2.2)
#lines(1:3,oncho_no_effect[10,4:6],type="b",lwd=2,col=2)
#lines(1:3,oncho_no_effect[14,1:3],type="b",lwd=2,col=3,lty=2)
#lines(1:3,oncho_no_effect[18,1:3],type="b",lwd=2,col=5,lty=2)
#mtext(side=4,"Global null hypothesis",cex=1,line=1)#,line=2.2)



plot(1:3,oncho_low_effect[7,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1,cex.axis=1.3)
legend("top",legend=c(expression(N[1]==80),expression(N[1]==120),"MA1","MA2"),lwd=2,col=c(1,2,3,4,"white","white"),lty=c(1,1,2,2),bty="n",cex=1.3)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
mtext(side=3,"Low Dose",cex=1.3)#,line=2.2)
lines(1:3,oncho_low_effect[7,4:6],type="b",lty=1,lwd=2,col=2)
lines(1:3,oncho_low_effect[11,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_low_effect[15,1:3],type="b",lwd=2,col=5,lty=2)

plot(1:3,oncho_low_effect[8,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
mtext(side=3,"Medium Dose",cex=1.3)#,line=2.2)
lines(1:3,oncho_low_effect[8,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_low_effect[12,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_low_effect[16,1:3],type="b",lwd=2,col=5,lty=2)

plot(1:3,oncho_low_effect[9,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
mtext(side=3,"High Dose",cex=1.3)#,line=2.2)
lines(1:3,oncho_low_effect[9,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_low_effect[13,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_low_effect[17,1:3],type="b",lwd=2,col=5,lty=2)

plot(1:3,oncho_low_effect[10,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_low_effect[10,4:6],type="b",lwd=2,col=2)
mtext(side=3,"Disjunctive",cex=1.3)#,line=2.2)lines(1:3,oncho_low_effect[10,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_low_effect[14,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_low_effect[18,1:3],type="b",lwd=2,col=5,lty=2)
mtext(side=4,"Scenario high dose eff.",cex=1,line=1)#,line=2.2)


plot(1:3,oncho_med_effect[7,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1,cex.axis=1.3)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
lines(1:3,oncho_med_effect[7,4:6],type="b",lty=1,lwd=2,col=2)
lines(1:3,oncho_med_effect[11,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_med_effect[15,1:3],type="b",lwd=2,col=5,lty=2)

plot(1:3,oncho_med_effect[8,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_med_effect[8,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_med_effect[12,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_med_effect[16,1:3],type="b",lwd=2,col=5,lty=2)

plot(1:3,oncho_med_effect[9,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_med_effect[9,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_med_effect[13,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_med_effect[17,1:3],type="b",lwd=2,col=5,lty=2)

plot(1:3,oncho_med_effect[10,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_med_effect[10,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_med_effect[14,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_med_effect[18,1:3],type="b",lwd=2,col=5,lty=2)
mtext(side=4,"Scenario Trend",cex=1,line=1)#,line=2.2)

plot(1:3,oncho_high_effect[7,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1,cex.axis=1.3)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
lines(1:3,oncho_high_effect[7,4:6],type="b",lty=1,lwd=2,col=2)
lines(1:3,oncho_high_effect[11,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_high_effect[15,1:3],type="b",lwd=2,col=5,lty=2)
axis(1,c(1:3),c("0.1","0.2","0.3"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,expression(alpha[1]),cex=1,line=2)#,line=2.2)

plot(1:3,oncho_high_effect[8,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_high_effect[8,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_high_effect[12,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_high_effect[16,1:3],type="b",lwd=2,col=5,lty=2)
axis(1,c(1:3),c("0.1","0.2","0.3"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,expression(alpha[1]),cex=1,line=2)#,line=2.2)

plot(1:3,oncho_high_effect[9,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_high_effect[9,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_high_effect[13,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_high_effect[17,1:3],type="b",lwd=2,col=5,lty=2)
axis(1,c(1:3),c("0.1","0.2","0.3"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,expression(alpha[1]),cex=1,line=2)#,line=2.2)

plot(1:3,oncho_high_effect[10,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_high_effect[10,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_high_effect[14,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_high_effect[18,1:3],type="b",lwd=2,col=5,lty=2)
mtext(side=4,"Scenario all dose eff.",cex=1,line=1)#,line=2.2)
axis(1,c(1:3),c("0.1","0.2","0.3"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,expression(alpha[1]),cex=1,line=2)#,line=2.2)

dev.off()



#Selection
pdf(file ="Oncho_selection.pdf", width = 9, height = 12, pointsize = 12, paper = "special")

par (mfrow=c(4,3), mar=c(1, 1, 1, 1), oma=c(3, 3, 1, 1))
plot(1:3,oncho_no_effect[1,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1)
legend("top",legend=c(expression(N[1]==80),expression(N[1]==120)),lwd=2,col=c(1,2),lty=c(1,1),bty="n",cex=1.3)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
mtext(side=3,"Low Dose",cex=1.3)#,line=2.2)
#legend("topright",legend=c(expression(N[1]==60),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3),expression(N[1]==90),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3)),cex=.8,lwd=2,col=c("white",1:3,"white",4:6),lwd=2.5,ncol=2,bty="n",pch=c(1,1:3,1,4:6))
lines(1:3,oncho_no_effect[1,4:6],type="b",lwd=2,col=2)

plot(1:3,oncho_no_effect[2,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
mtext(side=3,"Medium Dose",cex=1.3)#,line=2.2)
lines(1:3,oncho_no_effect[2,4:6],type="b",lwd=2,col=2)

plot(1:3,oncho_no_effect[3,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
mtext(side=3,"High Dose",cex=1.3)#,line=2.2)
lines(1:3,oncho_no_effect[3,4:6],type="b",lwd=2,col=2)

mtext(side=4,"Global null hypothesis",cex=1,line=1)#,line=2.2)



plot(1:3,oncho_low_effect[1,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1,cex.axis=1.3)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
lines(1:3,oncho_low_effect[1,4:6],type="b",lty=1,lwd=2,col=2)

plot(1:3,oncho_low_effect[2,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_low_effect[2,4:6],type="b",lwd=2,col=2)

plot(1:3,oncho_low_effect[3,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_low_effect[3,4:6],type="b",lwd=2,col=2)
mtext(side=4,"Scenario high dose eff.",cex=1,line=1)#,line=2.2)


plot(1:3,oncho_med_effect[1,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1,cex.axis=1.3)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
lines(1:3,oncho_med_effect[1,4:6],type="b",lty=1,lwd=2,col=2)

plot(1:3,oncho_med_effect[2,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_med_effect[2,4:6],type="b",lwd=2,col=2)

plot(1:3,oncho_med_effect[3,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_med_effect[3,4:6],type="b",lwd=2,col=2)
mtext(side=4,"Scenario Trend",cex=1,line=1)#,line=2.2)

plot(1:3,oncho_high_effect[1,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1,cex.axis=1.3)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
lines(1:3,oncho_high_effect[1,4:6],type="b",lty=1,lwd=2,col=2)
axis(1,c(1:3),c("0.1","0.2","0.3"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,expression(alpha[1]),cex=1,line=2)#,line=2.2)

plot(1:3,oncho_high_effect[2,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_high_effect[2,4:6],type="b",lwd=2,col=2)
axis(1,c(1:3),c("0.1","0.2","0.3"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,expression(alpha[1]),cex=1,line=2)#,line=2.2)

plot(1:3,oncho_high_effect[3,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_high_effect[3,4:6],type="b",lwd=2,col=2)
axis(1,c(1:3),c("0.1","0.2","0.3"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,expression(alpha[1]),cex=1,line=2)#,line=2.2)
mtext(side=4,"Scenario all dose eff.",cex=1,line=1)#,line=2.2)
dev.off()

#Conditional Power


pdf(file ="Oncho_condpower.pdf", width = 9, height = 9, pointsize = 12, paper = "special")
par (mfrow=c(3,3), mar=c(1, 1, 1, 1), oma=c(3, 3, 1, 1))
#plot(1:3,oncho_no_effect[4,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1)
#mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
#mtext(side=3,"Low Dose",cex=1.3)#,line=2.2)
#legend("topright",legend=c(expression(N[1]==60),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3),expression(N[1]==90),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3)),cex=.8,lwd=2,col=c("white",1:3,"white",4:6),lwd=2.5,ncol=2,bty="n",pch=c(1,1:3,1,4:6))
#lines(1:3,oncho_no_effect[4,4:6],type="b",lwd=2,col=2)

#plot(1:3,oncho_no_effect[5,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
#mtext(side=3,"Medium Dose",cex=1.3)#,line=2.2)
#lines(1:3,oncho_no_effect[5,4:6],type="b",lwd=2,col=2)

#plot(1:3,oncho_no_effect[6,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
#mtext(side=3,"High Dose",cex=1.3)#,line=2.2)
#lines(1:3,oncho_no_effect[6,4:6],type="b",lwd=2,col=2)

#mtext(side=4,"Global null hypothesis",cex=1,line=1)#,line=2.2)



plot(1:3,oncho_low_effect[4,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1,cex.axis=1.3)
legend("top",legend=c(expression(N[1]==80),expression(N[1]==120)),lwd=2,col=c(1,2),lty=c(1,1),bty="n",cex=1.3)
mtext(side=3,"Low Dose",cex=1.3)#,line=2.2)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
lines(1:3,oncho_low_effect[4,4:6],type="b",lty=1,lwd=2,col=2)

plot(1:3,oncho_low_effect[5,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
mtext(side=3,"Medium Dose",cex=1.3)#,line=2.2)
lines(1:3,oncho_low_effect[5,4:6],type="b",lwd=2,col=2)

plot(1:3,oncho_low_effect[6,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
mtext(side=3,"High Dose",cex=1.3)#,line=2.2)
lines(1:3,oncho_low_effect[6,4:6],type="b",lwd=2,col=2)
mtext(side=4,"Scenario high dose eff.",cex=1,line=1)#,line=2.2)


plot(1:3,oncho_med_effect[4,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1,cex.axis=1.3)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
lines(1:3,oncho_med_effect[4,4:6],type="b",lty=1,lwd=2,col=2)

plot(1:3,oncho_med_effect[5,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_med_effect[5,4:6],type="b",lwd=2,col=2)

plot(1:3,oncho_med_effect[6,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_med_effect[6,4:6],type="b",lwd=2,col=2)
mtext(side=4,"Scenario Trend",cex=1,line=1)#,line=2.2)

plot(1:3,oncho_high_effect[4,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1,cex.axis=1.3)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
lines(1:3,oncho_high_effect[4,4:6],type="b",lty=1,lwd=2,col=2)
axis(1,c(1:3),c("0.1","0.2","0.3"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,expression(alpha[1]),cex=1,line=2)#,line=2.2)

plot(1:3,oncho_high_effect[5,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_high_effect[5,4:6],type="b",lwd=2,col=2)
axis(1,c(1:3),c("0.1","0.2","0.3"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,expression(alpha[1]),cex=1,line=2)#,line=2.2)

plot(1:3,oncho_high_effect[6,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_high_effect[6,4:6],type="b",lwd=2,col=2)
axis(1,c(1:3),c("0.1","0.2","0.3"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,expression(alpha[1]),cex=1,line=2)#,line=2.2)
mtext(side=4,"Scenario all dose eff.",cex=1,line=1)#,line=2.2)

dev.off()



##################################
#
# Simulation 2
#




###############################################################################################
###############################################################################################

#
#
#     Compare  linear model - t-test Wilcoxon
#
#


N1<-120
N<-200-N1
alpha1<-.3
sel_scen=c(0,0,0,1,1,1)
test1=c(0,2,3,0,2,3)

#no_effect

oncho_no_effect_test<-mapply(simul_res,mu_raw_0, sd_raw_0 , r0_6no,r1_6no,r2_6no,r3_6no, r0_12no,r1_12no,r2_12no,r3_12no,  rho ,
                        n_trials,4,N1 , N, rmonth, alpha1 , alpha,
                        sim_out1,sel_scen, side1,test1,dropout,rr0,rr1,rr2,rr3)

#low_effect
oncho_low_effect_test<-mapply(simul_res,mu_raw_0, sd_raw_0 , r0_6lo,r1_6lo,r2_6lo,r3_6lo, r0_12lo,r1_12lo,r2_12lo,r3_12lo,  rho ,
                         n_trials,4,N1 , N, rmonth, alpha1 , alpha,
                         sim_out1,sel_scen, side1,test1,dropout,rr0,rr1,rr2,rr3)

#med_effect

oncho_med_effect_test<-mapply(simul_res,mu_raw_0, sd_raw_0 , r0_6me,r1_6me,r2_6me,r3_6me, r0_12me,r1_12me,r2_12me,r3_12me,  rho ,
                         n_trials,4,N1 , N, rmonth, alpha1 , alpha,
                         sim_out1,sel_scen, side1,test1,dropout,rr0,rr1,rr2,rr3)

#high_effect

oncho_high_effect_test<-mapply(simul_res,mu_raw_0, sd_raw_0 , r0_6hi,r1_6hi,r2_6hi,r3_6hi, r0_12hi,r1_12hi,r2_12hi,r3_12hi,  rho ,
                          n_trials,4,N1 , N, rmonth, alpha1 , alpha,
                          sim_out1,sel_scen, side1,test1,dropout,rr0,rr1,rr2,rr3)







#Power
pdf(file ="Oncho_Pow_test.pdf", width = 12, height = 9, pointsize = 12, paper = "special")

par (mfrow=c(3,4), mar=c(1, 1, 1, 1), oma=c(3, 3, 1, 1))
#plot(1:3,oncho_no_effect_test[7,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1)
#mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
#mtext(side=3,"Low Dose",cex=1.3)#,line=2.2)
##legend("topright",legend=c(expression(N[1]==60),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3),expression(N[1]==90),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3)),cex=.8,lwd=2,col=c("white",1:3,"white",4:6),lwd=2.5,ncol=2,bty="n",pch=c(1,1:3,1,4:6))
#lines(1:3,oncho_no_effect_test[7,4:6],type="b",lwd=2,col=2)
#lines(1:3,oncho_no_effect_test[11,1:3],type="b",lwd=2,col=3,lty=2)
##lines(1:3,oncho_no_effect_test[11,4:6],type="b",lwd=2,col=4) #sollte das selbe sein wie in voriger Zeile, daher weglassen
#lines(1:3,oncho_no_effect_test[15,1:3],type="b",lwd=2,col=5,lty=2)

#plot(1:3,oncho_no_effect_test[8,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
#mtext(side=3,"Medium Dose",cex=1.3)#,line=2.2)
#lines(1:3,oncho_no_effect_test[8,4:6],type="b",lwd=2,col=2)
#lines(1:3,oncho_no_effect_test[12,1:3],type="b",lwd=2,col=3,lty=2)
#lines(1:3,oncho_no_effect_test[16,1:3],type="b",lwd=2,col=5,lty=2)#

#plot(1:3,oncho_no_effect_test[9,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
#mtext(side=3,"High Dose",cex=1.3)#,line=2.2)
#lines(1:3,oncho_no_effect_test[9,4:6],type="b",lwd=2,col=2)
#lines(1:3,oncho_no_effect_test[13,1:3],type="b",lwd=2,col=3,lty=2)
#lines(1:3,oncho_no_effect_test[17,1:3],type="b",lwd=2,col=5,lty=2)

#plot(1:3,oncho_no_effect_test[10,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
#mtext(side=3,"Disjunctive",cex=1.3)#,line=2.2)
#lines(1:3,oncho_no_effect_test[10,4:6],type="b",lwd=2,col=2)
#lines(1:3,oncho_no_effect_test[14,1:3],type="b",lwd=2,col=3,lty=2)
#lines(1:3,oncho_no_effect_test[18,1:3],type="b",lwd=2,col=5,lty=2)
#mtext(side=4,"Global null hypothesis",cex=1,line=1)#,line=2.2)



plot(1:3,oncho_low_effect_test[7,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1,cex.axis=1.3)
legend("top",legend=c("Selection strategy B","Selection strategy A","MA1","MA2"),lwd=2,col=c(1,2,3,4),lty=c(1,1,2,2),bty="n",cex=1.3)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
mtext(side=3,"Low Dose",cex=1.3)#,line=2.2)
lines(1:3,oncho_low_effect_test[7,4:6],type="b",lty=1,lwd=2,col=2)
lines(1:3,oncho_low_effect_test[11,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_low_effect_test[15,1:3],type="b",lwd=2,col=5,lty=2)

plot(1:3,oncho_low_effect_test[8,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
mtext(side=3,"Medium Dose",cex=1.3)#,line=2.2)
lines(1:3,oncho_low_effect_test[8,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_low_effect_test[12,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_low_effect_test[16,1:3],type="b",lwd=2,col=5,lty=2)

plot(1:3,oncho_low_effect_test[9,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
mtext(side=3,"High Dose",cex=1.3)#,line=2.2)
lines(1:3,oncho_low_effect_test[9,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_low_effect_test[13,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_low_effect_test[17,1:3],type="b",lwd=2,col=5,lty=2)

plot(1:3,oncho_low_effect_test[10,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
mtext(side=3,"Disjunctive",cex=1.3)#,line=2.2)
lines(1:3,oncho_low_effect_test[10,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_low_effect_test[14,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_low_effect_test[18,1:3],type="b",lwd=2,col=5,lty=2)
mtext(side=4,"Scenario high dose eff.",cex=1,line=1)#,line=2.2)


plot(1:3,oncho_med_effect_test[7,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1,cex.axis=1.3)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
lines(1:3,oncho_med_effect_test[7,4:6],type="b",lty=1,lwd=2,col=2)
lines(1:3,oncho_med_effect_test[11,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_med_effect_test[15,1:3],type="b",lwd=2,col=5,lty=2)

plot(1:3,oncho_med_effect_test[8,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_med_effect_test[8,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_med_effect_test[12,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_med_effect_test[16,1:3],type="b",lwd=2,col=5,lty=2)

plot(1:3,oncho_med_effect_test[9,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_med_effect_test[9,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_med_effect_test[13,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_med_effect_test[17,1:3],type="b",lwd=2,col=5,lty=2)

plot(1:3,oncho_med_effect_test[10,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_med_effect_test[10,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_med_effect_test[14,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_med_effect_test[18,1:3],type="b",lwd=2,col=5,lty=2)
mtext(side=4,"Scenario Trend",cex=1,line=1)#,line=2.2)

plot(1:3,oncho_high_effect_test[7,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1,cex.axis=1.3)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
lines(1:3,oncho_high_effect_test[7,4:6],type="b",lty=1,lwd=2,col=2)
lines(1:3,oncho_high_effect_test[11,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_high_effect_test[15,1:3],type="b",lwd=2,col=5,lty=2)
axis(1,c(1:3),c("lm","t","cc"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,"Analysis",cex=1,line=2)#,line=2.2)

plot(1:3,oncho_high_effect_test[8,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_high_effect_test[8,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_high_effect_test[12,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_high_effect_test[16,1:3],type="b",lwd=2,col=5,lty=2)
axis(1,c(1:3),c("lm","t","cc"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,"Analysis",cex=1,line=2)#,line=2.2)

plot(1:3,oncho_high_effect_test[9,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_high_effect_test[9,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_high_effect_test[13,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_high_effect_test[17,1:3],type="b",lwd=2,col=5,lty=2)
axis(1,c(1:3),c("lm","t","cc"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,"Analysis",cex=1,line=2)#,line=2.2)

plot(1:3,oncho_high_effect_test[10,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_high_effect_test[10,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_high_effect_test[14,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_high_effect_test[18,1:3],type="b",lwd=2,col=5,lty=2)
mtext(side=4,"Scenario all dose eff.",cex=1,line=1)#,line=2.2)
axis(1,c(1:3),c("lm","t","cc"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,"Analysis",cex=1,line=2)#,line=2.2)

dev.off()



#Selection
pdf(file ="Oncho_selection_test.pdf", width = 9, height = 12, pointsize = 12, paper = "special")

par (mfrow=c(4,3), mar=c(1, 1, 1, 1), oma=c(3, 3, 1, 1))
plot(1:3,oncho_no_effect_test[1,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1,cex.axis=1.3)
legend("top",legend=c("Selection strategy B","Selection strategy A"),lwd=2,col=c(1,2,3,4),lty=c(1,1,2,2),bty="n",cex=1.3)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
mtext(side=3,"Low Dose",cex=1.3)#,line=2.2)
#legend("topright",legend=c(expression(N[1]==60),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3),expression(N[1]==90),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3)),cex=.8,lwd=2,col=c("white",1:3,"white",4:6),lwd=2.5,ncol=2,bty="n",pch=c(1,1:3,1,4:6))
lines(1:3,oncho_no_effect_test[1,4:6],type="b",lwd=2,col=2)

plot(1:3,oncho_no_effect_test[2,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
mtext(side=3,"Medium Dose",cex=1.3)#,line=2.2)
lines(1:3,oncho_no_effect_test[2,4:6],type="b",lwd=2,col=2)

plot(1:3,oncho_no_effect_test[3,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
mtext(side=3,"High Dose",cex=1.3)#,line=2.2)
lines(1:3,oncho_no_effect_test[3,4:6],type="b",lwd=2,col=2)

mtext(side=4,"Global null hypothesis",cex=1,line=1)#,line=2.2)



plot(1:3,oncho_low_effect_test[1,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1,cex.axis=1.3)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
lines(1:3,oncho_low_effect_test[1,4:6],type="b",lty=1,lwd=2,col=2)

plot(1:3,oncho_low_effect_test[2,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_low_effect_test[2,4:6],type="b",lwd=2,col=2)

plot(1:3,oncho_low_effect_test[3,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_low_effect_test[3,4:6],type="b",lwd=2,col=2)
mtext(side=4,"Scenario high dose eff.",cex=1,line=1)#,line=2.2)


plot(1:3,oncho_med_effect_test[1,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1,cex.axis=1.3)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
lines(1:3,oncho_med_effect_test[1,4:6],type="b",lty=1,lwd=2,col=2)

plot(1:3,oncho_med_effect_test[2,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_med_effect_test[2,4:6],type="b",lwd=2,col=2)

plot(1:3,oncho_med_effect_test[3,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_med_effect_test[3,4:6],type="b",lwd=2,col=2)
mtext(side=4,"Scenario Trend",cex=1,line=1)#,line=2.2)

plot(1:3,oncho_high_effect_test[1,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1,cex.axis=1.3)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
lines(1:3,oncho_high_effect_test[1,4:6],type="b",lty=1,lwd=2,col=2)
axis(1,c(1:3),c("lm","t","cc"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,"Analysis",cex=1,line=2)#,line=2.2)

plot(1:3,oncho_high_effect_test[2,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_high_effect_test[2,4:6],type="b",lwd=2,col=2)
axis(1,c(1:3),c("lm","t","cc"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,"Analysis",cex=1,line=2)#,line=2.2)

plot(1:3,oncho_high_effect_test[3,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_high_effect_test[3,4:6],type="b",lwd=2,col=2)
axis(1,c(1:3),c("lm","t","cc"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,"Analysis",cex=1,line=2)#,line=2.2)
mtext(side=4,"Scenario all dose eff.",cex=1,line=1)#,line=2.2)
dev.off()

#Conditional Power


pdf(file ="Oncho_condpower_test.pdf", width = 9, height = 9, pointsize = 12, paper = "special")
par (mfrow=c(3,3), mar=c(1, 1, 1, 1), oma=c(3, 3, 1, 1))
#plot(1:3,oncho_no_effect_test[4,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1)
#mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
#mtext(side=3,"Low Dose",cex=1.3)#,line=2.2)
##legend("topright",legend=c(expression(N[1]==60),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3),expression(N[1]==90),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3)),cex=.8,lwd=2,col=c("white",1:3,"white",4:6),lwd=2.5,ncol=2,bty="n",pch=c(1,1:3,1,4:6))
#lines(1:3,oncho_no_effect_test[4,4:6],type="b",lwd=2,col=2)

#plot(1:3,oncho_no_effect_test[5,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
#mtext(side=3,"Medium Dose",cex=1.3)#,line=2.2)
#lines(1:3,oncho_no_effect_test[5,4:6],type="b",lwd=2,col=2)

#plot(1:3,oncho_no_effect_test[6,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
#mtext(side=3,"High Dose",cex=1.3)#,line=2.2)
#lines(1:3,oncho_no_effect_test[6,4:6],type="b",lwd=2,col=2)
#mtext(side=4,"Global null hypothesis",cex=1,line=1)#,line=2.2)



plot(1:3,oncho_low_effect_test[4,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1,cex.axis=1.3)
legend("top",legend=c("Selection strategy B","Selection strategy A"),lwd=2,col=c(1,2,3,4),lty=c(1,1,2,2),bty="n",cex=1.3)
mtext(side=3,"Low Dose",cex=1.3)#,line=2.2)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
lines(1:3,oncho_low_effect_test[4,4:6],type="b",lty=1,lwd=2,col=2)

plot(1:3,oncho_low_effect_test[5,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
mtext(side=3,"Medium Dose",cex=1.3)#,line=2.2)
lines(1:3,oncho_low_effect_test[5,4:6],type="b",lwd=2,col=2)

plot(1:3,oncho_low_effect_test[6,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
mtext(side=3,"High Dose",cex=1.3)#,line=2.2)
lines(1:3,oncho_low_effect_test[6,4:6],type="b",lwd=2,col=2)
mtext(side=4,"Scenario high dose eff.",cex=1,line=1)#,line=2.2)


plot(1:3,oncho_med_effect_test[4,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1,cex.axis=1.3)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
lines(1:3,oncho_med_effect_test[4,4:6],type="b",lty=1,lwd=2,col=2)

plot(1:3,oncho_med_effect_test[5,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_med_effect_test[5,4:6],type="b",lwd=2,col=2)

plot(1:3,oncho_med_effect_test[6,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_med_effect_test[6,4:6],type="b",lwd=2,col=2)
mtext(side=4,"Scenario Trend",cex=1,line=1)#,line=2.2)

plot(1:3,oncho_high_effect_test[4,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1,cex.axis=1.3)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
lines(1:3,oncho_high_effect_test[4,4:6],type="b",lty=1,lwd=2,col=2)
axis(1,c(1:3),c("lm","t","cc"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,"Analysis",cex=1,line=2)#,line=2.2)

plot(1:3,oncho_high_effect_test[5,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_high_effect_test[5,4:6],type="b",lwd=2,col=2)
axis(1,c(1:3),c("lm","t","cc"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,"Analysis",cex=1,line=2)#,line=2.2)

plot(1:3,oncho_high_effect_test[6,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_high_effect_test[6,4:6],type="b",lwd=2,col=2)
axis(1,c(1:3),c("0.1","0.2","0.3"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
axis(1,c(1:3),c("lm","t","cc"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,"Analysis",cex=1,line=2)#,line=2.2)

dev.off()




#Error


pdf(file ="Oncho_Error.pdf", width = 12, height = 6, pointsize = 12, paper = "special")

par (mfrow=c(2,4), mar=c(1, 1, 1, 1), oma=c(3, 3, 1, 1))
plot(1:3,oncho_no_effect[7,1:3],type="b",xaxt="n",ylim=c(0,0.05),lty=1,lwd=2,col=1,cex.axis=1.3)
legend("top",legend=c(expression(N[1]==80),expression(N[1]==120),"MA1","MA2"),lwd=2,col=c(1,2,3,4),lty=c(1,1,2,2),bty="n",cex=1.3)
mtext(side=2,"FWER",cex=1,line=2.1)#,line=2.2)
mtext(side=3,"Low Dose",cex=1.3)#,line=2.2)
lines(1:3,oncho_no_effect[7,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_no_effect[11,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_no_effect[15,1:3],type="b",lwd=2,col=5,lty=2)
abline(h=0.025)
axis(1,c(1:3),c("0.1","0.2","0.3"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,expression(alpha[1]),cex=1,line=2)#,line=2.2)

plot(1:3,oncho_no_effect[8,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,0.05),lwd=2,col=1)
mtext(side=3,"Medium Dose",cex=1.3)#,line=2.2)
lines(1:3,oncho_no_effect[8,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_no_effect[12,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_no_effect[16,1:3],type="b",lwd=2,col=5,lty=2)
abline(h=0.025)
axis(1,c(1:3),c("0.1","0.2","0.3"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,expression(alpha[1]),cex=1,line=2)#,line=2.2)

plot(1:3,oncho_no_effect[9,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,0.05),lwd=2,col=1)
mtext(side=3,"High Dose",cex=1.3)#,line=2.2)
lines(1:3,oncho_no_effect[9,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_no_effect[13,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_no_effect[17,1:3],type="b",lwd=2,col=5,lty=2)
abline(h=0.025)
axis(1,c(1:3),c("0.1","0.2","0.3"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,expression(alpha[1]),cex=1,line=2)#,line=2.2)

plot(1:3,oncho_no_effect[10,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,0.05),lwd=2,col=1)
mtext(side=3,"Total Study FWER",cex=1.3)#,line=2.2)
lines(1:3,oncho_no_effect[10,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_no_effect[14,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_no_effect[18,1:3],type="b",lwd=2,col=5,lty=2)
mtext(side=4,"Global null hypothesis",cex=1,line=1)#,line=2.2)
abline(h=0.025)
axis(1,c(1:3),c("0.1","0.2","0.3"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,expression(alpha[1]),cex=1,line=2)#,line=2.2)



plot(1:3,oncho_no_effect_test[7,1:3],type="b",xaxt="n",ylim=c(0,0.05),lty=1,lwd=2,col=1,cex.axis=1.3)
legend("top",legend=c("Selection strategy B","Selection strategy A","MA1","MA2"),lwd=2,col=c(1,2,3,4),lty=c(1,1,2,2),bty="n",cex=1.3)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
##legend("topright",legend=c(expression(N[1]==60),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3),expression(N[1]==90),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3)),cex=.8,lwd=2,col=c("white",1:3,"white",4:6),lwd=2.5,ncol=2,bty="n",pch=c(1,1:3,1,4:6))
lines(1:3,oncho_no_effect_test[7,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_no_effect_test[11,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_no_effect_test[15,1:3],type="b",lwd=2,col=5,lty=2)
abline(h=0.025)
axis(1,c(1:3),c("lm","t","cc"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,"Analysis",cex=1,line=2)#,line=2.2)

plot(1:3,oncho_no_effect_test[8,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,0.05),lwd=2,col=1)
lines(1:3,oncho_no_effect_test[8,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_no_effect_test[12,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_no_effect_test[16,1:3],type="b",lwd=2,col=5,lty=2)#
abline(h=0.025)
axis(1,c(1:3),c("lm","t","cc"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,"Analysis",cex=1,line=2)#,line=2.2)

plot(1:3,oncho_no_effect_test[9,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,0.05),lwd=2,col=1)
lines(1:3,oncho_no_effect_test[9,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_no_effect_test[13,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_no_effect_test[17,1:3],type="b",lwd=2,col=5,lty=2)
abline(h=0.025)
axis(1,c(1:3),c("lm","t","cc"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,"Analysis",cex=1,line=2)#,line=2.2)

plot(1:3,oncho_no_effect_test[10,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,0.05),lwd=2,col=1)
lines(1:3,oncho_no_effect_test[10,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_no_effect_test[14,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_no_effect_test[18,1:3],type="b",lwd=2,col=5,lty=2)
mtext(side=4,"Global null hypothesis",cex=1,line=1)#,line=2.2)
abline(h=0.025)
axis(1,c(1:3),c("lm","t","cc"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,"Analysis",cex=1,line=2)#,line=2.2)

dev.off()















##################################
#
# Simulation 3
#




###############################################################################################
###############################################################################################

#
#
#     Add total responder; Compare  linear model - t-test Wilcoxon
#
#


sel_scen=0#c(0,0,0)
test1=c(0,2,3,4)
rr0<-.1
rr1<-.1
rr2<-.1
rr3<-.1
#no_effect

oncho_no_effect_resp<-mapply(simul_res,mu_raw_0, sd_raw_0 , r0_6no,r1_6no,r2_6no,r3_6no, r0_12no,r1_12no,r2_12no,r3_12no,  rho ,
                             n_trials,4,N1 , N, rmonth, alpha1 , alpha,
                             sim_out1,sel_scen, side1,test1,dropout,rr0,rr1,rr2,rr3)

rr0<-.1
rr1<-.1
rr2<-.1
rr3<-r3_12lo-.1

#low_effect
oncho_low_effect_resp<-mapply(simul_res,mu_raw_0, sd_raw_0 , r0_6lo,r1_6lo,r2_6lo,r3_6lo, r0_12lo,r1_12lo,r2_12lo,r3_12lo,  rho ,
                              n_trials,4,N1 , N, rmonth, alpha1 , alpha,
                              sim_out1,sel_scen, side1,test1,dropout,rr0,rr1,rr2,rr3)

#med_effect
rr0<-.1
rr1<-.1
rr2<-r2_12me-.1
rr3<-r3_12me-.1

oncho_med_effect_resp<-mapply(simul_res,mu_raw_0, sd_raw_0 , r0_6me,r1_6me,r2_6me,r3_6me, r0_12me,r1_12me,r2_12me,r3_12me,  rho ,
                              n_trials,4,N1 , N, rmonth, alpha1 , alpha,
                              sim_out1,sel_scen, side1,test1,dropout,rr0,rr1,rr2,rr3)

#high_effect

rr0<-.1
rr1<-r1_12hi-.1
rr2<-r2_12hi-.1
rr3<-r3_12hi-.1

oncho_high_effect_resp<-mapply(simul_res,mu_raw_0, sd_raw_0 , r0_6hi,r1_6hi,r2_6hi,r3_6hi, r0_12hi,r1_12hi,r2_12hi,r3_12hi,  rho ,
                               n_trials,4,N1 , N, rmonth, alpha1 , alpha,
                               sim_out1,sel_scen, side1,test1,dropout,rr0,rr1,rr2,rr3)







#Power



#Power
pdf(file ="Oncho_Pow_resp.pdf", width = 12, height = 9, pointsize = 12, paper = "special")

par (mfrow=c(3,4), mar=c(1, 1, 1, 1), oma=c(3, 3, 1, 1))
#plot(1:4,oncho_no_effect_resp[7,1:4],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1)
#mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
#mtext(side=3,"Low Dose",cex=1.3)#,line=2.2)
##legend("topright",legend=c(expression(N[1]==60),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3),expression(N[1]==90),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3)),cex=.8,lwd=2,col=c("white",1:4,"white",4:6),lwd=2.5,ncol=2,bty="n",pch=c(1,1:4,1,4:6))
#lines(1:4,oncho_no_effect_resp[7,4:6],type="b",lwd=2,col=2)
#lines(1:4,oncho_no_effect_resp[11,1:4],type="b",lwd=2,col=3,lty=2)
##lines(1:4,oncho_no_effect_resp[11,4:6],type="b",lwd=2,col=4) #sollte das selbe sein wie in voriger Zeile, daher weglassen
#lines(1:4,oncho_no_effect_resp[15,1:4],type="b",lwd=2,col=5,lty=2)

#plot(1:4,oncho_no_effect_resp[8,1:4],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
#mtext(side=3,"Medium Dose",cex=1.3)#,line=2.2)
#lines(1:4,oncho_no_effect_resp[8,4:6],type="b",lwd=2,col=2)
#lines(1:4,oncho_no_effect_resp[12,1:4],type="b",lwd=2,col=3,lty=2)
#lines(1:4,oncho_no_effect_resp[16,1:4],type="b",lwd=2,col=5,lty=2)#

#plot(1:4,oncho_no_effect_resp[9,1:4],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
#mtext(side=3,"High Dose",cex=1.3)#,line=2.2)
#lines(1:4,oncho_no_effect_resp[9,4:6],type="b",lwd=2,col=2)
#lines(1:4,oncho_no_effect_resp[13,1:4],type="b",lwd=2,col=3,lty=2)
#lines(1:4,oncho_no_effect_resp[17,1:4],type="b",lwd=2,col=5,lty=2)

#plot(1:4,oncho_no_effect_resp[10,1:4],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
#mtext(side=3,"Disjunctive",cex=1.3)#,line=2.2)
#lines(1:4,oncho_no_effect_resp[10,4:6],type="b",lwd=2,col=2)
#lines(1:4,oncho_no_effect_resp[14,1:4],type="b",lwd=2,col=3,lty=2)
#lines(1:4,oncho_no_effect_resp[18,1:4],type="b",lwd=2,col=5,lty=2)
#mtext(side=4,"Global null hypothesis",cex=1,line=1)#,line=2.2)



plot(1:4,oncho_low_effect_resp[7,1:4],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1,cex.axis=1.3)
#legend("top",legend=c(expression(N[1]==80),expression(N[1]==120),"MA1","MA2"),lwd=2,col=c(1,2,3,4),lty=c(1,1,2,2),bty="n")
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
mtext(side=3,"Low Dose",cex=1.3)#,line=2.2)
#lines(1:4,oncho_low_effect_resp[7,4:6],type="b",lty=1,lwd=2,col=2)
lines(1:4,oncho_low_effect_resp[11,1:4],type="b",lwd=2,col=3,lty=2)
lines(1:4,oncho_low_effect_resp[15,1:4],type="b",lwd=2,col=5,lty=2)
legend("top",legend=c("Adaptive design","MA1","MA2"),lwd=2,col=c(1,3,4),lty=c(1,2,2),bty="n",cex=1.3)

plot(1:4,oncho_low_effect_resp[8,1:4],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
mtext(side=3,"Medium Dose",cex=1.3)#,line=2.2)
#lines(1:4,oncho_low_effect_resp[8,4:6],type="b",lwd=2,col=2)
lines(1:4,oncho_low_effect_resp[12,1:4],type="b",lwd=2,col=3,lty=2)
lines(1:4,oncho_low_effect_resp[16,1:4],type="b",lwd=2,col=5,lty=2)

plot(1:4,oncho_low_effect_resp[9,1:4],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
mtext(side=3,"High Dose",cex=1.3)#,line=2.2)
#lines(1:4,oncho_low_effect_resp[9,4:6],type="b",lwd=2,col=2)
lines(1:4,oncho_low_effect_resp[13,1:4],type="b",lwd=2,col=3,lty=2)
lines(1:4,oncho_low_effect_resp[17,1:4],type="b",lwd=2,col=5,lty=2)

plot(1:4,oncho_low_effect_resp[10,1:4],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
mtext(side=3,"Disjunctive",cex=1.3)#,line=2.2)
#lines(1:4,oncho_low_effect_resp[10,4:6],type="b",lwd=2,col=2)
lines(1:4,oncho_low_effect_resp[14,1:4],type="b",lwd=2,col=3,lty=2)
lines(1:4,oncho_low_effect_resp[18,1:4],type="b",lwd=2,col=5,lty=2)
mtext(side=4,"Scenario high dose eff.",cex=1,line=1)#,line=2.2)


plot(1:4,oncho_med_effect_resp[7,1:4],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1,cex.axis=1.3)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
#lines(1:4,oncho_med_effect_resp[7,4:6],type="b",lty=1,lwd=2,col=2)
lines(1:4,oncho_med_effect_resp[11,1:4],type="b",lwd=2,col=3,lty=2)
lines(1:4,oncho_med_effect_resp[15,1:4],type="b",lwd=2,col=5,lty=2)

plot(1:4,oncho_med_effect_resp[8,1:4],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
#lines(1:4,oncho_med_effect_resp[8,4:6],type="b",lwd=2,col=2)
lines(1:4,oncho_med_effect_resp[12,1:4],type="b",lwd=2,col=3,lty=2)
lines(1:4,oncho_med_effect_resp[16,1:4],type="b",lwd=2,col=5,lty=2)

plot(1:4,oncho_med_effect_resp[9,1:4],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
#lines(1:4,oncho_med_effect_resp[9,4:6],type="b",lwd=2,col=2)
lines(1:4,oncho_med_effect_resp[13,1:4],type="b",lwd=2,col=3,lty=2)
lines(1:4,oncho_med_effect_resp[17,1:4],type="b",lwd=2,col=5,lty=2)

plot(1:4,oncho_med_effect_resp[10,1:4],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
#lines(1:4,oncho_med_effect_resp[10,4:6],type="b",lwd=2,col=2)
lines(1:4,oncho_med_effect_resp[14,1:4],type="b",lwd=2,col=3,lty=2)
lines(1:4,oncho_med_effect_resp[18,1:4],type="b",lwd=2,col=5,lty=2)
mtext(side=4,"Scenario Trend",cex=1,line=1)#,line=2.2)

plot(1:4,oncho_high_effect_resp[7,1:4],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1,cex.axis=1.3)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
#lines(1:4,oncho_high_effect_resp[7,4:6],type="b",lty=1,lwd=2,col=2)
lines(1:4,oncho_high_effect_resp[11,1:4],type="b",lwd=2,col=3,lty=2)
lines(1:4,oncho_high_effect_resp[15,1:4],type="b",lwd=2,col=5,lty=2)
axis(1,c(1:4),c("lm","t","cc","c"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,"Analysis",cex=1,line=2)#,line=2.2)

plot(1:4,oncho_high_effect_resp[8,1:4],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
#lines(1:4,oncho_high_effect_resp[8,4:6],type="b",lwd=2,col=2)
lines(1:4,oncho_high_effect_resp[12,1:4],type="b",lwd=2,col=3,lty=2)
lines(1:4,oncho_high_effect_resp[16,1:4],type="b",lwd=2,col=5,lty=2)
axis(1,c(1:4),c("lm","t","cc","c"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,"Analysis",cex=1,line=2)#,line=2.2)

plot(1:4,oncho_high_effect_resp[9,1:4],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
#lines(1:4,oncho_high_effect_resp[9,4:6],type="b",lwd=2,col=2)
lines(1:4,oncho_high_effect_resp[13,1:4],type="b",lwd=2,col=3,lty=2)
lines(1:4,oncho_high_effect_resp[17,1:4],type="b",lwd=2,col=5,lty=2)
axis(1,c(1:4),c("lm","t","cc","c"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,"Analysis",cex=1,line=2)#,line=2.2)

plot(1:4,oncho_high_effect_resp[10,1:4],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
#lines(1:4,oncho_high_effect_resp[10,4:6],type="b",lwd=2,col=2)
lines(1:4,oncho_high_effect_resp[14,1:4],type="b",lwd=2,col=3,lty=2)
lines(1:4,oncho_high_effect_resp[18,1:4],type="b",lwd=2,col=5,lty=2)
mtext(side=4,"Scenario all dose eff.",cex=1,line=1)#,line=2.2)
axis(1,c(1:4),c("lm","t","cc","c"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,"Analysis",cex=1,line=2)#,line=2.2)

dev.off()



#Selection
pdf(file ="Oncho_selection_resp.pdf", width = 9, height = 12, pointsize = 12, paper = "special")

par (mfrow=c(4,3), mar=c(1, 1, 1, 1), oma=c(3, 3, 1, 1))
plot(1:4,oncho_no_effect_resp[1,1:4],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1,cex.axis=1.3)
#legend("top",legend=c(expression(N[1]==80),expression(N[1]==120)),lwd=2,col=c(1,2,3,4),lty=c(1,1,2,2),bty="n")
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
mtext(side=3,"Low Dose",cex=1.3)#,line=2.2)
#legend("topright",legend=c(expression(N[1]==60),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3),expression(N[1]==90),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3)),cex=.8,lwd=2,col=c("white",1:4,"white",4:6),lwd=2.5,ncol=2,bty="n",pch=c(1,1:4,1,4:6))
#lines(1:4,oncho_no_effect_resp[1,4:6],type="b",lwd=2,col=2)
legend("top",legend=c("Adaptive design","MA1","MA2"),lwd=2,col=c(1,3,4),lty=c(1,2,2),bty="n",cex=1.3)

plot(1:4,oncho_no_effect_resp[2,1:4],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
mtext(side=3,"Medium Dose",cex=1.3)#,line=2.2)
#lines(1:4,oncho_no_effect_resp[2,4:6],type="b",lwd=2,col=2)

plot(1:4,oncho_no_effect_resp[3,1:4],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
mtext(side=3,"High Dose",cex=1.3)#,line=2.2)
#lines(1:4,oncho_no_effect_resp[3,4:6],type="b",lwd=2,col=2)

mtext(side=4,"Global null hypothesis",cex=1,line=1)#,line=2.2)



plot(1:4,oncho_low_effect_resp[1,1:4],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1,cex.axis=1.3)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
#lines(1:4,oncho_low_effect_resp[1,4:6],type="b",lty=1,lwd=2,col=2)

plot(1:4,oncho_low_effect_resp[2,1:4],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
#lines(1:4,oncho_low_effect_resp[2,4:6],type="b",lwd=2,col=2)

plot(1:4,oncho_low_effect_resp[3,1:4],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
#lines(1:4,oncho_low_effect_resp[3,4:6],type="b",lwd=2,col=2)
mtext(side=4,"Scenario high dose eff.",cex=1,line=1)#,line=2.2)


plot(1:4,oncho_med_effect_resp[1,1:4],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1,cex.axis=1.3)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
#lines(1:4,oncho_med_effect_resp[1,4:6],type="b",lty=1,lwd=2,col=2)

plot(1:4,oncho_med_effect_resp[2,1:4],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
#lines(1:4,oncho_med_effect_resp[2,4:6],type="b",lwd=2,col=2)

plot(1:4,oncho_med_effect_resp[3,1:4],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
#lines(1:4,oncho_med_effect_resp[3,4:6],type="b",lwd=2,col=2)
mtext(side=4,"Scenario Trend",cex=1,line=1)#,line=2.2)

plot(1:4,oncho_high_effect_resp[1,1:4],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1,cex.axis=1.3)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
#lines(1:4,oncho_high_effect_resp[1,4:6],type="b",lty=1,lwd=2,col=2)
axis(1,c(1:4),c("lm","t","cc","c"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,"Analysis",cex=1,line=2)#,line=2.2)

plot(1:4,oncho_high_effect_resp[2,1:4],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
#lines(1:4,oncho_high_effect_resp[2,4:6],type="b",lwd=2,col=2)
axis(1,c(1:4),c("lm","t","cc","c"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,"Analysis",cex=1,line=2)#,line=2.2)

plot(1:4,oncho_high_effect_resp[3,1:4],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
#lines(1:4,oncho_high_effect_resp[3,4:6],type="b",lwd=2,col=2)
axis(1,c(1:4),c("lm","t","cc","c"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,"Analysis",cex=1,line=2)#,line=2.2)
mtext(side=4,"Scenario all dose eff.",cex=1,line=1)#,line=2.2)
dev.off()

#Conditional Power


pdf(file ="Oncho_condpower_resp.pdf", width = 9, height = 9, pointsize = 12, paper = "special")
par (mfrow=c(3,3), mar=c(1, 1, 1, 1), oma=c(3, 3, 1, 1))
#plot(1:4,oncho_no_effect_resp[4,1:4],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1)
#mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
#mtext(side=3,"Low Dose",cex=1.3)#,line=2.2)
##legend("topright",legend=c(expression(N[1]==60),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3),expression(N[1]==90),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3)),cex=.8,lwd=2,col=c("white",1:4,"white",4:6),lwd=2.5,ncol=2,bty="n",pch=c(1,1:4,1,4:6))
#lines(1:4,oncho_no_effect_resp[4,4:6],type="b",lwd=2,col=2)

#plot(1:4,oncho_no_effect_resp[5,1:4],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
#mtext(side=3,"Medium Dose",cex=1.3)#,line=2.2)
#lines(1:4,oncho_no_effect_resp[5,4:6],type="b",lwd=2,col=2)

#plot(1:4,oncho_no_effect_resp[6,1:4],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
#mtext(side=3,"High Dose",cex=1.3)#,line=2.2)
#lines(1:4,oncho_no_effect_resp[6,4:6],type="b",lwd=2,col=2)
#mtext(side=4,"Global null hypothesis",cex=1,line=1)#,line=2.2)



plot(1:4,oncho_low_effect_resp[4,1:4],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1,cex.axis=1.3)
#legend("top",legend=c(expression(N[1]==80),expression(N[1]==120)),lwd=2,col=c(1,2,3,4),lty=c(1,1,2,2),bty="n")
mtext(side=3,"Low Dose",cex=1.3)#,line=2.2)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
#lines(1:4,oncho_low_effect_resp[4,4:6],type="b",lty=1,lwd=2,col=2)
legend("top",legend=c("Adaptive design","MA1","MA2"),lwd=2,col=c(1,3,4),lty=c(1,2,2),bty="n",cex=1.3)

plot(1:4,oncho_low_effect_resp[5,1:4],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
mtext(side=3,"Medium Dose",cex=1.3)#,line=2.2)
#lines(1:4,oncho_low_effect_resp[5,4:6],type="b",lwd=2,col=2)

plot(1:4,oncho_low_effect_resp[6,1:4],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
mtext(side=3,"High Dose",cex=1.3)#,line=2.2)
#lines(1:4,oncho_low_effect_resp[6,4:6],type="b",lwd=2,col=2)
mtext(side=4,"Scenario high dose eff.",cex=1,line=1)#,line=2.2)


plot(1:4,oncho_med_effect_resp[4,1:4],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1,cex.axis=1.3)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
#lines(1:4,oncho_med_effect_resp[4,4:6],type="b",lty=1,lwd=2,col=2)

plot(1:4,oncho_med_effect_resp[5,1:4],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
#lines(1:4,oncho_med_effect_resp[5,4:6],type="b",lwd=2,col=2)

plot(1:4,oncho_med_effect_resp[6,1:4],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
#lines(1:4,oncho_med_effect_resp[6,4:6],type="b",lwd=2,col=2)
mtext(side=4,"Scenario Trend",cex=1,line=1)#,line=2.2)

plot(1:4,oncho_high_effect_resp[4,1:4],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1,cex.axis=1.3)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
#lines(1:4,oncho_high_effect_resp[4,4:6],type="b",lty=1,lwd=2,col=2)
axis(1,c(1:4),c("lm","t","cc","c"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,"Analysis",cex=1,line=2)#,line=2.2)

plot(1:4,oncho_high_effect_resp[5,1:4],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
#lines(1:4,oncho_high_effect_resp[5,4:6],type="b",lwd=2,col=2)
axis(1,c(1:4),c("lm","t","cc","c"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,"Analysis",cex=1,line=2)#,line=2.2)

plot(1:4,oncho_high_effect_resp[6,1:4],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
#lines(1:4,oncho_high_effect_resp[6,4:6],type="b",lwd=2,col=2)
axis(1,c(1:4),c("lm","t","cc","c"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,"Analysis",cex=1,line=2)#,line=2.2)

dev.off()




#Error


pdf(file ="Oncho_Error_resp.pdf", width = 12, height = 3, pointsize = 12, paper = "special")

par (mfrow=c(1,4), mar=c(1, 1, 1, 1), oma=c(3, 3, 1, 1))



plot(1:4,oncho_no_effect_resp[7,1:4],type="b",xaxt="n",ylim=c(0,0.05),lty=1,lwd=2,col=1,cex.axis=1.3)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
##legend("topright",legend=c(expression(N[1]==60),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3),expression(N[1]==90),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3)),cex=.8,lwd=2,col=c("white",1:4,"white",4:6),lwd=2.5,ncol=2,bty="n",pch=c(1,1:4,1,4:6))
#lines(1:4,oncho_no_effect_resp[7,4:6],type="b",lwd=2,col=2)
lines(1:4,oncho_no_effect_resp[11,1:4],type="b",lwd=2,col=3,lty=2)
lines(1:4,oncho_no_effect_resp[15,1:4],type="b",lwd=2,col=5,lty=2)
abline(h=0.025)
axis(1,c(1:4),c("lm","t","cc","c"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,"Analysis",cex=1,line=2)#,line=2.2)
mtext(side=3,"Low Dose",cex=1.3)#,line=2.2)
legend("top",legend=c("Adaptive design","MA1","MA2"),lwd=2,col=c(1,3,4),lty=c(1,2,2),bty="n",cex=1.3)

plot(1:4,oncho_no_effect_resp[8,1:4],type="b",xaxt="n",yaxt="n",ylim=c(0,0.05),lwd=2,col=1)
#lines(1:4,oncho_no_effect_resp[8,4:6],type="b",lwd=2,col=2)
lines(1:4,oncho_no_effect_resp[12,1:4],type="b",lwd=2,col=3,lty=2)
lines(1:4,oncho_no_effect_resp[16,1:4],type="b",lwd=2,col=5,lty=2)#
abline(h=0.025)
axis(1,c(1:4),c("lm","t","cc","c"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,"Analysis",cex=1,line=2)#,line=2.2)
mtext(side=3,"Medium dose",cex=1.3)#,line=2.2)


plot(1:4,oncho_no_effect_resp[9,1:4],type="b",xaxt="n",yaxt="n",ylim=c(0,0.05),lwd=2,col=1)
#lines(1:4,oncho_no_effect_resp[9,4:6],type="b",lwd=2,col=2)
lines(1:4,oncho_no_effect_resp[13,1:4],type="b",lwd=2,col=3,lty=2)
lines(1:4,oncho_no_effect_resp[17,1:4],type="b",lwd=2,col=5,lty=2)
abline(h=0.025)
axis(1,c(1:4),c("lm","t","cc","c"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,"Analysis",cex=1,line=2)#,line=2.2)
mtext(side=3,"High dose",cex=1.3)#,line=2.2)

plot(1:4,oncho_no_effect_resp[10,1:4],type="b",xaxt="n",yaxt="n",ylim=c(0,0.05),lwd=2,col=1)
#lines(1:4,oncho_no_effect_resp[10,4:6],type="b",lwd=2,col=2)
lines(1:4,oncho_no_effect_resp[14,1:4],type="b",lwd=2,col=3,lty=2)
lines(1:4,oncho_no_effect_resp[18,1:4],type="b",lwd=2,col=5,lty=2)
mtext(side=4,"Global null hypothesis",cex=1,line=1)#,line=2.2)
abline(h=0.025)
axis(1,c(1:4),c("lm","t","cc","c"),padj=-0.3,cex.axis=1.3)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1.3)
mtext(side=1,"Analysis",cex=1,line=2)#,line=2.2)
mtext(side=3,"Total Study FWER",cex=1.3)#,line=2.2)

dev.off()