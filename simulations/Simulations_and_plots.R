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

pdf(file ="Scenarios.pdf", width = 8, height = 8, pointsize = 12, paper = "special")
plot(1:3,c(r1_6lo,r2_6lo,r3_6lo),type="b",lty=2,ylim=c(0,1),ylab="reduction rate",xlab="Doses",lwd=2)
legend("top",legend=c("","12months","6months","low","","","median","","","high","",""),lwd=2,col=c(0,0,1,0,0,1,2,3,0,1,2,3),bty="n",lty=c(1,1,1,1,1,1,1,1,1,2,2,2),ncol=4)
lines(1:3,c(r1_12lo,r2_12lo,r3_12lo),type="b",lty=1,lwd=2)
lines(1:3,c(r1_12me,r2_12me,r3_12me),type="b",lty=1,lwd=2,col=2)
lines(1:3,c(r1_6me,r2_6me,r3_6me),type="b",lty=2,lwd=2,col=2)            
lines(1:3,c(r1_12hi,r2_12hi,r3_12hi),type="b",lty=1,lwd=2,col=3)
lines(1:3,c(r1_6hi,r2_6hi,r3_6hi),type="b",lty=2,lwd=2,col=3)            
dev.off()


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
test1=0

#no_effect

#simulation.global.null<-c(mu_raw_0, sd_raw_0, r0_6,r1_6,r2_6,r3_6, r0_12,r1_12,r2_12,r3_12,  rho ,   #data.matrix(expand.grid
#            n_trials,n_arms,N1 , N, rmonth, alpha1 , alpha,  sim_out1,sel_scen, side1, test1)#)

oncho_no_effect<-mapply(simul_res,mu_raw_0, sd_raw_0 , r0_6no,r1_6no,r2_6no,r3_6no, r0_12no,r1_12no,r2_12no,r3_12no,  rho ,
       n_trials,4,N1 , N, rmonth, alpha1 , alpha,
       sim_out1,sel_scen, side1,test1)



#low_effect
oncho_low_effect<-mapply(simul_res,mu_raw_0, sd_raw_0 , r0_6lo,r1_6lo,r2_6lo,r3_6lo, r0_12lo,r1_12lo,r2_12lo,r3_12lo,  rho ,
                        n_trials,4,N1 , N, rmonth, alpha1 , alpha,
                        sim_out1,sel_scen, side1,test1)


#med_effect

oncho_med_effect<-mapply(simul_res,mu_raw_0, sd_raw_0 , r0_6me,r1_6me,r2_6me,r3_6me, r0_12me,r1_12me,r2_12me,r3_12me,  rho ,
                         n_trials,4,N1 , N, rmonth, alpha1 , alpha,
                         sim_out1,sel_scen, side1,test1)



#high_effect

oncho_high_effect<-mapply(simul_res,mu_raw_0, sd_raw_0 , r0_6hi,r1_6hi,r2_6hi,r3_6hi, r0_12hi,r1_12hi,r2_12hi,r3_12hi,  rho ,
                         n_trials,4,N1 , N, rmonth, alpha1 , alpha,
                         sim_out1,sel_scen, side1,test1)




#Power
pdf(file ="Oncho_Pow.pdf", width = 12, height = 12, pointsize = 12, paper = "special")

par (mfrow=c(4,4), mar=c(1, 1, 1, 1), oma=c(3, 3, 1, 1))
plot(1:3,oncho_no_effect[7,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
mtext(side=3,"Low Dose",cex=1.3)#,line=2.2)
#legend("topright",legend=c(expression(N[1]==60),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3),expression(N[1]==90),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3)),cex=.8,lwd=2,col=c("white",1:3,"white",4:6),lwd=2.5,ncol=2,bty="n",pch=c(1,1:3,1,4:6))
lines(1:3,oncho_no_effect[7,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_no_effect[11,1:3],type="b",lwd=2,col=3,lty=2)
#lines(1:3,oncho_no_effect[11,4:6],type="b",lwd=2,col=4) #sollte das selbe sein wie in voriger Zeile, daher weglassen
lines(1:3,oncho_no_effect[15,1:3],type="b",lwd=2,col=5,lty=2)

plot(1:3,oncho_no_effect[8,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
mtext(side=3,"Intermediate Dose",cex=1.3)#,line=2.2)
lines(1:3,oncho_no_effect[8,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_no_effect[12,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_no_effect[16,1:3],type="b",lwd=2,col=5,lty=2)

plot(1:3,oncho_no_effect[9,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
mtext(side=3,"High Dose",cex=1.3)#,line=2.2)
lines(1:3,oncho_no_effect[9,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_no_effect[13,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_no_effect[17,1:3],type="b",lwd=2,col=5,lty=2)

plot(1:3,oncho_no_effect[10,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
mtext(side=3,"Disjunctive",cex=1.3)#,line=2.2)
lines(1:3,oncho_no_effect[10,4:6],type="p",lwd=2,col=2)
lines(1:3,oncho_no_effect[14,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_no_effect[18,1:3],type="b",lwd=2,col=5,lty=2)
mtext(side=4,"Scenario 0",cex=1,line=1)#,line=2.2)



plot(1:3,oncho_low_effect[7,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
lines(1:3,oncho_low_effect[7,4:6],type="b",lty=1,lwd=2,col=2)
lines(1:3,oncho_low_effect[11,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_low_effect[15,1:3],type="b",lwd=2,col=5,lty=2)

plot(1:3,oncho_low_effect[8,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_low_effect[8,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_low_effect[12,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_low_effect[16,1:3],type="b",lwd=2,col=5,lty=2)

plot(1:3,oncho_low_effect[9,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_low_effect[9,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_low_effect[13,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_low_effect[17,1:3],type="b",lwd=2,col=5,lty=2)

plot(1:3,oncho_low_effect[10,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_low_effect[10,4:6],type="p",lwd=2,col=2)
lines(1:3,oncho_low_effect[14,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_low_effect[18,1:3],type="b",lwd=2,col=5,lty=2)
mtext(side=4,"Scenario low",cex=1,line=1)#,line=2.2)


plot(1:3,oncho_med_effect[7,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1)
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
lines(1:3,oncho_med_effect[10,4:6],type="p",lwd=2,col=2)
lines(1:3,oncho_med_effect[14,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_med_effect[18,1:3],type="b",lwd=2,col=5,lty=2)
mtext(side=4,"Scenario med",cex=1,line=1)#,line=2.2)

plot(1:3,oncho_high_effect[7,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
lines(1:3,oncho_high_effect[7,4:6],type="b",lty=1,lwd=2,col=2)
lines(1:3,oncho_high_effect[11,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_high_effect[15,1:3],type="b",lwd=2,col=5,lty=2)
axis(1,c(1:3),c("0.1","0.2","0.3"),padj=-0.3,cex.axis=1)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1)
mtext(side=1,expression(alpha[1]),cex=1,line=2)#,line=2.2)

plot(1:3,oncho_high_effect[8,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_high_effect[8,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_high_effect[12,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_high_effect[16,1:3],type="b",lwd=2,col=5,lty=2)
axis(1,c(1:3),c("0.1","0.2","0.3"),padj=-0.3,cex.axis=1)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1)
mtext(side=1,expression(alpha[1]),cex=1,line=2)#,line=2.2)

plot(1:3,oncho_high_effect[9,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_high_effect[9,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_high_effect[13,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_high_effect[17,1:3],type="b",lwd=2,col=5,lty=2)
axis(1,c(1:3),c("0.1","0.2","0.3"),padj=-0.3,cex.axis=1)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1)
mtext(side=1,expression(alpha[1]),cex=1,line=2)#,line=2.2)

plot(1:3,oncho_high_effect[10,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_high_effect[10,4:6],type="p",lwd=2,col=2)
lines(1:3,oncho_high_effect[14,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_high_effect[18,1:3],type="b",lwd=2,col=5,lty=2)
mtext(side=4,"Scenario high",cex=1,line=1)#,line=2.2)
axis(1,c(1:3),c("0.1","0.2","0.3"),padj=-0.3,cex.axis=1)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1)
mtext(side=1,expression(alpha[1]),cex=1,line=2)#,line=2.2)

dev.off()



#Selection
pdf(file ="Oncho_selection.pdf", width = 9, height = 12, pointsize = 12, paper = "special")

par (mfrow=c(4,3), mar=c(1, 1, 1, 1), oma=c(3, 3, 1, 1))
plot(1:3,oncho_no_effect[1,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
mtext(side=3,"Low Dose",cex=1.3)#,line=2.2)
#legend("topright",legend=c(expression(N[1]==60),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3),expression(N[1]==90),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3)),cex=.8,lwd=2,col=c("white",1:3,"white",4:6),lwd=2.5,ncol=2,bty="n",pch=c(1,1:3,1,4:6))
lines(1:3,oncho_no_effect[1,4:6],type="b",lwd=2,col=2)

plot(1:3,oncho_no_effect[2,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
mtext(side=3,"Intermediate Dose",cex=1.3)#,line=2.2)
lines(1:3,oncho_no_effect[2,4:6],type="b",lwd=2,col=2)

plot(1:3,oncho_no_effect[3,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
mtext(side=3,"High Dose",cex=1.3)#,line=2.2)
lines(1:3,oncho_no_effect[3,4:6],type="b",lwd=2,col=2)

mtext(side=4,"Scenario 0",cex=1,line=1)#,line=2.2)



plot(1:3,oncho_low_effect[1,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
lines(1:3,oncho_low_effect[1,4:6],type="b",lty=1,lwd=2,col=2)

plot(1:3,oncho_low_effect[2,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_low_effect[2,4:6],type="b",lwd=2,col=2)

plot(1:3,oncho_low_effect[3,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_low_effect[3,4:6],type="b",lwd=2,col=2)
mtext(side=4,"Scenario low",cex=1,line=1)#,line=2.2)


plot(1:3,oncho_med_effect[1,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
lines(1:3,oncho_med_effect[1,4:6],type="b",lty=1,lwd=2,col=2)

plot(1:3,oncho_med_effect[2,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_med_effect[2,4:6],type="b",lwd=2,col=2)

plot(1:3,oncho_med_effect[3,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_med_effect[3,4:6],type="b",lwd=2,col=2)
mtext(side=4,"Scenario med",cex=1,line=1)#,line=2.2)

plot(1:3,oncho_high_effect[1,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
lines(1:3,oncho_high_effect[1,4:6],type="b",lty=1,lwd=2,col=2)
axis(1,c(1:3),c("0.1","0.2","0.3"),padj=-0.3,cex.axis=1)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1)
mtext(side=1,expression(alpha[1]),cex=1,line=2)#,line=2.2)

plot(1:3,oncho_high_effect[2,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_high_effect[2,4:6],type="b",lwd=2,col=2)
axis(1,c(1:3),c("0.1","0.2","0.3"),padj=-0.3,cex.axis=1)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1)
mtext(side=1,expression(alpha[1]),cex=1,line=2)#,line=2.2)

plot(1:3,oncho_high_effect[3,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_high_effect[3,4:6],type="b",lwd=2,col=2)
axis(1,c(1:3),c("0.1","0.2","0.3"),padj=-0.3,cex.axis=1)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1)
mtext(side=1,expression(alpha[1]),cex=1,line=2)#,line=2.2)
mtext(side=4,"Scenario high",cex=1,line=1)#,line=2.2)
dev.off()

#Conditional Power


pdf(file ="Oncho_condpower.pdf", width = 9, height = 12, pointsize = 12, paper = "special")
par (mfrow=c(4,3), mar=c(1, 1, 1, 1), oma=c(3, 3, 1, 1))
plot(1:3,oncho_no_effect[4,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
mtext(side=3,"Low Dose",cex=1.3)#,line=2.2)
#legend("topright",legend=c(expression(N[1]==60),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3),expression(N[1]==90),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3)),cex=.8,lwd=2,col=c("white",1:3,"white",4:6),lwd=2.5,ncol=2,bty="n",pch=c(1,1:3,1,4:6))
lines(1:3,oncho_no_effect[4,4:6],type="b",lwd=2,col=2)

plot(1:3,oncho_no_effect[5,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
mtext(side=3,"Intermediate Dose",cex=1.3)#,line=2.2)
lines(1:3,oncho_no_effect[5,4:6],type="b",lwd=2,col=2)

plot(1:3,oncho_no_effect[6,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
mtext(side=3,"High Dose",cex=1.3)#,line=2.2)
lines(1:3,oncho_no_effect[6,4:6],type="b",lwd=2,col=2)

mtext(side=4,"Scenario 0",cex=1,line=1)#,line=2.2)



plot(1:3,oncho_low_effect[4,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
lines(1:3,oncho_low_effect[4,4:6],type="b",lty=1,lwd=2,col=2)

plot(1:3,oncho_low_effect[5,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_low_effect[5,4:6],type="b",lwd=2,col=2)

plot(1:3,oncho_low_effect[6,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_low_effect[6,4:6],type="b",lwd=2,col=2)
mtext(side=4,"Scenario low",cex=1,line=1)#,line=2.2)


plot(1:3,oncho_med_effect[4,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
lines(1:3,oncho_med_effect[4,4:6],type="b",lty=1,lwd=2,col=2)

plot(1:3,oncho_med_effect[5,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_med_effect[5,4:6],type="b",lwd=2,col=2)

plot(1:3,oncho_med_effect[6,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_med_effect[6,4:6],type="b",lwd=2,col=2)
mtext(side=4,"Scenario med",cex=1,line=1)#,line=2.2)

plot(1:3,oncho_high_effect[4,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
lines(1:3,oncho_high_effect[4,4:6],type="b",lty=1,lwd=2,col=2)
axis(1,c(1:3),c("0.1","0.2","0.3"),padj=-0.3,cex.axis=1)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1)
mtext(side=1,expression(alpha[1]),cex=1,line=2)#,line=2.2)

plot(1:3,oncho_high_effect[5,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_high_effect[5,4:6],type="b",lwd=2,col=2)
axis(1,c(1:3),c("0.1","0.2","0.3"),padj=-0.3,cex.axis=1)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1)
mtext(side=1,expression(alpha[1]),cex=1,line=2)#,line=2.2)

plot(1:3,oncho_high_effect[6,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_high_effect[6,4:6],type="b",lwd=2,col=2)
axis(1,c(1:3),c("0.1","0.2","0.3"),padj=-0.3,cex.axis=1)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1)
mtext(side=1,expression(alpha[1]),cex=1,line=2)#,line=2.2)
mtext(side=4,"Scenario high",cex=1,line=1)#,line=2.2)









###############################################################################################
###############################################################################################

#
#
#     Compare Wilcoxon linear model
#
#


N1<-90
alpha1<-.1
sel_scen=c(0,1)
test1=c(0,1,3)

#no_effect

oncho_no_effect_test<-mapply(simul_res,mu_raw_0, sd_raw_0 , r0_6no,r1_6no,r2_6no,r3_6no, r0_12no,r1_12no,r2_12no,r3_12no,  rho ,
                        n_trials,4,N1 , N, rmonth, alpha1 , alpha,
                        sim_out1,sel_scen, side1,test1)

#low_effect
oncho_low_effect_test<-mapply(simul_res,mu_raw_0, sd_raw_0 , r0_6lo,r1_6lo,r2_6lo,r3_6lo, r0_12lo,r1_12lo,r2_12lo,r3_12lo,  rho ,
                         n_trials,4,N1 , N, rmonth, alpha1 , alpha,
                         sim_out1,sel_scen, side1,test1)

#med_effect

oncho_med_effect_test<-mapply(simul_res,mu_raw_0, sd_raw_0 , r0_6me,r1_6me,r2_6me,r3_6me, r0_12me,r1_12me,r2_12me,r3_12me,  rho ,
                         n_trials,4,N1 , N, rmonth, alpha1 , alpha,
                         sim_out1,sel_scen, side1,test1)

#high_effect

oncho_high_effect_test<-mapply(simul_res,mu_raw_0, sd_raw_0 , r0_6hi,r1_6hi,r2_6hi,r3_6hi, r0_12hi,r1_12hi,r2_12hi,r3_12hi,  rho ,
                          n_trials,4,N1 , N, rmonth, alpha1 , alpha,
                          sim_out1,sel_scen, side1,test1)







#Power
pdf(file ="Oncho_Pow_test.pdf", width = 12, height = 12, pointsize = 12, paper = "special")

par (mfrow=c(4,4), mar=c(1, 1, 1, 1), oma=c(3, 3, 1, 1))
plot(1:3,oncho_no_effect_test[7,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
mtext(side=3,"Low Dose",cex=1.3)#,line=2.2)
#legend("topright",legend=c(expression(N[1]==60),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3),expression(N[1]==90),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3)),cex=.8,lwd=2,col=c("white",1:3,"white",4:6),lwd=2.5,ncol=2,bty="n",pch=c(1,1:3,1,4:6))
lines(1:3,oncho_no_effect_test[7,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_no_effect_test[11,1:3],type="b",lwd=2,col=3,lty=2)
#lines(1:3,oncho_no_effect_test[11,4:6],type="b",lwd=2,col=4) #sollte das selbe sein wie in voriger Zeile, daher weglassen
lines(1:3,oncho_no_effect_test[15,1:3],type="b",lwd=2,col=5,lty=2)

plot(1:3,oncho_no_effect_test[8,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
mtext(side=3,"Intermediate Dose",cex=1.3)#,line=2.2)
lines(1:3,oncho_no_effect_test[8,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_no_effect_test[12,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_no_effect_test[16,1:3],type="b",lwd=2,col=5,lty=2)

plot(1:3,oncho_no_effect_test[9,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
mtext(side=3,"High Dose",cex=1.3)#,line=2.2)
lines(1:3,oncho_no_effect_test[9,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_no_effect_test[13,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_no_effect_test[17,1:3],type="b",lwd=2,col=5,lty=2)

plot(1:3,oncho_no_effect_test[10,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
mtext(side=3,"Disjunctive",cex=1.3)#,line=2.2)
lines(1:3,oncho_no_effect_test[10,4:6],type="p",lwd=2,col=2)
lines(1:3,oncho_no_effect_test[14,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_no_effect_test[18,1:3],type="b",lwd=2,col=5,lty=2)
mtext(side=4,"Scenario 0",cex=1,line=1)#,line=2.2)



plot(1:3,oncho_low_effect_test[7,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
lines(1:3,oncho_low_effect_test[7,4:6],type="b",lty=1,lwd=2,col=2)
lines(1:3,oncho_low_effect_test[11,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_low_effect_test[15,1:3],type="b",lwd=2,col=5,lty=2)

plot(1:3,oncho_low_effect_test[8,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_low_effect_test[8,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_low_effect_test[12,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_low_effect_test[16,1:3],type="b",lwd=2,col=5,lty=2)

plot(1:3,oncho_low_effect_test[9,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_low_effect_test[9,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_low_effect_test[13,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_low_effect_test[17,1:3],type="b",lwd=2,col=5,lty=2)

plot(1:3,oncho_low_effect_test[10,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_low_effect_test[10,4:6],type="p",lwd=2,col=2)
lines(1:3,oncho_low_effect_test[14,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_low_effect_test[18,1:3],type="b",lwd=2,col=5,lty=2)
mtext(side=4,"Scenario low",cex=1,line=1)#,line=2.2)


plot(1:3,oncho_med_effect_test[7,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1)
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
lines(1:3,oncho_med_effect_test[10,4:6],type="p",lwd=2,col=2)
lines(1:3,oncho_med_effect_test[14,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_med_effect_test[18,1:3],type="b",lwd=2,col=5,lty=2)
mtext(side=4,"Scenario med",cex=1,line=1)#,line=2.2)

plot(1:3,oncho_high_effect_test[7,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
lines(1:3,oncho_high_effect_test[7,4:6],type="b",lty=1,lwd=2,col=2)
lines(1:3,oncho_high_effect_test[11,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_high_effect_test[15,1:3],type="b",lwd=2,col=5,lty=2)
axis(1,c(1:3),c("0.1","0.2","0.3"),padj=-0.3,cex.axis=1)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1)
mtext(side=1,expression(alpha[1]),cex=1,line=2)#,line=2.2)

plot(1:3,oncho_high_effect_test[8,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_high_effect_test[8,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_high_effect_test[12,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_high_effect_test[16,1:3],type="b",lwd=2,col=5,lty=2)
axis(1,c(1:3),c("0.1","0.2","0.3"),padj=-0.3,cex.axis=1)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1)
mtext(side=1,expression(alpha[1]),cex=1,line=2)#,line=2.2)

plot(1:3,oncho_high_effect_test[9,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_high_effect_test[9,4:6],type="b",lwd=2,col=2)
lines(1:3,oncho_high_effect_test[13,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_high_effect_test[17,1:3],type="b",lwd=2,col=5,lty=2)
axis(1,c(1:3),c("0.1","0.2","0.3"),padj=-0.3,cex.axis=1)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1)
mtext(side=1,expression(alpha[1]),cex=1,line=2)#,line=2.2)

plot(1:3,oncho_high_effect_test[10,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_high_effect_test[10,4:6],type="p",lwd=2,col=2)
lines(1:3,oncho_high_effect_test[14,1:3],type="b",lwd=2,col=3,lty=2)
lines(1:3,oncho_high_effect_test[18,1:3],type="b",lwd=2,col=5,lty=2)
mtext(side=4,"Scenario high",cex=1,line=1)#,line=2.2)
axis(1,c(1:3),c("0.1","0.2","0.3"),padj=-0.3,cex.axis=1)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1)
mtext(side=1,expression(alpha[1]),cex=1,line=2)#,line=2.2)

dev.off()



#Selection
pdf(file ="Oncho_selection_test.pdf", width = 9, height = 12, pointsize = 12, paper = "special")

par (mfrow=c(4,3), mar=c(1, 1, 1, 1), oma=c(3, 3, 1, 1))
plot(1:3,oncho_no_effect_test[1,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
mtext(side=3,"Low Dose",cex=1.3)#,line=2.2)
#legend("topright",legend=c(expression(N[1]==60),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3),expression(N[1]==90),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3)),cex=.8,lwd=2,col=c("white",1:3,"white",4:6),lwd=2.5,ncol=2,bty="n",pch=c(1,1:3,1,4:6))
lines(1:3,oncho_no_effect_test[1,4:6],type="b",lwd=2,col=2)

plot(1:3,oncho_no_effect_test[2,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
mtext(side=3,"Intermediate Dose",cex=1.3)#,line=2.2)
lines(1:3,oncho_no_effect_test[2,4:6],type="b",lwd=2,col=2)

plot(1:3,oncho_no_effect_test[3,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
mtext(side=3,"High Dose",cex=1.3)#,line=2.2)
lines(1:3,oncho_no_effect_test[3,4:6],type="b",lwd=2,col=2)

mtext(side=4,"Scenario 0",cex=1,line=1)#,line=2.2)



plot(1:3,oncho_low_effect_test[1,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
lines(1:3,oncho_low_effect_test[1,4:6],type="b",lty=1,lwd=2,col=2)

plot(1:3,oncho_low_effect_test[2,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_low_effect_test[2,4:6],type="b",lwd=2,col=2)

plot(1:3,oncho_low_effect_test[3,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_low_effect_test[3,4:6],type="b",lwd=2,col=2)
mtext(side=4,"Scenario low",cex=1,line=1)#,line=2.2)


plot(1:3,oncho_med_effect_test[1,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
lines(1:3,oncho_med_effect_test[1,4:6],type="b",lty=1,lwd=2,col=2)

plot(1:3,oncho_med_effect_test[2,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_med_effect_test[2,4:6],type="b",lwd=2,col=2)

plot(1:3,oncho_med_effect_test[3,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_med_effect_test[3,4:6],type="b",lwd=2,col=2)
mtext(side=4,"Scenario med",cex=1,line=1)#,line=2.2)

plot(1:3,oncho_high_effect_test[1,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
lines(1:3,oncho_high_effect_test[1,4:6],type="b",lty=1,lwd=2,col=2)
axis(1,c(1:3),c("0.1","0.2","0.3"),padj=-0.3,cex.axis=1)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1)
mtext(side=1,expression(alpha[1]),cex=1,line=2)#,line=2.2)

plot(1:3,oncho_high_effect_test[2,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_high_effect_test[2,4:6],type="b",lwd=2,col=2)
axis(1,c(1:3),c("0.1","0.2","0.3"),padj=-0.3,cex.axis=1)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1)
mtext(side=1,expression(alpha[1]),cex=1,line=2)#,line=2.2)

plot(1:3,oncho_high_effect_test[3,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_high_effect_test[3,4:6],type="b",lwd=2,col=2)
axis(1,c(1:3),c("0.1","0.2","0.3"),padj=-0.3,cex.axis=1)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1)
mtext(side=1,expression(alpha[1]),cex=1,line=2)#,line=2.2)
mtext(side=4,"Scenario high",cex=1,line=1)#,line=2.2)
dev.off()

#Conditional Power


pdf(file ="Oncho_condpower_test.pdf", width = 9, height = 12, pointsize = 12, paper = "special")
par (mfrow=c(4,3), mar=c(1, 1, 1, 1), oma=c(3, 3, 1, 1))
plot(1:3,oncho_no_effect_test[4,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
mtext(side=3,"Low Dose",cex=1.3)#,line=2.2)
#legend("topright",legend=c(expression(N[1]==60),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3),expression(N[1]==90),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3)),cex=.8,lwd=2,col=c("white",1:3,"white",4:6),lwd=2.5,ncol=2,bty="n",pch=c(1,1:3,1,4:6))
lines(1:3,oncho_no_effect_test[4,4:6],type="b",lwd=2,col=2)

plot(1:3,oncho_no_effect_test[5,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
mtext(side=3,"Intermediate Dose",cex=1.3)#,line=2.2)
lines(1:3,oncho_no_effect_test[5,4:6],type="b",lwd=2,col=2)

plot(1:3,oncho_no_effect_test[6,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
mtext(side=3,"High Dose",cex=1.3)#,line=2.2)
lines(1:3,oncho_no_effect_test[6,4:6],type="b",lwd=2,col=2)

mtext(side=4,"Scenario 0",cex=1,line=1)#,line=2.2)



plot(1:3,oncho_low_effect_test[4,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
lines(1:3,oncho_low_effect_test[4,4:6],type="b",lty=1,lwd=2,col=2)

plot(1:3,oncho_low_effect_test[5,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_low_effect_test[5,4:6],type="b",lwd=2,col=2)

plot(1:3,oncho_low_effect_test[6,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_low_effect_test[6,4:6],type="b",lwd=2,col=2)
mtext(side=4,"Scenario low",cex=1,line=1)#,line=2.2)


plot(1:3,oncho_med_effect_test[4,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
lines(1:3,oncho_med_effect_test[4,4:6],type="b",lty=1,lwd=2,col=2)

plot(1:3,oncho_med_effect_test[5,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_med_effect_test[5,4:6],type="b",lwd=2,col=2)

plot(1:3,oncho_med_effect_test[6,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_med_effect_test[6,4:6],type="b",lwd=2,col=2)
mtext(side=4,"Scenario med",cex=1,line=1)#,line=2.2)

plot(1:3,oncho_high_effect_test[4,1:3],type="b",xaxt="n",ylim=c(0,1),lty=1,lwd=2,col=1)
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
lines(1:3,oncho_high_effect_test[4,4:6],type="b",lty=1,lwd=2,col=2)
axis(1,c(1:3),c("0.1","0.2","0.3"),padj=-0.3,cex.axis=1)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1)
mtext(side=1,expression(alpha[1]),cex=1,line=2)#,line=2.2)

plot(1:3,oncho_high_effect_test[5,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_high_effect_test[5,4:6],type="b",lwd=2,col=2)
axis(1,c(1:3),c("0.1","0.2","0.3"),padj=-0.3,cex.axis=1)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1)
mtext(side=1,expression(alpha[1]),cex=1,line=2)#,line=2.2)

plot(1:3,oncho_high_effect_test[6,1:3],type="b",xaxt="n",yaxt="n",ylim=c(0,1),lwd=2,col=1)
lines(1:3,oncho_high_effect_test[6,4:6],type="b",lwd=2,col=2)
axis(1,c(1:3),c("0.1","0.2","0.3"),padj=-0.3,cex.axis=1)#,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1)
mtext(side=1,expression(alpha[1]),cex=1,line=2)#,line=2.2)
mtext(side=4,"Scenario high",cex=1,line=1)#,line=2.2)

dev.off()
















#OLD: Delete soon!

matplot(1:3,oncho_no_effect[c(1,4,7),],type="p",pch=1:6,xaxt="n",yaxt="n",ylim=c(0,1))
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
legend("topright",legend=c(expression(N[1]==60),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3),expression(N[1]==90),expression(alpha[1]==0.1),expression(alpha[1]==0.2),expression(alpha[1]==0.3)),cex=.8,lwd=2,col=c("white",1:3,"white",4:6),lwd=2.5,ncol=2,bty="n",pch=c(1,1:3,1,4:6))
mtext(side=3,"Low Dose",cex=1.3)#,line=2.2)
matplot(1:3,oncho_no_effect[c(2,5,8),],type="p",pch=1:6,xaxt="n",yaxt="n",yaxt="n",ylim=c(0,1))
mtext(side=3,"Intermediate Dose",cex=1.3)#,line=2.2)
matplot(1:3,oncho_no_effect[c(3,6,9),],type="p",pch=1:6,xaxt="n",yaxt="n",yaxt="n",ylim=c(0,1))
mtext(side=3,"High Dose",cex=1.3)#,line=2.2)
matplot(1:3,oncho_no_effect[c(10,13,17),],type="p",pch=1:6,xaxt="n",yaxt="n",yaxt="n",ylim=c(0,1))
mtext(side=3,"Disjunctive Power",cex=1.3)#,line=2.2)
mtext(side=4,"No effective",cex=1,line=1)#,line=2.2)


matplot(1:3,oncho_low_effect[c(1,4,7),],type="p",pch=1:6,xaxt="n",yaxt="n",ylim=c(0,1))
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
matplot(1:3,oncho_low_effect[c(2,5,8),],type="p",pch=1:6,xaxt="n",yaxt="n",yaxt="n",ylim=c(0,1))
matplot(1:3,oncho_low_effect[c(3,6,9),],type="p",pch=1:6,xaxt="n",yaxt="n",yaxt="n",ylim=c(0,1))
matplot(1:3,oncho_low_effect[c(10,13,17),],type="p",pch=1:6,xaxt="n",yaxt="n",yaxt="n",ylim=c(0,1))
mtext(side=4,"low effective",cex=1,line=1)#,line=2.2)

matplot(1:3,oncho_med_effect[c(1,4,7),],type="p",pch=1:6,xaxt="n",yaxt="n",ylim=c(0,1))
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
matplot(1:3,oncho_med_effect[c(2,5,8),],type="p",pch=1:6,xaxt="n",yaxt="n",yaxt="n",ylim=c(0,1))
matplot(1:3,oncho_med_effect[c(3,6,9),],type="p",pch=1:6,xaxt="n",yaxt="n",yaxt="n",ylim=c(0,1))
matplot(1:3,oncho_med_effect[c(10,13,17),],type="p",pch=1:6,xaxt="n",yaxt="n",yaxt="n",ylim=c(0,1))
mtext(side=4,"Median effective",cex=1,line=1)#,line=2.2)

matplot(1:3,oncho_high_effect[c(1,4,7),],type="p",pch=1:6,xaxt="n",yaxt="n",ylim=c(0,1))
mtext(side=2,"Proportion",cex=1,line=2.1)#,line=2.2)
axis(1,1:3,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1)
matplot(1:3,oncho_high_effect[c(2,5,8),],type="p",pch=1:6,xaxt="n",yaxt="n",yaxt="n",ylim=c(0,1))
axis(1,1:3,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1)
matplot(1:3,oncho_high_effect[c(3,6,9),],type="p",pch=1:6,xaxt="n",yaxt="n",yaxt="n",ylim=c(0,1))
axis(1,1:3,c("Selection","condPow","Pow"),padj=-0.3,cex.axis=1)
matplot(1:3,oncho_high_effect[c(10,13,17),],type="p",pch=1:6,xaxt="n",yaxt="n",yaxt="n",ylim=c(0,1))
mtext(side=4,"High effective",cex=1,line=1)
axis(1,1:3,c("Adaptive","MA1","MA2"),padj=-0.3,cex.axis=1)



###########################################
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
lines(4:6,unlist(oo.oncho1)[7:9],type="p",lwd=2,col="red")
lines(c(4,5,7),unlist(oo.oncho1)[7:9],type="p",lwd=2,col="red")
lines(4:6,unlist(oo.oncho1)[7:9],type="p",lwd=2,col="red")


title(expression(alpha[1]==0.1))
lines(unlist(oo.oncho1[,2]),ylim=c(0,1),type="p",lwd=2,col="red")
lines(unlist(oo.oncho1[,3]),ylim=c(0,1),type="p",lwd=2,col="blue")
legend("topleft",legend=c(expression(rho==0),expression(rho==.5),expression(rho==0.7)),cex=.8,lwd=2,col=c(1,"red","blue"),lwd=2.5,ncol=1,bty="n")
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
lines(unlist(oo.oncho1[,2]),ylim=c(0,1),type="p",lwd=2,col="red")
lines(unlist(oo.oncho1[,3]),ylim=c(0,1),type="p",lwd=2,col="blue")
legend("topleft",legend=c(expression(rho==0),expression(rho==.5),expression(rho==0.7)),cex=.8,lwd=2,col=c(1,"red","blue"),lwd=2.5,ncol=1,bty="n")
abline(v=3.5)

plot(unlist(oo.oncho1[,4]),ylim=c(0,1))
title(expression(alpha[1]==0.5))
lines(unlist(oo.oncho1[,5]),ylim=c(0,1),type="p",lwd=2,col="red")
lines(unlist(oo.oncho1[,6]),ylim=c(0,1),type="p",lwd=2,col="blue")
abline(v=3.5)

plot(unlist(oo.oncho2[,1]),ylim=c(0,1))
lines(unlist(oo.oncho2[,2]),ylim=c(0,1),type="p",lwd=2,col="red")
lines(unlist(oo.oncho2[,3]),ylim=c(0,1),type="p",lwd=2,col="blue")
abline(v=3.5)
plot(unlist(oo.oncho2[,4]),ylim=c(0,1))
lines(unlist(oo.oncho2[,5]),ylim=c(0,1),type="p",lwd=2,col="red")
lines(unlist(oo.oncho2[,6]),ylim=c(0,1),type="p",lwd=2,col="blue")
abline(v=3.5)

plot(unlist(oo.oncho3[,1]),ylim=c(0,1))
lines(unlist(oo.oncho3[,2]),ylim=c(0,1),type="p",lwd=2,col="red")
lines(unlist(oo.oncho3[,3]),ylim=c(0,1),type="p",lwd=2,col="blue")
abline(v=3.5)
plot(unlist(oo.oncho3[,4]),ylim=c(0,1))
lines(unlist(oo.oncho3[,5]),ylim=c(0,1),type="p",lwd=2,col="red")
lines(unlist(oo.oncho3[,6]),ylim=c(0,1),type="p",lwd=2,col="blue")
abline(v=3.5)

plot(unlist(oo.oncho4[,1]),ylim=c(0,1))
lines(unlist(oo.oncho4[,2]),ylim=c(0,1),type="p",lwd=2,col="red")
lines(unlist(oo.oncho4[,3]),ylim=c(0,1),type="p",lwd=2,col="blue")
abline(v=3.5)
plot(unlist(oo.oncho4[,4]),ylim=c(0,1))
lines(unlist(oo.oncho4[,5]),ylim=c(0,1),type="p",lwd=2,col="red")
lines(unlist(oo.oncho4[,6]),ylim=c(0,1),type="p",lwd=2,col="blue")
abline(v=3.5)



plot(unlist(oo.oncho5[,1]),ylim=c(0,1))
lines(unlist(oo.oncho5[,2]),ylim=c(0,1),type="p",lwd=2,col="red")
lines(unlist(oo.oncho5[,3]),ylim=c(0,1),type="p",lwd=2,col="blue")
abline(v=3.5)
plot(unlist(oo.oncho5[,4]),ylim=c(0,1))
lines(unlist(oo.oncho5[,5]),ylim=c(0,1),type="p",lwd=2,col="red")
lines(unlist(oo.oncho5[,6]),ylim=c(0,1),type="p",lwd=2,col="blue")
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
lines(unlist(oo.mans1[,2]),ylim=c(0,1),type="p",lwd=2,col="red")
lines(unlist(oo.mans1[,3]),ylim=c(0,1),type="p",lwd=2,col="blue")
legend("topleft",legend=c(expression(rho==0),expression(rho==.5),expression(rho==0.7)),cex=.8,lwd=2,col=c(1,"red","blue"),lwd=2.5,ncol=1,bty="n")
abline(v=3.5)

plot(unlist(oo.mans1[,4]),ylim=c(0,1))
title(expression(alpha[1]==0.5))
lines(unlist(oo.mans1[,5]),ylim=c(0,1),type="p",lwd=2,col="red")
lines(unlist(oo.mans1[,6]),ylim=c(0,1),type="p",lwd=2,col="blue")
abline(v=3.5)

plot(unlist(oo.mans2[,1]),ylim=c(0,1))
lines(unlist(oo.mans2[,2]),ylim=c(0,1),type="p",lwd=2,col="red")
lines(unlist(oo.mans2[,3]),ylim=c(0,1),type="p",lwd=2,col="blue")
abline(v=3.5)
plot(unlist(oo.mans2[,4]),ylim=c(0,1))
lines(unlist(oo.mans2[,5]),ylim=c(0,1),type="p",lwd=2,col="red")
lines(unlist(oo.mans2[,6]),ylim=c(0,1),type="p",lwd=2,col="blue")
abline(v=3.5)

plot(unlist(oo.mans3[,1]),ylim=c(0,1))
lines(unlist(oo.mans3[,2]),ylim=c(0,1),type="p",lwd=2,col="red")
lines(unlist(oo.mans3[,3]),ylim=c(0,1),type="p",lwd=2,col="blue")
abline(v=3.5)
plot(unlist(oo.mans3[,4]),ylim=c(0,1))
lines(unlist(oo.mans3[,5]),ylim=c(0,1),type="p",lwd=2,col="red")
lines(unlist(oo.mans3[,6]),ylim=c(0,1),type="p",lwd=2,col="blue")
abline(v=3.5)

plot(unlist(oo.mans4[,1]),ylim=c(0,1))
lines(unlist(oo.mans4[,2]),ylim=c(0,1),type="p",lwd=2,col="red")
lines(unlist(oo.mans4[,3]),ylim=c(0,1),type="p",lwd=2,col="blue")
abline(v=3.5)
plot(unlist(oo.mans4[,4]),ylim=c(0,1))
lines(unlist(oo.mans4[,5]),ylim=c(0,1),type="p",lwd=2,col="red")
lines(unlist(oo.mans4[,6]),ylim=c(0,1),type="p",lwd=2,col="blue")
abline(v=3.5)

plot(unlist(oo.mans5[,1]),ylim=c(0,1))
lines(unlist(oo.mans5[,2]),ylim=c(0,1),type="p",lwd=2,col="red")
lines(unlist(oo.mans5[,3]),ylim=c(0,1),type="p",lwd=2,col="blue")
abline(v=3.5)
plot(unlist(oo.mans5[,4]),ylim=c(0,1))
lines(unlist(oo.mans5[,5]),ylim=c(0,1),type="p",lwd=2,col="red")
lines(unlist(oo.mans5[,6]),ylim=c(0,1),type="p",lwd=2,col="blue")
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
lines(unlist(oo.loa1[,2]),ylim=c(0,1),type="p",lwd=2,col="red")
lines(unlist(oo.loa1[,3]),ylim=c(0,1),type="p",lwd=2,col="blue")
legend("topleft",legend=c(expression(rho==0),expression(rho==.5),expression(rho==0.7)),cex=.8,lwd=2,col=c(1,"red","blue"),lwd=2.5,ncol=1,bty="n")
abline(v=3.5)

plot(unlist(oo.loa1[,4]),ylim=c(0,1))
title(expression(alpha[1]==0.5))
lines(unlist(oo.loa1[,5]),ylim=c(0,1),type="p",lwd=2,col="red")
lines(unlist(oo.loa1[,6]),ylim=c(0,1),type="p",lwd=2,col="blue")
abline(v=3.5)

plot(unlist(oo.loa2[,1]),ylim=c(0,1))
lines(unlist(oo.loa2[,2]),ylim=c(0,1),type="p",lwd=2,col="red")
lines(unlist(oo.loa2[,3]),ylim=c(0,1),type="p",lwd=2,col="blue")
abline(v=3.5)
plot(unlist(oo.loa2[,4]),ylim=c(0,1))
lines(unlist(oo.loa2[,5]),ylim=c(0,1),type="p",lwd=2,col="red")
lines(unlist(oo.loa2[,6]),ylim=c(0,1),type="p",lwd=2,col="blue")
abline(v=3.5)

plot(unlist(oo.loa3[,1]),ylim=c(0,1))
lines(unlist(oo.loa3[,2]),ylim=c(0,1),type="p",lwd=2,col="red")
lines(unlist(oo.loa3[,3]),ylim=c(0,1),type="p",lwd=2,col="blue")
abline(v=3.5)
plot(unlist(oo.loa3[,4]),ylim=c(0,1))
lines(unlist(oo.loa3[,5]),ylim=c(0,1),type="p",lwd=2,col="red")
lines(unlist(oo.loa3[,6]),ylim=c(0,1),type="p",lwd=2,col="blue")
abline(v=3.5)

plot(unlist(oo.loa4[,1]),ylim=c(0,1))
lines(unlist(oo.loa4[,2]),ylim=c(0,1),type="p",lwd=2,col="red")
lines(unlist(oo.loa4[,3]),ylim=c(0,1),type="p",lwd=2,col="blue")
abline(v=3.5)
plot(unlist(oo.loa4[,4]),ylim=c(0,1))
lines(unlist(oo.loa4[,5]),ylim=c(0,1),type="p",lwd=2,col="red")
lines(unlist(oo.loa4[,6]),ylim=c(0,1),type="p",lwd=2,col="blue")
abline(v=3.5)


plot(unlist(oo.loa5[,1]),ylim=c(0,1))
lines(unlist(oo.loa5[,2]),ylim=c(0,1),type="p",lwd=2,col="red")
lines(unlist(oo.loa5[,3]),ylim=c(0,1),type="p",lwd=2,col="blue")
abline(v=3.5)
plot(unlist(oo.loa5[,4]),ylim=c(0,1))
lines(unlist(oo.loa5[,5]),ylim=c(0,1),type="p",lwd=2,col="red")
lines(unlist(oo.loa5[,6]),ylim=c(0,1),type="p",lwd=2,col="blue")
abline(v=3.5)

plot(unlist(oo.loa6[,1]),ylim=c(0,1))
lines(unlist(oo.loa6[,2]),ylim=c(0,1),type="p",lwd=2,col="red")
lines(unlist(oo.loa6[,3]),ylim=c(0,1),type="p",lwd=2,col="blue")
abline(v=3.5)
plot(unlist(oo.loa6[,4]),ylim=c(0,1))
lines(unlist(oo.loa6[,5]),ylim=c(0,1),type="p",lwd=2,col="red")
lines(unlist(oo.loa6[,6]),ylim=c(0,1),type="p",lwd=2,col="blue")
abline(v=3.5)

dev.off()