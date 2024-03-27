#mu_raw_0, sd_raw_0 , r0_6lo,r1_6lo,r2_6lo,r3_6lo, r0_12lo,r1_12lo,r2_12lo,r3_12lo,  rho ,
#n_trials,4,N1 , N, rmonth, alpha1 , alpha,
#sim_out1,sel_scen, side1,test,dropout,rr0,rr1,rr2,rr3


scen.d<-data.frame(disease=c("onchocerciasis","mansonellosis","loiasis"),mu_raw_0=c(19,1838,5000), sd_raw_0=c(30,2565,4000),bound=c(-2.05,0,0))

scen.e<-data.frame(r0_6=c(0,0,0,0,0,0,0),r1_6=c(0,0,0,0,0,0.5,0),r2_6=c(0,0,.2,.3,.4,.5,0),r3_6=c(.5,.6,.5,.5,.5,.5,0),
                              r0_12=c(0,0,0,0,0,0,0),r1_12=c(0,0,0,0,0,0.6,0),r2_12=c(0,0,.3,.4,.5,.6,0),r3_12=c(.6,.7,.6,.6,.6,.6,0))



f.rr<-function(r)
{
  max(r-0.2,.1)
}

simulation.scenario<-function(d,i,n_trials,alpha1,test)
{
  rho<-.5
  n_arms<-4
  N1<-120
  N2<-80
  rmonth<-1
  alpha<-.025
  sim_out1<-1
  sel_scen<-1
  side1<-1
  dropout<-.1
  
  
  
  
    res<-mapply(simul_res,scen.d$mu_raw_0[d], scen.d$sd_raw_0[d] , scen.e$r0_6[i],scen.e$r1_6[i],scen.e$r2_6[i],scen.e$r3_6[i], scen.e$r0_12[i],
                                            scen.e$r1_12[i],scen.e$r2_12[i],scen.e$r3_12[i],  rho ,
                                            n_trials,4,N1 , N2, rmonth, alpha1 , alpha,
                                            sim_out1,sel_scen, side1,test,dropout,
                                            f.rr(scen.e$r0_12[i]),f.rr(scen.e$r1_12[i]),f.rr(scen.e$r2_12[i]),f.rr(scen.e$r3_12[i]),scen.d$bound[d])
  
  
    unlist(c(d,scen.e[i,],round(c(res[7],res[8],res[9],res[6],res[10],res[1],res[2],res[3]),2)*100))
    
  
}
  

wrap<-function(d,n_trials=10000,alpha1=.3,test=4)
{
  d.1<-matrix(nrow=7,ncol=17)

  for (i in 1:(dim(scen.e)[1]))
      d.1[i,]<-simulation.scenario(d,i,n_trials,alpha1=alpha1,test=test)

  d.1<-as.data.frame(d.1)
  names(d.1)<-c("disease","r0_6","r1_6","r2_6","r3_6","r0_12","r1_12","r2_12","r3_12","Pow.Low","Pow.Me","Pow.Hi","Pow.cond","Pow.Disj","SelP.Low","SelP.Me","SelP.Hi")

  d.1
}

w1<-wrap(1)
w1
w2<-wrap(2)
w2
w3<-wrap(3)
w3

w1.1<-wrap(1,alpha1=1)
w1.1
w2.1<-wrap(2,alpha1=1)
w2.1
w3.1<-wrap(3,alpha1=1)
w3.1

w1f<-wrap(1,test=5)
w1f
w2f<-wrap(2,test=5)
w2f
w3f<-wrap(3,test=5)
w3f

w1.1f<-wrap(1,alpha1=1,test=5)
w1.1f
w2.1f<-wrap(2,alpha1=1,test=5)
w2.1f
w3.1f<-wrap(3,alpha1=1,test=5)
w3.1f









d.2<-matrix(nrow=7,ncol=17)

for (i in 1:(dim(scen.e)[1]))
  d.2[i,]<-simulation.scenario(2,i)

d.2<-as.data.frame(d.2)
names(d.2)<-c("disease","r0_6","r1_6","r2_6","r3_6","r0_12","r1_12","r2_12","r3_12","Pow.Low","Pow.Me","Pow.Hi","Pow.cond","Pow.Disj","SelP.Low","SelP.Me","SelP.Hi")

d.2     


d.3<-matrix(nrow=7,ncol=17)

for (i in 1:(dim(scen.e)[1]))
  d.3[i,]<-simulation.scenario(3,i)

d.3<-as.data.frame(d.3)
names(d.3)<-c("disease","r0_6","r1_6","r2_6","r3_6","r0_12","r1_12","r2_12","r3_12","Pow.Low","Pow.Me","Pow.Hi","Pow.cond","Pow.Disj","SelP.Low","SelP.Me","SelP.Hi")

d.3    




#alpha1=1


d.1<-matrix(nrow=7,ncol=17)

for (i in 1:(dim(scen.e)[1]))
  d.1[i,]<-simulation.scenario(1,i,alpha1=1)

d.1<-as.data.frame(d.1)
names(d.1)<-c("disease","r0_6","r1_6","r2_6","r3_6","r0_12","r1_12","r2_12","r3_12","Pow.Low","Pow.Me","Pow.Hi","Pow.cond","Pow.Disj","SelP.Low","SelP.Me","SelP.Hi")

d.1

d.2<-matrix(nrow=7,ncol=17)

for (i in 1:(dim(scen.e)[1]))
  d.2[i,]<-simulation.scenario(2,i,alpha1=1)

d.2<-as.data.frame(d.2)
names(d.2)<-c("disease","r0_6","r1_6","r2_6","r3_6","r0_12","r1_12","r2_12","r3_12","Pow.Low","Pow.Me","Pow.Hi","Pow.cond","Pow.Disj","SelP.Low","SelP.Me","SelP.Hi")

d.2     


d.3<-matrix(nrow=7,ncol=17)

for (i in 1:(dim(scen.e)[1]))
  d.3[i,]<-simulation.scenario(3,i,alpha1=1)

d.3<-as.data.frame(d.3)
names(d.3)<-c("disease","r0_6","r1_6","r2_6","r3_6","r0_12","r1_12","r2_12","r3_12","Pow.Low","Pow.Me","Pow.Hi","Pow.cond","Pow.Disj","SelP.Low","SelP.Me","SelP.Hi")

d.3                            