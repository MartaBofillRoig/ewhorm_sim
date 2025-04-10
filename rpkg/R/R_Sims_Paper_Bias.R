#define baseline-parameters

scen.d<-data.frame(disease=c("onchocerciasis","mansonellosis","loiasis"),mu_raw_0=c(19,1838,5000), sd_raw_0=c(30,2565,4000),bound=c(-2.05,0,0))


#define scenarios of reduction rates for 6 and 12 months
scen.eI<-data.frame(r0_6=c(0,0,0,0,0,0,0,0),
                    r1_6=c(0,0,0,0,0,0,.4,0.5),
                    r2_6=c(0,0,0,.2,.3,.4,.4,.5),
                    r3_6=c(0,.5,.6,.5,.5,.5,.4,.5),
                    r0_12=c(0,0,0,0,0,0,0,0),
                    r1_12=c(0,0,0,0,0,0,.5,0.6),
                    r2_12=c(0,0,0,.3,.4,.5,.5,.6),
                    r3_12=c(0,.6,.7,.6,.6,.6,.5,.6))

scen.eI<-scen.eI[c(-3,-4,-8),]


scen.eII<-data.frame(r0_6=c(0,0,0,0,0,0,0,0),r1_6=c(0,0,0,0,0,0,.3,0.4),r2_6=c(0,0,0,.1,.2,.3,.3,.4),r3_6=c(0,.4,.5,.4,.4,.4,.3,.4),
                    r0_12=c(0,0,0,0,0,0,0,0),r1_12=c(0,0,0,0,0,0,.5,0.6),r2_12=c(0,0,0,.3,.4,.5,.5,.6),r3_12=c(0,.6,.7,.6,.6,.6,.5,.6))

scen.eII<-scen.eII[c(-3,-4,-8),]




f.rr<-function(r) #function for total responder
{
  max(r-0.2,.1)
}


simulation.scenario<-function(d,i,n_trials,j,t,l,k,scen)
{
  n_arms<-4
  #rmonth<-1
  alpha<-.025
  #sim_out1<-1
  sel_scen<-1
  side1<-1
  dropout<-.1
  #rho<-.5
  #N1<-120
  #N2<-80
  scen.e<-scen
  N2<-200-N1.v[k]
  
    res<-mapply(simul_res,scen.d$mu_raw_0[d], scen.d$sd_raw_0[d] , scen.e$r0_6[i],scen.e$r1_6[i],scen.e$r2_6[i],scen.e$r3_6[i], scen.e$r0_12[i],
                                            scen.e$r1_12[i],scen.e$r2_12[i],scen.e$r3_12[i],  rho.v[l] ,
                                            n_trials,4,N1.v[k] , N2, alpha1.v[j] , alpha,
                                            sel_scen, side1,test.v[t],dropout,
                                            f.rr(scen.e$r0_12[i]),f.rr(scen.e$r1_12[i]),f.rr(scen.e$r2_12[i]),f.rr(scen.e$r3_12[i]),scen.d$bound[d])
  
    unlist(c(d,i,scen.e[i,],alpha1.v[j],N1.v[k],rho.v[l],test.v[t],#round(
             c(res[7],res[8],res[9],res[6],res[10],res[1],res[2],res[3],
               res[c(19:22)],NA,res[(15:18)],NA,res[c(25:72)])#,2)*100
             ))
}
  

wrap<-function(d,scen,n_trials=50000)
{
  d.1<-matrix(nrow=75,ncol=(17+15+6+4+1+10+12+3+3+3+6))#(#))

  h<-0
  
  for (i in 1:(dim(scen)[1]))
    for (j in 1:length(alpha1.v))
      for (k in 1:length(N1.v))
        for(l in 1:length(rho.v))
          for (t in 1:length(test.v))
        {
        h<-h+1
        print(h)
        d.1[h,]<-simulation.scenario(d,i,n_trials,j,t,l,k,scen)
        }
        
  d.1<-as.data.frame(d.1)
  names(d.1)<-c("disease","scenario","r0_6","r1_6","r2_6","r3_6","r0_12","r1_12","r2_12","r3_12","alpha1","N1","rho","test",
                "Pow.Low","Pow.Me","Pow.Hi","Pow.cond",
                "Pow.Disj","SelP.Low","SelP.Me","SelP.Hi",
                "M1Pow.Low","M1Pow.Me","M1Pow.Hi","M1Pow.Disj","M1Pow.cond","M2Pow.Low","M2Pow.Me","M2Pow.Hi","M2Pow.Disj","M2Pow.cond",
                "TrueC.Low","TrueC.Me","TrueC.Hi","C.Low","C.Me","C.Hi","sdC.Low","sdC.Me","sdC.Hi","Bias.Low","Bias.Me","Bias.Hi",
                "M1C.Low","M1C.Me","M1C.Hi","CIupC.Low","CIupC.Me","CIupC.Hi","CICov.Low","CICov.Me","CICov.Hi","Ccond.Low","Ccond.Me","Ccond.Hi",
                "biasCcond.Low","biasCcond.Me","biasCcond.Hi",#"CICovBH.Low","CICovBH.Me","CICovBH.Hi",
                "CICovcond.Low","CICovcond.Me","CICovcond.Hi",
                "CImeaninvnorm.Low","CImeaninvnorm.Me","CImeaninvnorm.Hi","CICovinvnorm.Low","CICovinvnorm.Me","CICovinvnorm.Hi",
                "CImeancond.Low","CImeancond.Me","CImeancond.Hi","meanCinv.Low","meanCinv.Me","meanCinv.Hi","biasinv.Low","biasinv.Me","biasinv.H")
  d.1
}



#####Simulation Start



alpha1.v<-seq(0.1,0.5,.1)
N1.v<-120#c(80,120)
rho.v<-.5#c(0.4,0.5)
test.v<-c(0,3,4)


#test and alpha1

set.seed(1234)
w1I<-wrap(1,scen.eI) #disease 1
w1I
save(w1I,file="w1ICI.RData")
w2I<-wrap(2,scen.eI) #disease 2
w2I
w3I<-wrap(3,scen.eI) #disease 3
w3I

#
save(w1I,file="w1I.RData")
save(w2I,file="w2I.RData")
save(w3I,file="w3I.RData")
load("w1I.RData")
load("w2I.RData")
load("w3I.RData")


#####nur Wilcoxon wegen concordance;
#alpha1.v<-seq(0.1,0.5,.1)
#N1.v<-120
#rho.v<-.5
#test.v<-4

#test and alpha1

#w1IW<-wrap(1,scen.eI) #disease 1
#w1IW
#w2IW<-wrap(2,scen.eI) #disease 2
#w2IW
#w3IW<-wrap(3,scen.eI) #disease 3
#w3IW

#save(w1IW,file="w1IW.RData")
#save(w2IW,file="w2IW.RData")
#save(w3IW,file="w3IW.RData")
#load("w1IW.RData")
#load("w2IW.RData")
#load("w3IW.RData")


w1II<-wrap(1,scen.eII) #disease 1
w1II
w2II<-wrap(2,scen.eII) #disease 2
w2II
w3II<-wrap(3,scen.eII) #disease 3
w3II

save(w1I,file="w1II.RData")
save(w2I,file="w2II.RData")
save(w3I,file="w3II.RData")
load("w1II.RData")
load("w2II.RData")
load("w3II.RData")


#rho and test
alpha1.v<-.3#seq(0.1,0.5,.1)
N1.v<-120#c(80,120)
rho.v<-c(0.4,0.5,0.6)
test.v<-c(0,3,4)

w1I.rho<-wrap(1,scen.eI) #disease 1
w1I.rho
w2I.rho<-wrap(2,scen.eI) #disease 2
w2I.rho
w3I.rho<-wrap(3,scen.eI) #disease 3
w3I.rho

save(w1I.rho,file="w1I.rho.RData")
save(w2I.rho,file="w2I.rho.RData")
save(w3I.rho,file="w3I.rho.RData")
load("w1I.rho.RData")
load("w2I.rho.RData")
load("w3I.rho.RData")

#N1 and test
alpha1.v<-.3#seq(0.1,0.5,.1)
N1.v<-c(80,100,120)
rho.v<-.5#c(0.4,0.5,0.6)
test.v<-c(0,3,4)

w1I.N1<-wrap(1,scen.eI) #disease 1
w1I.N1
w2I.N1<-wrap(2,scen.eI) #disease 2
w2I.N1
w3I.N1<-wrap(3,scen.eI) #disease 3
w3I.N1

save(w1I.N1,file="w1I.N1.RData")
save(w2I.N1,file="w2I.N1.RData")
save(w3I.N1,file="w3I.N1.RData")
load("w1I.N1.RData")
load("w2I.N1.RData")
load("w3I.N1.RData")


scen.d<-data.frame(disease=c("onchocerciasis","mansonellosis","loiasis"),mu_raw_0=c(15,1000,4000), sd_raw_0=c(40,3500,5000),bound=c(-2.05,0,0))


alpha1.v<-seq(0.1,0.5,.1)
N1.v<-120#c(80,120)
rho.v<-.5#c(0.4,0.5)
test.v<-c(0,3,4)

#test and alpha1

w1I.worst<-wrap(1,scen.eI) #disease 1
w1I.worst
w2I.worst<-wrap(2,scen.eI) #disease 2
w2I.worst
w3I.worst<-wrap(3,scen.eI) #disease 3
w3I.worst

save(w1I.worst,file="w1I.worst.RData")
save(w2I.worst,file="w2I.worst.RData")
save(w3I.worst,file="w3I.worst.RData")
load("w1I.worst.RData")
load("w2I.worst.RData")
load("w3I.worst.RData")


#####Simulation End
##############################################################################################

