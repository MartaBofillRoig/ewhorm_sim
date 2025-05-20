

#Plot of all Figures from Manuscript and Supplementary material ()


#functions to make long data sets for depiction in Figures

dat_reshape<-function(dat)
{
  dat<-dat[!is.na(dat$disease),]
  
  dat.a<-reshape(data = dat,
                 #idvar= "id",
                 varying = matrix(c(15:19,23,24,25,27,26,28,29,30,32,31),nrow=3,byrow=TRUE), #We need to specify here the columns to be reshaped
                 sep= "",
                 timevar= "H",
                 times = c(1:5),
                 v.names=c("Pow","PowMa1","PowMa2"),
                 #new.row.names= 1:10,
                 direction = "long")
  
  dat.a<-dat.a[,-ncol(dat.a)] #clear id
  
  dat<-reshape(data = dat.a,
               #idvar= "id",
               varying = matrix(c("Pow","PowMa1","PowMa2"),nrow=1),#We need to specify here the columns to be reshaped
               sep= "",
               timevar= "testMA",
               times = c(1:3),
               v.names=c("PowAll"),
               direction = "long")
  
  
  dat$H<-as.factor(dat$H)
  levels(dat$H)<-c("Low dose","Medium dose","High dose","Cond. Pow High dose","Disjunctive")
  dat$test<-as.factor(dat$test)
  levels(dat$test)<-c("LM","WilcoxCC","WilcoxC")
  levels(dat$testMA)<-c("AD","MA1","MA2")
  
  dat$Int_testMA<-interaction(dat$testMA,dat$test)
  
  dat$Int_testMA<-as.factor(dat$Int_testMA)
  levels(dat$Int_testMA)<-c("LM AD","LM MA1","LM MA2","WilcoxCC AD","WilcoxCC MA1","WilcoxCC MA2","WilcoxC AD","WilcoxC MA1",
                            "WilcoxC MA2")
  dat
}


#make long only for selection prob.
dat_reshape_selP<-function(dat)
{
  dat<-dat[!is.na(dat$disease),]
  
  dat<-reshape(data = dat,
               #idvar= "id",
               varying = c(20:22,18), #We need to specify here the columns to be reshaped
               sep= "",
               timevar= "SelP.dose",
               times = c(1:4),
               v.names=c("SelP"),
               #new.row.names= 1:10,
               direction = "long")
  
  
  dat$SelP.dose<-as.factor(dat$SelP.dose)
  levels(dat$SelP.dose)<-c("Sel.Pr Low dose","Sel.Pr Medium dose","Sel.Pr High dose","Cond.Pow High dose")
  dat$test<-as.factor(dat$test)
  levels(dat$test)<-c("LM","WilcoxCC","WilcoxC")
  
  dat
}


#generate long data sets from each simulation data set
w1I.long.selP<-dat_reshape_selP(w1I)
w1II.long.selP<-dat_reshape_selP(w1II)
w1I.rho.long.selP<-dat_reshape_selP(w1I.rho)
w1I.N1.long.selP<-dat_reshape_selP(w1I.N1)
w1I.worst.long.selP<-dat_reshape_selP(w1I.worst)
w2I.long.selP<-dat_reshape_selP(w2I)
w2II.long.selP<-dat_reshape_selP(w2II)
w2I.rho.long.selP<-dat_reshape_selP(w2I.rho)
w2I.N1.long.selP<-dat_reshape_selP(w2I.N1)
w2I.worst.long.selP<-dat_reshape_selP(w2I.worst)
w3I.long.selP<-dat_reshape_selP(w3I)
w3II.long.selP<-dat_reshape_selP(w3II)
w3I.rho.long.selP<-dat_reshape_selP(w3I.rho)
w3I.N1.long.selP<-dat_reshape_selP(w3I.N1)
w3I.worst.long.selP<-dat_reshape_selP(w3I.worst)

w1I.long<-dat_reshape(w1I)
w1II.long<-dat_reshape(w1II)
w1I.rho.long<-dat_reshape(w1I.rho)
w1I.N1.long<-dat_reshape(w1I.N1)
w1I.worst.long<-dat_reshape(w1I.worst)
w2I.long<-dat_reshape(w2I)
w2II.long<-dat_reshape(w2II)
w2I.rho.long<-dat_reshape(w2I.rho)
w2I.N1.long<-dat_reshape(w2I.N1)
w2I.worst.long<-dat_reshape(w2I.worst)
w3I.long<-dat_reshape(w3I)
w3II.long<-dat_reshape(w3II)
w3I.rho.long<-dat_reshape(w3I.rho)
w3I.N1.long<-dat_reshape(w3I.N1)
w3I.worst.long<-dat_reshape(w3I.worst)


w1I.long.selP$scenario1<-factor(w1I.long.selP$scenario,labels=c("No effect","Effic. in high dose","Trend (a)", "Trend (b)", "All doses eff."))
w1II.long.selP$scenario1<-factor(w1II.long.selP$scenario,labels=c("No effect","Effic. in high dose","Trend (a)", "Trend (b)", "All doses eff."))
w1I.rho.long.selP$scenario1<-factor(w1I.rho.long.selP$scenario,labels=c("No effect","Effic. in high dose","Trend (a)", "Trend (b)", "All doses eff."))
w1I.N1.long.selP$scenario1<-factor(w1I.N1.long.selP$scenario,labels=c("No effect","Effic. in high dose","Trend (a)", "Trend (b)", "All doses eff."))
w1I.worst.long.selP$scenario1<-factor(w1I.worst.long.selP$scenario,labels=c("No effect","Effic. in high dose","Trend (a)", "Trend (b)", "All doses eff."))
w2I.long.selP$scenario1<-factor(w2I.long.selP$scenario,labels=c("No effect","Effic. in high dose","Trend (a)", "Trend (b)", "All doses eff."))
w2II.long.selP$scenario1<-factor(w2II.long.selP$scenario,labels=c("No effect","Effic. in high dose","Trend (a)", "Trend (b)", "All doses eff."))
w2I.rho.long.selP$scenario1<-factor(w2I.rho.long.selP$scenario,labels=c("No effect","Effic. in high dose","Trend (a)", "Trend (b)", "All doses eff."))
w2I.N1.long.selP$scenario1<-factor(w2I.N1.long.selP$scenario,labels=c("No effect","Effic. in high dose","Trend (a)", "Trend (b)", "All doses eff."))
w2I.worst.long.selP$scenario1<-factor(w2I.worst.long.selP$scenario,labels=c("No effect","Effic. in high dose","Trend (a)", "Trend (b)", "All doses eff."))
w3I.long.selP$scenario1<-factor(w3I.long.selP$scenario,labels=c("No effect","Effic. in high dose","Trend (a)", "Trend (b)", "All doses eff."))
w3II.long.selP$scenario1<-factor(w3II.long.selP$scenario,labels=c("No effect","Effic. in high dose","Trend (a)", "Trend (b)", "All doses eff."))
w3I.rho.long.selP$scenario1<-factor(w3I.rho.long.selP$scenario,labels=c("No effect","Effic. in high dose","Trend (a)", "Trend (b)", "All doses eff."))
w3I.N1.long.selP$scenario1<-factor(w3I.N1.long.selP$scenario,labels=c("No effect","Effic. in high dose","Trend (a)", "Trend (b)", "All doses eff."))
w3I.worst.long.selP$scenario1<-factor(w3I.worst.long.selP$scenario,labels=c("No effect","Effic. in high dose","Trend (a)", "Trend (b)", "All doses eff."))


w1I.long$scenario1<-factor(w1I.long$scenario,labels=c("No effect","Effic. in high dose","Trend (a)", "Trend (b)", "All doses eff."))
w1II.long$scenario1<-factor(w1II.long$scenario,labels=c("No effect","Effic. in high dose","Trend (a)", "Trend (b)", "All doses eff."))
w1I.rho.long$scenario1<-factor(w1I.rho.long$scenario,labels=c("No effect","Effic. in high dose","Trend (a)", "Trend (b)", "All doses eff."))
w1I.N1.long$scenario1<-factor(w1I.N1.long$scenario,labels=c("No effect","Effic. in high dose","Trend (a)", "Trend (b)", "All doses eff."))
w1I.worst.long$scenario1<-factor(w1I.worst.long$scenario,labels=c("No effect","Effic. in high dose","Trend (a)", "Trend (b)", "All doses eff."))
w2I.long$scenario1<-factor(w2I.long$scenario,labels=c("No effect","Effic. in high dose","Trend (a)", "Trend (b)", "All doses eff."))
w2II.long$scenario1<-factor(w2II.long$scenario,labels=c("No effect","Effic. in high dose","Trend (a)", "Trend (b)", "All doses eff."))
w2I.rho.long$scenario1<-factor(w2I.rho.long$scenario,labels=c("No effect","Effic. in high dose","Trend (a)", "Trend (b)", "All doses eff."))
w2I.N1.long$scenario1<-factor(w2I.N1.long$scenario,labels=c("No effect","Effic. in high dose","Trend (a)", "Trend (b)", "All doses eff."))
w2I.worst.long$scenario1<-factor(w2I.worst.long$scenario,labels=c("No effect","Effic. in high dose","Trend (a)", "Trend (b)", "All doses eff."))
w3I.long$scenario1<-factor(w3I.long$scenario,labels=c("No effect","Effic. in high dose","Trend (a)", "Trend (b)", "All doses eff."))
w3II.long$scenario1<-factor(w3II.long$scenario,labels=c("No effect","Effic. in high dose","Trend (a)", "Trend (b)", "All doses eff."))
w3I.rho.long$scenario1<-factor(w3I.rho.long$scenario,labels=c("No effect","Effic. in high dose","Trend (a)", "Trend (b)", "All doses eff."))
w3I.N1.long$scenario1<-factor(w3I.N1.long$scenario,labels=c("No effect","Effic. in high dose","Trend (a)", "Trend (b)", "All doses eff."))
w3I.worst.long$scenario1<-factor(w3I.worst.long$scenario,labels=c("No effect","Effic. in high dose","Trend (a)", "Trend (b)", "All doses eff."))




#Functions to generate ggplots from generated data (long format)
#plot power as function of alpha1;
library(ggplot2)
plotpow.alpha1<-function(dat)
{
  plotpow<-ggplot(dat[which((dat$scenario>1)&(dat$H %in% c("Low dose","Medium dose","High dose","Disjunctive"))),]##[which(w1I.long$scenario %in% c(2,5,6,7,8)),]#[which(w1I$xx %in% c("r", "f")),]
                  , aes(x=alpha1)) +                                     
    geom_line(aes(y=PowAll,group=Int_testMA,color=Int_testMA,linetype=Int_testMA), linewidth = .5) + 
    #geom_line(aes(y=Power.Nix.,group = test), size = 1.5,linetype="twodash") +# Increase line size to 2
    scale_color_manual(values = c("black","red","green","black","red","green","black","red","green"))+#,guide=guide_legend(title=NA)) +
    scale_linetype_manual(values=c(1,1,1,2,2,2,6,6,6))+
    facet_grid(scenario1~H)+ylim(0, 1)+theme(legend.position="bottom")+
    guides(color=guide_legend(title="Design",ncol=3),
           linetype=guide_legend(title="Design",ncol=3))
  
  #facet_matrix
  plot1<-plotpow+labs(x=expression(alpha[1]),y="Power")
  plot1
}
plotpow.alpha1(w1I.long)
dat<-w1I.long

#plot power as function of rho;
plotpow.rho<-function(dat)
{
  plotpow<-ggplot(dat[which((dat$scenario>1)&(dat$H %in% c("Low dose","Medium dose","High dose","Disjunctive"))),]# %in% c(2,5,6,7,8)),]#[which(w1I$xx %in% c("r", "f")),]
                  , aes(x=rho)) +                                     
    geom_line(aes(y=PowAll,group=Int_testMA,color=Int_testMA,linetype=Int_testMA), linewidth = .5) + 
    #geom_line(aes(y=Power.Nix.,group = test), size = 1.5,linetype="twodash") +# Increase line size to 2
    scale_color_manual(values = c("black","red","green","black","red","green","black","red","green"))+#,guide=guide_legend(title=NA)) +
    scale_linetype_manual(values=c(1,1,1,2,2,2,6,6,6))+
    facet_grid(scenario1~H)+ylim(0, 1)+theme(legend.position="bottom")+
    guides(color=guide_legend(title="Design",ncol=3),
           linetype=guide_legend(title="Design",ncol=3))
  
  #facet_matrix
  plot1<-plotpow+labs(x=expression(rho),y="Power")
  plot1
}

#plot power as function of N1;
plotpow.N1<-function(dat)
{
  plotpow<-ggplot(dat[which((dat$scenario>1)&(dat$H %in% c("Low dose","Medium dose","High dose","Disjunctive"))),]##[which(w1I.long$scenario %in% c(2,5,6,7,8)),]#[which(w1I$xx %in% c("r", "f")),]
                  , aes(x=N1)) +                                     
    geom_line(aes(y=PowAll,group=Int_testMA,color=Int_testMA,linetype=Int_testMA), linewidth = .5) + 
    #geom_line(aes(y=Power.Nix.,group = test), size = 1.5,linetype="twodash") +# Increase line size to 2
    scale_color_manual(values = c("black","red","green","black","red","green","black","red","green"))+#,guide=guide_legend(title=NA)) +
    scale_linetype_manual(values=c(1,1,1,2,2,2,6,6,6))+
    facet_grid(scenario1~H)+ylim(0, 1)+theme(legend.position="bottom")+
    guides(color=guide_legend(title="Design",ncol=3),
           linetype=guide_legend(title="Design",ncol=3))
  
  #facet_matrix
  plot1<-plotpow+labs(x=expression(N[1]),y="Power")
  plot1
}

#Plot type I error as function of alpha1
plotpow.alpha1.error<-function(dat)
{
  plotpow<-ggplot(dat[which((dat$scenario==1)&(dat$H %in% c("Low dose","Medium dose","High dose","Disjunctive"))),]##[which(w1I.long$scenario %in% c(2,5,6,7,8)),]#[which(w1I$xx %in% c("r", "f")),]
                  , aes(x=alpha1)) +                                     
    geom_line(aes(y=PowAll,group=Int_testMA,color=Int_testMA,linetype=Int_testMA), linewidth = .5) + 
    #geom_line(aes(y=Power.Nix.,group = test), size = 1.5,linetype="twodash") +# Increase line size to 2
    scale_color_manual(values = c("black","red","green","black","red","green","black","red","green"))+#,guide=guide_legend(title=NA)) +
    scale_linetype_manual(values=c(1,1,1,2,2,2,6,6,6))+
    facet_grid(scenario1~H)+ylim(0, 0.05)+theme(legend.position="bottom")+
    guides(color=guide_legend(title="Design",ncol=3),
           linetype=guide_legend(title="Design",ncol=3))
  
  #facet_matrix
  plot1<-plotpow+labs(x=expression(alpha[1]),y="Power")
  plot1
}


#Plot type I error as function of rho
plotpow.rho.error<-function(dat)
{
  plotpow<-ggplot(dat[which((dat$scenario==1)&(dat$H %in% c("Low dose","Medium dose","High dose","Disjunctive"))),]##[which(w1I.long$scenario %in% c(2,5,6,7,8)),]#[which(w1I$xx %in% c("r", "f")),]
                  , aes(x=rho)) +                                     
    geom_line(aes(y=PowAll,group=Int_testMA,color=Int_testMA,linetype=Int_testMA), linewidth = .5) + 
    #geom_line(aes(y=Power.Nix.,group = test), size = 1.5,linetype="twodash") +# Increase line size to 2
    scale_color_manual(values = c("black","red","green","black","red","green","black","red","green"))+#,guide=guide_legend(title=NA)) +
    scale_linetype_manual(values=c(1,1,1,2,2,2,6,6,6))+
    facet_grid(scenario1~H)+ylim(0, 0.05)+theme(legend.position="bottom")+
    guides(color=guide_legend(title="Design",ncol=3),
           linetype=guide_legend(title="Design",ncol=3))
  
  #facet_matrix
  plot1<-plotpow+labs(x=expression(rho),y="Power")
  plot1
}
plotpow.rho.error(w1I.rho.long)


#Plot type I error as function of N1
plotpow.N1.error<-function(dat)
{
  plotpow<-ggplot(dat[which((dat$scenario==1)&(dat$H %in% c("Low dose","Medium dose","High dose","Disjunctive"))),]##[which(w1I.long$scenario %in% c(2,5,6,7,8)),]#[which(w1I$xx %in% c("r", "f")),]
                  , aes(x=N1)) +                                     
    geom_line(aes(y=PowAll,group=Int_testMA,color=Int_testMA,linetype=Int_testMA), linewidth = .5) + 
    #geom_line(aes(y=Power.Nix.,group = test), size = 1.5,linetype="twodash") +# Increase line size to 2
    scale_color_manual(values = c("black","red","green","black","red","green","black","red","green"))+#,guide=guide_legend(title=NA)) +
    scale_linetype_manual(values=c(1,1,1,2,2,2,6,6,6))+
    facet_grid(scenario1~H)+ylim(0, 0.05)+theme(legend.position="bottom")+
    guides(color=guide_legend(title="Design",ncol=3),
           linetype=guide_legend(title="Design",ncol=3))
  
  #facet_matrix
  plot1<-plotpow+labs(x=expression(N[1]),y="Power")
  plot1
}
#



#Plot selection prob. as function of alpha1
plotpow.alpha1.selP<-function(dat)
{
  plotpow<-ggplot(dat#[which(dat$scenario>1),]##[which(w1I.long$scenario %in% c(2,5,6,7,8)),]#[which(w1I$xx %in% c("r", "f")),]
                  , aes(x=alpha1)) +                                     
    geom_line(aes(y=SelP,group=test,color=test,linetype=test), linewidth = .5) + 
    scale_color_manual(values = c("black","red","green"))+#,guide=guide_legend(title=NA)) +
    scale_linetype_manual(values=c(1,1,1))+
    facet_grid(scenario1~SelP.dose)+ylim(0, 1)+theme(legend.position="bottom")+
    guides(color=guide_legend(title="Design",ncol=3),
           linetype=guide_legend(title="Design",ncol=3))
  
  #facet_matrix
  plot1<-plotpow+labs(x=expression(alpha[1]),y="Power")
  plot1
}

plotpow.alpha1.selP(w3I.long.selP)
#w3I.worst.long.selP$SelP.dose

#Plot selection prob. as function of rho
plotpow.rho.selP<-function(dat)
{
  plotpow<-ggplot(dat[which(dat$scenario>1),]# %in% c(2,5,6,7,8)),]#[which(w1I$xx %in% c("r", "f")),]
                  , aes(x=rho)) +                                     
    geom_line(aes(y=SelP,group=test,color=test,linetype=test), linewidth = .5) + 
    scale_color_manual(values = c("black","red","green"))+#,guide=guide_legend(title=NA)) +
    scale_linetype_manual(values=c(1,1,1))+
    facet_grid(scenario1~SelP.dose)+ylim(0, 1)+theme(legend.position="bottom")+
    guides(color=guide_legend(title="Design",ncol=3),
           linetype=guide_legend(title="Design",ncol=3))
  
  #facet_matrix
  plot1<-plotpow+labs(x=expression(rho),y="Power")
  plot1
}

#Plot selection prob. as function of N1
plotpow.N1.selP<-function(dat)
{
  plotpow<-ggplot(dat[which(dat$scenario>1),]##[which(w1I.long$scenario %in% c(2,5,6,7,8)),]#[which(w1I$xx %in% c("r", "f")),]
                  , aes(x=N1)) +                                     
    geom_line(aes(y=SelP,group=test,color=test,linetype=test), linewidth = .5) + 
    scale_color_manual(values = c("black","red","green"))+#,guide=guide_legend(title=NA)) +
    scale_linetype_manual(values=c(1,1,1))+
    facet_grid(scenario1~SelP.dose)+ylim(0, 1)+theme(legend.position="bottom")+
    guides(color=guide_legend(title="Design",ncol=3),
           linetype=guide_legend(title="Design",ncol=3))
  
  #facet_matrix
  plot1<-plotpow+labs(x=expression(N[1]),y="Power")
  plot1
}

#Generate the individual plots

pdf("plot_onchoR.pdf")
print(plotpow.alpha1(w1I.long)) 
print(plotpow.alpha1(w1II.long))
print(plotpow.alpha1(w1I.worst.long))
print(plotpow.rho(w1I.rho.long))
print(plotpow.N1(w1I.N1.long))
print(plotpow.alpha1.error(w1I.long))
print(plotpow.alpha1.error(w1II.long))
print(plotpow.alpha1.error(w1I.worst.long))
print(plotpow.rho.error(w1I.rho.long))
print(plotpow.N1.error(w1I.N1.long))
print(plotpow.alpha1.selP(w1I.long.selP)) 
print(plotpow.alpha1.selP(w1II.long.selP))
print(plotpow.alpha1.selP(w1I.worst.long.selP))
print(plotpow.rho.selP(w1I.rho.long.selP))
print(plotpow.N1.selP(w1I.N1.long.selP))
dev.off() 


pdf("plot_mansR.pdf")
print(plotpow.alpha1(w2I.long))     
print(plotpow.alpha1(w2II.long))
print(plotpow.alpha1(w2I.worst.long))
print(plotpow.rho(w2I.rho.long))
print(plotpow.N1(w2I.N1.long))
print(plotpow.alpha1.error(w2I.long))   
print(plotpow.alpha1.error(w2II.long))
print(plotpow.alpha1.error(w2I.worst.long))
print(plotpow.rho.error(w2I.rho.long))
print(plotpow.N1.error(w2I.N1.long))
print(plotpow.alpha1.selP(w2I.long.selP))  
print(plotpow.alpha1.selP(w2II.long.selP))
print(plotpow.alpha1.selP(w2I.worst.long.selP))
print(plotpow.rho.selP(w2I.rho.long.selP))
print(plotpow.N1.selP(w2I.N1.long.selP))
dev.off() 


pdf("plot_loaR.pdf")
print(plotpow.alpha1(w3I.long))
print(plotpow.alpha1(w3II.long))
print(plotpow.alpha1(w3I.worst.long))
print(plotpow.rho(w3I.rho.long))
print(plotpow.N1(w3I.N1.long))
print(plotpow.alpha1.error(w3I.long))   
print(plotpow.alpha1.error(w3II.long))
print(plotpow.alpha1.error(w3I.worst.long))
print(plotpow.rho.error(w3I.rho.long))
print(plotpow.N1.error(w3I.N1.long))
print(plotpow.alpha1.selP(w3I.long.selP))  
print(plotpow.alpha1.selP(w3II.long.selP))
print(plotpow.alpha1.selP(w3I.worst.long.selP))
print(plotpow.rho.selP(w3I.rho.long.selP))
print(plotpow.N1.selP(w3I.N1.long.selP))
dev.off() 









###############################################



#Reproduce tables of Bias and confidence limit of concordance in Manuscript and supplementary material

w1I.Bias<-w1I[w1I$test==4,]



alph1<-rep(1:5,5)
o11<-cbind(w1I.Bias,alph1)
o13<-round(o11[o11$alph1==3,],3)

#True Concordance
o13$TrueC.Low
o13$TrueC.Me
o13$TrueC.Hi

#Table 4
#Bias inverse normal;
cbind(o13$biasinv.Low,o13$biasinv.Me,o13$Bias.Hi)
#CI inverse normal
cbind(o13$CImeaninvnorm.Low,o13$CImeaninvnorm.Me,o13$CImeaninvnorm.Hi)
#Coverage inverse normal;
cbind(o13$CICovinvnorm.Low,o13$CICovinvnorm.Me,o13$CICovinvnorm.Hi)*cbind(o13$SelP.Lo,o13$SelP.Me,c(1,1,1,1,1))+(1-cbind(o13$SelP.Lo,o13$SelP.Me,c(1,1,1,1,1)))



#SUPPL
#Bias unconditional;
cbind(o13$Bias.Low,o13$Bias.Me,o13$Bias.Hi)
#CI unconditional
cbind(o13$CIupC.Low,o13$CIupC.Me,o13$CIupC.Hi)
#Coverage unconditional;
cbind(o13$CICov.Low,o13$CICov.Me,o13$CICov.Hi)*cbind(o13$SelP.Lo,o13$SelP.Me,c(1,1,1,1,1))+(1-cbind(o13$SelP.Lo,o13$SelP.Me,c(1,1,1,1,1)))



#SUPPL
#Bias conditional;
cbind(o13$biasCcond.Low,o13$biasCcond.Me,o13$biasCcond.Hi)
#CI conditional
cbind(o13$CImeancond.Low,o13$CImeancond.Me,o13$CImeancond.Hi)
#Coverage conditional;
cbind(o13$CICovcond.Low,o13$CICovcond.Me,o13$CICovcond.Hi)*cbind(o13$SelP.Lo,o13$SelP.Me,c(1,1,1,1,1))+(1-cbind(o13$SelP.Lo,o13$SelP.Me,c(1,1,1,1,1)))