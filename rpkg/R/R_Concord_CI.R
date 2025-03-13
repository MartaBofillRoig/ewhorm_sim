#vorbereitende Programme für die Sims.


########
#########DAS hier brauch ich nicht fuers package!!!!!!!!!!!!!!!!! Steht in einzelenen Funktionen;
########


#m<-60; n<-60; mu=1

#true<-wilcox.test(rnorm(m*10000,mean=mu),rnorm(n*10000),alternative="greater")$statistic/(m*10000*n*10000)
#true

#x<-rnorm(m,mean=mu); y<-rnorm(n)

#theta.data<-wilcox.test(x,y,alternative="greater")$statistic/(m*n)


var2<-function(theta,m,n)  #modified Hanley–McNeil approach for variance  estimation  (Newcombe, 2006)
{
  mstar<-1/2*(m+n)-1
  nstar<-1/2*(m+n)-1
  theta*(1-theta)*(1+nstar*(1-theta)/(2-theta)+mstar*theta/(1+theta))/(m*n)
}

#var2(.5,m,n)

#confidence limits are given by thetahat +- zV2(theta):

CI.help<-function(theta,thetahat,m,n,alpha) #pooled variance
{
  #theta: confidence limit - to be found by uniroot function
  #thetahat: observed concordance value
  #m: group size placebo
  #n: group size treatment
  mstar<-1/2*(m+n)-1
  nstar<-1/2*(m+n)-1
  variance<-theta*(1-theta)*(1+nstar*(1-theta)/(2-theta)+mstar*theta/(1+theta))/(m*n)
  
  sqrt(variance)*qnorm(alpha)-theta+thetahat
}
CI(.6,.5,30,30,0.025)

CI2<-function(theta,thetahat,alpha,N1,N2) #gewichtete Varianz der stageweisen Varianz
{
  #theta: confidence limit
  #thetahat: observed concordance value
  #m: group size placebo
  #n: group size treatment
  N<-N1+N2
  var<-sqrt(N1/N)*var2(theta,N1/3,N1/3)+sqrt(N2/N)*var2(theta,N2/3,N2/3)
  #sqrt(var2(theta,m,m))*qnorm(alpha)-theta+thetahat
  sqrt(var)*qnorm(alpha)-theta+thetahat
}


#with uniroot, this CI is resolved for theta (CI limit)

#uniroot(CI, lower = .01, upper = .99, tol = .Machine$double.eps^0.25, maxiter = 10000,thetahat=theta.data,m=m,n=n,alpha=alpha)
#uniroot(CI, lower = .01, upper = .99, tol = .Machine$double.eps^0.25, maxiter = 10000,thetahat=.5,m=m,n=n,alpha=alpha)
#uniroot(CI, lower = .1, upper = .9, tol = .Machine$double.eps^0.25, maxiter = 10000,thetahat=.5,m=1000,n=1000)
#true

#ueber inverse normal
p<-c(.51,.51)
1-pnorm((sqrt(5/10)*qnorm(1-p[1])+sqrt(5/10)*qnorm(1-p[2]))/sqrt(2))
#1-pnorm(sum(qnorm(1-p))/sqrt(length(p)))

#CIinvn<-function(theta,thetahat1,thetahat2,m1,m2,alpha,N1,N2)
#{
#  N<-N1+N2
#  (sqrt(N1/N)*(thetahat1-theta)/sqrt(var2(theta,m1,m1))+sqrt((N2/N))*(thetahat2-theta)/sqrt(var2(theta,m2,m2)))-qnorm(1-alpha)
#  #(1-pnorm(sqrt(N1/N)*(thetahat1-theta)/sqrt(var2(theta,m1,m1))+sqrt((N2/N))*(thetahat2-theta)/sqrt(var2(theta,m2,m2))))-alpha#qnorm(1-alpha)
#}
CIinvn<-function(theta,thetahat1,thetahat2,m1,m2,value,N1,N2) #gleich wie CIinvn, hab hier statt alpha val eingegeben, kann Median berechnen.
{
  N<-N1+N2
  (sqrt(N1/N)*(thetahat1-theta)/sqrt(var2(theta,m1,m1))+sqrt((N2/N))*(thetahat2-theta)/sqrt(var2(theta,m2,m2)))-qnorm(1-value)
  #(1-pnorm(sqrt(N1/N)*(thetahat1-theta)/sqrt(var2(theta,m1,m1))+sqrt((N2/N))*(thetahat2-theta)/sqrt(var2(theta,m2,m2))))-alpha#qnorm(1-alpha)
}
CIinvn(.6,.5,.5,15,15,0.025,15,30)
CIinvn(.5,.5,.5,15,15,0.025,15,30)

uniroot(CI, lower = .01, upper = .99, tol = .Machine$double.eps^0.25, maxiter = 10000,thetahat=.6,m=30,n=30,alpha=0.025)
#uniroot(CIinvn, lower = .01, upper = .99, tol = .Machine$double.eps^0.25, maxiter = 10000,thetahat1=.7,thetahat2=.7,m1=15,m2=15,alpha=0.025,N1=15,N2=30)
uniroot(CIinvn, lower = .01, upper = .99, tol = .Machine$double.eps^0.25, maxiter = 10000,thetahat1=.7,thetahat2=.7,m1=15,m2=15,value=0.025,N1=15,N2=30)
