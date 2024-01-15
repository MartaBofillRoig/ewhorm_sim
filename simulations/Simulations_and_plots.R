
#oncho

mu<-20
sigma<-28



oo.oncho1<-mapply(simul_res, mu, sigma, c(0, 0, 0, 0),c(0,0,0,0),  c(0,0.5,.7,0,0.5,.7), 10000,4, 60 , 90, 1, c(.1,.1,.1,.5,.5,.5), .025, T,0, T,"t")


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