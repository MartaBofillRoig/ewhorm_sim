# install.packages("gMCP")
library(gMCP)
# graphGUI() 

graph_bh <- BonferroniHolm(3)

p1=c(.1,.12)
z1 <- c(qnorm(1-p1),0)
v <- c(1/2,1/2,0)
A_matrix <- doInterim(graph=graph_bh,z1=z1,v=v,alpha=0.025)
A_matrix@Aj

z_gamma <- qnorm(1-.025)
N=(N1+N2)/3
n1=N/2
w1=sqrt(n1/N)
w2=sqrt((N-n1)/N)

1-pnorm((z_gamma-w1*z1[2])/(w2))