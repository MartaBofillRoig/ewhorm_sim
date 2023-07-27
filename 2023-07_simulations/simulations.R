##########################
# eWHORM simulations
# July 2023
# Marta Bofill Roig
##########################

setwd("C:/Users/mbofi/Dropbox/CeMSIIS/GitHub/ewhorm_sim/2023-07_simulations")
source("aux.R")

# Settings
set.seed(123)
mu=c(0,1,2,5)
N1 = 30*4

# # Example
# db = sim_data(n_arms=4, N=30*4, mu_6m=mu, mu_12m=mu+c(0,1,1,2), sd_y=0.1)
# model_6m = lm(y_6m~treat,data=db)
# model_12m = lm(y_12m~treat,data=db)
# summary(model_6m)
# summary(model_12m)

# Load the multcomp package
# library(multcomp)
# 
# g<-glht(mmm(model_6m,model_12m),
#      mlf(mcp(treat = "Dunnett")), alternative = "greater")

# var_cov_matrix <- vcov(g)

# Function to simulate trial data (2-stages, with dose selection)
sim_trial <- function(n_arms=4, N1=30*4, N2=30*2, mu=c(0,1,2,5), sd_y=0.1, alpha1=0.5){
  
  # stage1
  db_stage1 = sim_data(n_arms=n_arms, N=N1, mu_6m=mu, mu_12m=mu+c(0,1,1,2), sd_y=sd_y)
  
  
  res_stage1 = DunnettTest(x=db_stage1$y_6m, g=db_stage1$treat) 
  
  # selection
  #--- Worst case: No trend is seen in any of the doses (e.g., all p> alpha1): select highest dose
  if(sum(res_stage1$Placebo[,4]>alpha1)==3){ 
    sel=4
  }
  #--- Intermediate case: some doses show a trend: select the (highest) dose, no new recruitment for the other doses
  if(sum(res_stage1$Placebo[,4]<alpha1)<3){
    sel=which.min(res_stage1$Placebo[,4])+1
  }
  #--- Intermediate case 2: some doses show a trend: select the lowest effective dose, no new recruitment for the other doses
  if(sum(res_stage1$Placebo[,4]<alpha1)<3){
    sel=which.min(res_stage1$Placebo[,4])+1
  }
  #--- Best case
  if(sum(res_stage1$Placebo[,4]<alpha1)==3){
    sel=2
  }
  
  hyp <- get_hyp_mat(3,sel) 
  hyp <- hyp + (hyp != 0)
  
  sset_hyp1 <- subset(db_stage1,db_stage1$treat==levels(db_stage1$treat)[c(1,hyp[1,][hyp[1,] != 0])])
  sset_hyp2 <- subset(db_stage1,db_stage1$treat==levels(db_stage1$treat)[c(1,hyp[2,][hyp[2,] != 0])])
  sset_hyp3 <- subset(db_stage1,db_stage1$treat==levels(db_stage1$treat)[c(1,hyp[3,][hyp[3,] != 0])])
  sset_hyp4 <- subset(db_stage1,db_stage1$treat==levels(db_stage1$treat)[c(1,hyp[4,][hyp[4,] != 0])])
  
  # model_12m = lm(y_12m~treat,data=sset_hyp1)
  # mm=summary(model_12m)
  # mm$coefficients
  # result_anova1[[1]][[5]][[1]]
  
  pvalue_anova1 <- summary(aov(y_12m ~ treat, data = sset_hyp1))[[1]][[5]][[1]]
  pvalue_anova2 <- summary(aov(y_12m ~ treat, data = sset_hyp2))[[1]][[5]][[1]]
  pvalue_anova3 <- summary(aov(y_12m ~ treat, data = sset_hyp3))[[1]][[5]][[1]]
  pvalue_anova4 <- summary(aov(y_12m ~ treat, data = sset_hyp4))[[1]][[5]][[1]]
  
  # stage2
  db_stage2 = sim_data(n_arms=2, N=N2, mu_6m=mu[c(1,sel)], mu_12m=mu[c(1,sel)]+c(0,1,1,2)[c(1,sel)], sd_y=sd_y)
  levels(db_stage2$treat) = levels(db_stage1$treat)[c(1,sel)]
  
  list_res=list(db_stage1,db_stage2)
  
  return(list_res)
}

# Example
# n_arms=4; N1=30*4; N2=30*2; mu=c(0,1,2,5); sd_y=0.1
res = sim_trial(n_arms=4, N1=30*4, N2=30*2, mu=c(0,1,2,5), sd_y=0.1)
res



