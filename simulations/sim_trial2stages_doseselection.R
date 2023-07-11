####################################################################
# Simulation function
# simulates a trial with three treatment arms (doses) and a control
# two stages
####################################################################
library(DescTools)
# ss_treat=100; ss_contr=100; ss_2=100; p_resp=c(0.1,0.2,0.4,0.6); alpha=0.05

sim_trial_disease <- function(ss_treat, ss_contr, ss_2, p_resp,alpha){
  
  n_arms <- length(p_resp)  

  r_hat <- c(sum(runif(ss_contr)<p_resp[1]))
  
  # 
  m <- matrix(runif((n_arms-1)*ss_treat), ncol = n_arms-1)
  mbin <- matrix(NA, ncol = n_arms-1, nrow = ss_treat)
  for(i in 1:(n_arms-1)){
    mbin[,i] = m[,i]<p_resp[i+1]  
  }
  
  r_hat <- c(r_hat,colSums(mbin))
  
  m_res <- matrix(
    c(c(ss_contr,rep(ss_treat,3)) - r_hat, r_hat),
    byrow = T,
    nrow = 2,
    dimnames = list(resp = 0:1, dose = 0:(n_arms - 1))
  )

  # trend test
  pval1 = CochranArmitageTest(m_res, alternative = "one.sided")$p.value
  
  if(pval1<alpha){
    
    # selecting the most promising arm
    sel = which.max(r_hat[2:n_arms]) 
    
    # CLOSED-TEST FIRST STAGE
    ########################## 
    db_list2 <- as.data.frame(combn(0:3,2, simplify = F))
    db_list2_tests <- which(db_list2[2,]==sel) 
    db_list2[,db_list2_tests]
    
    list3 <- combn(0:3,3, simplify = F)
    
    # ANALYSIS SECOND STAGE
    ##########################
    # evaluating the efficacy of the selected dose in stage 2
    r_hat2 <- c(sum(runif(ss_2/2)<p_resp[1]))
    r_hat2 <- c(r_hat2,sum(runif(ss_2/2)<p_resp[sel+1]))
    
    qval2 <- prop.test(r_hat2, n = rep(ss_2/2, 2), alternative = "less")$p.value
    
  }else{
    
  }
  
  
  
  
  # inverse normal combination test
  pval = 1 - pnorm(qnorm(1 - pval1) / sqrt(2) + qnorm(1 - pval2) / sqrt(2))
  
  return(pval)
}



sim_trial_disease(ss_treat=100,ss_contr=100,ss_2=100,p_resp=c(0.1,0.2,0.4,0.6))


