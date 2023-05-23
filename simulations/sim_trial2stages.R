# Simulation function
# simulates a trial with three treatment arms (doses) and a control
# two stages

# ss_treat=100; ss_contr=100; ss_2=100; p_resp=c(0.1,0.2,0.4,0.6)

sim_trial_disease <- function(ss_treat, ss_contr, ss_2, p_resp){
  
  n_arms <- length(p_resp)  

  # p_hat <- c(sum(runif(ss_contr)<p_resp[1])/ss_contr)
  r_hat <- c(sum(runif(ss_contr)<p_resp[1]))
  
  # 
  m <- matrix(runif((n_arms-1)*ss_treat), ncol = n_arms-1)
  mbin <- matrix(NA, ncol = n_arms-1, nrow = ss_treat)
  for(i in 1:(n_arms-1)){
    mbin[,i] = m[,i]<p_resp[i+1]  
  }
  
  # p_hat <- c(p_hat,colSums(mbin)/ss_treat) 
  r_hat <- c(r_hat,colSums(mbin))
  
  m_res <- matrix(
    c(c(ss_contr,rep(ss_treat,3)) - r_hat, r_hat),
    byrow = T,
    nrow = 2,
    dimnames = list(resp = 0:1, dose = 0:(n_arms - 1))
  )

  pval1 = CochranArmitageTest(m_res, alternative = "one.sided")$p.value
  
  # selecting the most promising arm
  sel = which.max(r_hat[2:n_arms]) 
  
  # second stage
  r_hat2 <- c(sum(runif(ss_2/2)<p_resp[1]))
  r_hat2 <- c(r_hat2,sum(runif(ss_2/2)<p_resp[sel+1]))
  
  prop.test(x, n = rep(m2, 2), alternative = "less")$p.value
  # inverse normal combination test
  pval = 1 - pnorm(qnorm(1 - p1) / sqrt(2) + qnorm(1 - p2) / sqrt(2))
  
}



# sim_trial_disease(ss_treat=100,ss_contr=100,ss_2=100,p_resp=c(0.1,0.2,0.4,0.6))


