#' Simulate data from a multi-arm multi-stage trial with shared control and dose selection
#' @description Function to simulate trial data (2-stages, with dose selection)
#'
#' @param n_arms number of arms (including control)
#' @param N1 sample size stage 1
#' @param N2 sample size stage 2
#' @param mu_6m 6-month mean response per arm (vector of length `n_arm`)
#' @param mu_12m 12-month mean response per arm (vector of length `n_arm`)
#' @param sigma covariance matrix between 6- and 12-month responses assumed equal across arms (matrix of dim 2x2)
#' @param alpha1 significance level for dose selection
#' @param alpha significance level for selected dose vs control comparison
#' @param rmonth recruitment per month (recruitment speed assumed constant over time)
#' @param p_safety probability of each dose to be safe
#' @param safety indicator - if true, it simulates safety according to  `p_safety`
#' @returns Combined p-value, selected dose and safety for each dose (if argument `safety=TRUE`)
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats runif
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @importFrom stats aov
#' @importFrom stats t.test
#' @importFrom multcomp glht
#' @importFrom multcomp mcp
#' @export
#' @details eWHORM simulations
#' @author Marta Bofill Roig


# Function to simulate trial data (2-stages, with dose selection)
sim_trial <- function(n_arms=4, N1=30*4, N2=30*2, mu_6m, mu_12m, sigma, rmonth, alpha1=0.5, alpha=0.05, p_safety=c(0.9,0.8,0.7), safety=T){

  # stage1
  db_stage1 = sim_data(n_arms=n_arms, N=N1, mu_6m=mu_6m, mu_12m=mu_12m, sigma=sigma, rmonth=rmonth)

  recruit_time1 = max(db_stage1$recruit_time)

  model_aov = aov(y_6m ~ treat, db_stage1)
  model_dunnet = summary(glht(model = model_aov, linfct=mcp(treat="Dunnett"), alternative = "less"))
  pval_dunnet = model_dunnet$test$pvalues

  # Safety
  if(safety==T){
    safety_dose1 <- (runif(1)<p_safety[1])
    safety_dose2 <- (runif(1)<p_safety[2])
    safety_dose3 <- (runif(1)<p_safety[3])
  }
  # if safety_dose(i) is true, then consider safe
  safety=c(safety_dose1,safety_dose2,safety_dose3)

  # Selection
  #--- Arms indicators: 1: Low; 2: Medium, 3:High

  #--- Worst case: No trend is seen in any of the doses (e.g., all p> alpha1): select highest dose
  if(sum(pval_dunnet>alpha1)==3){
    sel=3
  }

  #--- Intermediate case: some doses show a trend: select the lowest effective dose, no new recruitment for the other doses
  if(sum(pval_dunnet<alpha1)<3 & sum(pval_dunnet<alpha1)>=1){
    sel=which(pval_dunnet<alpha1)[1]
  }

  #--- Intermediate case 2: some doses show a trend: select the (highest) effective dose, no new recruitment for the other doses
  # if(sum(pval_dunnet<alpha1)<3){
  #   sel=which.min(pval_dunnet)+1
  # }

  #--- Best case: All doses show a trend: select the lowest effective dose, no new recruitment for the other doses
  if(sum(pval_dunnet<alpha1)==3){
    sel=1
  }

  sel <- sel+1
  hyp <- get_hyp_mat(3,sel-1)
  hyp <- hyp + (hyp != 0)
  # (hyp)

  sset_hyp1 <- subset(db_stage1,db_stage1$treat==levels(db_stage1$treat)[c(1,hyp[1,][hyp[1,] != 0])])
  sset_hyp2 <- subset(db_stage1,db_stage1$treat==levels(db_stage1$treat)[c(1,hyp[2,][hyp[2,] != 0])])
  sset_hyp3 <- subset(db_stage1,db_stage1$treat==levels(db_stage1$treat)[c(1,hyp[3,][hyp[3,] != 0])])
  sset_hyp4 <- subset(db_stage1,db_stage1$treat==levels(db_stage1$treat)[c(1,hyp[4,][hyp[4,] != 0])])


  #
  # pvalue_Dunnett1 <- min(DunnettTest(x=sset_hyp1$y_12m, g=sset_hyp1$treat)$Placebo[,4])
  # pvalue_Dunnett2 <- min(DunnettTest(x=sset_hyp2$y_12m, g=sset_hyp2$treat)$Placebo[,4])
  # pvalue_Dunnett3 <- min(DunnettTest(x=sset_hyp3$y_12m, g=sset_hyp3$treat)$Placebo[,4])
  # pvalue_Dunnett4 <- min(DunnettTest(x=sset_hyp4$y_12m, g=sset_hyp4$treat)$Placebo[,4])
  #
  pvalue_Dunnett1 = min(summary(glht(model = aov(y_12m ~ treat, sset_hyp1), linfct=mcp(treat="Dunnett"), alternative = "less"))$test$pvalues)
  pvalue_Dunnett2 = min(summary(glht(model = aov(y_12m ~ treat, sset_hyp2), linfct=mcp(treat="Dunnett"), alternative = "less"))$test$pvalues)
  pvalue_Dunnett3 = min(summary(glht(model = aov(y_12m ~ treat, sset_hyp3), linfct=mcp(treat="Dunnett"), alternative = "less"))$test$pvalues)
  pvalue_Dunnett4 = min(summary(glht(model = aov(y_12m ~ treat, sset_hyp4), linfct=mcp(treat="Dunnett"), alternative = "less"))$test$pvalues)

  #
  pvalue_stage1 <- max(pvalue_Dunnett1,pvalue_Dunnett2,pvalue_Dunnett3,pvalue_Dunnett4)
  wmax=which.max(c(pvalue_Dunnett1,pvalue_Dunnett2,pvalue_Dunnett3,pvalue_Dunnett4))
  if(wmax==4){
    N1s=4*30
    N2s=2*30
  }
  if(wmax==2 || wmax==3){
    N1s=3*30
    N2s=2*30
  }
  if(wmax==1){
    N1s=2*30
    N2s=2*30
  }

  # stage2
  db_stage2 = sim_data(n_arms=2, N=N2, mu_6m=mu_6m[c(1,sel)], mu_12m=mu_12m[c(1,sel)], sigma=sigma, rmonth=rmonth)
  levels(db_stage2$treat) = levels(db_stage1$treat)[c(1,sel)]
  pvalue_stage2 <- t.test(y_12m ~ treat, data = db_stage2, alternative = c("less"))$p.value

  recruit_time2 = max(db_stage2$recruit_time)

  # Inverse normal combination test
  combined_pvalue = 1 - pnorm(qnorm(1 - pvalue_stage1) / sqrt(N1s/(N1s+N2s)) + qnorm(1 - pvalue_stage2) / sqrt(N2s/(N1s+N2s)))

  # Fisher combination test
  # combined_test = -2*(log(pvalue_stage1)+log(pvalue_stage2))
  # combined_pvalue = 1- pchisq(combined_test, 2)

  return(list(combined_pvalue=combined_pvalue, selected_dose=sel, safety=safety,
              pvalue_stage1=pvalue_stage1, pvalue_stage2=pvalue_stage2,
              recruit_time1=recruit_time1, recruit_time2=recruit_time2))
}
