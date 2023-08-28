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
#' @importFrom gMCP doInterim
#' @importFrom gMCP BonferroniHolm
#' @export
#' @details eWHORM simulations
#' @author Marta Bofill Roig


# Function to simulate trial data (2-stages, with dose selection)
sim_trial_pce <- function(n_arms=4, N1=30*4, N2=30*2, mu_6m, mu_12m, sigma, rmonth, alpha1=0.1, alpha=0.05, p_safety=c(0.9,0.8,0.7), safety=T, promising=F){


  #######################################
  # stage1
  db_stage1 = sim_data(n_arms=n_arms-1, N=N1, mu_6m=mu_6m[1:n_arms-1], mu_12m=mu_12m[1:n_arms-1], sigma=sigma, rmonth=rmonth)
  recruit_time1 = max(db_stage1$recruit_time)

  model_aov = aov(y_6m ~ treat, db_stage1)
  model_dunnett = summary(glht(model = model_aov, linfct=mcp(treat="Dunnett"), alternative = "less"))
  pval_dunnett = model_dunnett$test$pvalues

  #######################################
  # decisions based on pvalues from Dunnett test at 6 month

  if(sum(pval_dunnett<alpha1)==2){sc=2}
  if(sum(pval_dunnett<alpha1)==1){sc=1;sel=3}
  if(sum(pval_dunnett<alpha1)==0){sc=0}

  #######################################
  # pvalues ttest 12 months

  pval <- c()

  for(j in 1:(n_arms-2)){
    sub1 = subset(db_stage1,(db_stage1$treat==levels(db_stage1$treat)[1])+(db_stage1$treat==levels(db_stage1$treat)[j+1])==1)
    mod1 = lm(y_12m ~ treat, sub1) #are we using this model or should we use individual models?
    res1 = summary(mod1)
    pval[j] <- pt(coef(res1)[2,3], mod1$df, lower.tail = FALSE)
  }
  z = qnorm(1-pval)

  #######################################
  # preplanning adaptive conditional error

  N=N1+N2
  graph_bh <- BonferroniHolm(3)

  # the package assumes that wj are equal for all j
  v <- c(1/2,1/2,0)
  z1 <- c(z,0)
  preplan <- doInterim(graph=graph_bh,z1=z1,v=v,alpha=0.025)
  # preplan@Aj
  # preplan@BJ

  #######################################
  # stage2
  # sc=1 --> Arm A or B continue to stage 2
  # sc=0 --> Arm A and B stop and do not continue to stage 2

  if(sc==2){
    db_stage2 = sim_data(n_arms=3, N=N2, mu_6m=mu_6m[c(1,2,3)], mu_12m=mu_12m[c(1,2,3)], sigma=sigma, rmonth=rmonth)
    levels(db_stage2$treat) = levels(db_stage1$treat)[c(1,2,3)]
    recruit_time2 = max(db_stage2$recruit_time)

    pval2 <- c()

    for(j in 1:2){
      sub2 = subset(db_stage2,(db_stage2$treat==levels(db_stage2$treat)[1])+(db_stage2$treat==levels(db_stage2$treat)[j+1])==1)
      mod2 = lm(y_12m ~ treat, sub2) #are we using this model or should we use individual models?
      res2 = summary(mod2)
      pval2[j] <- pt(coef(res2)[2,3], mod2$df, lower.tail = FALSE)
    }

    Avalues <- c(preplan@BJ[7]/2, #H123
                 preplan@BJ[6], #H12
                 preplan@BJ[5], #H13
                 preplan@BJ[3], #H23
                 preplan@BJ[2], #H2
                 preplan@BJ[4]  #H1
    )

    # pval2[1] <= Avalues[c(1,2,3,6)] #p1
    # pval2[2] <= Avalues[c(1,2,4,5)] #p2

    decision_stage2 = matrix(c(pval2[1] <= Avalues[c(1,2,3,6)],
                               pval2[2] <= Avalues[c(1,2,4,5)],
                               byrow = F, ncol = 2)

  }
  if(sc==1){
    db_stage2 = sim_data(n_arms=3, N=N2, mu_6m=mu_6m[c(1,sel,4)], mu_12m=mu_12m[c(1,sel,4)], sigma=sigma, rmonth=rmonth)
    levels(db_stage2$treat) = c(levels(db_stage1$treat)[c(1,sel)],"High")
    recruit_time2 = max(db_stage2$recruit_time)

    pval2 <- c()

    for(j in 1:2){
      sub2 = subset(db_stage2,(db_stage2$treat==levels(db_stage2$treat)[1])+(db_stage2$treat==levels(db_stage2$treat)[j+1])==1)
      mod2 = lm(y_12m ~ treat, sub2) #are we using this model or should we use individual models?
      res2 = summary(mod2)
      pval2[j] <- pt(coef(res2)[2,3], mod2$df, lower.tail = FALSE)
    }

    Avalues <- c(preplan@BJ[7]/2, #H123
                 preplan@BJ[6], #H12
                 preplan@BJ[5], #H13
                 preplan@BJ[3], #H23
                 preplan@BJ[2], #H2
                 preplan@BJ[1]  #H3
    )

    # pval2[1] <= Avalues[c(1,2,4,5)] #p2
    # pval2[2] <= Avalues[c(1,3,4,6)] #p3

    decision_stage2 = matrix(c(pval2[1] <= Avalues[c(1,2,4,5)],
                               pval2[2] <= Avalues[c(1,3,4,6)]),
                               byrow = F, ncol = 2)

  }
  if(sc==0){
    db_stage2 = sim_data(n_arms=2, N=N2, mu_6m=mu_6m[c(1,4)], mu_12m=mu_12m[c(1,4)], sigma=sigma, rmonth=rmonth)
    levels(db_stage2$treat) = c(levels(db_stage1$treat)[c(1)],"High")
    recruit_time2 = max(db_stage2$recruit_time)

    mod2 = lm(y_12m ~ treat, db_stage2) #are we using this model or should we use individual models?
    res2 = summary(mod2)
    pval2 <- pt(coef(res2)[2,3], mod2$df, lower.tail = FALSE)

    Avalues <- c(preplan@BJ[7], #H123
                 preplan@BJ[5], #H13
                 preplan@BJ[3], #H23
                 preplan@BJ[1]  #H3
                )

    decision_stage2 = matrix(pval2 <= Avalues, ncol = 1) #p3

  }

  #######################################
  # TO BE UPDATED --> what shall we report?
  return(list(pvalue_stage1=pval,
              pvalue_stage2=pval2,
              decision_stage1=(pval<alpha1),
              decision_stage2=decision_stage2,
              selected_dose=sel,
              recruit_time1=recruit_time1,
              recruit_time2=recruit_time2))
}

# library(multcomp);library(ewhorm);library(gMCP)
# mu = c(0,0,0,0)
# sg_m=matrix(c(1,.9,.9,1),nrow=2,byrow = T)
# N1=30*4; N2=30*2
# mu_6m <- mu; mu_12m <- mu
# rmonth=10
# n_arms=4
