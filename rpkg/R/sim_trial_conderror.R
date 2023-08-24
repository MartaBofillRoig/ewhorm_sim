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
  db_stage1 = sim_data(n_arms=n_arms-1, N=N1, mu_6m=mu_6m[1:n_arms-1], mu_12m=mu_12m[1:n_arms-1], sigma=sg_m, rmonth=rmonth)
  recruit_time1 = max(db_stage1$recruit_time)

  model_aov = aov(y_6m ~ treat, db_stage1)
  model_dunnet = summary(glht(model = model_aov, linfct=mcp(treat="Dunnett"), alternative = "less"))
  pval_dunnet = model_dunnet$test$pvalues

  #######################################
  pval <- c()

  for(j in 1:(n_arms-1)){
    sub1 = subset(db_stage1,(db_stage1$treat==levels(db_stage1$treat)[1])+(db_stage1$treat==levels(db_stage1$treat)[j+1])==1)
    mod1 = lm(y_6m ~ treat, sub1) #are we using this model or should we use individual models?
    res1 = summary(mod1)
    pval[j] <- pt(coef(res1)[2,3], mod1$df, lower.tail = FALSE)
  }
  z = qnorm(1-pval)


# --> PENDING: check if for works
  # sub1 = subset(db_stage1,(db_stage1$treat==levels(db_stage1$treat)[1])+(db_stage1$treat==levels(db_stage1$treat)[2])==1)
  # sub2 = subset(db_stage1,(db_stage1$treat==levels(db_stage1$treat)[1])+(db_stage1$treat==levels(db_stage1$treat)[3])==1)
  #
  # mod1 = lm(y_6m ~ treat, sub1) #are we using this model or should we use individual models?
  # res1 = summary(mod1)
  # pval[1] <- pt(coef(res1)[2,3], mod$df, lower.tail = FALSE)
  #
  # mod = lm(y_6m ~ treat, db_stage1)



  #######################################

  # preplan
  N=N1+N2

  db_stage1

  model=summary(lm(y_6m~treat,data=db_stage1))
  model$coefficients[2,3]

  graph_bh <- BonferroniHolm(3)


  # the package assumes that wj are equal for all j
  p1=c(.1,.12)
  z1 <- c(qnorm(1-p1),0)
  v <- c(1/2,1/2,0)
  A_matrix <- doInterim(graph=graph_bh,z1=z1,v=v,alpha=0.025)
  A_matrix@Aj

  # z_gamma <- qnorm(1-.025)
  # N=(N1+N2)/3
  # n1=N/2
  # w1=sqrt(n1/N)
  # w2=sqrt((N-n1)/N)
  #
  # 1-pnorm((z_gamma-w1*z1[2])/(w2))

  # define A --> pending

  #######################################

  # decisions based on pvalues ttest1 and ttest2 -->pending

  if(sum(pval<alpha1)==2){sc=3}
  if(sum(pval<alpha1)==1){sc=1;sel=3}
  if(sum(pval<alpha1)==0){sc=0}

  #######################################
  # stage2
  # sc=1 --> Arm A or B continue to stage 2
  # sc=0 --> Arm A and B stop and do not continue to stage 2

  if(sc==2){
    db_stage2 = sim_data(n_arms=3, N=N2, mu_6m=mu_6m[c(1,2,3)], mu_12m=mu_12m[c(1,2,3)], sigma=sg_m, rmonth=rmonth)
    levels(db_stage2$treat) = levels(db_stage1$treat)[c(1,2,3)]
  }
  if(sc==1){
    db_stage2 = sim_data(n_arms=3, N=N2, mu_6m=mu_6m[c(1,sel,4)], mu_12m=mu_12m[c(1,sel,4)], sigma=sg_m, rmonth=rmonth)
    levels(db_stage2$treat) = c(levels(db_stage1$treat)[c(1,sel)],"High")
  }
  if(sc==0){
    db_stage2 = sim_data(n_arms=2, N=N2, mu_6m=mu_6m[c(1,4)], mu_12m=mu_12m[c(1,4)], sigma=sigma, rmonth=rmonth)
    levels(db_stage2$treat) = levels(db_stage1$treat)[c(1,4)]

  }
  recruit_time2 = max(db_stage2$recruit_time)


  #######################################
  return(list(combined_pvalue=combined_pvalue, selected_dose=sel, safety=safety,
              pvalue_stage1=pvalue_stage1, pvalue_stage2=pvalue_stage2,
              recruit_time1=recruit_time1, recruit_time2=recruit_time2))
}

# library(multcomp);library(ewhorm)
# mu = c(0,0,0,0)
# sg_m=matrix(c(1,.9,.9,1),nrow=2,byrow = T)
# N1=30*4; N2=30*2
# mu_6m <- mu; mu_12m <- mu
# rmonth=10
# n_arms=4
