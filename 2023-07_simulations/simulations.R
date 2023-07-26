##########################
# eWHORM simulations
# July 2023
# Marta Bofill Roig
##########################

# Settings
set.seed(123)
mu=c(0,1,2,5)
N1 = 30*4

# install.packages("DescTools")
library(DescTools)

# Prelude - Auxiliary functions
# Function to get the column index of the maximum value in a row
get_max_col_index <- function(row) {
  return(which.max(row))
}

# Function to compute the hypotheses to test (closed test)
get_hyp_mat <- function(n_hypothesis = 3, selected_hypothesis = 1){
  elements <- c(rep(0, n_hypothesis), 1:n_hypothesis)
  
  hyp_mat <- unique(combinations(n=length(elements), r=n_hypothesis, v=elements, set = F, repeats.allowed = F))
  
  selected_rows <- hyp_mat[which(apply(hyp_mat == selected_hypothesis, 1, any)), ]
  
  selected_rows
}
# get_hyp_mat(3,2)



# Function to simulate trial data (1-stage, multiple arms)
sim_data <- function(n_arms, N, mu_6m, mu_12m, sd_y=0.1){ 

  treatments <- factor(sample(1:n_arms, N, replace = TRUE), 
                       levels = 1:n_arms, 
                       labels = c("Placebo", "Low", "Medium", "High")[1:n_arms])
  X <- model.matrix(~ treatments - 1)
  
  # simulate trial data
  y_6m = X %*% mu_6m + rnorm(N, mean = 0, sd = sd_y)
  y_12m = X %*% mu_12m + rnorm(N, mean = 0, sd = sd_y)
  
  # Treatment indicator for dataframe 
  max_col_indices <- apply(X, 1, get_max_col_index)
  treat=unname(max_col_indices)
  
  treat=factor(treat, levels = 1:n_arms,
               labels = c("Placebo", "Low", "Medium", "High")[1:n_arms]) 
  
  # Output
  data=data.frame(y_6m=y_6m, y_12m=y_12m, treat=treat)

  return(data)
}

# Example
db = sim_data(n_arms=4, N=30*4, mu_6m=mu, mu_12m=mu+c(0,1,1,2), sd_y=0.1)
model_6m = lm(y_6m~treat,data=db)
model_12 = lm(y_12m~treat,data=db)
summary(model_6m)
summary(model_12)


# Function to simulate trial data (2-stages, with dose selection)
sim_trial <- function(n_arms=4, N1=30*4, N2=30*2, mu=c(0,1,2,5), sd_y=0.1, alpha1=0.5){
  
  # stage1
  db_stage1 = sim_data(n_arms=n_arms, N=N1, mu_6m=mu, mu_12m=mu+c(0,1,1,2), sd_y=sd_y)
  
  
  res_stage1 = DunnettTest(x=db_stage1$y_6m, g=db_stage1$treat) 
  
  # selection
  
  #Worst case: No trend is seen in any of the doses (e.g., all p> alpha1): select highest dose
  if(sum(res_stage1$Placebo[,4]>alpha1)==3){ 
    sel=4
  }
  #Intermediate case: some doses show a trend: select the (highest) dose, no new recruitment for the other doses
  if(sum(res_stage1$Placebo[,4]<alpha1)<3){
    sel=which.min(res_stage1$Placebo[,4])+1
  }
  #Intermediate case 2: some doses show a trend: select the lowest effective dose, no new recruitment for the other doses
  if(sum(res_stage1$Placebo[,4]<alpha1)<3){
    sel=which.min(res_stage1$Placebo[,4])+1
  }
  # Best case
  if(sum(res_stage1$Placebo[,4]<alpha1)==3){
    sel=2
  }
  
  db_hyp= subset(db_stage1,db_stage1$treat==levels(db_stage1$treat)[c(1,2,4)])
  DunnettTest(x=db_hyp$y_6m, g=db_hyp$treat) 
  
  t.test(db_stage1)
  hyp_closedtest(sel=3)
  
  
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

# Function to decide the selected dose
select_dose <- function(){
  
}

select_dose_borrow <- function(){
  
}



