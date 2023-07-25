##########################
# eWHORM simulations
# July 2023
# Marta Bofill Roig
##########################

# Prelude - Auxiliar functions
# Function to get the column index of the maximum value in a row
get_max_col_index <- function(row) {
  return(which.max(row))
}

# Settings
set.seed(123)
mu=c(0,1,2,5)
N1 = 30*4

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
sim_trial <- function(n_arms=4, N1=30*4, N2=30*2, mu=c(0,1,2,5), sd_y=0.1){
  
  # stage1
  db_stage1 = sim_data(n_arms=n_arms, N=N1, mu_6m=mu, mu_12m=mu+c(0,1,1,2), sd_y=sd_y)
  
  # stage2
  sel=sample(2:4,1)
  db_stage2 = sim_data(n_arms=2, N=N2, mu_6m=mu[c(1,sel)], mu_12m=mu[c(1,sel)]+c(0,1,1,2)[c(1,sel)], sd_y=sd_y)
  
  
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



