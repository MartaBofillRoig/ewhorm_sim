##########################
# eWHORM simulations
# Auxiliary functions
# July 2023
# Marta Bofill Roig
##########################

# install.packages("DescTools")
library(DescTools)
library(gtools)

# Auxiliary functions 
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
# Example
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
# db = sim_data(n_arms=4, N=30*4, mu_6m=mu, mu_12m=mu+c(0,1,1,2), sd_y=0.1) 
# model_6m = lm(y_6m~treat,data=db)
# model_12m = lm(y_12m~treat,data=db)
# summary(model_6m)
# summary(model_12m) 