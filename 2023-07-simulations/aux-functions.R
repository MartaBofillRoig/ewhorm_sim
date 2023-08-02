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

