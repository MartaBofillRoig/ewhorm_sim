library(gtools)
elements_1 <- c(0,1,2,3)

hyp_mat <- combinations(n=4, r=3, v=elements_1, set = T, repeats.allowed = T)

selected_rows <- hyp_mat[which(apply(hyp_mat, 1, function(x) sum(x==1)==1 & sum(x==2)<=1 & sum(x==3)<=1)), ]

selected_rows


library(utils)
hyp_closedtest <- function(sel){
  # Define the elements
  elements <- c(1, 2, 3)
  
  # Compute all combinations of 1, 2, and 3 elements 
  combinations_2 <- combn(elements, 2) 
  
  hyp=matrix(combinations_2,ncol=2,byrow=T)
  hyp=cbind(hyp,c(0,0,0))
  hyp=rbind(c(1,0,0),c(2,0,0),c(3,0,0),hyp,c(1,2,3))
  hyp
  
  # Find rows that include the number sel
  rows_with_sel <- which(apply(hyp == sel, 1, any))
  
  # Select rows with the number sel
  selected_rows <- hyp[rows_with_sel, ]
  
  return(selected_rows)
}

get_hyp_mat <- function(n_hypothesis = 3, selected_hypothesis = 1){
  elements <- c(rep(0, n_hypothesis), 1:n_hypothesis)
  
  hyp_mat <- unique(combinations(n=length(elements), r=n_hypothesis, v=elements, set = F, repeats.allowed = F))
  
  selected_rows <- hyp_mat[which(apply(hyp_mat == selected_hypothesis, 1, any)), ]
  
  selected_rows
}
# get_hyp_mat(3,2)