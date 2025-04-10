#' Function to compute the hypotheses to test (closed test)
#' @description Function to compute the hypotheses to test (closed test)
#'
#' @param n_hypothesis num elementary hypotheses
#' @param selected_hypothesis selected hypothesis for closed test
#' @keywords internal
#' @returns maximum value in a row
#' @export
#' @details eWHORM simulations
#' @importFrom gtools combinations
#' @author Marta Bofill Roig
#'
get_hyp_mat <- function(n_hypothesis = 3, selected_hypothesis = 1){
  elements <- c(rep(0, n_hypothesis), 1:n_hypothesis)

  hyp_mat <- unique(combinations(n=length(elements), r=n_hypothesis, v=elements, set = F, repeats.allowed = F))

  selected_rows <- hyp_mat[which(apply(hyp_mat == selected_hypothesis, 1, any)), ]

  selected_rows
}
