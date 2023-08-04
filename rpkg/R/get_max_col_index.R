#' Function to get the column index of the maximum value in a row
#' @description Function to get the column index of the maximum value in a row
#'
#' @param row selected row 
#' @returns maximum value in a row
#' @export
#' @details eWHORM simulations
#' @author Marta Bofill Roig
#' 
# Function to get the column index of the maximum value in a row
get_max_col_index <- function(row) {
  return(which.max(row))
}