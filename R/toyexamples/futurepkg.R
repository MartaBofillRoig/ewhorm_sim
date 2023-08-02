# install.packages("future")
library(future)
library(future.apply)
plan(multisession, workers = 1)
fun <- function(x, y) {
  return(x + y)
}

replicate_fun <- function(scenarios, nsim) {
  results <- vector("list", length = nrow(scenarios))
  for (i in 1:nrow(scenarios)) {
    x <- scenarios[i, 1]
    y <- scenarios[i, 2]
    results[[i]] <- replicate(nsim, fun(x, y))
  }
  return(results)
}

# Replace "set" with your actual matrix of scenarios
set <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2, byrow = TRUE)

# Define the number of replications you want
nsim <- 5

# Use future_lapply for parallel execution
future::plan(multisession,workers=1)
results <- future_lapply(set, replicate_fun, nsim = nsim)

scenarios=set
nsim=5
