
# Simulating efficacy and safety endpoints
# Test code

library(mvtnorm) 

n_arms = 4
ss_arm = 100
N=ss_arm*n_arms
mu=c(0,0,0,0)
mu1=c(0,0.5,1,2)
mu_12m=mu_6m=mu1

(sigma=matrix(c(1,.5,.5,
                .5,1,.9,
                .5,.9,1), nrow=3,byrow = T))
eigen(sigma)

treatments <- factor(c(sample(rep(1:n_arms, floor(N/n_arms))), sample(1:n_arms, N-floor(N/n_arms)*n_arms, replace=T)),
                     # sample(1:n_arms, N, replace = TRUE),
                     levels = 1:n_arms,
                     labels = c("Placebo", "Low", "Medium", "High")[1:n_arms])
X <- model.matrix(~ treatments - 1)
y <- X %*% matrix(c(mu_6m, mu_12m, mu), nrow=n_arms, byrow = F) + rmvnorm(n=N, mean = c(0,0, 0), sigma = sigma )

# Treatment indicator for dataframe
max_col_indices <- apply(X, 1, get_max_col_index)
treat = unname(max_col_indices)

treat = factor(treat, levels = 1:n_arms,
               labels = c("Placebo", "Low", "Medium", "High")[1:n_arms])

data = data.frame(y_6m=y[,1], y_12m=y[,2], safety=y[,3], treat=treat)

data$safety = ifelse(y[,3]<pnorm(.8),1,0)

sum(y[,3]<pnorm(.8))/N

sum(data$safety)/N
sum(subset(data,data$treat=="Placebo")$safety/ss_arm)
sum(subset(data,data$treat=="Low")$safety/ss_arm)
sum(subset(data,data$treat=="Medium")$safety/ss_arm)
sum(subset(data,data$treat=="High")$safety/ss_arm)

