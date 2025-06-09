#Sys.setenv("_R_CHECK_SYSTEM_CLOCK_" = 0)

# Create r package folder
# setwd("C:/Users/Marta/Dropbox/C5/GitHub/ewhorm_sim/rpkg)
# usethis::create_package("C:/Users/Marta/Dropbox/C5/GitHub/ewhorm_sim/rpkg")

# Copy in R folder the functions of the r package 
setwd("C:/Users/Marta/Dropbox/C5/GitHub/ewhorm_sim/rpkg")
devtools::document()
devtools::load_all()

# Build & check the package
devtools::build(pkg = "C:/Users/Marta/Dropbox/C5/GitHub/ewhorm_sim/rpkg", path = NULL, binary = FALSE, manual = TRUE, vignettes = TRUE)
devtools::check_built(path = "C:/Users/Marta/Dropbox/C5/GitHub/ewhorm_sim/rpkg", cran = TRUE, manual = TRUE, incoming = TRUE)
devtools::build_manual(pkg = "C:/Users/Marta/Dropbox/C5/GitHub/ewhorm_sim/rpkg", path = NULL)

#create vignette
# usethis::use_vignette("my-vignette")

# pkgdown::build_site(pkg = "C:/Users/mbofi/Dropbox/CeMSIIS/GitHub/ewhorm_sim/rpkg")

# https://www.r-bloggers.com/2017/08/building-a-website-with-pkgdown-a-short-guide/
# https://r-pkgs.org/vignettes.html

library(ewhorm)
db <- ewhorm::sim_data(n_arms = 4,
                 N = 30 * 4,
                 mu_6m = c(0,0,0,0),
                 mu_12m= c(0,0,0,0),
                 sigma=diag(1,2),
                 rmonth =12)

summary(db)

sim_trial(n_arms=4, 
          N1=120, N2=60, 
          mu_6m = c(0,0,0,0),
          mu_12m= c(0,0,0,0), 
          sigma=diag(1,2), 
          rmonth=12, 
          alpha1=0.5, 
          alpha=0.05, 
          p_safety=c(0.9,0.8,0.7), 
          safety=T)

# sim_trial(n_arms=4, N1=120, N2=60, mu_6m, mu_12m, sigma, rmonth, alpha1=0.5, alpha=0.05, p_safety=c(0.9,0.8,0.7), safety=T)
#' @examples
#' mu=c(0,0,0,0); sigma=matrix(c(0.1,0,0,0.1), nrow = 2, byrow = T)
#' y=sim_data(n_arms=4, N=30*4, mu_6m=mu, mu_12m=mu+c(0,1,1,2), sigma=sigma, rmonth=10)