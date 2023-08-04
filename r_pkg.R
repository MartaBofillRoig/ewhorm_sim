#Sys.setenv("_R_CHECK_SYSTEM_CLOCK_" = 0)

setwd("C:/Users/mbofi/Dropbox/CeMSIIS/GitHub/ewhorm_sim/rpkg")

# Create r package folder
usethis::create_package("C:/Users/mbofi/Dropbox/CeMSIIS/GitHub/ewhorm_sim/rpkg")

# Copy in R folder the functions of the r package
setwd("C:/Users/mbofi/Dropbox/CeMSIIS/GitHub/ewhorm_sim/rpkg")
devtools::document()
devtools::load_all()

# Build & check the package
devtools::build(pkg = "C:/Users/mbofi/Dropbox/CeMSIIS/GitHub/ewhorm_sim/rpkg", path = NULL, binary = FALSE, manual = TRUE, vignettes = TRUE)
devtools::check_built(path = "C:/Users/mbofi/Dropbox/CeMSIIS/GitHub/ewhorm_sim/rpkg", cran = TRUE, manual = TRUE, incoming = TRUE)
devtools::build_manual(pkg = "C:/Users/mbofi/Dropbox/CeMSIIS/GitHub/ewhorm_sim/rpkg", path = NULL)

#create vignette
# usethis::use_vignette("my-vignette")

# pkgdown::build_site(pkg = "C:/Users/mbofi/Dropbox/CeMSIIS/GitHub/ewhorm_sim/rpkg")

# https://www.r-bloggers.com/2017/08/building-a-website-with-pkgdown-a-short-guide/
# https://r-pkgs.org/vignettes.html
