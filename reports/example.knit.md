---
title: "Examples"
format: pdf
editor: visual
---


## Examples R package

Install R Package


::: {.cell}

```{.r .cell-code}
# setwd("C:/Users/mbofi/Dropbox/CeMSIIS/GitHub/ewhorm_sim/rpkg")
install.packages('C:/Users/mbofi/Dropbox/CeMSIIS/GitHub/ewhorm_sim/ewhorm_0.1.tar.gz',repos=NULL)
library(ewhorm)
```
:::


## Simulate data


::: {.cell}

```{.r .cell-code}
simulated_data <- ewhorm::sim_data(n_arms = 4,
                                   N = 30 * 4,
                                   mu_6m = c(0,0,0,0),
                                   mu_12m= c(0,0,0,0),
                                   sigma=diag(1,2),
                                   rmonth =12)
head(simulated_data,10)
```

::: {.cell-output .cell-output-stdout}
```
          y_6m       y_12m   treat recruit_time
1  -0.53139934  0.35453665 Placebo            3
2   0.05556132 -0.07613965    High            4
3  -2.16353049 -0.46860407  Medium            9
4  -0.42822239  0.78068245     Low            4
5  -1.35028541 -0.43646163     Low            8
6   0.06814993  0.28519301     Low            9
7   0.27439431 -0.73811365     Low            1
8   0.16952060 -0.09384885 Placebo            1
9   1.05579968 -2.01543372     Low            5
10  1.74841645 -1.17738823  Medium            6
```
:::
:::


## Simulate trial


::: {.cell}

```{.r .cell-code}
result_trial <- sim_trial(n_arms=4, 
                          N1=120, N2=60, 
                          mu_6m = c(0,0,0,0),
                          mu_12m= c(0,0,0,0), 
                          sigma=diag(1,2), 
                          rmonth=12, 
                          alpha1=0.5, 
                          alpha=0.05, 
                          p_safety=c(0.9,0.8,0.7), 
                          safety=T)

(result_trial)
```

::: {.cell-output .cell-output-stdout}
```
$result1
[1] FALSE

$result2
Low-Placebo 
          2 

$safety
[1] TRUE TRUE TRUE
```
:::
:::

