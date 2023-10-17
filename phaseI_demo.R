library(BOIN)
source("CFO_tox_utils.R")
source("ACFO_tox_utils.R")

phi <- 0.20
ncohort <- 20
cohortsize <- 3
init.level <- 1
init.dose <- 1

tau <- 3
accrual <- 6
tite.dist <- 2
accrual.dist <- 1

ndose <- 7

alp.prior=phi; bet.prior=1-phi
add.args <- list(alp.prior=phi, bet.prior=1-phi)

p.trues <- list()
p.trues[[1]] <- c(0.05, 0.20, 0.46, 0.50, 0.60, 0.70, 0.80)
p.trues[[2]] <- c(0.02, 0.05, 0.20, 0.28, 0.34, 0.40, 0.44)
p.trues[[3]] <- c(0.01, 0.05, 0.10, 0.20, 0.32, 0.50, 0.70) 
p.trues[[4]] <- c(0.01, 0.04, 0.07, 0.10, 0.50, 0.70, 0.90) 
p.trues[[5]] <- c(0.01, 0.05, 0.10, 0.14, 0.20, 0.26, 0.34) 
p.trues[[6]] <- c(0.01, 0.02, 0.03, 0.05, 0.20, 0.40, 0.50) 
p.trues[[7]] <- c(0.01, 0.04, 0.07, 0.10, 0.15, 0.20, 0.25)
p.trues[[8]] <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.20, 0.45) 


idx <- 5
p.true <- p.trues[[idx]]
ps <- p.trues[[idx]]

p.true <- c(0.01, 0.05, 0.10, 0.20, 0.32, 0.50, 0.70) 

