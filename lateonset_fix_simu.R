library(magrittr)
library(parallel)
library(survival)

source("./setting/Sce_5dose_phi33_mod.R")
source("utilities.R")
source("CFO_tox_utils.R")
source("ACFO_tox_utils.R")
source("lateonset_utils.R")

run.fn <- function(i){
  set.seed(seeds[i])

  if (i %% 100 == 0){
    print(i)
  }
  fcfo.res <- Simu.Fn(target, p.true, tau, cohortsize, ncohort, 
                      accrual, tite.dist, accrual.dist, design=1, impute.method="frac", add.args=add.args)
  titecfo.res <- Simu.Fn(target, p.true, tau, cohortsize, ncohort, 
                         accrual, tite.dist, accrual.dist, design=1, impute.method="TITE", add.args=add.args)
  facfo.res <- Simu.Fn(target, p.true, tau, cohortsize, ncohort, 
                       accrual, tite.dist, accrual.dist, design=2, impute.method="frac", add.args=add.args)
  titeacfo.res <- Simu.Fn(target, p.true, tau, cohortsize, ncohort, 
                          accrual, tite.dist, accrual.dist, design=2, impute.method="TITE", add.args=add.args)
  benchacfo.res <- Simu.Fn(target, p.true, tau, cohortsize, ncohort, 
                           accrual, tite.dist, accrual.dist, design=3, impute.method="No", add.args=add.args)
  ress <- list(
    fcfo=fcfo.res,
    titecfo=titecfo.res,
    facfo=facfo.res,
    titeacfo=titeacfo.res,  
    benchacfo=benchacfo.res,
    paras=list(p.true=p.true, 
               mtd=tmtd, 
               add.args=add.args,
               target=target,
               ncohort=ncohort,
               cohortsize=cohortsize)
  )
  ress
}

for (i in 1:length(p.trues)){
  tau <- 3
  accrual <- 6
  tite.dist <- 2
  accrual.dist <- 1
  init.dose <- 1
  
  add.args <- list(alp.prior=target, bet.prior=1-target, CV=0.95, suspend=F, crmCI.CV=0.80)
  
  nsimu <- 5000
  
  idx <- i
  p.true <- p.trues[[idx]]
  print(p.true)
  ndose <- length(p.true)
  tmtd <- MTD.level(target, p.true)
  
  seeds <- c(1:5000)
  
  t <- system.time({
    results <- mclapply(1:nsimu, run.fn, mc.cores=40)
  })
  print(t)
  file.name <- paste0("./results/","TITE_MTDSimu_",ndose, "dose_phi_", target*100, "_fix_",  idx, ".Rdata")
  save(results, file=file.name)
  
  print(post.process.random(results))
}
