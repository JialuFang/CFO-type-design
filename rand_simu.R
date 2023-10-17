library(magrittr)
library(parallel)
source("utilities.R")
source("CFO_tox_utils.R")
source("CRM_tox_utils.R")
source("BOIN_tox_utils.R")
source("ACFO_tox_utils.R")
source("./setting/Sce_7dose_phi20.R")

ndose <- length(p.trues[[1]])
init.level <- 1

add.args <- list(alp.prior=target, bet.prior=1-target)

#mu <- 0.27
#phi <- 0.3
#ress <- lapply(1:5000, function(i)gen.rand.doses(8, phi, mu1=mu, mu2=mu))
#sapply(ress, prob.diff.fn, target=phi) %>% mean

## Taget = 0.2
# dose 7, mu1=mu2=0.36, 0.05
# dose 7, mu1=mu2=0.50, 0.07
# dose 7, mu1=mu2=0.65, 0.1
# dose 7, mu1=mu2=0.87, 0.15

# dose 10, mu1=mu2=0.39, 0.05
# dose 10, mu1=mu2=0.53, 0.07
# dose 10, mu1=mu2=0.68, 0.1
# dose 10, mu1=mu2=0.91, 0.15

## Target = 0.3
# dose 5, mu1=mu2=0.23, 0.05
# dose 5, mu1=mu2=0.38, 0.07 # no this one
# dose 5, mu1=mu2=0.53, 0.1
# dose 5, mu1=mu2=0.71, 0.15

# dose 7, mu1=mu2=0.27, 0.05
# dose 7, mu1=mu2=0.42, 0.07
# dose 7, mu1=mu2=0.56, 0.1
# dose 7, mu1=mu2=0.74, 0.15

# dose 8, mu1=mu2=0.28, 0.05
# dose 8, mu1=mu2=0.43, 0.07
# dose 8, mu1=mu2=0.57, 0.1
# dose 8, mu1=mu2=0.75, 0.15

# dose 10, mu1=mu2=0.30, 0.05
# dose 10, mu1=mu2=0.44, 0.07
# dose 10, mu1=mu2=0.58, 0.1
# dose 10, mu1=mu2=0.77, 0.15

## Target = 0.33
# dose 5, mu1=mu2=0.21, 0.05
# dose 5, mu1=mu2=0.36, 0.07 # no this one
# dose 5, mu1=mu2=0.52, 0.1
# dose 5, mu1=mu2=0.69, 0.15

diff_list <- c(0.05, 0.07, 0.1, 0.15)
mu_list <- c(0.36, 0.50, 0.65, 0.87)
for (i in 1:4){
  mu <- mu_list[i]
  diff <- diff_list[i]
  run.fn <- function(k){
    set.seed(seeds[k])
    if (k%%100 == 0){
      print(k)
    }
    
    p.true.all <- gen.rand.doses(ndose, target, mu1=mu, mu2=mu)
    p.true <- p.true.all$p.true
    tmtd <- p.true.all$mtd.level
    
    cfo.res <- CFO.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize,
                           init.level=init.level,  add.args=add.args)
    crm.res <- crm.simu.fn(target=target, p.true=p.true, 
                           init.level=init.level, cohortsize=cohortsize, ncohort=ncohort)
    boin.res <- boin.simu.fn(target=target, p.true=p.true, ncohort=ncohort, 
                             init.level=init.level, cohortsize=cohortsize)
    acfo.res <- ACFO.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize, 
                             init.level=init.level, add.args=add.args)
    
    ress <- list(
      cfo=cfo.res,
      crm=crm.res,
      boin=boin.res,
      acfo=acfo.res, 
      paras=list(p.true=p.true, 
                 mtd=tmtd, 
                 add.args=add.args,
                 target=target,
                 ncohort=ncohort,
                 cohortsize=cohortsize)
    )
    ress
  }
  
  nsimu <- 5000
  seeds <- 1:nsimu
  
  file.name <- paste0("./results/","MTDSimu_",ndose, "dose_phi_", target*100, "_random_",  diff*100, ".Rdata")
  
  t <- system.time({
    results <- mclapply(1:nsimu, run.fn, mc.cores=40)
  })
  print(t)
  
  post.process.random(results)
  #save(results, file=file.name)
}

